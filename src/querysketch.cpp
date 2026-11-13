#include <getopt.h>

#include <algorithm>
#include <deque>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>

#include "ntHashIterator.hpp"
#include "ntcard.hpp"
#include "QuerySketchUtil.hpp"
#include "TreeSketchUtil.hpp"
#include "GroupSketchUtil.hpp"
#include "NodeCarrierUtil.hpp"
#include "FileUtil.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "querysketch"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 2.0.0\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTIONS] [ARGS]\n"
    "Identify references and their abundance in FILE(S).\n\n"
    "Acceptable file formats: fastq in compressed formats gz, bz, zip, xz.\n"
    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N\tuse N parallel threads [1]\n"
    "  -b, --bit=N\t\tuse N bits per element in Bloom filter [64]\n"
    "  -o, --out=STRING\tthe output file name\n"
    "  -l, --r1=STRING\tinput read 1\n"
    "  -r, --r2=STRING\tinput read 2\n"
    "  -g, --ref=STRING\tinput uniqsketch reference\n"
    "  -h, --hit=N\t\tnumber of uniqsketch hits to call a reference [10]\n"
    "  -a, --acutoff=N\tabundance cutoff to report [0.0]\n"
    "  -s, --rcutoff=N\tread cutoff to report [2; 0 with --tree on a strict-margin index]\n"
    "      --sensitive\tsets sensitivity parameter h=10\n"
    "      --very-sensitive\tsets sensitivity parameter h=5\n"
    "      --solid\t\tonly use solid k-mers in reads\n"
    "      --breadth=F\tmin breadth of a call: distinct signatures hit, as a fraction of\n"
    "\t\t\twhat its hit count would spread over if the genome were present;\n"
    "\t\t\tbelow it the hits sit in one locus and the call is withheld\n"
    "\t\t\t(flat query; --tree-breadth governs --tree) [0.5; 0=off]\n"
    "      --tree=DIR\thierarchical query using a --tree index directory\n"
    "      --tree-prune=N\tmin coarse-level hits to keep a sub-tree in play (with --tree) [1]\n"
    "      --tree-breadth=F\tmin breadth of a call: distinct signatures hit, as a fraction of\n"
    "\t\t\twhat its hit count would spread over if the genome were present;\n"
    "\t\t\tbelow it the hits sit in one locus and the call is withheld (with --tree) [0.5]\n"
    "      --no-group-signatures\tleaf level uses globally unique signatures only (with --tree)\n"
    "      --help\t\tdisplay this help and exit\n"
    "      --version\t\tversion information and exit\n"
    "\n";

static const char shortopts[] = "t:g:l:r:o:h:a:s:b:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_BREADTH, OPT_TREE, OPT_TREEPRUNE, OPT_TREEBREADTH, OPT_NOGROUPSIGS };

// --sensitive / --very-sensitive / --solid are getopt flag entries, which store
// through a pointer at parse time and so cannot target a struct member. They write
// here and the values are copied into the config once parsing is complete. -h also
// routes through g_sketchHit so the flags and the explicit option share one source.
static int g_sketchHit = 10;
static int g_solid = 0;

static const struct option longopts[] = {
    {"threads",        required_argument, nullptr, 't'},
    {"out",            required_argument, nullptr, 'o'},
    {"r1",             required_argument, nullptr, 'l'},
    {"r2",             required_argument, nullptr, 'r'},
    {"ref",            required_argument, nullptr, 'g'},
    {"hit",            required_argument, nullptr, 'h'},
    {"acutoff",        required_argument, nullptr, 'a'},
    {"rcutoff",        required_argument, nullptr, 's'},
    {"bit",            required_argument, nullptr, 'b'},
    {"sensitive",      no_argument, &g_sketchHit, 10},
    {"very-sensitive", no_argument, &g_sketchHit, 5},
    {"solid",          no_argument, &g_solid, 1},
    {"breadth",        required_argument, nullptr, OPT_BREADTH},
    {"tree",           required_argument, nullptr, OPT_TREE},
    {"tree-prune",     required_argument, nullptr, OPT_TREEPRUNE},
    {"tree-breadth",   required_argument, nullptr, OPT_TREEBREADTH},
    {"no-group-signatures", no_argument,  nullptr, OPT_NOGROUPSIGS},
    {"help",           no_argument, nullptr, OPT_HELP},
    {"version",        no_argument, nullptr, OPT_VERSION},
    {nullptr, 0, nullptr, 0}
};

/**
 * Hierarchical (tree) query: narrow candidate references from coarse to fine.
 *
 * The coarsest level sketch is queried first; cluster representatives that clear the
 * hit/read thresholds mark their descendant leaves as candidates. Each finer level is
 * then queried with its sketch restricted to signatures of references that are still
 * candidates, so the search narrows down the tree. At the leaf level a leaf is in
 * play if the narrowing kept it or if its own signatures, global or within-group,
 * carry evidence of their own (a member need not resemble its representative; see
 * the leaf-level block), and its result is written.
 */
static void runTreeQuery(QuerySketchConfig& cfg, const std::string& treeDir,
                         const std::vector<std::string>& sampleFiles, bool readCutoffGiven) {
    std::string treeFile = treeDir + "/tree_index.tsv";
    if (!validateFile(treeFile, "tree index file")) exit(EXIT_FAILURE);
    if (!validateFile(cfg.ref, "leaf sketch (--ref)")) exit(EXIT_FAILURE);

    // The read cutoff exists to keep single-read noise out of the report. Whether a
    // single read is noise depends on how the index was built: margin-2 strict
    // signatures need two coincident sequencing errors to be hit by a read from
    // another reference, so one supporting read is real evidence and, at 0.01%
    // abundance, often all the evidence a sample holds. Margin-1 signatures are
    // hit by single-error reads of any abundant relative, and there the cutoff
    // must stay. The index records which case applies; an explicit --rcutoff wins.
    TreeBuildParams bp = loadTreeBuildParams(treeBuildParamsPath(treeDir));
    // Per leaf and ancestor node, the node signatures the leaf's genome carries
    // (ancestor_signatures.tsv). Lets a family row stand beside a resolved member
    // when the member cannot account for the node's hits; see
    // generateTreeQueryResult. Empty for an index without the table.
    const auto ancestorCarried = loadAncestorSignatures(ancestorSignaturesPath(treeDir));
    if (ancestorCarried.empty())
        std::cerr << "Tree query: index has no ancestor_signatures.tsv; a resolved member always "
                     "suppresses its family row (rebuild with the current uniqsketch to lift that)\n";
    // The two filter layouts select different signatures from the same pool, so
    // an index is meant to be queried by a binary built with the layout that
    // built it. The sketch is read verbatim either way; this is a warning only.
    if (bp.loaded && !bp.bloomFilter.empty() && bp.bloomFilter != bloomFilterName()) {
        std::cerr << "Tree query: index was built with the " << bp.bloomFilter
                  << " Bloom filter, this binary uses the " << bloomFilterName()
                  << "; results are valid but will not match an index built by this binary\n";
    }
    if (!readCutoffGiven) {
        if (bp.loaded && bp.minMargin >= 2 && bp.strictMargin) {
            cfg.readCutoff = 0;
            std::cerr << "Tree query: index built with --min-margin=2 --strict-margin; "
                         "single-read calls are reported (--rcutoff=0)\n";
        } else {
            std::cerr << "Tree query: index "
                      << (bp.loaded ? "was built without a strict margin" : "predates build_params.tsv")
                      << "; keeping --rcutoff=" << cfg.readCutoff
                      << " (rebuild with the current uniqsketch to enable single-read calls)\n";
        }
    }

    std::vector<TreeNode> tree = loadTree(treeFile);

    // level -> (representative basename -> covered leaf names), from the same
    // partition the index was built with (treeLevelPartition), so every leaf is
    // covered by exactly one representative at every level. The earlier version
    // derived this from node levels alone and lost every leaf that first appears
    // at a finer level (a never-clustered leaf, or a node whose parent is more
    // than one level up): such a leaf was never scored at the coarser level, so it
    // was absent from the candidate set and filtered out before it could be seen.
    // On a 120-genome E. coli tree that pruned every leaf before the leaf level.
    std::unordered_map<unsigned,
        std::unordered_map<std::string, std::unordered_set<std::string>>> levelMap;
    // Node name per (level, representative) for family rows, only for entries that
    // stand for two or more leaves: a leaf standing for itself at a coarse level is
    // the same genome the leaf level scores, not a family.
    std::unordered_map<unsigned,
        std::unordered_map<std::string, std::string>> levelNodeName;
    for (const auto& [level, entries] : treeLevelPartition(tree)) {
        for (const LevelEntry& e : entries) {
            std::unordered_set<std::string> leaves;
            for (const std::string& p : e.leafPaths) leaves.insert(getBaseId(p));
            levelMap[level][getBaseId(e.repPath)] = std::move(leaves);
            if (e.isNode && e.leafPaths.size() >= 2)
                levelNodeName[level][getBaseId(e.repPath)] = e.name;
        }
    }

    // Level sketches coarsest-first, then the leaf/full-resolution sketch (cfg.ref).
    std::vector<std::pair<unsigned, std::string>> levelSketches;
    for (unsigned lvl = 0; lvl < 100; lvl++) {
        std::string p = treeDir + "/sketch_level_" + std::to_string(lvl) + ".tsv";
        if (fileExists(p)) levelSketches.emplace_back(lvl, p);
    }
    unsigned leafLevel = 0;
    for (const TreeNode& node : tree)
        if (node.isLeaf) leafLevel = std::max(leafLevel, node.level);
    levelSketches.emplace_back(leafLevel, cfg.ref);
    std::sort(levelSketches.begin(), levelSketches.end());

    // Every level's sketch, coarsest first, into one index. Targets are (level,
    // reference) pairs; slots point at the per-signature counters in each level's
    // maps, so the pass writes straight into what the logs print. The sequence
    // maps are dropped once the index holds the sequences.
    struct Level {
        unsigned level;
        std::vector<std::string> ref;                   // refId -> name
        std::vector<SketchHash> sig;                    // per ref: signature -> hits; owns the index's sequences
        std::deque<std::string> extra;                  // sequences the maps do not hold (stable addresses)
        unsigned targetBase = 0;
        unsigned slotBegin = 0, slotEnd = 0;            // this level's slot range
    };
    std::vector<unsigned> slotRef;                      // per slot: the reference it belongs to
    auto now = [] { return std::chrono::steady_clock::now(); };
    auto secs = [](std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };
    auto t0 = now();
    std::vector<Level> levels(levelSketches.size());
    QueryIndex index;
    std::vector<unsigned*> slotCounter;

    // The level files are independent; parse them side by side (each on a copy
    // of the config, since loadSketch records k in it).
    std::vector<SketchHash> levelHash(levelSketches.size());
    // Carriers of the coarse-level signatures (NodeCarrierUtil.hpp): the seventh
    // column of a level sketch, absent from an index built without node carriers.
    std::vector<std::unordered_map<std::string, std::vector<unsigned>>> levelCarriers(levelSketches.size());
    #pragma omp parallel for schedule(dynamic) num_threads(static_cast<int>(std::min<size_t>(levelSketches.size(), cfg.threads)))
    for (size_t li = 0; li < levelSketches.size(); li++) {
        QuerySketchConfig local = cfg;
        levels[li].level = levelSketches[li].first;
        loadSketch(local, levelSketches[li].second, levelHash[li], levels[li].ref, levels[li].sig);
        if (li + 1 < levelSketches.size()) levelCarriers[li] = loadLevelCarriers(levelSketches[li].second);
    }
    {
        // k, as loadSketch would have recorded it: the length of the first signature.
        std::ifstream in(levelSketches.front().second);
        std::string first;
        std::getline(in, first);
        cfg.kmerLen = static_cast<unsigned>(first.find_first_of(" \t"));
    }

    // Node signatures with carriers: a coarse-level signature that genomes outside
    // its node also carry is its own target, so that after the pass it is added to
    // its node only when no carrier is present (see the admission before the
    // walk). An index without carriers has none of these and behaves as before.
    struct CarriedSig {
        size_t level;                 // index into `levels`
        unsigned ref;                 // the node's position at that level
        const std::string* seq;       // key in L.sig[ref], owned by the map
        std::vector<unsigned> carriers;
    };
    std::vector<CarriedSig> carriedSigs;
    for (size_t li = 0; li < levelSketches.size(); li++) {
        Level& L = levels[li];
        SketchHash& hash = levelHash[li];
        L.targetBase = index.numTargets;
        index.numTargets += static_cast<unsigned>(L.ref.size());
        L.slotBegin = static_cast<unsigned>(slotCounter.size());
        auto& carriers = levelCarriers[li];
        auto carriedBy = [&](const std::string& seq, unsigned ref, const std::string* key) {
            auto c = carriers.find(seq);
            if (c == carriers.end()) return false;
            carriedSigs.push_back({li, ref, key, std::move(c->second)});
            return true;
        };
        // Each signature's slot is its counter in its reference's map, whose key
        // also gives the index its sequence. When every line of the sketch is a
        // distinct k-mer (always, unless a file was edited by hand) the maps hold
        // exactly the signatures and are walked directly; otherwise the seq -> ref
        // map decides, as loadSketch's first-line-wins rule has it, and a
        // signature a map does not hold still counts for its reference.
        size_t mapped = 0;
        for (const SketchHash& m : L.sig) mapped += m.size();
        if (mapped == hash.size()) {
            for (unsigned i = 0; i < L.ref.size(); i++) {
                for (auto& kv : L.sig[i]) {
                    if (carriedBy(kv.first, i, &kv.first)) continue;
                    index.add(kv.first, L.targetBase + i, static_cast<unsigned>(slotCounter.size()));
                    slotCounter.push_back(&kv.second);
                    slotRef.push_back(i);
                }
            }
        } else {
            for (const auto& [seq, refId] : hash) {
                auto it = L.sig[refId].find(seq);
                if (it != L.sig[refId].end()) {
                    if (carriedBy(it->first, refId, &it->first)) continue;
                    index.add(it->first, L.targetBase + refId, static_cast<unsigned>(slotCounter.size()));
                    slotCounter.push_back(&it->second);
                    slotRef.push_back(refId);
                } else {
                    L.extra.push_back(seq);
                    index.add(L.extra.back(), L.targetBase + refId, NO_SLOT);
                }
            }
        }
        L.slotEnd = static_cast<unsigned>(slotCounter.size());
        carriers = {};
        // The index and the level's maps hold the sequences now; the leaf level's
        // map stays a little longer, to screen the group signatures against.
        if (li + 1 < levelSketches.size()) hash = SketchHash();
    }
    index.numSlots = static_cast<unsigned>(slotCounter.size());
    const SketchHash& leafHash = levelHash.back();

    // Within-group signatures (GroupSketchUtil.hpp) are each their own target,
    // so that after the pass the ones whose leaf is a candidate and none of whose
    // carriers is can be added to that leaf; the rest are discarded unseen. A
    // group signature that duplicates a leaf signature is dropped, as it was when
    // the leaf sketch was a map the group signatures were merged into.
    std::vector<GroupSignature> groupSigs;
    std::vector<char> groupDropped;
    unsigned groupBase = index.numTargets;
    if (!cfg.noGroupSignatures) {
        groupSigs = loadGroupSketch(groupSketchPath(treeDir));
        if (groupSigs.empty()) {
            std::cerr << "Tree query: index has no within-group signatures (built with "
                         "--no-group-signatures or by an older uniqsketch); leaf level uses "
                         "globally unique signatures only\n";
        }
        groupDropped.assign(groupSigs.size(), 0);
        for (unsigned g = 0; g < groupSigs.size(); g++) {
            if (leafHash.count(groupSigs[g].seq)) { groupDropped[g] = 1; continue; }
            index.add(groupSigs[g].seq, groupBase + g, NO_SLOT);
        }
        index.numTargets += static_cast<unsigned>(groupSigs.size());
    }
    const unsigned carriedBase = index.numTargets;
    for (unsigned j = 0; j < carriedSigs.size(); j++) index.add(*carriedSigs[j].seq, carriedBase + j, NO_SLOT);
    index.numTargets += static_cast<unsigned>(carriedSigs.size());
    levelHash.clear();
    index.finalize();

    auto t1 = now();

    // Size the --solid read filter from the sample. Without --solid nothing uses
    // the estimate, and the pass over the reads it costs is as long as the query.
    if (cfg.solid) {
        size_t totalKmers = 0;
        getCardinality(totalKmers, cfg.dbfSize, cfg.sbfSize, cfg.kmerLen, cfg.threads, sampleFiles);
        if (!resolveBloomSizes(totalKmers, cfg.dbfSize, cfg.sbfSize)) {
            std::cerr << "Error: insufficient k-mer content in the sample for query.\n";
            exit(EXIT_FAILURE);
        }
    }
    auto t2 = now();

    // The one pass over the reads.
    QueryTally tally(index);
    queryReads(cfg, sampleFiles, index, tally);
    auto t3 = now();

    // Candidate leaves. Every leaf starts in play; a leaf drops out only when a
    // representative covering it was scored at some level and fell below the
    // pruning floor. A leaf that no scored representative covers at a level (a
    // level file from an older index that omits self-represented leaves, say) is
    // simply carried forward, so a gap in the index can cost sensitivity in the
    // narrowing but can never delete a reference from the result.
    std::unordered_set<std::string> candidates;
    for (const TreeNode& node : tree)
        if (node.isLeaf) candidates.insert(getBaseId(node.filePath));
    // Leaves covered, at any level, by a representative with evidence of its
    // own, whether or not the narrowing had already dropped them. A member
    // carries its ancestors' signatures at some levels and not others (on the
    // 5 493-genome E. coli tree a leaf's level-0 node could be hit on 83 % of
    // its signatures and its level-1 node on none), so the narrowing, which
    // needs survival at every level, loses genomes the coarse levels did see.
    // The leaf level uses this as the mask's notion of "may be present".
    std::unordered_set<std::string> everInPlay;
    // Leaves covered, at some level, by a representative whose evidence stands
    // on its own: the reporting hit threshold, on two or more reads. This is
    // what lets a leaf join on a single read (see the leaf level). everInPlay
    // is too wide for that: one hit on one of a representative's thousands of
    // signatures is a stray that nearly every lineage collects from a few
    // million error-carrying reads of the dominant genome (on the 1 000-genome
    // sweep all four false single-read rows had such an ancestor; on the
    // 5 493-genome sweep 41 of 43), so gating on it changed nothing.
    std::unordered_set<std::string> lineageAttested;
    const size_t attestReads = std::max<size_t>(static_cast<size_t>(cfg.readCutoff), 1);
    std::vector<FamilyCall> familyCalls;         // coarse detections, for level-appropriate reporting

    // Carried node signatures (NodeCarrierUtil.hpp) are settled now, before the
    // walk. A hit on one is evidence that the node or one of the signature's
    // carriers is present. When every carrier is absent the signature joins its
    // node: its hits go into the node's tally and its signature map, so the walk,
    // the rejoin rule and the family row see them as ordinary evidence. When a
    // carrier is present, or may be, the signature is set aside, hits or not, and
    // leaves the node's signature set for this sample, so the node's breadth is
    // judged on the signatures that can speak for it (kept in the total, the
    // set-aside hits of a node whose set is mostly carried would read as misses
    // at high coverage). A node whose carried signatures are set aside is still
    // reached through its own node-unique ones; its parent, which contains the
    // node and its sibling carriers alike, is where the ambiguous evidence lands
    // (the parent's own signatures see the same reads), and that is the honest
    // row: on the 5 493-genome E. coli index a strain at 1 % that is a sibling of a
    // 21-leaf node and carries 1 715 of its 1 812 carried signatures put a family
    // row on that node, which it is not in, when the rule below judged it absent
    // on its own signatures alone.
    //
    // A carrier is *present* when anything of its own was hit: one of its leaf
    // signatures, or one of its within-group signatures none of whose own carriers
    // may be present (such a hit is evidence for it or for an unseen genome; the
    // strain above had 42 on two reads while its 72 leaf signatures, packed into
    // three sites, had none); or when a representative covering it at some level
    // that does not also cover the node clears the reporting threshold on two or
    // more reads, on everything hit under it, node-unique and carried alike (a
    // common ancestor says nothing about which of the two is present; the K-12
    // node's carriers are the K-12 relatives outside it, all under one parent).
    // An unseen carrier is *excluded* when, at the coverage the node's hits imply
    // for it (its hits per carried signature it carries), its own signature loci
    // would have collected cfg.carrierExcludeReads reads or more: the number of
    // loci, not of signatures, since a strain a few substitutions from its
    // relatives has all its signatures at those sites and is seen only when a read
    // covers one of them without error. An unseen carrier that cannot be excluded
    // *may be present* and blocks what it carries. One with no signature of any
    // kind cannot be seen or excluded: it is absent when a lineage of its own
    // (a representative covering it and not the node) was silent, and may be
    // present when that lineage has any trace or when it has no lineage apart
    // from the node's, a signature-less leaf beside the node under one parent,
    // for which nothing in the index can speak. All of this is read from the
    // tallies before any admission, so the outcome does not depend on the order
    // the nodes are taken in. Measured on the 5 493-genome index: pure
    // K-12 MG1655 or 4A at 10x, an empty table before, is the 345-leaf level-1
    // family row (the K-12 level-3 node's carried candidates all fail the strict
    // margin, and the level-2 node's carriers are K-12-like genomes in sibling
    // nodes that the same reads light); the 1 %-sibling case above gives the
    // 29-leaf level-2 row the strain belongs to.
    if (!carriedSigs.empty()) {
        const size_t leafIdx = levels.size() - 1;
        const Level& leafL = levels[leafIdx];
        const size_t nLeaves = leafL.ref.size();
        const unsigned hitFloor = static_cast<unsigned>(cfg.sketchHit);
        std::vector<std::string> leafName = loadTreeLeaves(treeLeavesPath(treeDir));
        if (leafName.size() != nLeaves) leafName = leafL.ref;
        // Per node, its carried signatures. A node's evidence for the tests below
        // is everything hit under it, node-unique and carried alike, before any
        // admission (the wide reading: a carrier whose own node is lit by any of
        // its signatures counts as attested, which errs towards setting aside).
        std::map<std::pair<size_t, unsigned>, std::vector<unsigned>> byNode;
        for (unsigned j = 0; j < carriedSigs.size(); j++)
            byNode[{carriedSigs[j].level, carriedSigs[j].ref}].push_back(j);
        // Per coarse level: leaf name -> covering representative; which
        // representatives have evidence at the reporting standard; and the leaves
        // under a representative with any evidence at all (the everInPlay notion
        // of the leaf level, for a carrier with no signatures of its own).
        std::vector<std::unordered_map<std::string, unsigned>> cover(leafIdx);
        std::vector<std::vector<char>> nodeEvidence(leafIdx), nodeTrace(leafIdx);
        std::unordered_set<std::string> lineageTrace;
        for (size_t li = 0; li < leafIdx; li++) {
            const Level& L = levels[li];
            const auto& lm = levelMap[L.level];
            nodeEvidence[li].assign(L.ref.size(), 0);
            nodeTrace[li].assign(L.ref.size(), 0);
            for (unsigned i = 0; i < L.ref.size(); i++) {
                auto it = lm.find(L.ref[i]);
                if (it != lm.end()) for (const auto& lf : it->second) cover[li][lf] = i;
                unsigned c = tally.count[L.targetBase + i];
                std::unordered_set<std::string> reads = tally.reads[L.targetBase + i];
                auto bn = byNode.find({li, i});
                if (bn != byNode.end())
                    for (unsigned j : bn->second) {
                        c += tally.count[carriedBase + j];
                        reads.insert(tally.reads[carriedBase + j].begin(), tally.reads[carriedBase + j].end());
                    }
                if (cfg.solid && c) ++c;
                nodeEvidence[li][i] = (c >= hitFloor && reads.size() > attestReads);
                nodeTrace[li][i] = (c >= static_cast<unsigned>(cfg.treePruneHit) && !reads.empty());
                if (nodeTrace[li][i] && it != lm.end())
                    for (const auto& lf : it->second) lineageTrace.insert(lf);
            }
        }
        // Seen: a hit on a leaf signature; then a hit on a within-group signature
        // none of whose carriers may be present, where a carrier may be present
        // when it has a leaf hit or, having no signatures of its own to be hit,
        // when its lineage has any trace (a 99 % main with no signature of any
        // kind, 266 such leaves on the 5 493-genome index, carries the within-group
        // signatures of a thousand leaves and would otherwise make them all seen).
        // Own signatures, leaf and within-group, with their loci, for the
        // exclusion test.
        std::vector<char> leafSeen(nLeaves, 0), seen(nLeaves, 0), mayBe(nLeaves, 0);
        std::vector<unsigned> ownSigs(nLeaves, 0);
        for (unsigned c = 0; c < nLeaves; c++) {
            leafSeen[c] = seen[c] = tally.count[leafL.targetBase + c] > 0;
            ownSigs[c] = static_cast<unsigned>(leafL.sig[c].size());
        }
        for (unsigned g = 0; g < groupSigs.size(); g++)
            if (!groupDropped[g] && groupSigs[g].refId < nLeaves) ownSigs[groupSigs[g].refId]++;
        for (unsigned c = 0; c < nLeaves; c++)
            mayBe[c] = leafSeen[c] || (ownSigs[c] == 0 && lineageTrace.count(leafName[c]) > 0);
        std::vector<std::tuple<unsigned, std::string, unsigned>> groupPos;
        groupPos.reserve(groupSigs.size());
        for (unsigned g = 0; g < groupSigs.size(); g++) {
            const GroupSignature& s = groupSigs[g];
            if (groupDropped[g] || s.refId >= nLeaves) continue;
            groupPos.emplace_back(s.refId, s.contig, s.pos);
            if (seen[s.refId] || tally.count[groupBase + g] == 0) continue;
            bool ambiguous = false;
            for (unsigned c : s.carriers) if (c < nLeaves && mayBe[c]) { ambiguous = true; break; }
            if (!ambiguous) seen[s.refId] = 1;
        }
        const std::vector<unsigned> loci = countSignatureLoci(cfg.ref, nLeaves, groupPos, cfg.locusGap);
        groupPos = {};
        size_t leavesLeafSeen = 0, leavesSeen = 0;
        for (unsigned c = 0; c < nLeaves; c++) { leavesLeafSeen += leafSeen[c]; leavesSeen += seen[c]; }
        // What a carrier's own lineage says about it, apart from the ancestors it
        // shares with the node: ATTESTED when a representative covering it, and not
        // the node, clears the reporting threshold on two or more reads; TRACE when
        // one has any hit at the pruning floor; SILENT when it has such
        // representatives and none was hit; NONE when every representative covering
        // it also covers the node, so nothing in the index can speak for it alone.
        enum Lineage : char { NONE, SILENT, TRACE, ATTESTED };
        auto lineage = [&](unsigned c, const std::string& nodeRep) {
            Lineage best = NONE;
            for (size_t li = 0; li < leafIdx; li++) {
                auto it = cover[li].find(leafName[c]);
                if (it == cover[li].end()) continue;
                auto nt = cover[li].find(nodeRep);
                if (nt != cover[li].end() && nt->second == it->second) continue;   // common ancestor
                if (nodeEvidence[li][it->second]) return ATTESTED;
                if (nodeTrace[li][it->second]) best = TRACE;
                else if (best == NONE) best = SILENT;
            }
            return best;
        };
        // Per node, per carrier: the signatures it carries and the hits on them,
        // which give the coverage the carrier would have were it the one present.
        // A carrier with no signature of any kind cannot be seen or excluded; it is
        // absent when its own lineage was silent, and may be present when its
        // lineage has a trace or when it has no lineage apart from the node's. The
        // first version treated it as absent in both cases, and on the 928-genome
        // mock community a signature-less strain at 10x, a direct child of the
        // level-2 node its 20-leaf sibling node hangs off, put a 78-read family row
        // on that sibling node, which it is not in, through the 45 filled
        // signatures it carries; the parent's row, which covers it, was then
        // suppressed as the coarser of two nested calls.
        enum Status : char { ABSENT, PRESENT, POSSIBLE };
        size_t used = 0, hit = 0, asidePresent = 0, asidePossible = 0, asideHit = 0;
        size_t carriersSeen = 0, carriersExcluded = 0, carriersPossible = 0, carriersUntestable = 0,
               carriersUntestablePossible = 0;
        std::set<std::pair<size_t, unsigned>> nodesUsing;
        for (const auto& [node, sigs] : byNode) {
            Level& L = levels[node.first];
            const std::string& nodeRep = L.ref[node.second];
            std::unordered_map<unsigned, std::pair<unsigned, unsigned>> perCarrier;   // carrier -> (carried, hits)
            bool anyHit = false;
            for (unsigned j : sigs) {
                const unsigned h = tally.count[carriedBase + j];
                anyHit |= h > 0;
                for (unsigned c : carriedSigs[j].carriers) {
                    if (c >= nLeaves) continue;
                    auto& e = perCarrier[c];
                    e.first++;
                    e.second += h;
                }
            }
            std::unordered_map<unsigned, Status> status;
            for (const auto& [c, e] : perCarrier) {
                Status st;
                const Lineage ln = seen[c] ? ATTESTED : lineage(c, nodeRep);
                if (ln == ATTESTED) { st = PRESENT; carriersSeen += anyHit; }
                else if (!anyHit) { st = ABSENT; }     // nothing to explain; the node's signatures are idle
                else if (ownSigs[c] == 0) {
                    st = (ln == SILENT) ? ABSENT : POSSIBLE;
                    ++carriersUntestable;
                    carriersUntestablePossible += st == POSSIBLE;
                } else {
                    const double perSig = static_cast<double>(e.second) / static_cast<double>(e.first);
                    const double expectedReads = perSig * static_cast<double>(loci[c]);
                    if (expectedReads >= cfg.carrierExcludeReads) { st = ABSENT; ++carriersExcluded; }
                    else { st = POSSIBLE; ++carriersPossible; }
                }
                status[c] = st;
            }
            for (unsigned j : sigs) {
                const CarriedSig& s = carriedSigs[j];
                const unsigned t = carriedBase + j;
                Status worst = ABSENT;
                for (unsigned c : s.carriers) {
                    if (c >= nLeaves) continue;
                    const Status st = status[c];
                    if (st == PRESENT) { worst = PRESENT; break; }
                    if (st == POSSIBLE) worst = POSSIBLE;
                }
                auto it = L.sig[s.ref].find(*s.seq);
                if (worst != ABSENT) {
                    if (worst == PRESENT) ++asidePresent; else ++asidePossible;
                    asideHit += tally.count[t] > 0;
                    if (it != L.sig[s.ref].end()) L.sig[s.ref].erase(it);   // s.seq is not read again
                    continue;
                }
                ++used;
                if (tally.count[t] == 0) continue;
                ++hit;
                const unsigned nt = L.targetBase + s.ref;
                tally.count[nt] += tally.count[t];
                tally.reads[nt].insert(tally.reads[t].begin(), tally.reads[t].end());
                if (it != L.sig[s.ref].end()) it->second += tally.count[t];
                nodesUsing.insert(node);
            }
        }
        std::cerr << "Tree query: carried node signatures: " << carriedSigs.size() << " in the index, "
                  << used << " in use (" << hit << " with hits, for " << nodesUsing.size() << " node(s)), "
                  << asidePresent + asidePossible << " set aside with their hits (" << asideHit
                  << " had hits): " << asidePresent << " for a carrier seen or attested by its lineage, "
                  << asidePossible << " for a carrier that may be present; carriers of nodes with hits: "
                  << carriersSeen << " seen or attested, " << carriersExcluded << " excluded on their own loci, "
                  << carriersPossible << " possible, " << carriersUntestable << " without signatures ("
                  << carriersUntestablePossible << " of them possible on their lineage); "
                  << leavesSeen << " of " << nLeaves << " leaves seen (" << leavesLeafSeen
                  << " on their leaf signatures)\n";
    }

    // Walk the levels coarse to fine over the tallies. At a coarse level only
    // references that still cover a candidate leaf are scored; the others'
    // tallies are discarded, as their signatures used to be left out of that
    // level's pass. The leaf level adds back any leaf with evidence of its own.
    for (size_t li = 0; li < levels.size(); li++) {
        Level& L = levels[li];
        const unsigned level = L.level;
        const bool isLeaf = (li + 1 == levels.size());
        const size_t n = L.ref.size();

        std::vector<bool> scored(n, false);
        std::vector<bool> present;                      // leaf level: what the within-group mask treats as present
        size_t rejoined = 0;                            // coarse level: nodes scored on their own evidence after being pruned above
        if (!isLeaf) {
            const auto& lm = levelMap[level];
            const auto& nodeName = levelNodeName[level];
            for (unsigned i = 0; i < n; i++) {
                auto it = lm.find(L.ref[i]);
                if (it == lm.end()) continue;
                for (const auto& lf : it->second)
                    if (candidates.count(lf)) { scored[i] = true; break; }
                // Evidence at this level, narrowed or not (the tallies cover every
                // representative; `scored` only decides which are reported from).
                unsigned c = tally.count[L.targetBase + i];
                if (cfg.solid && c) ++c;
                const size_t r = tally.reads[L.targetBase + i].size();
                if (c >= static_cast<unsigned>(cfg.treePruneHit) && r > 0)
                    for (const auto& lf : it->second) everInPlay.insert(lf);
                if (c >= static_cast<unsigned>(cfg.sketchHit) && r > attestReads)
                    for (const auto& lf : it->second) lineageAttested.insert(lf);
                // A node whose leaves were all pruned at a coarser level is scored
                // anyway when its own representative clears the reporting threshold
                // on two or more reads. The narrowing assumes a member carries its
                // ancestors' signatures, and a node's signatures are its
                // representative's k-mers absent outside the node, in a dense set
                // its private accessory sequence; on the 5 493-genome E. coli index
                // 10x of K-12 MG1655, alone in the sample, hit none of the 2 413
                // level-0 signatures of the 896-leaf node that holds it, the whole
                // lineage was gone at the top level, and the query reported nothing
                // at all. The leaf level's rescue by own evidence (07eed42) cannot
                // reach a leaf that has no signatures, and a cluster of K-12
                // derivatives one to three substitutions apart has none at margin 2
                // (928-genome mock communities: five such strains invisible in all
                // 18 of their occurrences, including at 10x). The evidence is in
                // hand: the tallies cover every representative at every level, and
                // everInPlay already reads them past the narrowing. Such a node is
                // reported as a family row, under the same suppression as any other
                // (no row when a member is resolved; the finest nested call wins).
                // The candidate set is not touched, so strain rows and the carrier
                // mask are exactly as before; only family rows can be added. Two
                // reads, not the reporting read floor: a node pruned above has no
                // corroboration from its ancestors, and a single read with two
                // coincident errors reproduces runs of 10 to 38 consecutive
                // signature k-mers of an unrelated genome (the 5 493-genome sweep's
                // false strain rows were all of that kind), so the standard is the
                // one lineageAttested uses for a single-read leaf.
                if (!scored[i] && c >= static_cast<unsigned>(cfg.sketchHit) && r > attestReads &&
                    nodeName.count(L.ref[i]) > 0) {
                    scored[i] = true;
                    ++rejoined;
                }
            }
        } else {
            // A leaf is in play if the narrowing kept it, or if its own signatures,
            // global or within-group, clear the floor. The set serves twice: it is
            // what the leaf level scores and reports, and it is what the
            // within-group mask below checks a signature's carriers against.
            //
            // The narrowing alone assumed two things that hold on a hundred
            // genomes and fail on thousands (docs/v2.0-validation.md,
            // "Scalability"). First, that a member carries its ancestors'
            // signatures: a node's signatures are its representative's k-mers
            // absent outside the node, and in a dense set those are the
            // representative's private accessory sequence, so a member's reads
            // hit none of them and the member was pruned with its own signatures
            // never consulted (44 % of 5 493 E. coli leaves had such an ancestor;
            // a 90 % main vanished from the output). Hence the rescue by own
            // evidence. Second, that a genome present in the sample is a scored
            // candidate: a leaf with no global signatures has no line in the leaf
            // sketch, so it could not be named here, was never scored, and never
            // counted as the carrier of a within-group signature; its reads went
            // into its relatives' signatures (222 802 of 223 003 false hits in
            // that benchmark). Hence the names from leaves.tsv: such a leaf is in
            // play when its lineage is, which is the only evidence it can leave,
            // and then blocks the signatures it carries.
            //
            // A carrier in play by lineage alone, with no hit on its own
            // signatures, still blocks. Treating it as absent (tried first) let a
            // 0.01 % contaminant's few reads, too few to hit its own signatures,
            // land on the within-group signatures of its siblings and produce 49
            // false single-read rows on the 120-genome sweep that had none. And
            // for the mask, lineage means evidence at any level (everInPlay), not
            // survival at all of them: a genome with no global signatures that was
            // pruned at one level is otherwise invisible, and on the 5 493-genome
            // index a 99.99 % main in that position put 976 hits into 20
            // within-group signatures of the 0.01 % contaminant, 10 of which it
            // carries verbatim, inflating that contaminant's count fifty-fold.
            std::vector<std::string> leafName = loadTreeLeaves(treeLeavesPath(treeDir));
            if (leafName.size() != n) {
                if (!leafName.empty())
                    std::cerr << "Tree query: " << treeLeavesPath(treeDir) << " lists " << leafName.size()
                              << " leaves, the leaf sketch " << n << "; ignoring it\n";
                leafName = L.ref;
                for (const GroupSignature& s : groupSigs)
                    if (s.refId < n && leafName[s.refId].empty()) leafName[s.refId] = s.refName;
            }
            // The reporter is handed L.ref, which loadSketch fills from the leaf
            // sketch, so a leaf with no global signatures has no name there. Such a
            // leaf is still reported when its within-group signatures carry the
            // sample (928 E. coli, a strain with 0 global and 72 within-group
            // signatures, 8 reads): the row came out with an empty name, and the
            // family-row suppression, which looks a resolved leaf up by name,
            // missed it and reported the leaf's family alongside. Name it here.
            for (unsigned i = 0; i < n; i++)
                if (L.ref[i].empty()) L.ref[i] = leafName[i];
            // The floor for joining by own evidence is the reporting threshold, not
            // the pruning floor. A leaf below it cannot be reported, so bringing it
            // in gains nothing, and as a carrier it would withhold its relatives'
            // within-group signatures: with the pruning floor (tried first), 268 of
            // 300 rescues on the 120-genome sweep were one read on one signature,
            // the single-error noise the strict margin keeps below threshold, and
            // each cost a present genome up to 7 % of its hits. Every genuine
            // rescue seen on the 5 493-genome index had 20 to 14 006 hits.
            // And a single read counts only for a leaf whose lineage is attested.
            // One read is evidence on a strict-margin index because a read from
            // another organism needs two coincident errors to hit a signature;
            // across thousands of close genomes the pairs of k-mers two
            // substitutions apart are numerous enough that such reads do occur,
            // in runs of 10 to 38 consecutive k-mers from one read of the
            // dominant genome. On the 5 493-genome sweep they made 43 false
            // strain rows, every one a single read, 35 of them leaves joining by
            // own evidence and 39 with no ancestor in common with the sample's
            // main; every true single-read row belonged to a family whose node
            // had evidence of its own (48 to 52 858 reads). So a leaf covered at
            // some level by a representative that cleared the reporting hit
            // threshold on two or more reads (lineageAttested) joins at
            // --rcutoff, and a leaf whose lineage has no such support needs two
            // reads of its own. Two variants were tried first. Two reads of
            // every rescue removed the 35 but cost real resolution, because the
            // rescue is judged on the leaf's global signatures and the
            // within-group ones nobody else carries, not on everything it is
            // later scored with: a strain whose global signatures all sit on a
            // 65-kb plasmid had one read there and 29 on its within-group
            // signatures, and the 1 000-genome sweep lost a close contaminant's
            // strain rows at 0.01 % (29 to 23) and 0.1 % (46 to 38), the
            // 5 493-genome sweep at 1 % (40 to 28). Gating on everInPlay instead
            // (one hit on a covering representative) removed 2 of the 43: a
            // stray hit on an ancestor's thousands of signatures is as common as
            // the stray run on the leaf's.
            const unsigned floor = static_cast<unsigned>(cfg.sketchHit);
            const size_t readFloor = static_cast<size_t>(cfg.readCutoff);
            const size_t unattestedFloor = std::max<size_t>(readFloor, 1);
            auto leafReadFloor = [&](unsigned i) {
                return (!leafName[i].empty() && lineageAttested.count(leafName[i]) > 0) ? readFloor : unattestedFloor;
            };
            size_t narrowed = 0, unseen = 0, rescuedGlobal = 0, rescuedGroup = 0, unnamed = 0;
            size_t rescuedOneRead = 0;                   // rescues on a single read, lineage attested
            for (unsigned i = 0; i < n; i++) {
                const unsigned t = L.targetBase + i;
                const bool inLineage = !leafName[i].empty() && candidates.count(leafName[i]) > 0;
                const bool ownGlobal = tally.count[t] >= floor && tally.reads[t].size() > leafReadFloor(i);
                scored[i] = inLineage || ownGlobal;
                if (inLineage) {
                    ++narrowed;
                    if (L.sig[i].empty()) ++unseen;         // cannot be reported; blocks what it carries
                } else if (ownGlobal) {
                    ++rescuedGlobal;
                    if (tally.reads[t].size() == 1) ++rescuedOneRead;
                } else if (leafName[i].empty()) {
                    // Nothing is known about this leaf, not even its name; it may be
                    // present. Keep it in play so that the signatures it carries
                    // are withheld.
                    scored[i] = true;
                    ++unnamed;
                }
            }
            // Rescue by within-group evidence: only the signatures with no carrier
            // at all count. Those are k-mers unique in the whole index that the
            // global selection did not draw (a quarter of the within-group set on
            // the 120-genome index), so a hit on one is evidence for its leaf
            // alone. A hit on a carried signature is evidence that the leaf *or a
            // carrier* is present, and the carrier is exactly the genome the
            // narrowing may have missed: counting such hits (tried first) put a
            // 0.01 % contaminant's reads on the within-group signatures of eight
            // unrelated leaves it carried, and reported them as strains; counting
            // them when no carrier was in play (tried second) did the same, since
            // a contaminant at that abundance leaves nothing else to be in play
            // by. The mask below judges the same hits the same way once the leaf
            // is in play, so this keeps the two decisions consistent.
            std::vector<unsigned> grpCount(n, 0);
            std::vector<std::unordered_set<std::string>> grpReads(n);
            for (unsigned g = 0; g < groupSigs.size(); g++) {
                const GroupSignature& s = groupSigs[g];
                if (groupDropped[g] || s.refId >= n || scored[s.refId] || !s.carriers.empty()) continue;
                unsigned t = groupBase + g;
                if (tally.count[t] == 0) continue;
                grpCount[s.refId] += tally.count[t];
                grpReads[s.refId].insert(tally.reads[t].begin(), tally.reads[t].end());
            }
            std::vector<size_t> grpRead(n, 0);
            for (unsigned i = 0; i < n; i++) grpRead[i] = grpReads[i].size();
            for (unsigned i = 0; i < n; i++) {
                if (!scored[i] && grpCount[i] >= floor && grpRead[i] > leafReadFloor(i)) {
                    scored[i] = true;
                    ++rescuedGroup;
                    if (grpRead[i] == 1) ++rescuedOneRead;
                }
            }
            // The mask's set: in play, or, for a leaf with no global signatures,
            // covered at some level by a representative with evidence. A leaf that
            // has global signatures and is present hits them (they are spread
            // along its genome), so its absence from `scored` means absence;
            // only the signature-less leaf leaves no trace of its own, and that is
            // the leak this closes. Extending the lineage rule to every leaf (tried
            // first) made all 20 members of a present main's level-0 node block,
            // and cost a 0.01 % sibling 175-388 of its within-group signatures on
            // the 120-genome sweep.
            present.assign(n, false);
            size_t presentByLineage = 0;
            for (unsigned i = 0; i < n; i++) {
                present[i] = scored[i] ||
                             (L.sig[i].empty() && !leafName[i].empty() && everInPlay.count(leafName[i]) > 0);
                if (present[i] && !scored[i]) ++presentByLineage;
            }
            std::cerr << "Tree query: leaf level: " << narrowed << " leaves in play after narrowing ("
                      << unseen << " with no global signatures), " << rescuedGlobal
                      << " more by hits on their own global signatures, " << rescuedGroup
                      << " by hits on their within-group signatures (" << rescuedOneRead
                      << " of these on one read, their lineage attested at a coarse level); "
                      << presentByLineage
                      << " signature-less leaves treated as present carriers on coarse-level evidence";
            if (unnamed) {
                std::cerr << ", " << unnamed << " unnamed and kept in play (no "
                          << treeLeavesPath(treeDir) << "; rebuild the index to name them)";
            }
            std::cerr << "\n";
        }

        std::vector<unsigned> levelCount(n, 0);
        std::vector<std::unordered_set<std::string>> levelRead(n);
        for (unsigned i = 0; i < n; i++) {
            if (!scored[i]) continue;
            unsigned t = L.targetBase + i;
            levelCount[i] = tally.count[t];
            levelRead[i].swap(tally.reads[t]);
        }
        for (unsigned s = L.slotBegin; s < L.slotEnd; s++)
            if (scored[slotRef[s]]) *slotCounter[s] += tally.slotCount[s];

        // Within-group signatures join the leaf level for leaves in play, each
        // admitted only if none of the genomes carrying it may be present: in
        // play itself, or covered at some level by a representative with
        // evidence. A carrier with any trace of evidence blocks the signature;
        // the cost is a smaller supplementary set in samples that contain the
        // carrier, never a hit attributed to the wrong genome.
        if (isLeaf && !groupSigs.empty()) {
            size_t used = 0, blocked = 0, idle = 0;
            std::unordered_set<unsigned> leavesUsing;
            std::unordered_set<std::string> admitted;
            for (unsigned g = 0; g < groupSigs.size(); g++) {
                const GroupSignature& s = groupSigs[g];
                if (s.refId >= n || !scored[s.refId]) { ++idle; continue; }
                bool carrierPresent = false;
                for (unsigned c : s.carriers)
                    if (c < n && present[c]) { carrierPresent = true; break; }
                if (carrierPresent) { ++blocked; continue; }
                if (groupDropped[g] || !admitted.insert(s.seq).second) continue;
                unsigned t = groupBase + g;
                L.sig[s.refId].emplace(s.seq, tally.count[t]);
                levelCount[s.refId] += tally.count[t];
                if (!tally.reads[t].empty()) {
                    if (levelRead[s.refId].empty()) levelRead[s.refId].swap(tally.reads[t]);
                    else levelRead[s.refId].insert(tally.reads[t].begin(), tally.reads[t].end());
                }
                ++used;
                leavesUsing.insert(s.refId);
            }
            std::cerr << "Tree query: within-group signatures: " << used << " in use for "
                      << leavesUsing.size() << " leaves, " << blocked
                      << " withheld (a carrier may be present), " << idle
                      << " belong to leaves not in play\n";
        }

        // Adjust for the first occurrence the --solid filter swallowed.
        if (cfg.solid)
            for (unsigned& c : levelCount) if (c) ++c;

        if (isLeaf) {
            generateTreeQueryResult(cfg, levelCount, L.ref, L.sig, levelRead,
                                    familyCalls, cfg.out);
            std::cerr << "Tree query: leaf level " << level << " ("
                      << n << " leaves) -> " << cfg.out << "\n";
            std::cerr << "Tree query time (s): load index " << std::setprecision(1) << std::fixed
                      << secs(t0, t1) << ", sample cardinality " << secs(t1, t2)
                      << ", read pass " << secs(t2, t3) << ", report " << secs(t3, now()) << "\n";
            break;
        }

        const auto& lm = levelMap[level];

        // Two independent decisions are made from the same coarse-level evidence.
        //
        // Pruning decides which sub-trees are still worth querying, and is
        // deliberately permissive: anything with a trace of evidence stays in play,
        // because a sub-tree dropped here can never be recovered at a finer level.
        // Requiring the full reporting thresholds to survive pruning made the tree
        // strictly less sensitive than a flat index, since a representative that
        // fell just short took all of its leaves down with it.
        //
        // Reporting decides whether this level's evidence is strong enough to stand
        // on its own as a family-level call, and uses the ordinary thresholds. Both
        // levels emit the same number of signatures per reference, so the counts are
        // on a comparable scale and no separate reporting threshold is needed.
        const auto& nameMap = levelNodeName[level];
        for (unsigned i = 0; i < n; i++) {
            if (!scored[i]) continue;
            auto it = lm.find(L.ref[i]);
            if (it == lm.end()) continue;

            bool survivesPruning =
                levelCount[i] >= static_cast<unsigned>(cfg.treePruneHit) &&
                !levelRead[i].empty();
            if (!survivesPruning) {
                for (const auto& lf : it->second) candidates.erase(lf);
                continue;
            }

            // A family call needs a real group behind it (levelNodeName holds only
            // nodes standing for two or more leaves); a leaf representing itself at a
            // coarse level is scored again, as a strain, at the leaf level.
            auto nit = nameMap.find(L.ref[i]);
            bool worthReporting =
                nit != nameMap.end() &&
                levelCount[i] >= static_cast<unsigned>(cfg.sketchHit) &&
                levelRead[i].size() > static_cast<size_t>(cfg.readCutoff);
            if (worthReporting) {
                // Held as a candidate family call; only emitted if the leaf level
                // fails to resolve any member of this group.
                FamilyCall fc;
                fc.name = nit->second;
                fc.level = level;
                fc.count = levelCount[i];
                fc.reads = levelRead[i].size();
                fc.sigsHit = countSigsHit(L.sig[i]);
                fc.sigsTotal = L.sig[i].size();
                fc.coveredLeaves = it->second;
                if (!ancestorCarried.empty()) {
                    const std::string key = ancestorKey(fc.name, level);
                    for (const std::string& lf : fc.coveredLeaves) {
                        auto li = ancestorCarried.find(lf);
                        if (li == ancestorCarried.end()) continue;
                        auto ni = li->second.find(key);
                        if (ni != li->second.end()) fc.memberCarried[lf] = ni->second;
                    }
                }
                familyCalls.push_back(std::move(fc));
            }
        }
        std::cerr << "Tree query: level " << level << " -> "
                  << candidates.size() << " candidate leaves";
        if (rejoined)
            std::cerr << "; " << rejoined << " node(s) pruned at a coarser level scored on their own evidence";
        std::cerr << "\n";
    }
}


int main(int argc, char* const argv[]) {
#ifdef _OPENMP
    double sTime = omp_get_wtime();
#endif

    QuerySketchConfig cfg;
    std::string treeDir;  // --tree: hierarchical query index directory (empty = flat)
    bool readCutoffGiven = false;  // --rcutoff passed explicitly (tree mode picks a default otherwise)

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) {
        std::istringstream arg(optarg != nullptr ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 't': arg >> cfg.threads;         break;
        case 'o': arg >> cfg.out;             break;
        case 'l': arg >> cfg.r1;              break;
        case 'r': arg >> cfg.r2;              break;
        case 'g': arg >> cfg.ref;             break;
        case 'h': arg >> g_sketchHit;          break;
        case 'a': arg >> cfg.abundanceCutoff; break;
        case 's': arg >> cfg.readCutoff; readCutoffGiven = true; break;
        case 'b': arg >> cfg.bits;            break;
        case OPT_BREADTH: arg >> cfg.breadth; break;
        case OPT_TREE: arg >> treeDir;         break;
        case OPT_TREEPRUNE: arg >> cfg.treePruneHit; break;
        case OPT_TREEBREADTH: arg >> cfg.treeBreadth; break;
        case OPT_NOGROUPSIGS: cfg.noGroupSignatures = true; break;
        case OPT_HELP:
            std::cerr << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cerr << VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        }
        if (optarg != nullptr && !arg.eof()) {
            std::cerr << PROGRAM ": invalid option: `-"
                      << static_cast<char>(c) << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }
    cfg.sketchHit = g_sketchHit;
    cfg.solid = g_solid;

    if (cfg.r1.empty()) {
        std::cerr << "Error: no input file r1 provided.\n";
        die = true;
    }
    if (cfg.ref.empty()) {
        std::cerr << "Error: no sketch file provided.\n";
        die = true;
    }
    if (die) {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

    // Validate input files exist
    if (!validateFile(cfg.r1, "read file r1")) die = true;
    if (!cfg.r2.empty() && !validateFile(cfg.r2, "read file r2")) die = true;
    if (!validateFile(cfg.ref, "sketch reference file")) die = true;
    if (die) {
        exit(EXIT_FAILURE);
    }

#ifdef _OPENMP
    omp_set_num_threads(cfg.threads);
#endif

    std::vector<std::string> sampleFiles;
    sampleFiles.push_back(cfg.r1);
    if (!cfg.r2.empty()) {
        sampleFiles.push_back(cfg.r2);
    }

    if (!treeDir.empty()) {
        // Hierarchical query: narrow candidates coarse -> fine down the tree.
        runTreeQuery(cfg, treeDir, sampleFiles, readCutoffGiven);
    } else {
        // Flat query against a single sketch.
        SketchHash sketchHash;
        std::vector<std::string> sketchRef;
        std::vector<SketchHash> refSigCount;

        loadSketch(cfg, cfg.ref, sketchHash, sketchRef, refSigCount);

        std::cout << "k: " << cfg.kmerLen << "\n";
        std::cout << "Number of references: " << sketchRef.size() << "\n";
        std::cout << "Size of sketch: " << sketchHash.size() << "\n";

        // The sample's k-mer cardinality sizes the --solid read filter. Without
        // --solid nothing uses it, and estimating it is a full pass over the reads
        // (and a 512 MB counter table), about as long as the query itself.
        if (cfg.solid) {
            size_t totalKmers = 0;
            getCardinality(totalKmers, cfg.dbfSize, cfg.sbfSize, cfg.kmerLen, cfg.threads, sampleFiles);

            if (!resolveBloomSizes(totalKmers, cfg.dbfSize, cfg.sbfSize)) {
                std::cerr << "Error: insufficient k-mer content for Bloom filter construction.\n"
                          << "  Total k-mers found: " << totalKmers << "\n"
                          << "  Distinct k-mers (F0): " << cfg.dbfSize << "\n"
                          << "  Solid k-mers: " << cfg.sbfSize << "\n"
                          << "  k-mer length: " << cfg.kmerLen << "\n"
                          << "Possible causes:\n"
                          << "  - Input reads are shorter than k (" << cfg.kmerLen << ")\n"
                          << "  - Input files are empty or contain too few reads\n"
                          << "  - Try a smaller k value in the sketch\n";
                exit(EXIT_FAILURE);
            }

            std::cout << "Number of weak k-mers: " << cfg.dbfSize - cfg.sbfSize << "\n";
            std::cout << "Number of solid k-mers: " << cfg.sbfSize << "\n";
        }

        std::vector<unsigned> sketchCount(sketchRef.size());
        std::vector<std::unordered_set<std::string>> refRead(sketchRef.size());

        querySampleBatch(cfg, sampleFiles, sketchHash, sketchCount, refSigCount, refRead);
        generateQueryResult(cfg, sketchCount, sketchRef, refSigCount, refRead, cfg.out);
    }

#ifdef _OPENMP
    std::cout << "Runtime(sec): " << std::setprecision(4) << std::fixed
              << omp_get_wtime() - sTime << "\n";
#endif
    return 0;
}