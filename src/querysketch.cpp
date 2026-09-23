#include <getopt.h>

#include <algorithm>
#include <deque>
#include <numeric>
#include <unordered_map>

#include "ntHashIterator.hpp"
#include "ntcard.hpp"
#include "QuerySketchUtil.hpp"
#include "TreeSketchUtil.hpp"
#include "GroupSketchUtil.hpp"
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
    "      --tree=DIR\thierarchical query using a --tree index directory\n"
    "      --tree-prune=N\tmin coarse-level hits to keep a sub-tree in play (with --tree) [1]\n"
    "      --no-group-signatures\tleaf level uses globally unique signatures only (with --tree)\n"
    "      --help\t\tdisplay this help and exit\n"
    "      --version\t\tversion information and exit\n"
    "\n";

static const char shortopts[] = "t:g:l:r:o:h:a:s:b:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_TREE, OPT_TREEPRUNE, OPT_NOGROUPSIGS };

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
    {"tree",           required_argument, nullptr, OPT_TREE},
    {"tree-prune",     required_argument, nullptr, OPT_TREEPRUNE},
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
 * candidates, so the search narrows down the tree. The final (leaf) level sketch
 * (cfg.ref) is queried restricted to the surviving candidate leaves and its result is
 * written. This applies the candidate narrowing that the flat port computed but never
 * used, so the tree actually focuses the query on promising sub-trees.
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
    if (!readCutoffGiven) {
        TreeBuildParams bp = loadTreeBuildParams(treeBuildParamsPath(treeDir));
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
    #pragma omp parallel for schedule(dynamic) num_threads(static_cast<int>(std::min<size_t>(levelSketches.size(), cfg.threads)))
    for (size_t li = 0; li < levelSketches.size(); li++) {
        QuerySketchConfig local = cfg;
        levels[li].level = levelSketches[li].first;
        loadSketch(local, levelSketches[li].second, levelHash[li], levels[li].ref, levels[li].sig);
    }
    {
        // k, as loadSketch would have recorded it: the length of the first signature.
        std::ifstream in(levelSketches.front().second);
        std::string first;
        std::getline(in, first);
        cfg.kmerLen = static_cast<unsigned>(first.find_first_of(" \t"));
    }

    for (size_t li = 0; li < levelSketches.size(); li++) {
        Level& L = levels[li];
        SketchHash& hash = levelHash[li];
        L.targetBase = index.numTargets;
        index.numTargets += static_cast<unsigned>(L.ref.size());
        L.slotBegin = static_cast<unsigned>(slotCounter.size());
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
                    index.add(kv.first, L.targetBase + i, static_cast<unsigned>(slotCounter.size()));
                    slotCounter.push_back(&kv.second);
                    slotRef.push_back(i);
                }
            }
        } else {
            for (const auto& [seq, refId] : hash) {
                auto it = L.sig[refId].find(seq);
                if (it != L.sig[refId].end()) {
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
    std::vector<FamilyCall> familyCalls;         // coarse detections, for level-appropriate reporting

    // Walk the levels coarse to fine over the tallies. Only references that still
    // cover a candidate leaf are scored at a level; the others' tallies are
    // discarded, exactly as their signatures used to be left out of that level's
    // pass, so a leaf whose family was pruned cannot resurface on its own hits.
    for (size_t li = 0; li < levels.size(); li++) {
        Level& L = levels[li];
        const unsigned level = L.level;
        const bool isLeaf = (li + 1 == levels.size());
        const size_t n = L.ref.size();

        std::vector<bool> scored(n, false);
        {
            const auto& lm = levelMap[level];
            for (unsigned i = 0; i < n; i++) {
                bool inPlay = false;
                if (isLeaf) {
                    inPlay = candidates.count(L.ref[i]) > 0;
                } else {
                    auto it = lm.find(L.ref[i]);
                    if (it != lm.end())
                        for (const auto& lf : it->second)
                            if (candidates.count(lf)) { inPlay = true; break; }
                }
                scored[i] = inPlay;
            }
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

        // Within-group signatures join the leaf level for candidate leaves, each
        // admitted only if none of the genomes carrying it is itself a candidate.
        // The candidate set is the permissive one the coarse levels left, so a
        // carrier with any trace of evidence blocks the signature; the cost is a
        // smaller supplementary set in samples that contain the carrier, never a
        // hit attributed to the wrong genome.
        if (isLeaf && !groupSigs.empty()) {
            size_t used = 0, blocked = 0, idle = 0;
            std::unordered_set<unsigned> leavesUsing;
            std::unordered_set<std::string> admitted;
            for (unsigned g = 0; g < groupSigs.size(); g++) {
                const GroupSignature& s = groupSigs[g];
                if (s.refId >= n || !scored[s.refId]) { ++idle; continue; }
                bool carrierPresent = false;
                for (unsigned c : s.carriers)
                    if (c < n && scored[c]) { carrierPresent = true; break; }
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
                      << leavesUsing.size() << " candidate leaves, " << blocked
                      << " withheld (a carrier is a candidate), " << idle
                      << " belong to non-candidates\n";
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
                familyCalls.push_back(std::move(fc));
            }
        }
        std::cerr << "Tree query: level " << level << " -> "
                  << candidates.size() << " candidate leaves\n";
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
        case OPT_TREE: arg >> treeDir;         break;
        case OPT_TREEPRUNE: arg >> cfg.treePruneHit; break;
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