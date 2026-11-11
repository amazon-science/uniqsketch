#ifndef GROUPSKETCHUTIL_HPP_
#define GROUPSKETCHUTIL_HPP_

#include <sys/stat.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "UniqSketchUtil.hpp"
#include "QuerySketchUtil.hpp"
#include "TreeSketchUtil.hpp"
#include "FileUtil.hpp"
#include "ntcard.hpp"

/**
 * Within-group leaf signatures.
 *
 * A leaf's ordinary signatures are its k-mers that occur in no other genome of
 * the index. For a strain with close relatives that pool is small, and it is
 * small for a reason that has nothing to do with the relatives: most of what the
 * strain does not share with its siblings it does share with some unrelated
 * genome elsewhere in the index, so global uniqueness throws it away. On 120
 * E. coli genomes the k-mers of a strain absent from its siblings but present
 * somewhere else outnumbered its globally unique k-mers by a median 3.2x, and by
 * 26x for the strain with the smallest global pool.
 *
 * Those k-mers can still separate the strain from its siblings, provided the
 * genomes that carry them are not in the sample. This header builds them as a
 * second signature set: for each leaf, k-mers unique within its group (the
 * leaves under its direct parent node), selected like ordinary signatures, each
 * recorded with the list of other genomes in the index that contain it (its
 * carriers). The query uses a within-group signature only when none of its
 * carriers is a candidate for the sample; see runTreeQuery. Global signatures
 * stay in use throughout, so the extra set can add evidence but never remove it.
 *
 * The residual risk is a carrier present below the candidate floor: its reads
 * then land on the shared signatures at that level. The floor is one hit at any
 * coarse level, so this needs a genome too rare for a single read to reach one
 * of its group's several thousand coarse signatures.
 */

/** A leaf's within-group signature with its provenance and carriers. */
struct GroupSignature {
    std::string seq;                 // canonical k-mer
    unsigned refId;                  // leaf it belongs to (leaf sketch numbering)
    std::string refName;
    std::string contig;
    unsigned pos;
    std::string node;                // group node the uniqueness was computed within
    std::vector<unsigned> carriers;  // other leaves (refIds) containing this k-mer
};

/**
 * Groups for within-group selection: one per internal node that has at least
 * one leaf as a direct child. The node's members are all its descendant leaves;
 * candidates are computed for its direct leaf children only (deeper leaves get
 * theirs from the finer node they hang off).
 */
struct SignatureGroup {
    std::string nodeName;
    std::vector<unsigned> members;      // all descendant leaves (refIds)
    std::vector<unsigned> directLeaves; // members whose parent is this node
};

inline std::vector<SignatureGroup> signatureGroups(const std::vector<TreeNode>& tree,
                                                   const std::vector<std::string>& refFiles) {
    std::unordered_map<std::string, unsigned> idOf;
    for (unsigned i = 0; i < refFiles.size(); i++) idOf[refFiles[i]] = i;
    std::unordered_map<unsigned, size_t> byId;
    for (size_t i = 0; i < tree.size(); i++) byId[tree[i].id] = i;

    // Descendant leaves per node.
    std::function<void(size_t, std::vector<unsigned>&)> collect =
        [&](size_t idx, std::vector<unsigned>& out) {
            const TreeNode& n = tree[idx];
            if (n.isLeaf) {
                auto it = idOf.find(n.filePath);
                if (it != idOf.end()) out.push_back(it->second);
                return;
            }
            for (unsigned c : n.children) {
                auto ci = byId.find(c);
                if (ci != byId.end()) collect(ci->second, out);
            }
        };

    std::map<unsigned, SignatureGroup> groups;   // by node id, for a stable order
    for (size_t i = 0; i < tree.size(); i++) {
        const TreeNode& n = tree[i];
        if (!n.isLeaf || n.parentId == -1) continue;
        auto pit = byId.find(static_cast<unsigned>(n.parentId));
        auto lit = idOf.find(n.filePath);
        if (pit == byId.end() || lit == idOf.end()) continue;
        SignatureGroup& g = groups[static_cast<unsigned>(n.parentId)];
        if (g.members.empty()) {
            g.nodeName = tree[pit->second].name;
            collect(pit->second, g.members);
            std::sort(g.members.begin(), g.members.end());
        }
        g.directLeaves.push_back(lit->second);
    }
    std::vector<SignatureGroup> out;
    for (auto& kv : groups)
        if (kv.second.members.size() >= 2) out.push_back(std::move(kv.second));
    return out;
}

/**
 * Fill in the carriers of every signature: each genome of the index other than
 * the signature's own leaf that contains the k-mer. One pass over each genome,
 * probing a hash-keyed view of the signatures and confirming the sequence on a
 * probe hit, exactly as the query does with reads.
 */
inline void findCarriers(const UniqSketchConfig& cfg, std::vector<GroupSignature>& sigs,
                         const std::vector<std::string>& refFiles) {
    if (sigs.empty()) return;
    std::unordered_map<uint64_t, std::vector<unsigned>> index;
    index.reserve(sigs.size());
    for (unsigned i = 0; i < sigs.size(); i++)
        index[NTC64(sigs[i].seq.c_str(), cfg.kmerLen)].push_back(i);

    std::vector<std::vector<bool>> hit(refFiles.size());
    #pragma omp parallel for schedule(dynamic)
    for (unsigned g = 0; g < refFiles.size(); g++) {
        std::vector<bool> mark(sigs.size(), false);
        scanFasta(refFiles[g], [&](const std::string&, const std::string& seq, size_t) {
            if (seq.size() < cfg.kmerLen) return;
            ntHashIterator itr(seq, 1, cfg.kmerLen);
            while (itr != itr.end()) {
                auto it = index.find((*itr)[0]);
                if (it != index.end()) {
                    const char* kp = seq.c_str() + itr.get_pos();
                    for (unsigned idx : it->second)
                        if (sigs[idx].refId != g && kmerMatches(kp, sigs[idx].seq)) mark[idx] = true;
                }
                ++itr;
            }
        }, wrapPad(cfg));
        hit[g] = std::move(mark);
    }
    for (unsigned g = 0; g < refFiles.size(); g++)
        for (unsigned i = 0; i < sigs.size(); i++)
            if (hit[g][i]) sigs[i].carriers.push_back(g);
}

inline std::string groupSketchPath(const std::string& outdir) {
    return outdir + "/sketch_group.tsv";
}

inline void writeGroupSketch(const std::string& path, const std::vector<GroupSignature>& sigs,
                             unsigned numRef) {
    std::ofstream out(path);
    for (const GroupSignature& s : sigs) {
        out << s.seq << "\t" << s.refId << "\t" << s.refName << "\t" << numRef
            << "\t" << s.contig << "\t" << s.pos << "\t";
        if (s.carriers.empty()) out << "-";
        for (size_t i = 0; i < s.carriers.size(); i++) out << (i ? "," : "") << s.carriers[i];
        out << "\t" << s.node << "\n";
    }
}

/**
 * Build the within-group signature set for a tree index.
 *
 * Must run after the leaf sketch, with the same `refFiles` order (it defines the
 * refIds) and the universe distinct filter from that build (for the margin
 * check). Signatures already present in the leaf sketch are not repeated.
 *
 * The work is laid out so that a many-core machine stays busy on many small
 * groups. Filter sizes come from ntCard one group at a time (it runs its own
 * parallel region over the group's files and sets the OpenMP thread count, so
 * it cannot be called from inside a parallel loop). Candidate identification is
 * the ordinary distinct/solid filter pass: each group's members go into its
 * filters with every thread on that group (GroupScanner, so the result does not
 * depend on the thread count), then the direct leaves of a whole batch of groups
 * are scanned for candidates as one flat list of (group, leaf) jobs with the
 * batch's filters live at once, so a group of two genomes does not leave the
 * rest of the machine idle; a batch takes as many groups as fit in the memory
 * the main build's own filters already used, so peak memory does not grow with
 * the number of groups. Candidates are written only for the direct leaves, the
 * only ones selected from here.
 * Selection is one flat list of (group, leaf) jobs, each writing to its own
 * buffer, so the output order is fixed regardless of scheduling. The sequential
 * form of this (one group after another, one leaf after another) took 23
 * minutes on 120 E. coli genomes, most of it in the single-threaded margin
 * checks over 21 million candidates.
 *
 * @return Number of signatures written.
 */
inline size_t buildGroupSketch(const UniqSketchConfig& cfg, const std::vector<TreeNode>& tree,
                               const std::vector<std::string>& refFiles,
                               const BloomFilter& universeDbf) {
    using Clock = std::chrono::steady_clock;
    auto since = [](Clock::time_point t) {
        return std::chrono::duration<double>(Clock::now() - t).count();
    };

    std::vector<SignatureGroup> groups = signatureGroups(tree, refFiles);
    if (groups.empty()) return 0;

    // Global signatures per leaf, so the within-group set holds only the extra ones.
    std::vector<std::unordered_set<std::string>> globalSigs(refFiles.size());
    {
        std::ifstream in(cfg.outfile);
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream ls(line);
            std::string seq; unsigned refId;
            if (ls >> seq >> refId && refId < globalSigs.size()) globalSigs[refId].insert(seq);
        }
    }

    std::string groupRoot = cfg.outdir + "/groups";
    mkdir(groupRoot.c_str(), 0777);

    // Per-group configuration: the ordinary pipeline's settings over the group's
    // genomes alone, output under groups/g<i>. The leaf numbering stays the
    // index-wide one; only the filter sizes are the group's.
    Clock::time_point t = Clock::now();
    std::vector<UniqSketchConfig> gcfg(groups.size(), cfg);
    std::vector<bool> usable(groups.size(), false);
    for (size_t gi = 0; gi < groups.size(); gi++) {
        UniqSketchConfig& c = gcfg[gi];
        c.outdir = groupRoot + "/g" + std::to_string(gi);
        c.refstat = c.outdir + "/db_uniq_count.tsv";
        c.numRef = groups[gi].members.size();
        mkdir(c.outdir.c_str(), 0777);
        std::vector<std::string> files;
        for (unsigned m : groups[gi].members) files.push_back(refFiles[m]);
        size_t total = 0;
        getCardinality(total, c.dbfSize, c.sbfSize, c.kmerLen, c.threads, files);
        usable[gi] = resolveBloomSizes(total, c.dbfSize, c.sbfSize);
    }
    const double secCardinality = since(t);

    // Candidates unique within each group. Every member genome goes into the
    // group's distinct/solid filter pair; every direct leaf is then scanned
    // against the solid filter (a k-mer distinct but not solid within the group
    // is in exactly one member). Same filters and sizes as identifyUniqKmers
    // would build for the group, so the candidate files are the ones it would
    // write, minus the members that are not direct leaves.
    t = Clock::now();
    std::vector<std::vector<int>> candCount(groups.size());   // per direct leaf
    const size_t budgetBits = cfg.m1() + cfg.m2();
    GroupScanner scanner(cfg);                                 // scratch shared by every group
    for (size_t batchStart = 0; batchStart < groups.size();) {
        size_t batchEnd = batchStart;
        for (size_t bits = 0; batchEnd < groups.size(); batchEnd++) {
            size_t need = usable[batchEnd] ? gcfg[batchEnd].m1() + gcfg[batchEnd].m2() : 0;
            if (batchEnd > batchStart && bits + need > budgetBits) break;
            bits += need;
        }

        std::vector<std::unique_ptr<BloomFilter>> dbf(groups.size()), sbf(groups.size());
        std::vector<std::pair<unsigned, unsigned>> scanJobs;   // (group, refId)
        for (size_t gi = batchStart; gi < batchEnd; gi++) {
            if (!usable[gi]) continue;
            dbf[gi].reset(new BloomFilter(gcfg[gi].m1(), gcfg[gi].nhash1, gcfg[gi].kmerLen));
            sbf[gi].reset(new BloomFilter(gcfg[gi].m2(), gcfg[gi].nhash2, gcfg[gi].kmerLen));
            for (unsigned l : groups[gi].directLeaves) scanJobs.emplace_back(gi, l);
            candCount[gi].assign(groups[gi].directLeaves.size(), 0);
        }

        // Each group's members into its filters, one group at a time with every
        // thread on it (GroupScanner): the same filters on any thread count, and
        // a group of twenty genomes no longer waits on one thread. Members went
        // into the filters one per thread before, over a flat list of (group,
        // member) jobs, which was order-dependent in the same way as the level
        // scans.
        for (size_t gi = batchStart; gi < batchEnd; gi++) {
            if (!usable[gi]) continue;
            std::vector<std::vector<std::string>> singles;
            for (unsigned m : groups[gi].members) singles.push_back({refFiles[m]});
            scanner.load(singles, *dbf[gi], *sbf[gi]);
        }
        #pragma omp parallel for schedule(dynamic)
        for (size_t j = 0; j < scanJobs.size(); j++) {
            unsigned gi = scanJobs[j].first, leaf = scanJobs[j].second;
            const auto& dl = groups[gi].directLeaves;
            int slot = static_cast<int>(std::find(dl.begin(), dl.end(), leaf) - dl.begin());
            checkRef(gcfg[gi], refFiles[leaf], slot, *sbf[gi], candCount[gi]);
        }

        // The group's stat file, as identifyUniqKmers writes it: candidate path,
        // count, and the leaf's index-wide refId.
        for (size_t gi = batchStart; gi < batchEnd; gi++) {
            if (!usable[gi]) continue;
            std::ofstream stat(gcfg[gi].refstat);
            const auto& dl = groups[gi].directLeaves;
            for (size_t i = 0; i < dl.size(); i++)
                stat << candidatePath(gcfg[gi].outdir, refFiles[dl[i]]) << "\t"
                     << candCount[gi][i] << "\t" << dl[i] << "\n";
        }
        batchStart = batchEnd;
    }
    const double secCandidates = since(t);

    // Select for every direct leaf, with the leaf sketch's numbering and the
    // universe margin. Candidates are unique within their group by construction
    // (and across groups by the leaf-sketch exclusion below), so no claim filter
    // is needed between leaves; see getUniqSet.
    t = Clock::now();
    struct SelectJob { unsigned group, leaf; };
    std::vector<SelectJob> selectJobs;
    for (size_t gi = 0; gi < groups.size(); gi++) {
        if (!usable[gi]) continue;
        for (unsigned l : groups[gi].directLeaves)
            selectJobs.push_back({static_cast<unsigned>(gi), l});
    }
    std::vector<std::ostringstream> picked(selectJobs.size());
    size_t totalSelected = 0;
    #pragma omp parallel for schedule(dynamic) reduction(+:totalSelected)
    for (size_t j = 0; j < selectJobs.size(); j++) {
        const SelectJob& sj = selectJobs[j];
        std::string candPath = candidatePath(gcfg[sj.group].outdir, refFiles[sj.leaf]);
        totalSelected += getUniqSet(cfg, candPath, sj.leaf, picked[j], universeDbf).selected;
    }

    std::vector<GroupSignature> all;
    for (size_t j = 0; j < selectJobs.size(); j++) {
        std::istringstream lines(picked[j].str());
        std::string line;
        while (std::getline(lines, line)) {
            std::istringstream ls(line);
            GroupSignature s;
            unsigned numRef;
            ls >> s.seq >> s.refId >> s.refName >> numRef >> s.contig >> s.pos;
            if (ls.fail()) continue;
            if (globalSigs[s.refId].count(s.seq)) continue;
            s.node = groups[selectJobs[j].group].nodeName;
            all.push_back(std::move(s));
        }
    }
    const double secSelection = since(t);

    t = Clock::now();
    findCarriers(cfg, all, refFiles);
    writeGroupSketch(groupSketchPath(cfg.outdir), all, static_cast<unsigned>(refFiles.size()));
    const double secCarriers = since(t);

    size_t totalCandidates = 0, withCarriers = 0;
    for (const auto& counts : candCount)
        for (int n : counts) totalCandidates += static_cast<size_t>(n);
    for (const auto& s : all) withCarriers += !s.carriers.empty();
    std::cerr << "Within-group signatures: " << groups.size() << " groups, "
              << totalCandidates << " candidates, " << totalSelected << " selected, "
              << all.size() << " written (" << totalSelected - all.size()
              << " already in the leaf sketch), " << withCarriers << " with at least one carrier\n"
              << "Within-group build time (s): cardinality " << std::setprecision(1) << std::fixed
              << secCardinality << ", candidates " << secCandidates << ", selection "
              << secSelection << ", carriers " << secCarriers << "\n";
    return all.size();
}

/**
 * Load sketch_group.tsv. Returns an empty vector when the file is absent. Lines
 * with fewer than the eight fields are skipped, as the stream-based parser this
 * replaced skipped them (an istringstream per line cost 0.4 s on 200 k lines).
 */
inline std::vector<GroupSignature> loadGroupSketch(const std::string& path) {
    std::vector<GroupSignature> out;
    std::ifstream in(path);
    if (!in) return out;
    std::string line;
    std::vector<std::pair<size_t, size_t>> f(8);
    while (std::getline(in, line)) {
        size_t n = 0, i = 0;
        while (n < 8) {
            while (i < line.size() && (line[i] == '\t' || line[i] == ' ')) ++i;
            if (i >= line.size()) break;
            size_t b = i;
            while (i < line.size() && line[i] != '\t' && line[i] != ' ') ++i;
            f[n++] = {b, i - b};
        }
        if (n < 8) continue;
        GroupSignature s;
        s.seq.assign(line, f[0].first, f[0].second);
        s.refId = static_cast<unsigned>(std::strtoul(line.c_str() + f[1].first, nullptr, 10));
        s.refName.assign(line, f[2].first, f[2].second);
        s.contig.assign(line, f[4].first, f[4].second);
        s.pos = static_cast<unsigned>(std::strtoul(line.c_str() + f[5].first, nullptr, 10));
        if (line[f[6].first] != '-') {
            const char* p = line.c_str() + f[6].first;
            const char* end = p + f[6].second;
            while (p < end) {
                char* q = nullptr;
                unsigned long v = std::strtoul(p, &q, 10);
                if (q == p) break;
                s.carriers.push_back(static_cast<unsigned>(v));
                p = (q < end && *q == ',') ? q + 1 : q;
            }
        }
        s.node.assign(line, f[7].first, f[7].second);
        out.push_back(std::move(s));
    }
    return out;
}

/**
 * What each leaf of a tree index can be reported through, and which leaves no
 * query can report at all.
 *
 * A leaf is reported as a strain through its own signatures (the leaf sketch) or
 * its within-group signatures, and as a member of a family through the signatures
 * of an ancestor node that its genome carries. A node's signatures are its
 * representative's k-mers absent from every genome outside the node, which in a
 * dense set is often the representative's private sequence, so a member may carry
 * none of them (on 5 493 E. coli genomes, 2 408 leaves had at least one such
 * ancestor). A leaf with no signatures of its own, none within its group and none
 * of any ancestor's carried cannot be reported by any query at any threshold: a
 * sample that is only that genome yields an empty table. 56 of the 5 493 were in
 * that position, 42 of them the K-12 core (MG1655, BW25113 and their derivatives),
 * whose lineage the clustering had split across nodes so that nothing of the core
 * was unique to any node.
 *
 * One pass over every genome against a hash-keyed view of all level signatures,
 * confirming the sequence on a probe hit as findCarriers does. Writes
 * <outdir>/ancestor_signatures.tsv, one row per (leaf, ancestor node) with the
 * node's signature count and how many the leaf carries, and
 * <outdir>/reachability.tsv, one row per leaf with its own and within-group
 * counts, its ancestors, and a status: own, group, family or none. The stderr
 * summary names the leaves with status none (the first 20).
 *
 * @param levels  Tree levels whose sketch files exist (sketch_level_<L>.tsv).
 */
inline std::string ancestorSignaturesPath(const std::string& outdir) {
    return outdir + "/ancestor_signatures.tsv";
}
/** Key of a node at a level in the ancestor-signature table: "<node>@<level>". A node
 *  spanning several levels has a different signature draw at each. */
inline std::string ancestorKey(const std::string& node, unsigned level) {
    return node + "@" + std::to_string(level);
}
/**
 * Load ancestor_signatures.tsv: leaf -> (node@level -> signatures of that node the leaf
 * genome carries). Empty when the file is absent (an index built before the table).
 */
inline std::unordered_map<std::string, std::unordered_map<std::string, size_t>>
loadAncestorSignatures(const std::string& path) {
    std::unordered_map<std::string, std::unordered_map<std::string, size_t>> out;
    std::ifstream in(path);
    if (!in.good()) return out;
    std::string line;
    std::getline(in, line);   // header
    while (std::getline(in, line)) {
        std::vector<std::string> f;
        size_t b = 0;
        for (size_t e; (e = line.find('\t', b)) != std::string::npos; b = e + 1) f.push_back(line.substr(b, e - b));
        f.push_back(line.substr(b));
        if (f.size() < 5) continue;
        out[f[0]][ancestorKey(f[2], static_cast<unsigned>(std::strtoul(f[1].c_str(), nullptr, 10)))] =
            std::strtoul(f[4].c_str(), nullptr, 10);
    }
    return out;
}
inline std::string reachabilityPath(const std::string& outdir) {
    return outdir + "/reachability.tsv";
}
inline void reportReachability(const UniqSketchConfig& cfg, const std::vector<TreeNode>& tree,
                               const std::vector<std::string>& refFiles,
                               const std::vector<unsigned>& levels) {
    using Clock = std::chrono::steady_clock;
    const auto t0 = Clock::now();
    const size_t numRef = refFiles.size();

    // Own and within-group signatures per leaf, counted from the files as written.
    auto countByRef = [&](const std::string& path) {
        std::vector<size_t> n(numRef, 0);
        std::ifstream in(path);
        std::string line;
        while (std::getline(in, line)) {
            size_t a = line.find('\t');
            if (a == std::string::npos) continue;
            size_t b = line.find('\t', a + 1);
            unsigned long id = std::strtoul(line.c_str() + a + 1, nullptr, 10);
            if (b != std::string::npos && id < numRef) n[id]++;
        }
        return n;
    };
    const std::vector<size_t> leafSigs = countByRef(cfg.outfile);
    const std::vector<size_t> groupSigs = countByRef(groupSketchPath(cfg.outdir));

    // The partition of the leaves at every level, and each leaf's entry there.
    std::vector<std::vector<LevelEntry>> entries;                      // per level, in `levels` order
    std::vector<std::unordered_map<std::string, unsigned>> repEntry;   // per level: rep basename -> entry
    std::vector<std::vector<int>> leafEntry;                           // per level: refId -> entry, -1 if none
    {
        std::unordered_map<std::string, unsigned> leafId;
        for (unsigned g = 0; g < numRef; g++) leafId[getBaseId(refFiles[g])] = g;
        std::map<unsigned, std::vector<LevelEntry>> byLevel;
        for (auto& [lvl, ents] : treeLevelPartition(tree)) byLevel[lvl] = std::move(ents);
        for (unsigned lvl : levels) {
            entries.push_back(byLevel.count(lvl) ? byLevel[lvl] : std::vector<LevelEntry>());
            std::unordered_map<std::string, unsigned> re;
            std::vector<int> le(numRef, -1);
            for (unsigned e = 0; e < entries.back().size(); e++) {
                re[getBaseId(entries.back()[e].repPath)] = e;
                for (const std::string& p : entries.back()[e].leafPaths) {
                    auto it = leafId.find(getBaseId(p));
                    if (it != leafId.end()) le[it->second] = static_cast<int>(e);
                }
            }
            repEntry.push_back(std::move(re));
            leafEntry.push_back(std::move(le));
        }
    }

    // Every level signature: sequence, level index, entry index; and per (level,
    // entry) the number of signatures.
    std::vector<std::string> sigSeq;
    std::vector<unsigned> sigLevel, sigEntry;
    std::vector<std::vector<size_t>> nodeTotal(levels.size());
    for (size_t li = 0; li < levels.size(); li++) {
        nodeTotal[li].assign(entries[li].size(), 0);
        std::ifstream in(cfg.outdir + "/sketch_level_" + std::to_string(levels[li]) + ".tsv");
        std::string line;
        while (std::getline(in, line)) {
            size_t a = line.find('\t');
            if (a == std::string::npos) continue;
            size_t b = line.find('\t', a + 1);
            if (b == std::string::npos) continue;
            size_t c = line.find('\t', b + 1);
            if (c == std::string::npos) continue;
            auto it = repEntry[li].find(line.substr(b + 1, c - b - 1));
            if (it == repEntry[li].end()) continue;
            sigSeq.push_back(line.substr(0, a));
            sigLevel.push_back(static_cast<unsigned>(li));
            sigEntry.push_back(it->second);
            nodeTotal[li][it->second]++;
        }
    }
    std::unordered_map<uint64_t, std::vector<unsigned>> index;
    index.reserve(sigSeq.size());
    for (unsigned i = 0; i < sigSeq.size(); i++)
        index[NTC64(sigSeq[i].c_str(), cfg.kmerLen)].push_back(i);

    // One pass per genome: distinct signatures of each ancestor that the genome carries.
    std::vector<std::vector<size_t>> carried(numRef, std::vector<size_t>(levels.size(), 0));
    #pragma omp parallel for schedule(dynamic)
    for (unsigned g = 0; g < numRef; g++) {
        std::vector<bool> seen(sigSeq.size(), false);
        scanFasta(refFiles[g], [&](const std::string&, const std::string& seq, size_t) {
            if (seq.size() < cfg.kmerLen) return;
            ntHashIterator itr(seq, 1, cfg.kmerLen);
            while (itr != itr.end()) {
                auto it = index.find((*itr)[0]);
                if (it != index.end()) {
                    const char* kp = seq.c_str() + itr.get_pos();
                    for (unsigned idx : it->second) {
                        if (seen[idx] || !kmerMatches(kp, sigSeq[idx])) continue;
                        seen[idx] = true;
                        const unsigned li = sigLevel[idx];
                        if (leafEntry[li][g] == static_cast<int>(sigEntry[idx])) carried[g][li]++;
                    }
                }
                ++itr;
            }
        }, wrapPad(cfg));
    }

    // Tables and summary.
    std::ofstream anc(ancestorSignaturesPath(cfg.outdir));
    anc << "leaf\tlevel\tnode\tnode_signatures\tcarried\n";
    std::ofstream reach(reachabilityPath(cfg.outdir));
    reach << "ref_id\tname\tleaf_signatures\tgroup_signatures\tancestors\tstatus\n";
    size_t own = 0, group = 0, family = 0, none = 0;
    std::vector<std::string> unreachable;
    for (unsigned g = 0; g < numRef; g++) {
        const std::string name = getBaseId(refFiles[g]);
        std::string ancestors;
        size_t familySigs = 0;
        for (size_t li = 0; li < levels.size(); li++) {
            const int e = leafEntry[li][g];
            if (!ancestors.empty()) ancestors += ";";
            ancestors += "L" + std::to_string(levels[li]) + ":";
            if (e < 0 || !entries[li][e].isNode) { ancestors += "-"; continue; }
            const LevelEntry& ent = entries[li][e];
            ancestors += ent.name + ":" + std::to_string(carried[g][li]) + "/" + std::to_string(nodeTotal[li][e]);
            familySigs += carried[g][li];
            anc << name << "\t" << levels[li] << "\t" << ent.name << "\t" << nodeTotal[li][e]
                << "\t" << carried[g][li] << "\n";
        }
        const char* status;
        if (leafSigs[g] > 0) { status = "own"; ++own; }
        else if (groupSigs[g] > 0) { status = "group"; ++group; }
        else if (familySigs > 0) { status = "family"; ++family; }
        else { status = "none"; ++none; unreachable.push_back(name); }
        reach << g << "\t" << name << "\t" << leafSigs[g] << "\t" << groupSigs[g] << "\t"
              << ancestors << "\t" << status << "\n";
    }
    std::cerr << "Reachability: " << own << " references have leaf signatures, " << group
              << " only within-group signatures, " << family
              << " only a family row (an ancestor's signatures carried), " << none
              << " nothing at any level";
    if (none) {
        std::cerr << " and cannot be reported by any query:";
        const size_t show = std::min<size_t>(unreachable.size(), 20);
        for (size_t i = 0; i < show; i++) std::cerr << (i ? ", " : " ") << unreachable[i];
        if (unreachable.size() > show) std::cerr << " and " << unreachable.size() - show << " more";
    }
    std::cerr << " (" << reachabilityPath(cfg.outdir) << ")\n";
    std::cerr << "Reachability scan time (s): " << std::setprecision(1) << std::fixed
              << std::chrono::duration<double>(Clock::now() - t0).count() << "\n";
}

#endif // GROUPSKETCHUTIL_HPP_
