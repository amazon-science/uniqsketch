#ifndef TREESKETCHUTIL_HPP_
#define TREESKETCHUTIL_HPP_

#include <algorithm>
#include <climits>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "SequenceUtil.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"
#include "ClusterUtil.hpp"
#include "kseq_gz.h"

/**
 * A node in the hierarchical sketch tree.
 * Internal nodes represent clusters; leaves represent individual references.
 */
struct TreeNode {
    unsigned id;                            // node id
    std::string name;                       // reference or cluster name
    std::string filePath;                   // file path (leaf nodes only)
    unsigned level;                         // depth in tree (0 = root level)
    int parentId;                           // parent node id (-1 for roots)
    std::vector<unsigned> children;         // child node ids
    bool isLeaf;                            // true if this is an individual reference
};

/**
 * Build a hierarchical clustering tree from reference files.
 * Uses multiple thresholds to create levels of increasing resolution.
 *
 * @param refFiles  All reference file paths.
 * @param thresholds  Descending list of unique k-mer thresholds for each level.
 * @param kmerLen  K-mer length.
 * @param nhash  Number of hash functions.
 * @param bits  Bits per element in Bloom filter.
 * @param dbfSize  Approximate genome size for BF sizing.
 * @param treeFile  Output TSV path for the tree structure.
 * @return Vector of TreeNode describing the full hierarchy.
 */
inline std::vector<TreeNode> buildHierarchicalTree(
        const std::vector<std::string>& refFiles,
        const std::vector<unsigned>& thresholds,
        unsigned kmerLen, unsigned nhash, unsigned bits, size_t dbfSize,
        const std::string& treeFile) {

    size_t n = refFiles.size();
    (void) nhash; (void) bits; (void) dbfSize;  // sizing is internal to the MinHash sketch

    // MinHash sketch (bottom-s) + exact k-mer count per reference (near-linear;
    // replaces the previous all-pairs O(n^2 * genome) unique-k-mer computation).
    const unsigned sketchSize = MINHASH_SKETCH_SIZE;
    std::vector<std::vector<uint64_t>> sketches(n);
    std::vector<size_t> kmerCount(n, 0);
    std::vector<std::string> names(n);

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < n; i++) {
        sketches[i] = clusterMinhashSketch(refFiles[i], kmerLen, sketchSize, kmerCount[i]);
        names[i] = getBaseId(refFiles[i]);
    }

    // Candidate distinguishing-k-mer estimates: only pairs that share sketch hashes
    // are scored (via an inverted index), rather than all O(n^2) pairs.
    std::unordered_map<uint64_t, std::vector<unsigned>> hashIndex;
    for (unsigned i = 0; i < n; i++)
        for (uint64_t h : sketches[i]) hashIndex[h].push_back(i);

    const size_t POPULAR_CAP = 20000;  // skip hashes shared by too many refs (stay linear)
    std::unordered_map<uint64_t, uint32_t> sharedCnt;
    for (auto& kv : hashIndex) {
        auto& refs = kv.second;
        if (refs.size() < 2 || refs.size() > POPULAR_CAP) continue;
        for (size_t a = 0; a < refs.size(); a++)
            for (size_t b = a + 1; b < refs.size(); b++) {
                unsigned i = refs[a], j = refs[b];
                if (i > j) std::swap(i, j);
                sharedCnt[(static_cast<uint64_t>(i) << 32) | j]++;
            }
    }

    // Estimated distinguishing k-mers between i and j (huge if they share no sketch
    // hashes, i.e. too distant to be a merge candidate).
    auto uniqueEst = [&](unsigned i, unsigned j) -> double {
        if (i > j) std::swap(i, j);
        auto it = sharedCnt.find((static_cast<uint64_t>(i) << 32) | j);
        if (it == sharedCnt.end()) return 1e18;
        size_t denom = std::min(sketches[i].size(), sketches[j].size());
        if (denom == 0) return 1e18;
        double jac = static_cast<double>(it->second) / static_cast<double>(denom);
        return static_cast<double>(std::min(kmerCount[i], kmerCount[j])) * (1.0 - jac);
    };

    // Build tree nodes: start with leaf nodes for each reference
    std::vector<TreeNode> tree;
    tree.reserve(n * 2);

    for (unsigned i = 0; i < n; i++) {
        TreeNode leaf;
        leaf.id = i;
        leaf.name = names[i];
        leaf.filePath = refFiles[i];
        leaf.level = thresholds.size(); // leaves are at the deepest level
        leaf.parentId = -1;
        leaf.isLeaf = true;
        tree.push_back(std::move(leaf));
    }

    // Bottom-up agglomerative nesting. Thresholds are processed finest-first so that
    // coarser levels are unions of finer clusters (a proper nested hierarchy, rather
    // than re-clustering from scratch per level). Clusters merge by representative
    // distance (representative-linkage) via union-find over the sparse candidate
    // graph: near-linear, and it avoids single-linkage chaining. Level numbering is
    // unchanged (0 = coarsest, thresholds.size() = leaves).
    struct ActiveCluster { unsigned rep; unsigned node; };
    std::vector<ActiveCluster> active;
    active.reserve(n);
    for (unsigned i = 0; i < n; i++) active.push_back({i, i});

    for (int lvl = static_cast<int>(thresholds.size()) - 1; lvl >= 0; --lvl) {
        double threshold = static_cast<double>(thresholds[lvl]);
        unsigned m = static_cast<unsigned>(active.size());

        std::vector<unsigned> parent(m);
        std::iota(parent.begin(), parent.end(), 0);
        auto find = [&](unsigned x) -> unsigned {
            while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
            return x;
        };

        // Map each current cluster's representative ref to its active-cluster index.
        std::unordered_map<unsigned, unsigned> repToIdx;
        for (unsigned c = 0; c < m; c++) repToIdx[active[c].rep] = c;

        // Merge clusters whose representatives form a candidate pair below threshold.
        for (auto& kv : sharedCnt) {
            unsigned i = static_cast<unsigned>(kv.first >> 32);
            unsigned j = static_cast<unsigned>(kv.first & 0xffffffffu);
            auto ii = repToIdx.find(i);
            auto jj = repToIdx.find(j);
            if (ii == repToIdx.end() || jj == repToIdx.end()) continue;
            if (uniqueEst(i, j) < threshold) {
                unsigned a = find(ii->second), b = find(jj->second);
                if (a != b) parent[a] = b;
            }
        }

        std::unordered_map<unsigned, std::vector<unsigned>> groups;
        for (unsigned c = 0; c < m; c++) groups[find(c)].push_back(c);

        std::vector<ActiveCluster> next;
        next.reserve(groups.size());
        for (auto& kv : groups) {
            auto& members = kv.second;
            if (members.size() == 1) {
                next.push_back(active[members[0]]);   // no merge at this level; carry up
                continue;
            }
            // Representative = largest genome among the merged clusters' representatives.
            unsigned repIdx = active[members[0]].rep;
            for (unsigned c : members)
                if (kmerCount[active[c].rep] > kmerCount[repIdx]) repIdx = active[c].rep;

            TreeNode internal;
            internal.id = static_cast<unsigned>(tree.size());
            internal.name = "cluster_L" + std::to_string(lvl) + "_" + names[repIdx];
            internal.filePath = refFiles[repIdx];
            internal.level = static_cast<unsigned>(lvl);
            internal.parentId = -1;
            internal.isLeaf = false;
            for (unsigned c : members) {
                unsigned childNode = active[c].node;
                tree[childNode].parentId = static_cast<int>(internal.id);
                internal.children.push_back(childNode);
            }
            unsigned newNode = static_cast<unsigned>(tree.size());
            tree.push_back(std::move(internal));
            next.push_back({repIdx, newNode});
        }
        active = std::move(next);
    }

    // Write tree structure to file
    std::ofstream treeOut(treeFile);
    treeOut << "node_id\tname\tlevel\tparent_id\tis_leaf\tchildren\tfile_path\n";
    for (const auto& node : tree) {
        treeOut << node.id << "\t" << node.name << "\t" << node.level << "\t"
                << node.parentId << "\t" << (node.isLeaf ? "yes" : "no") << "\t";
        for (size_t c = 0; c < node.children.size(); c++) {
            if (c > 0) treeOut << ",";
            treeOut << node.children[c];
        }
        treeOut << "\t" << node.filePath << "\n";
    }
    treeOut.close();

    // Report tree stats
    unsigned numLeaves = 0, numInternal = 0;
    for (const auto& node : tree) {
        if (node.isLeaf) ++numLeaves;
        else ++numInternal;
    }
    std::cerr << "Hierarchical tree built: " << numLeaves << " leaves, "
              << numInternal << " internal nodes, " << thresholds.size() << " levels\n";

    return tree;
}

/**
 * One entry of a level's partition: the reference whose genome stands for a set
 * of leaves at that level.
 */
struct LevelEntry {
    std::string repPath;                 // genome file used as the representative
    std::string name;                    // node name (cluster_L..) or the leaf's name
    std::vector<std::string> leafPaths;  // leaves represented (one path for a leaf standing for itself)
    bool isNode;                         // true for an internal node, false for a self-represented leaf
};

/**
 * Partition of the leaves at every internal level of the tree.
 *
 * buildHierarchicalTree creates a node only at the level where a merge happens,
 * so a node whose parent sits two levels up is implicitly present at the level in
 * between, an orphan node (no parent) is present at every coarser level, and a
 * leaf that never merged is present, as itself, at every level. A node therefore
 * spans the levels from just below its parent down to its own; a leaf is
 * represented at level L by the ancestor spanning L, or by itself if none does.
 *
 * Every leaf appears in exactly one entry at every level. That is what both the
 * index build and the query need: the build makes each level's signatures unique
 * against the complete set of genomes present at that level rather than against
 * whichever nodes happen to be numbered with it, and the query can narrow
 * candidates level by level knowing that a leaf absent from a level's entries
 * cannot exist. Both sides call this function so they cannot disagree.
 *
 * @return (level, entries) pairs, coarsest first, for levels 0..leafLevel-1.
 */
inline std::vector<std::pair<unsigned, std::vector<LevelEntry>>>
treeLevelPartition(const std::vector<TreeNode>& tree) {
    std::unordered_map<unsigned, size_t> byId;
    unsigned leafLevel = 0;
    for (size_t i = 0; i < tree.size(); i++) {
        byId[tree[i].id] = i;
        if (tree[i].isLeaf) leafLevel = std::max(leafLevel, tree[i].level);
    }

    std::vector<std::pair<unsigned, std::vector<LevelEntry>>> result;
    for (unsigned L = 0; L < leafLevel; L++) {
        // representative node id (or leaf id) -> entry
        std::map<unsigned, LevelEntry> entries;
        for (const TreeNode& leaf : tree) {
            if (!leaf.isLeaf) continue;
            // Walk up: the ancestor spanning L has level >= L and a parent with
            // level < L (or no parent). Ancestors are visited finest-first, so the
            // first one with level >= L whose parent is above L is the answer.
            const TreeNode* rep = nullptr;
            int pid = leaf.parentId;
            while (pid != -1) {
                const TreeNode& anc = tree[byId.at(static_cast<unsigned>(pid))];
                if (anc.level >= L) {
                    bool parentAbove = (anc.parentId == -1) ||
                                       tree[byId.at(static_cast<unsigned>(anc.parentId))].level < L;
                    if (parentAbove) { rep = &anc; break; }
                }
                pid = anc.parentId;
            }
            unsigned key = rep ? rep->id : leaf.id;
            auto it = entries.find(key);
            if (it == entries.end()) {
                LevelEntry e;
                e.repPath = rep ? rep->filePath : leaf.filePath;
                e.name = rep ? rep->name : leaf.name;
                e.isNode = (rep != nullptr);
                it = entries.emplace(key, std::move(e)).first;
            }
            it->second.leafPaths.push_back(leaf.filePath);
        }
        std::vector<LevelEntry> list;
        for (auto& kv : entries) list.push_back(std::move(kv.second));
        result.emplace_back(L, std::move(list));
    }
    return result;
}

/**
 * True when two levels partition the leaves identically (same representatives
 * standing for the same leaves). Such levels would produce the same sketch, so
 * the build writes only one and the query scores it once.
 */
inline bool samePartition(const std::vector<LevelEntry>& a, const std::vector<LevelEntry>& b) {
    if (a.size() != b.size()) return false;
    auto key = [](const std::vector<LevelEntry>& v) {
        std::vector<std::pair<std::string, std::vector<std::string>>> k;
        for (const auto& e : v) {
            std::vector<std::string> leaves = e.leafPaths;
            std::sort(leaves.begin(), leaves.end());
            k.emplace_back(e.repPath, std::move(leaves));
        }
        std::sort(k.begin(), k.end());
        return k;
    };
    return key(a) == key(b);
}

/**
 * Representative genome files per level, coarsest first, for hierarchical sketch
 * building. Thin wrapper over treeLevelPartition kept for callers that only need
 * the file lists.
 */
inline std::vector<std::pair<unsigned, std::vector<std::string>>>
getTreeLevelRefs(const std::vector<TreeNode>& tree) {
    std::vector<std::pair<unsigned, std::vector<std::string>>> result;
    for (const auto& [lvl, entries] : treeLevelPartition(tree)) {
        std::vector<std::string> refs;
        for (const auto& e : entries) refs.push_back(e.repPath);
        result.emplace_back(lvl, std::move(refs));
    }
    return result;
}

/**
 * Signature-selection parameters recorded next to a tree index
 * (<outdir>/build_params.tsv, one `key<TAB>value` per line).
 *
 * The query tool reads these to choose reporting defaults that match how the
 * signatures were selected. The one that matters today is the margin: with
 * --min-margin=2 --strict-margin every signature is at Hamming distance >= 2 from
 * every k-mer of every other reference, so a read from an abundant relative needs
 * two coincident sequencing errors to hit it and a single supporting read is
 * credible evidence. With margin-1 signatures one error suffices, and at high
 * coverage of a relative such reads outnumber the true ones, so the read cutoff
 * has to stay up. querysketch cannot tell the two indexes apart from the sketch
 * file alone; this record is how it knows.
 */
struct TreeBuildParams {
    unsigned kmerLen = 0;
    int sketchnum = 0;
    unsigned minMargin = 1;
    bool strictMargin = false;
    unsigned maxHomopolymer = 0;
    bool groupSignatures = false; // within-group leaf signatures built (sketch_group.tsv)
    std::string thresholds;      // comma-separated, coarsest first
    std::string bloomFilter;     // filter layout the index was built with (see bloomFilterName)
    bool loaded = false;         // false when the file is absent (pre-2.0 index)
};

/**
 * Name of the Bloom filter layout this binary was compiled with. The blocked and
 * standard filters draw different signatures from the same candidate pool (their
 * false positives fall on different k-mers), so an index should be queried by a
 * binary using the layout that built it. The build records this name and the
 * query compares it; a mismatch is a warning, not an error, since the signatures
 * in the sketch file are what they are and the query reads them verbatim. The
 * only query-side filter is the --solid read filter, which is sample-local.
 */
inline const char* bloomFilterName() {
#ifdef UNIQSKETCH_STANDARD_BF
    return "standard";
#else
    return "blocked";
#endif
}

inline std::string treeBuildParamsPath(const std::string& treeDir) {
    return treeDir + "/build_params.tsv";
}

inline void writeTreeBuildParams(const std::string& path, const TreeBuildParams& p) {
    std::ofstream out(path);
    out << "kmer\t" << p.kmerLen << "\n"
        << "cov\t" << p.sketchnum << "\n"
        << "min_margin\t" << p.minMargin << "\n"
        << "strict_margin\t" << (p.strictMargin ? 1 : 0) << "\n"
        << "max_homopolymer\t" << p.maxHomopolymer << "\n"
        << "group_signatures\t" << (p.groupSignatures ? 1 : 0) << "\n"
        << "tree\t" << p.thresholds << "\n"
        << "bloom_filter\t" << p.bloomFilter << "\n";
}

/**
 * The leaf numbering of a tree index: `<outdir>/leaves.tsv`, one row per leaf
 * reference in refId order (`ref_id`, `name`, `file_path`).
 *
 * The leaf sketch, the within-group sketch and its carriers column all refer to
 * leaves by this refId, but a leaf with no signatures of its own has no line in
 * either sketch, so until this file existed the query could not name it: it
 * could neither report it nor recognise it as the carrier of a within-group
 * signature. On a 5 493-genome E. coli index 416 leaves were in that position,
 * and one of them (a K-12 derivative at 90 % abundance) put 5 667 of its reads
 * into a relative's within-group signatures because the carrier mask could not
 * see it. Written by uniqsketch; loadTreeLeaves returns an empty vector for an
 * index that predates it, and the query then falls back to the names the
 * sketches carry.
 */
inline std::string treeLeavesPath(const std::string& treeDir) {
    return treeDir + "/leaves.tsv";
}

inline void writeTreeLeaves(const std::string& path, const std::vector<std::string>& refFiles) {
    std::ofstream out(path);
    out << "ref_id\tname\tfile_path\n";
    for (size_t i = 0; i < refFiles.size(); i++)
        out << i << "\t" << getBaseId(refFiles[i]) << "\t" << refFiles[i] << "\n";
}

/** refId -> leaf name, or an empty vector when the file is absent. */
inline std::vector<std::string> loadTreeLeaves(const std::string& path) {
    std::vector<std::string> names;
    std::ifstream in(path);
    if (!in) return names;
    std::string line;
    std::getline(in, line);                       // header
    while (std::getline(in, line)) {
        size_t a = line.find('\t');
        if (a == std::string::npos) continue;
        size_t b = line.find('\t', a + 1);
        unsigned id = 0;
        try { id = static_cast<unsigned>(std::stoul(line.substr(0, a))); } catch (const std::exception&) { continue; }
        if (id >= names.size()) names.resize(id + 1);
        names[id] = line.substr(a + 1, (b == std::string::npos ? line.size() : b) - a - 1);
    }
    return names;
}

inline TreeBuildParams loadTreeBuildParams(const std::string& path) {
    TreeBuildParams p;
    std::ifstream in(path);
    if (!in) return p;
    std::string line;
    while (std::getline(in, line)) {
        size_t tab = line.find('\t');
        if (tab == std::string::npos) continue;
        std::string key = line.substr(0, tab), val = line.substr(tab + 1);
        try {
            if (key == "kmer") p.kmerLen = static_cast<unsigned>(std::stoul(val));
            else if (key == "cov") p.sketchnum = std::stoi(val);
            else if (key == "min_margin") p.minMargin = static_cast<unsigned>(std::stoul(val));
            else if (key == "strict_margin") p.strictMargin = (val == "1");
            else if (key == "max_homopolymer") p.maxHomopolymer = static_cast<unsigned>(std::stoul(val));
            else if (key == "group_signatures") p.groupSignatures = (val == "1");
            else if (key == "tree") p.thresholds = val;
            else if (key == "bloom_filter") p.bloomFilter = val;
        } catch (const std::exception&) {
            // A malformed value leaves that field at its default; the record is
            // advisory, so this must never take the query down.
        }
    }
    p.loaded = true;
    return p;
}

/**
 * Load a hierarchical tree from a tree TSV file.
 *
 * @param treeFile  Path to the tree TSV file.
 * @return Vector of TreeNode.
 */
inline std::vector<TreeNode> loadTree(const std::string& treeFile) {
    std::vector<TreeNode> tree;
    std::ifstream in(treeFile);
    std::string line;
    std::getline(in, line); // skip header

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        // Split on tabs explicitly: the children column (index 5) is empty for leaf
        // nodes, and stream (>>) extraction would silently skip that empty field and
        // shift every later column left (reading the file path as the children list).
        std::vector<std::string> col;
        std::string field;
        std::istringstream ls(line);
        while (std::getline(ls, field, '\t')) {
            col.push_back(field);
        }
        if (col.size() < 5) continue;  // malformed row

        TreeNode node;
        node.id = static_cast<unsigned>(std::stoul(col[0]));
        node.name = col[1];
        node.level = static_cast<unsigned>(std::stoul(col[2]));
        node.parentId = std::stoi(col[3]);
        node.isLeaf = (col[4] == "yes");
        std::string childrenStr = (col.size() > 5) ? col[5] : "";
        node.filePath = (col.size() > 6) ? col[6] : "";

        if (!childrenStr.empty() && childrenStr != "-") {
            std::istringstream cs(childrenStr);
            std::string token;
            while (std::getline(cs, token, ',')) {
                if (!token.empty()) {
                    node.children.push_back(std::stoi(token));
                }
            }
        }
        tree.push_back(std::move(node));
    }
    return tree;
}

#endif // TREESKETCHUTIL_HPP_