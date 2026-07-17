#ifndef TREESKETCHUTIL_HPP_
#define TREESKETCHUTIL_HPP_

#include <algorithm>
#include <climits>
#include <cstring>
#include <fstream>
#include <iostream>
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
 * Get the set of references at each level of the tree for hierarchical sketch building.
 * Returns a vector of (level, refFiles) pairs, ordered from coarsest to finest.
 *
 * @param tree  The hierarchical tree.
 * @return Vector of pairs: (level, vector of file paths for that level's representatives).
 */
inline std::vector<std::pair<unsigned, std::vector<std::string>>>
getTreeLevelRefs(const std::vector<TreeNode>& tree) {
    std::unordered_map<unsigned, std::vector<std::string>> levelRefs;

    for (const auto& node : tree) {
        // Internal nodes at each level serve as representatives for that level
        if (!node.isLeaf && !node.filePath.empty()) {
            levelRefs[node.level].push_back(node.filePath);
        }
    }

    // Also add leaf-only references (those not part of any cluster)
    for (const auto& node : tree) {
        if (node.isLeaf && node.parentId == -1) {
            // Find the deepest level
            unsigned maxLevel = 0;
            for (const auto& [lvl, refs] : levelRefs) {
                maxLevel = std::max(maxLevel, lvl);
            }
            levelRefs[maxLevel].push_back(node.filePath);
        }
    }

    // Sort by level
    std::vector<std::pair<unsigned, std::vector<std::string>>> result;
    for (auto& [lvl, refs] : levelRefs) {
        result.emplace_back(lvl, std::move(refs));
    }
    std::sort(result.begin(), result.end());
    return result;
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