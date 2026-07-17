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

    // Pre-compute k-mer hashes
    std::vector<std::vector<uint64_t>> hashSets(n);
    std::vector<std::string> names(n);

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < n; i++) {
        gzFile fp = gzopen(refFiles[i].c_str(), "r");
        kseq_t* ks = kseq_init(fp);
        while (kseq_read(ks) >= 0) {
            std::string seq(ks->seq.s, ks->seq.l);
            ntHashIterator itr(seq, nhash, kmerLen);
            while (itr != itr.end()) {
                const uint64_t* hVec = *itr;
                for (unsigned h = 0; h < nhash; h++) {
                    hashSets[i].push_back(hVec[h]);
                }
                ++itr;
            }
        }
        kseq_destroy(ks);
        gzclose(fp);
        names[i] = getBaseId(refFiles[i]);
    }

    // Compute pairwise unique k-mer counts (symmetric: store both directions)
    // uniqueKmers[i][j] = number of k-mers in j not found in i's BF
    std::vector<std::vector<size_t>> uniqueKmers(n, std::vector<size_t>(n, 0));
    size_t bfBits = bits * dbfSize;

    #pragma omp parallel
    {
        BloomFilter dbFilter(bfBits, nhash, kmerLen);
        unsigned char* filterPtr = dbFilter.getFilter();
        size_t filterBytes = (bfBits + CHAR_BIT - 1) / CHAR_BIT;

        #pragma omp for schedule(dynamic)
        for (unsigned i = 0; i < n; i++) {
            std::memset(filterPtr, 0, filterBytes);
            for (size_t h = 0; h < hashSets[i].size(); h += nhash) {
                dbFilter.insert(&hashSets[i][h]);
            }
            for (unsigned j = 0; j < n; j++) {
                if (i == j) continue;
                size_t countUnique = 0;
                for (size_t h = 0; h < hashSets[j].size(); h += nhash) {
                    if (!dbFilter.contains(&hashSets[j][h])) {
                        ++countUnique;
                    }
                }
                uniqueKmers[i][j] = countUnique;
            }
        }
    }

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

    // For each threshold level (coarsest to finest), cluster and create internal nodes
    // Track which tree node each original reference currently belongs to
    std::vector<unsigned> currentNodeId(n);
    std::iota(currentNodeId.begin(), currentNodeId.end(), 0);

    for (unsigned lvl = 0; lvl < thresholds.size(); lvl++) {
        unsigned threshold = thresholds[lvl];

        // Complete-linkage agglomerative clustering:
        // Two clusters merge only if the MAXIMUM unique k-mer count
        // between any cross-cluster pair is below the threshold.
        // This prevents chaining artifacts from single-linkage.

        // Initialize: each reference is its own cluster
        std::vector<std::vector<unsigned>> clusters(n);
        for (unsigned i = 0; i < n; i++) {
            clusters[i] = {i};
        }

        // Iteratively merge the closest pair of clusters (complete-linkage distance)
        bool merged = true;
        while (merged) {
            merged = false;
            size_t bestDist = threshold; // only merge if below threshold
            int bestA = -1, bestB = -1;

            for (size_t a = 0; a < clusters.size(); a++) {
                if (clusters[a].empty()) continue;
                for (size_t b = a + 1; b < clusters.size(); b++) {
                    if (clusters[b].empty()) continue;

                    // Complete-linkage: max distance between any pair across clusters
                    size_t maxDist = 0;
                    for (unsigned i : clusters[a]) {
                        for (unsigned j : clusters[b]) {
                            size_t d = std::min(uniqueKmers[i][j], uniqueKmers[j][i]);
                            maxDist = std::max(maxDist, d);
                        }
                    }

                    if (maxDist < bestDist) {
                        bestDist = maxDist;
                        bestA = a;
                        bestB = b;
                    }
                }
            }

            if (bestA >= 0) {
                // Merge bestB into bestA
                for (unsigned m : clusters[bestB]) {
                    clusters[bestA].push_back(m);
                }
                clusters[bestB].clear();
                merged = true;
            }
        }

        // Collect non-empty clusters into groups
        std::unordered_map<unsigned, std::vector<unsigned>> groups;
        unsigned gid = 0;
        for (size_t c = 0; c < clusters.size(); c++) {
            if (!clusters[c].empty()) {
                groups[gid++] = std::move(clusters[c]);
            }
        }

        // Create internal nodes for clusters with >1 member
        for (auto& [root, members] : groups) {
            if (members.size() <= 1) continue;

            // Pick representative (most k-mers)
            unsigned repIdx = members[0];
            size_t maxK = hashSets[repIdx].size();
            for (unsigned m : members) {
                if (hashSets[m].size() > maxK) {
                    maxK = hashSets[m].size();
                    repIdx = m;
                }
            }

            TreeNode internal;
            internal.id = tree.size();
            internal.name = "cluster_L" + std::to_string(lvl) + "_" + names[repIdx];
            internal.filePath = refFiles[repIdx];
            internal.level = lvl;
            internal.parentId = -1;
            internal.isLeaf = false;

            // Set children: the current nodes of each member
            for (unsigned m : members) {
                unsigned childId = currentNodeId[m];
                if (tree[childId].parentId == -1) {
                    tree[childId].parentId = internal.id;
                    internal.children.push_back(childId);
                }
            }

            tree.push_back(std::move(internal));

            // Update current node mapping
            unsigned newId = tree.size() - 1;
            for (unsigned m : members) {
                currentNodeId[m] = newId;
            }
        }
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