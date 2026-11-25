#ifndef TREEFILEUTIL_HPP_
#define TREEFILEUTIL_HPP_

/**
 * The tree of a `uniqsketch --tree-file` build.
 *
 * A tree build clusters the references itself (buildHierarchicalTree) and then
 * builds a sketch for every level of the tree. With --tree-file the first step is
 * replaced: the tree is read from a file in the format of tree_index.tsv, so that
 * another clustering can be taken through the level builds, the within-group
 * signatures and a query as they are. The build takes from the file which
 * references are grouped at which level, and nothing else: the paths are those of
 * the reference list, and the member that stands for a node is decided as the
 * merges decide it (treeFileRepresentatives).
 *
 * Of the three tools only uniqsketch includes this header. TreeSketchUtil.hpp is
 * compiled into querysketch as well, and with these functions in it, called by
 * nothing there, g++ 7.3.1 gave another querysketch binary (the same section
 * sizes and the same 633 symbols, 24 of its function clones numbered differently
 * and 24 726 of the 275 074 bytes of .text changed). Kept here, the option changes
 * no file that querysketch or comparesketch compiles.
 */

#include <algorithm>
#include <climits>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "SequenceUtil.hpp"
#include "ClusterUtil.hpp"
#include "TreeSketchUtil.hpp"

/**
 * The k-mer count clusterMinhashSketch returns, without the sketch: the number of
 * windows of k bases, within one sequence, in which every character has an ntHash
 * seed. Those are the windows ntHashIterator yields (it restarts after a character
 * without a seed), read through kseq as the sketch pass reads them, so the count is
 * the same number, and nothing is hashed. A tree read from a file has no sketches,
 * and the count is what decides which member stands for a node. On the 120 E. coli
 * genomes of the validation index (617 M positions) it takes 0.6 s on one thread
 * where the sketch pass takes 5.0 s, and agrees with it on every genome.
 */
inline size_t countKmerPositions(const std::string& fPath, unsigned kmerLen) {
    size_t counted = 0;
    gzFile fp = gzopen(fPath.c_str(), "r");
    kseq_t* ks = kseq_init(fp);
    while (kseq_read(ks) >= 0) {
        size_t run = 0;                    // characters with a seed since the last one without
        for (size_t i = 0; i < ks->seq.l; i++) {
            run = seedTab[static_cast<unsigned char>(ks->seq.s[i])] != seedN ? run + 1 : 0;
            if (run >= kmerLen) ++counted;
        }
    }
    kseq_destroy(ks);
    gzclose(fp);
    return counted;
}

/**
 * What keeps a tree read from a file from standing in for the tree
 * buildHierarchicalTree hands the build, or an empty string. `tree` is the file as
 * loadTreeFile read it and has passed treeDefect (ids used once, parents in the
 * file, no cycle); `refFiles` are the references of the build in list order;
 * `numThresholds` is the number of --tree thresholds.
 *
 * The file must describe the references of this build in the shape the tool's own
 * trees have:
 *  - its leaves are the references, each once, leaf i being reference i of the
 *    list under the name the sketches give it (getBaseId). The build lists its
 *    references in the order of the tree's leaves and the numbering of the leaf
 *    sketch follows from that order (uniqsketch.cpp), so a file that numbered the
 *    leaves differently would make another index of the same grouping;
 *  - the leaves are at level `numThresholds` and every node above them at a
 *    smaller one: the index has the levels the thresholds name;
 *  - a node has two or more children, each at a deeper level and each naming it
 *    as its parent, and a node that names a parent is among that parent's
 *    children. A group that is the same at several levels is one node, at the
 *    finest of them (treeLevelPartition lets it span the others), as the merges
 *    make it; a chain of single children would make it a node per level, each
 *    with a representative of its own.
 * Two or more children at deeper levels, linked in both directions, also mean
 * that every node has two or more leaves under it.
 */
inline std::string treeFileDefect(const std::vector<TreeNode>& tree,
                                  const std::vector<std::string>& refFiles, size_t numThresholds) {
    const size_t n = refFiles.size();
    std::unordered_map<unsigned, size_t> byId;
    for (size_t i = 0; i < tree.size(); i++) byId[tree[i].id] = i;
    auto label = [](const TreeNode& x) { return std::to_string(x.id) + " (" + x.name + ")"; };

    std::unordered_map<std::string, size_t> listed;        // reference name -> position in the list
    for (size_t i = 0; i < n; i++) listed[getBaseId(refFiles[i])] = i;
    size_t leaves = 0;
    for (const TreeNode& x : tree) leaves += x.isLeaf;
    if (leaves != n)
        return "it has " + std::to_string(leaves) + " leaves, the list has " + std::to_string(n) + " references";
    for (const TreeNode& x : tree) {
        if (!x.isLeaf) continue;
        auto it = listed.find(x.name);
        if (it == listed.end()) return "leaf " + label(x) + " is not a reference of the list";
        if (it->second != x.id)
            return "leaf " + label(x) + " is reference " + std::to_string(it->second) +
                   " of the list; a leaf's id is its position in the list, from 0";
        if (x.level != numThresholds)
            return "leaf " + label(x) + " is at level " + std::to_string(x.level) + "; with " +
                   std::to_string(numThresholds) + " threshold(s) in --tree the leaves are at level " +
                   std::to_string(numThresholds);
        if (!x.children.empty()) return "leaf " + label(x) + " has children";
    }
    // As many leaves as references, each under the name of the reference at its id:
    // every reference once.

    for (const TreeNode& x : tree) {
        if (x.parentId != -1) {
            const TreeNode& p = tree[byId.at(static_cast<unsigned>(x.parentId))];
            if (p.isLeaf) return "node " + label(x) + " has the leaf " + label(p) + " as its parent";
            if (std::find(p.children.begin(), p.children.end(), x.id) == p.children.end())
                return "node " + label(x) + " names " + label(p) + " as its parent and is not among its children";
        }
        if (x.isLeaf) continue;
        if (x.level >= numThresholds)
            return "node " + label(x) + " is at level " + std::to_string(x.level) + "; with " +
                   std::to_string(numThresholds) + " threshold(s) in --tree the nodes are at levels 0 to " +
                   std::to_string(numThresholds - 1);
        if (x.children.size() < 2)
            return "node " + label(x) + " has " + (x.children.empty() ? "no child" : "one child") +
                   "; a node groups two or more (a group that is the same at several levels is one node, "
                   "at the finest of them)";
        std::unordered_set<unsigned> seen;
        for (unsigned c : x.children) {
            auto it = byId.find(c);
            if (it == byId.end())
                return "node " + label(x) + " has the child " + std::to_string(c) + ", which is not in the file";
            const TreeNode& child = tree[it->second];
            if (!seen.insert(c).second) return "node " + label(x) + " has the child " + label(child) + " twice";
            if (child.parentId != static_cast<int>(x.id))
                return "node " + label(x) + " has the child " + label(child) + ", whose parent is " +
                       (child.parentId == -1 ? std::string("not given") : std::to_string(child.parentId));
            if (child.level <= x.level)
                return "node " + label(x) + ", at level " + std::to_string(x.level) + ", has the child " +
                       label(child) + " at level " + std::to_string(child.level) + "; a child is at a deeper level";
        }
    }
    return std::string();
}

/**
 * A column of a tree file as a plain integer from 0 to `most`: digits and nothing
 * else, or false. loadTree reads the numbers of an index the tool wrote itself
 * with std::stoul and std::stoi, which take the leading digits and stop, accept
 * a sign and leading blanks, and wrap a value above the type. A file written by
 * hand or by a script is held to more: read with those two, the 8-phage tree
 * with the node id 5x, the levels 3.9 and 1e2 and the children 7;9,6,0,4 built
 * the stored index with exit status 0, and the id 4294967301 was leaf 5.
 */
inline bool treeFileNumber(const std::string& field, unsigned long long most, unsigned long long& value) {
    if (field.empty() || field.size() > 10 || field.find_first_not_of("0123456789") != std::string::npos)
        return false;
    value = std::stoull(field);
    return value <= most;
}

/**
 * Read the tree of a --tree-file build from `path`, checked (treeDefect,
 * treeFileDefect) and laid out as buildHierarchicalTree lays out its own, the
 * leaves first in list order, then the nodes by id. Returns an empty string, or
 * what is wrong with the file, the line named where it is one line.
 *
 * The file is read here, once, so a pipe is read like a file, and not through
 * loadTree, which reads an index the tool wrote: it drops the first line unread,
 * passes over a row of fewer than five columns, takes any word but yes for no,
 * keeps a carriage return in the last column and stops the run with "damaged
 * index?" on a number it cannot convert. Read that way, a file saved with CRLF
 * line endings is a tree without a leaf when its rows end after the fifth
 * column and a row that cannot be parsed when they end after the sixth, and one
 * with is_leaf written Yes or 1 is a tree without a leaf. Here a carriage
 * return before the line feed is part of the line ending, as in a list file,
 * and a row is told what is wrong with it. A node id, a parent id and a child
 * are at most 2147483647: loadTree, which reads the tree_index.tsv this build
 * writes, takes the last two as int.
 *
 * Of a leaf the file's id, name, level and parent are used, of a node its id,
 * level, parent and children in the order given. A leaf's path is that of the
 * list, whatever the file says (the file may come from another machine). A node's
 * name and path say which member stands for it; they are not read, and the nodes
 * come back without them (treeFileRepresentatives).
 */
inline std::string loadTreeFile(const std::string& path, const std::vector<std::string>& refFiles,
                                size_t numThresholds, std::vector<TreeNode>& out) {
    std::ifstream in(path);
    std::string line;
    if (!std::getline(in, line)) return "cannot be read, or is empty";
    // Without the header the first node would be taken for it, and the complaint
    // would be about a missing leaf.
    if (line.compare(0, 8, "node_id\t") != 0)
        return "its first line is not the header of a tree_index.tsv (node_id, name, level, parent_id, "
               "is_leaf, children, file_path; tab-separated)";
    std::vector<TreeNode> nodes;
    for (size_t lineNo = 2; std::getline(in, line); lineNo++) {
        // A carriage return before the line feed is part of the line ending, as in
        // a list file.
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::vector<std::string> col;
        size_t from = 0;
        for (size_t tab; (tab = line.find('\t', from)) != std::string::npos; from = tab + 1)
            col.push_back(line.substr(from, tab - from));
        col.push_back(line.substr(from));
        const std::string at = "line " + std::to_string(lineNo) + ": ";
        if (col.size() < 5)
            return at + std::to_string(col.size()) + " tab-separated column(s); a row has node_id, name, level, "
                   "parent_id and is_leaf, then children and file_path";
        TreeNode node;
        unsigned long long v = 0;
        if (!treeFileNumber(col[0], INT_MAX, v))
            return at + "node_id `" + col[0] + "' is not a plain integer (at most " + std::to_string(INT_MAX) + ")";
        node.id = static_cast<unsigned>(v);
        node.name = col[1];
        if (!treeFileNumber(col[2], INT_MAX, v)) return at + "level `" + col[2] + "' is not a plain integer";
        node.level = static_cast<unsigned>(v);
        if (col[3] == "-1") node.parentId = -1;
        else if (treeFileNumber(col[3], INT_MAX, v)) node.parentId = static_cast<int>(v);
        else return at + "parent_id `" + col[3] + "' is neither a node id nor -1";
        if (col[4] != "yes" && col[4] != "no") return at + "is_leaf `" + col[4] + "' is neither yes nor no";
        node.isLeaf = col[4] == "yes";
        if (col.size() > 5 && !col[5].empty() && col[5] != "-") {
            std::istringstream cs(col[5]);
            std::string token;
            while (std::getline(cs, token, ',')) {
                if (!treeFileNumber(token, INT_MAX, v))
                    return at + "children `" + col[5] + "' is not a comma-separated list of node ids";
                node.children.push_back(static_cast<unsigned>(v));
            }
        }
        if (col.size() > 6) node.filePath = col[6];
        nodes.push_back(std::move(node));
    }
    if (in.bad()) return "could not be read to its end";
    std::string defect = treeDefect(nodes);
    if (defect.empty()) defect = treeFileDefect(nodes, refFiles, numThresholds);
    if (!defect.empty()) return defect;

    out.assign(refFiles.size(), TreeNode());
    std::vector<const TreeNode*> internal;
    for (const TreeNode& x : nodes) {
        if (!x.isLeaf) { internal.push_back(&x); continue; }
        TreeNode& leaf = out[x.id];
        leaf.id = x.id;
        leaf.name = x.name;
        leaf.filePath = refFiles[x.id];
        leaf.level = x.level;
        leaf.parentId = x.parentId;
        leaf.isLeaf = true;
    }
    std::sort(internal.begin(), internal.end(), [](const TreeNode* a, const TreeNode* b) { return a->id < b->id; });
    for (const TreeNode* x : internal) {
        out.push_back(*x);
        out.back().name.clear();
        out.back().filePath.clear();
    }
    return std::string();
}

/**
 * Name the nodes of a tree read from a file (loadTreeFile) after their
 * representatives, by the rule of buildHierarchicalTree: a node stands for its
 * children through the largest of their representatives, a leaf being its own, and
 * of several largest ones through the first in the order of its children column
 * (the merges compare with a strict >, in that order). `kmerCount[i]` is the k-mer
 * count of leaf i.
 *
 * The file's own name and path of a node are not used. The tree_index.tsv of a
 * finished index names a node after the member its signatures come from, which the
 * level builds choose (chooseNodeRepresentative) starting from the genome the
 * merges named, and that starting point is part of the build: it is the genome
 * whose repeats the level's filter scan excludes (loadBFfaGroup, applied before
 * the choice), it heads the group the choice is made in, and
 * node_representatives.tsv records it as `tree_genome`. Started from the member an
 * earlier build chose, the build gives another index than the one the file came
 * from. Measured on the 120-genome E. coli index (27 nodes, 21 renamed), rebuilt
 * from its own tree_index.tsv: with the genome the file names as the starting
 * point the same 27 members are chosen, but 44 of the 653 files differ: the
 * level-0 and level-1 sketches (12 866 of 234 030 and 3 094 of 328 437 lines), 37
 * candidate files of those levels with their two stat files, the two reachability
 * tables and node_representatives.tsv. With the rule of the merges every file is
 * the index's, and so it is from the tree as first written, before the level
 * builds rename its nodes.
 */
inline void treeFileRepresentatives(std::vector<TreeNode>& tree, const std::vector<size_t>& kmerCount) {
    std::unordered_map<unsigned, size_t> byId;
    std::vector<size_t> order;                             // the nodes, deepest level first: children before parents
    for (size_t i = 0; i < tree.size(); i++) {
        byId[tree[i].id] = i;
        if (!tree[i].isLeaf) order.push_back(i);
    }
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return tree[a].level > tree[b].level; });
    std::unordered_map<unsigned, unsigned> repOf;          // node id -> leaf id of its representative
    for (size_t i : order) {
        TreeNode& node = tree[i];
        auto rep = [&](unsigned child) {
            const TreeNode& c = tree[byId.at(child)];
            return c.isLeaf ? c.id : repOf.at(c.id);
        };
        unsigned best = rep(node.children[0]);
        for (unsigned c : node.children)
            if (kmerCount[rep(c)] > kmerCount[best]) best = rep(c);
        repOf[node.id] = best;
        const TreeNode& leaf = tree[byId.at(best)];
        node.filePath = leaf.filePath;
        node.name = "cluster_L" + std::to_string(node.level) + "_" + leaf.name;
    }
}

/**
 * treeFileRepresentatives with the counts taken from the genomes
 * (countKmerPositions), one genome per thread. A leaf under no node is not read.
 */
inline void treeFileRepresentatives(std::vector<TreeNode>& tree, unsigned kmerLen) {
    size_t n = 0;
    while (n < tree.size() && tree[n].isLeaf) ++n;         // the leaves come first, leaf i at i
    std::vector<size_t> kmerCount(n, 0);
    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < n; i++)
        if (tree[i].parentId != -1) kmerCount[i] = countKmerPositions(tree[i].filePath, kmerLen);
    treeFileRepresentatives(tree, kmerCount);
}

/**
 * Record in build_params.tsv that the tree of the index was read from `treeFile`:
 * one line, `tree_file<TAB>path as given`, after those writeTreeBuildParams wrote.
 * Without it the `tree` line would read as if its thresholds had made the tree.
 * Only a --tree-file build calls this, so the file of any other build is the one
 * it always was. The line has no reader: loadTreeBuildParams passes over a key it
 * does not know, as it has since the file exists (1b3c01b), so it changes nothing
 * for querysketch, this one or an earlier one.
 */
inline void recordTreeFile(const std::string& buildParamsPath, const std::string& treeFile) {
    std::ofstream out(buildParamsPath, std::ios::app);
    out << "tree_file\t" << treeFile << "\n";
}

#endif // TREEFILEUTIL_HPP_
