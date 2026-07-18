#include <cassert>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TreeSketchUtil.hpp"

// Deterministic pseudo-random DNA (varied, so its k-mers are mostly distinct).
static std::string randDNA(size_t n, uint64_t seed) {
    static const char B[4] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.reserve(n);
    uint64_t x = seed;
    for (size_t i = 0; i < n; i++) {
        x = x * 6364136223846793005ULL + 1442695040888963407ULL;
        s.push_back(B[(x >> 33) & 3]);
    }
    return s;
}

static void writeFasta(const std::string& path, const std::string& seq) {
    std::ofstream o(path);
    o << ">seq\n";
    for (size_t i = 0; i < seq.size(); i += 70) o << seq.substr(i, 70) << "\n";
}

// loadTree must parse a leaf row whose children column is EMPTY without shifting the
// file-path field left (the tab-parse bug that motivated the explicit split).
void testLoadTree() {
    std::cerr << "START: Test loadTree... \t";
    const char* tf = "tree_test_index.tsv";
    {
        std::ofstream o(tf);
        o << "node_id\tname\tlevel\tparent_id\tis_leaf\tchildren\tfile_path\n";
        o << "0\tleafA\t2\t2\tyes\t\t/data/A.fa\n";
        o << "1\tleafB\t2\t2\tyes\t\t/data/B.fa\n";
        o << "2\tcluster_L0_A\t0\t-1\tno\t0,1\t/data/A.fa\n";
    }
    std::vector<TreeNode> tree = loadTree(tf);
    assert(tree.size() == 3);
    assert(tree[0].isLeaf && tree[0].parentId == 2 && tree[0].children.empty());
    assert(tree[0].filePath == "/data/A.fa");   // not shifted by the empty children column
    assert(!tree[2].isLeaf && tree[2].children.size() == 2);
    assert(tree[2].children[0] == 0 && tree[2].children[1] == 1);
    assert(tree[2].filePath == "/data/A.fa");
    std::remove(tf);
    std::cerr << "PASSED: Test loadTree\n";
}

// getTreeLevelRefs must group internal-node representatives by their level.
void testGetTreeLevelRefs() {
    std::cerr << "START: Test getTreeLevelRefs... \t";
    std::vector<TreeNode> tree;
    TreeNode a; a.id = 0; a.name = "A"; a.filePath = "/d/A.fa"; a.level = 1; a.parentId = 2; a.isLeaf = true; tree.push_back(a);
    TreeNode b; b.id = 1; b.name = "B"; b.filePath = "/d/B.fa"; b.level = 1; b.parentId = 2; b.isLeaf = true; tree.push_back(b);
    TreeNode c; c.id = 2; c.name = "cl"; c.filePath = "/d/A.fa"; c.level = 0; c.parentId = -1; c.isLeaf = false; c.children = {0, 1}; tree.push_back(c);

    auto levels = getTreeLevelRefs(tree);
    bool foundL0 = false;
    for (const auto& lr : levels)
        if (lr.first == 0) { assert(lr.second.size() == 1 && lr.second[0] == "/d/A.fa"); foundL0 = true; }
    assert(foundL0);
    std::cerr << "PASSED: Test getTreeLevelRefs\n";
}

// buildHierarchicalTree must group near-identical references and keep distinct ones apart.
void testBuildHierarchicalTree() {
    std::cerr << "START: Test buildHierarchicalTree... \t";
    std::string ab = randDNA(1000, 12345);   // A and B share this sequence exactly
    std::string cc = randDNA(1000, 99999);   // C is unrelated
    writeFasta("tree_A.fa", ab);
    writeFasta("tree_B.fa", ab);
    writeFasta("tree_C.fa", cc);

    std::vector<std::string> refs = {"tree_A.fa", "tree_B.fa", "tree_C.fa"};
    std::vector<unsigned> thresholds = {100};   // merge when < 100 distinguishing k-mers
    auto tree = buildHierarchicalTree(refs, thresholds, /*k*/ 21, /*nhash*/ 1,
                                      /*bits*/ 16, /*dbfSize*/ 100000, "tree_out.tsv");

    int leaves = 0, internal = 0;
    for (const auto& n : tree) { if (n.isLeaf) ++leaves; else ++internal; }
    assert(leaves == 3);
    assert(internal >= 1);

    // Some internal node groups A and B but not C.
    bool abGrouped = false;
    for (const auto& n : tree) {
        if (n.isLeaf) continue;
        bool hasA = false, hasB = false, hasC = false;
        for (unsigned c : n.children) {
            if (c >= tree.size() || !tree[c].isLeaf) continue;
            if (tree[c].name == "tree_A") hasA = true;
            if (tree[c].name == "tree_B") hasB = true;
            if (tree[c].name == "tree_C") hasC = true;
        }
        if (hasA && hasB) { abGrouped = true; assert(!hasC); }
    }
    assert(abGrouped);

    std::remove("tree_A.fa"); std::remove("tree_B.fa"); std::remove("tree_C.fa");
    std::remove("tree_out.tsv");
    std::cerr << "PASSED: Test buildHierarchicalTree\n";
}

int main() {
    testLoadTree();
    testGetTreeLevelRefs();
    testBuildHierarchicalTree();
    std::cerr << "TreeSketch: All tests PASSED!\n\n";
    return 0;
}
