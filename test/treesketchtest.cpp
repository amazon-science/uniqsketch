#include <cassert>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
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

// build_params.tsv round-trips, and its absence is reported rather than guessed:
// querysketch lowers the read cutoff only for a record that says strict margin 2.
void testBuildParamsRoundTrip() {
    std::cerr << "START: Test build params round trip... \t";
    const char* dir = ".";
    std::string path = treeBuildParamsPath(dir);
    TreeBuildParams p;
    p.kmerLen = 81; p.sketchnum = 500; p.minMargin = 2; p.strictMargin = true;
    p.maxHomopolymer = 6; p.thresholds = "10000,1000,100";
    writeTreeBuildParams(path, p);
    TreeBuildParams q = loadTreeBuildParams(path);
    assert(q.loaded);
    assert(q.kmerLen == 81 && q.sketchnum == 500);
    assert(q.minMargin == 2 && q.strictMargin);
    assert(q.maxHomopolymer == 6 && q.thresholds == "10000,1000,100");
    std::remove(path.c_str());

    TreeBuildParams missing = loadTreeBuildParams(path);
    assert(!missing.loaded);
    assert(missing.minMargin == 1 && !missing.strictMargin);  // pre-2.0 index: no margin assumed
    std::cerr << "PASSED: Test build params round trip\n";
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

// treeLevelPartition must cover every leaf exactly once at every internal level:
// a node stands for its leaves from just below its parent down to its own level,
// an orphan node from the top down, and a leaf that never merged stands for
// itself everywhere. This is the shape that broke the query on the E. coli set
// (singletons and orphans were never scored at the coarse levels and vanished).
void testTreeLevelPartition() {
    std::cerr << "START: Test treeLevelPartition... \t";
    // Thresholds T0 > T1 > T2 -> internal levels 0,1,2; leaves at level 3.
    //   node N (level 2) = {A, B}, parent P (level 0) = {A, B, C}   -> N spans 1..2
    //   node O (level 1) = {D, E}, orphan                            -> O spans 0..1
    //   leaf F never merges                                          -> itself at 0,1,2
    std::vector<TreeNode> t;
    auto leaf = [&](unsigned id, const char* n, int parent) {
        TreeNode x; x.id = id; x.name = n; x.filePath = std::string("/d/") + n + ".fa";
        x.level = 3; x.parentId = parent; x.isLeaf = true; t.push_back(x);
    };
    leaf(0, "A", 6); leaf(1, "B", 6); leaf(2, "C", 7); leaf(3, "D", 8); leaf(4, "E", 8); leaf(5, "F", -1);
    TreeNode N; N.id = 6; N.name = "N"; N.filePath = "/d/A.fa"; N.level = 2; N.parentId = 7; N.isLeaf = false; N.children = {0, 1}; t.push_back(N);
    TreeNode P; P.id = 7; P.name = "P"; P.filePath = "/d/A.fa"; P.level = 0; P.parentId = -1; P.isLeaf = false; P.children = {6, 2}; t.push_back(P);
    TreeNode O; O.id = 8; O.name = "O"; O.filePath = "/d/D.fa"; O.level = 1; O.parentId = -1; O.isLeaf = false; O.children = {3, 4}; t.push_back(O);

    auto part = treeLevelPartition(t);
    assert(part.size() == 3);
    for (const auto& [lvl, entries] : part) {
        std::set<std::string> covered;
        size_t n = 0;
        for (const auto& e : entries) for (const auto& p : e.leafPaths) { covered.insert(p); ++n; }
        assert(n == 6 && covered.size() == 6);          // every leaf exactly once
    }
    auto find = [&](unsigned lvl, const std::string& rep) -> const LevelEntry* {
        for (const auto& e : part[lvl].second) if (e.repPath == rep) return &e;
        return nullptr;
    };
    // level 0: P{A,B,C}, O{D,E} (orphan reaches the top), F
    assert(part[0].second.size() == 3);
    assert(find(0, "/d/A.fa") && find(0, "/d/A.fa")->name == "P" && find(0, "/d/A.fa")->leafPaths.size() == 3);
    assert(find(0, "/d/D.fa") && find(0, "/d/D.fa")->name == "O" && find(0, "/d/D.fa")->isNode);
    assert(find(0, "/d/F.fa") && !find(0, "/d/F.fa")->isNode);
    // level 1: N{A,B} (spans the gap below P), C itself, O{D,E}, F
    assert(part[1].second.size() == 4);
    assert(find(1, "/d/A.fa") && find(1, "/d/A.fa")->name == "N" && find(1, "/d/A.fa")->leafPaths.size() == 2);
    assert(find(1, "/d/C.fa") && !find(1, "/d/C.fa")->isNode);
    assert(find(1, "/d/D.fa") && find(1, "/d/D.fa")->name == "O");
    // level 2: N{A,B}, C, D, E, F
    assert(part[2].second.size() == 5);
    assert(find(2, "/d/A.fa") && find(2, "/d/A.fa")->name == "N");
    assert(find(2, "/d/D.fa") && !find(2, "/d/D.fa")->isNode);
    // levels 1 and 2 differ; a level identical to the next is detectable
    assert(!samePartition(part[1].second, part[2].second));
    assert(samePartition(part[2].second, part[2].second));
    std::cerr << "PASSED: Test treeLevelPartition\n";
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
    testBuildParamsRoundTrip();
    testGetTreeLevelRefs();
    testTreeLevelPartition();
    testBuildHierarchicalTree();
    std::cerr << "TreeSketch: All tests PASSED!\n\n";
    return 0;
}
