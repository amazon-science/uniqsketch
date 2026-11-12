#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <vector>

#include "TreeSketchUtil.hpp"
#include "GroupSketchUtil.hpp"
#include "NodeCarrierUtil.hpp"

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
    p.maxHomopolymer = 6; p.groupSignatures = true; p.thresholds = "10000,1000,100";
    p.bloomFilter = bloomFilterName();
    writeTreeBuildParams(path, p);
    TreeBuildParams q = loadTreeBuildParams(path);
    assert(q.loaded);
    assert(q.kmerLen == 81 && q.sketchnum == 500);
    assert(q.minMargin == 2 && q.strictMargin);
    assert(q.maxHomopolymer == 6 && q.thresholds == "10000,1000,100");
    assert(q.groupSignatures);
    assert(q.bloomFilter == bloomFilterName());               // the layout that built it, named
    assert(q.bloomFilter == "blocked" || q.bloomFilter == "standard");
    std::remove(path.c_str());

    TreeBuildParams missing = loadTreeBuildParams(path);
    assert(!missing.loaded);
    assert(missing.bloomFilter.empty());                      // no layout recorded
    assert(missing.minMargin == 1 && !missing.strictMargin);  // pre-2.0 index: no margin assumed
    assert(!missing.groupSignatures);                         // and no within-group signatures
    // A file from before the field existed loads with the field empty, so the
    // query's mismatch check stays silent on older indexes.
    { std::ofstream old(path); old << "kmer\t81\ncov\t500\nmin_margin\t2\nstrict_margin\t1\n"; }
    TreeBuildParams older = loadTreeBuildParams(path);
    assert(older.loaded && older.bloomFilter.empty());
    std::remove(path.c_str());
    std::cerr << "PASSED: Test build params round trip\n";
}

// leaves.tsv names every leaf by refId, including the ones no sketch mentions.
// The query relies on it to treat a signature-less leaf as the carrier of its
// relatives' within-group signatures; an index without it must load as empty so
// the query can fall back rather than fail.
void testTreeLeavesRoundTrip() {
    std::cerr << "START: Test leaves.tsv round trip... \t";
    const char* dir = ".";
    std::string path = treeLeavesPath(dir);
    std::vector<std::string> files = {"refs/GCF_1_A.fna", "refs/GCF_2_B.fna.gz", "GCF_3_C with space.fa"};
    writeTreeLeaves(path, files);
    std::vector<std::string> names = loadTreeLeaves(path);
    assert(names.size() == 3);
    assert(names[0] == "GCF_1_A" && names[1] == "GCF_2_B" && names[2] == "GCF_3_C with space");
    // The name column is what getBaseId gives the sketches, so the two agree.
    for (size_t i = 0; i < files.size(); i++) assert(names[i] == getBaseId(files[i]));
    std::remove(path.c_str());
    assert(loadTreeLeaves(path).empty());                       // pre-leaves.tsv index
    // A damaged line is skipped, not fatal; ids need not be dense.
    { std::ofstream bad(path); bad << "ref_id\tname\tfile_path\n0\tX\tx.fa\nnot-a-number\tY\ty.fa\n2\tZ\tz.fa\n"; }
    names = loadTreeLeaves(path);
    assert(names.size() == 3 && names[0] == "X" && names[1].empty() && names[2] == "Z");
    std::remove(path.c_str());
    std::cerr << "PASSED: Test leaves.tsv round trip\n";
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

// Within-group signatures: the groups come from direct parents, and every
// signature's carriers are exactly the other genomes that contain it.
void testGroupSignatures() {
    std::cerr << "START: Test within-group signature groups and carriers... \t";
    const unsigned k = 21;
    // A and B differ at two sites 300 bp apart; C is unrelated; D is C plus a copy
    // of A's neighbourhood around the first site, so A's within-group k-mers
    // spanning site 1 are carried by D and those spanning site 2 are not.
    std::string a = randDNA(1200, 4242);
    std::string b = a;
    b[400] = (a[400] == 'A') ? 'C' : 'A';
    b[700] = (a[700] == 'G') ? 'T' : 'G';
    std::string c = randDNA(1200, 9797);
    std::string d = c + a.substr(340, 120);
    writeFasta("grp_A.fa", a); writeFasta("grp_B.fa", b);
    writeFasta("grp_C.fa", c); writeFasta("grp_D.fa", d);
    std::vector<std::string> refs = {"grp_A.fa", "grp_B.fa", "grp_C.fa", "grp_D.fa"};

    // Tree: node N (level 0) = {A, B}; C and D never merged.
    std::vector<TreeNode> t;
    auto leaf = [&](unsigned id, const char* n, const std::string& p, int parent) {
        TreeNode x; x.id = id; x.name = n; x.filePath = p; x.level = 1; x.parentId = parent; x.isLeaf = true; t.push_back(x);
    };
    leaf(0, "A", "grp_A.fa", 4); leaf(1, "B", "grp_B.fa", 4); leaf(2, "C", "grp_C.fa", -1); leaf(3, "D", "grp_D.fa", -1);
    TreeNode N; N.id = 4; N.name = "N"; N.filePath = "grp_A.fa"; N.level = 0; N.parentId = -1; N.isLeaf = false; N.children = {0, 1}; t.push_back(N);

    auto groups = signatureGroups(t, refs);
    assert(groups.size() == 1);
    assert(groups[0].nodeName == "N");
    assert(groups[0].members == std::vector<unsigned>({0, 1}));
    assert(groups[0].directLeaves.size() == 2);

    UniqSketchConfig cfg;
    cfg.kmerLen = k;
    std::vector<GroupSignature> sigs;
    auto add = [&](size_t pos) {   // A's k-mer starting at pos, as A's signature
        GroupSignature s; s.seq = get_canonical(a.substr(pos, k)); s.refId = 0; s.refName = "A";
        s.contig = "seq"; s.pos = static_cast<unsigned>(pos); s.node = "N"; sigs.push_back(s);
    };
    add(390);   // spans site 1 (400): absent from B, present in D
    add(690);   // spans site 2 (700): absent from B, absent from D
    add(100);   // shared by A and B, nowhere else
    findCarriers(cfg, sigs, refs);
    assert(sigs[0].carriers == std::vector<unsigned>({3}));
    assert(sigs[1].carriers.empty());
    assert(sigs[2].carriers == std::vector<unsigned>({1}));

    // Round trip through the index file.
    writeGroupSketch("grp_sketch.tsv", sigs, 4);
    auto back = loadGroupSketch("grp_sketch.tsv");
    assert(back.size() == 3);
    assert(back[0].carriers == std::vector<unsigned>({3}) && back[0].refId == 0 && back[0].node == "N");
    assert(back[1].carriers.empty() && back[1].seq == sigs[1].seq);
    assert(back[2].carriers == std::vector<unsigned>({1}) && back[2].pos == 100);
    assert(loadGroupSketch("grp_sketch_missing.tsv").empty());

    for (const char* f : {"grp_A.fa", "grp_B.fa", "grp_C.fa", "grp_D.fa", "grp_sketch.tsv"}) std::remove(f);
    std::cerr << "PASSED: Test within-group signature groups and carriers\n";
}

// buildGroupSketch end to end, on the same fixture: A and B differ at two sites,
// D carries A's neighbourhood around the first. A's k-mers spanning site 1 are
// therefore unique within the group {A, B} but not globally, and are the only
// within-group signatures the build should add: A's site-2 k-mers and all of
// B's are globally unique and already in the leaf sketch. The build runs its
// candidate and selection phases as parallel job lists; the output must not
// depend on that, so two builds must be byte-identical.
void testBuildGroupSketch() {
    std::cerr << "START: Test buildGroupSketch... \t";
    const unsigned k = 21;
    std::string a = randDNA(1200, 4242);
    std::string b = a;
    b[400] = (a[400] == 'A') ? 'C' : 'A';
    b[700] = (a[700] == 'G') ? 'T' : 'G';
    std::string c = randDNA(1200, 9797);
    std::string d = c + a.substr(340, 120);
    writeFasta("gb_A.fa", a); writeFasta("gb_B.fa", b);
    writeFasta("gb_C.fa", c); writeFasta("gb_D.fa", d);
    std::vector<std::string> refs = {"gb_A.fa", "gb_B.fa", "gb_C.fa", "gb_D.fa"};

    std::vector<TreeNode> t;
    auto leaf = [&](unsigned id, const char* n, const std::string& p, int parent) {
        TreeNode x; x.id = id; x.name = n; x.filePath = p; x.level = 1; x.parentId = parent; x.isLeaf = true; t.push_back(x);
    };
    leaf(0, "A", "gb_A.fa", 4); leaf(1, "B", "gb_B.fa", 4); leaf(2, "C", "gb_C.fa", -1); leaf(3, "D", "gb_D.fa", -1);
    TreeNode N; N.id = 4; N.name = "N"; N.filePath = "gb_A.fa"; N.level = 0; N.parentId = -1; N.isLeaf = false; N.children = {0, 1}; t.push_back(N);

    // The leaf sketch, as uniqsketch's main builds it, then the group sketch.
    auto build = [&](const std::string& dir) {
        mkdir(dir.c_str(), 0777);
        UniqSketchConfig cfg;
        cfg.kmerLen = k;
        cfg.threads = 4;
        cfg.outdir = dir;
        cfg.outfile = dir + "/sketch_uniq.tsv";
        cfg.refstat = dir + "/db_uniq_count.tsv";
        cfg.numRef = refs.size();
        size_t total = 0;
        getCardinality(total, cfg.dbfSize, cfg.sbfSize, cfg.kmerLen, cfg.threads, refs);
        assert(resolveBloomSizes(total, cfg.dbfSize, cfg.sbfSize));
        BloomFilter dbf(cfg.m1(), cfg.nhash1, cfg.kmerLen);
        identifyUniqKmers(cfg, refs, dbf);
        buildSketch(cfg, cfg.refstat, dbf);
        return buildGroupSketch(cfg, t, refs, dbf);
    };
    size_t n1 = build("gb_idx1");
    size_t n2 = build("gb_idx2");
    assert(n1 == n2 && n1 > 0);

    std::vector<GroupSignature> sigs = loadGroupSketch(groupSketchPath("gb_idx1"));
    assert(sigs.size() == n1);
    std::set<std::string> leafSigs;
    {
        std::ifstream in("gb_idx1/sketch_uniq.tsv");
        std::string seq, rest;
        while (in >> seq && std::getline(in, rest)) leafSigs.insert(seq);
    }
    for (const GroupSignature& s : sigs) {
        assert(s.refId == 0 && s.refName == "gb_A" && s.node == "N");
        assert(s.pos + k > 400 && s.pos <= 400);                    // spans site 1
        assert(s.seq == get_canonical(a.substr(s.pos, k)));         // and is A's k-mer there
        assert(s.carriers == std::vector<unsigned>({3}));           // carried by D only
        assert(leafSigs.count(s.seq) == 0);                         // not repeated from the leaf sketch
    }

    // Byte-identical across builds: the parallel job lists must not leak scheduling.
    auto contents = [](const std::string& p) { std::ifstream in(p); return std::string((std::istreambuf_iterator<char>(in)), {}); };
    assert(contents(groupSketchPath("gb_idx1")) == contents(groupSketchPath("gb_idx2")));
    assert(contents("gb_idx1/sketch_uniq.tsv") == contents("gb_idx2/sketch_uniq.tsv"));

    for (const char* f : {"gb_A.fa", "gb_B.fa", "gb_C.fa", "gb_D.fa"}) std::remove(f);
    system("rm -rf gb_idx1 gb_idx2");
    std::cerr << "PASSED: Test buildGroupSketch (" << n1 << " within-group signatures)\n";
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

// Node signatures with carriers (NodeCarrierUtil.hpp), on a node whose core is
// carried outside it. A is the representative, with a 60-bp repeat of [500, 560)
// appended; the members are A with one substitution (A1), A without [100, 250)
// (A2), and A without [100, 250), [300, 330) and [350, 380) (A3, A4). B is
// unrelated plus A's [100, 400), so A's k-mers starting in [100, 380) are carried
// outside the node and are never plain candidates; C is unrelated. Of the carried
// k-mers, the carried pool holds those at least half the members (3 of 5) carry:
// [250, 280) (all five) and [280, 330) (A, A1, A2), and not [100, 250) or
// [330, 380) (A and A1 only). The repeat's k-mers occur twice in A and are
// excluded whatever their carriage; the junction k-mers it creates are plain.
void testNodeCarriers() {
    std::cerr << "START: Test node signatures with carriers... \t";
    const unsigned k = 21;
    std::string base = randDNA(1200, 2026);
    // A deletion's junction must not recreate the k-mer that overlapped the deleted
    // region by one base (a[349] a[380..400) is a[379..400) when a[349] == a[379]),
    // so the bases on either side of each deletion are made to differ.
    for (auto [i, j] : {std::pair<int, int>{99, 249}, {100, 250}, {299, 329}, {300, 330}, {349, 379}, {350, 380}})
        if (base[i] == base[j]) base[j] = (base[j] == 'A') ? 'C' : 'A';
    std::string a = base + base.substr(500, 60);
    std::string a1 = a; a1[600] = (a[600] == 'A') ? 'C' : 'A';
    std::string a2 = a.substr(0, 100) + a.substr(250, 100) + a.substr(380);
    std::string a3 = a.substr(0, 100) + a.substr(250, 50) + a.substr(330, 20) + a.substr(380);
    std::string b = randDNA(1200, 7171) + a.substr(100, 300);
    std::string c = randDNA(1200, 3131);
    writeFasta("nc_A.fa", a); writeFasta("nc_A1.fa", a1); writeFasta("nc_A2.fa", a2);
    writeFasta("nc_A3.fa", a3); writeFasta("nc_A4.fa", a3);
    writeFasta("nc_B.fa", b); writeFasta("nc_C.fa", c);
    std::vector<std::string> refs = {"nc_A.fa", "nc_A1.fa", "nc_A2.fa", "nc_A3.fa", "nc_A4.fa", "nc_B.fa", "nc_C.fa"};
    std::vector<std::vector<std::string>> groups = {{refs[0], refs[1], refs[2], refs[3], refs[4]}, {refs[5]}, {refs[6]}};

    UniqSketchConfig cfg;
    cfg.kmerLen = k;
    cfg.threads = 4;
    cfg.minMargin = 2;
    cfg.strictMargin = true;
    cfg.sketchnum = 5000;
    size_t total = 0;
    getCardinality(total, cfg.dbfSize, cfg.sbfSize, cfg.kmerLen, cfg.threads, refs);
    assert(resolveBloomSizes(total, cfg.dbfSize, cfg.sbfSize));

    // Carriage: the number of genomes containing a k-mer, exact here (no counter
    // is shared at this load), and once per genome however often it occurs there.
    assert(defaultNodeCarriers(8) == 0 && defaultNodeCarriers(120) == 12 && defaultNodeCarriers(5493) == 50);
    CarriageCounter cc(cfg.dbfSize);
    countCarriage(cfg, refs, cc);
    auto carriage = [&](const std::string& kmer) {
        ntHashIterator itr(kmer, CarriageCounter::HASHES, k);
        return cc.count(*itr);
    };
    assert(carriage(a.substr(260, k)) == 6);    // all five members and B
    assert(carriage(a.substr(300, k)) == 4);    // A, A1, A2 and B
    assert(carriage(a.substr(150, k)) == 3);    // A, A1 and B
    assert(carriage(a.substr(700, k)) == 5);    // the node only
    assert(carriage(a.substr(590, k)) == 4);    // spans A1's substitution
    assert(carriage(a.substr(510, k)) == 5);    // the repeat: once per genome
    assert(carriage(c.substr(10, k)) == 1);
    assert(carriage(randDNA(k, 4444)) == 0);

    // Candidates. A has 1240 k-mers: 80 positions hold the 40 repeated k-mers,
    // 280 (starting in [100, 380)) are carried by B. Plain: 1240 - 80 - 280 = 880,
    // in candidate files and a stat file byte for byte those of the plain scan,
    // whatever the cap. The carried pool, a second file, holds with a cap of 1
    // the 80 starting in [250, 330); with a cap of 0 there is no second file.
    auto level = [&](const std::string& dir) {
        UniqSketchConfig lvl = cfg;
        lvl.outdir = dir;
        lvl.outfile = dir + "/sketch_level_0.tsv";
        lvl.refstat = dir + "/db_uniq_count.tsv";
        lvl.numRef = groups.size();
        mkdir(dir.c_str(), 0777);
        return lvl;
    };
    auto contents = [](const std::string& p) { std::ifstream in(p); return std::string((std::istreambuf_iterator<char>(in)), {}); };
    auto candidates = [](const std::string& p) { CandidateSet c; return readCandidates(p, c) ? c.total : 0; };
    auto counts = [](const std::string& statFile) {   // the stat file's counts, in its order
        std::vector<std::pair<int, int>> n;
        std::ifstream in(statFile);
        std::string path; int cnt, id;
        while (in >> path >> cnt >> id) n.push_back({cnt, id});
        return n;
    };
    UniqSketchConfig plain = level("nc_plain"), cap0 = level("nc_cap0"), cap1 = level("nc_cap1");
    BloomFilter dbfPlain(plain.m1(), plain.nhash1, plain.kmerLen);
    identifyGroupUniqKmers(plain, groups, dbfPlain);
    {
        BloomFilter dbf(cap0.m1(), cap0.nhash1, cap0.kmerLen);
        identifyGroupUniqKmersCarried(cap0, groups, dbf, cc, 0);
    }
    BloomFilter dbf1(cap1.m1(), cap1.nhash1, cap1.kmerLen);
    identifyGroupUniqKmersCarried(cap1, groups, dbf1, cc, 1);
    assert(candidates(candidatePath(plain.outdir, refs[0])) == 880);
    for (const std::string& r : {refs[0], refs[5], refs[6]}) {
        assert(contents(candidatePath(plain.outdir, r)) == contents(candidatePath(cap0.outdir, r)));
        assert(contents(candidatePath(plain.outdir, r)) == contents(candidatePath(cap1.outdir, r)));
        assert(!fileExists(carriedCandidatePath(candidatePath(cap0.outdir, r))));
    }
    assert(counts(plain.refstat).size() == 3 && counts(cap0.refstat) == counts(plain.refstat) &&
           counts(cap1.refstat) == counts(plain.refstat));
    assert(candidates(carriedCandidatePath(candidatePath(cap1.outdir, refs[0]))) == 80);
    assert(!fileExists(carriedCandidatePath(candidatePath(cap1.outdir, refs[5]))));   // singles: no pool

    // Selection: the node-unique signatures are those of the plain build, byte for
    // byte; the carried pool only fills what they leave of the quota. Each
    // signature's position is in the sketch, so the sets are compared by it.
    auto positions = [&](const std::string& sketch, const std::string& ref, bool carriedOnly) {
        std::set<unsigned> pos;
        std::ifstream in(sketch);
        std::string line;
        while (std::getline(in, line)) {
            std::vector<std::string> f;
            size_t p = 0, q;
            while ((q = line.find('\t', p)) != std::string::npos) { f.push_back(line.substr(p, q - p)); p = q + 1; }
            f.push_back(line.substr(p));
            assert(f.size() >= 6);
            if (f[2] != ref) continue;
            if (carriedOnly && (f.size() < 7 || f[6] == "-")) continue;
            pos.insert(static_cast<unsigned>(std::stoul(f[5])));
        }
        return pos;
    };
    const std::string A = getBaseId(refs[0]);
    for (int quota : {5000, 900, 700, 600}) {
        plain.sketchnum = cap1.sketchnum = quota;
        buildSketch(plain, plain.refstat, dbfPlain);
        buildSketch(cap1, cap1.refstat, dbf1);
        std::set<unsigned> pPlain = positions(plain.outfile, A, false), pCap1 = positions(cap1.outfile, A, false);
        assert(pPlain.size() <= static_cast<size_t>(quota));
        for (unsigned p : pPlain) assert(pCap1.count(p) > 0);                     // the plain selection is kept
        size_t fill = 0;
        for (unsigned p : pCap1) if (!pPlain.count(p)) { ++fill; assert(p >= 250 && p < 330); }
        assert(pCap1.size() == pPlain.size() + fill);
        // Of the 80 carried candidates two fail the strict margin: the k-mers at 280
        // and 329 are one substitution from the junction k-mers A3's deletions
        // create (a[280..300) a[330] and a[299] a[330..350)). A shortfall of 78 or
        // more takes the whole usable pool; a smaller one is filled to the quota, or
        // one short of it when a selection run holds only an unusable candidate.
        const size_t shortfall = static_cast<size_t>(quota) - pPlain.size();
        if (shortfall == 0) assert(fill == 0 && pCap1 == pPlain);
        else if (shortfall >= 78) assert(fill == 78);
        else assert(fill <= shortfall && fill + 1 >= shortfall);
        // The other references have no carried pool and are unchanged.
        for (const std::string& r : {refs[5], refs[6]})
            assert(positions(plain.outfile, getBaseId(r), false) == positions(cap1.outfile, getBaseId(r), false));
    }

    // The carriers of the selected signatures as the level sketch's seventh column
    // (the sketch from the last build, quota 800, has none; rebuild at 5000): B's
    // leaf id for every filled signature, `-` for every node-unique one.
    plain.sketchnum = cap1.sketchnum = 5000;
    buildSketch(cap1, cap1.refstat, dbf1);
    std::vector<TreeNode> t;
    auto leaf = [&](unsigned id, const std::string& p, int parent) {
        TreeNode x; x.id = id; x.name = getBaseId(p); x.filePath = p; x.level = 1; x.parentId = parent; x.isLeaf = true; t.push_back(x);
    };
    for (unsigned i = 0; i < 5; i++) leaf(i, refs[i], 7);
    leaf(5, refs[5], -1); leaf(6, refs[6], -1);
    TreeNode N; N.id = 7; N.name = "N"; N.filePath = refs[0]; N.level = 0; N.parentId = -1; N.isLeaf = false; N.children = {0, 1, 2, 3, 4}; t.push_back(N);
    UniqSketchConfig top = cfg;
    top.outdir = "nc_cap1";
    top.outfile = "nc_cap1/sketch_uniq.tsv";       // no leaf sketch here: the report counts none
    reportReachability(top, t, refs, {0}, true);
    size_t rows = 0, carried = 0, inWindow = 0;
    {
        std::ifstream in(cap1.outfile);
        std::string line;
        while (std::getline(in, line)) {
            std::vector<std::string> f;
            size_t p = 0, q;
            while ((q = line.find('\t', p)) != std::string::npos) { f.push_back(line.substr(p, q - p)); p = q + 1; }
            f.push_back(line.substr(p));
            assert(f.size() == 7);
            ++rows;
            if (f[2] != A) { assert(f[6] == "-"); continue; }
            const unsigned pos = static_cast<unsigned>(std::stoul(f[5]));
            const bool window = pos >= 250 && pos < 330;
            inWindow += window;
            if (f[6] != "-") { ++carried; assert(window && f[6] == "5"); }
            else assert(!window);
        }
    }
    assert(rows > 0 && carried == 78 && carried == inWindow);
    auto back = loadLevelCarriers(cap1.outfile);
    assert(back.size() == carried);
    for (const auto& kv : back) assert(kv.second == std::vector<unsigned>({5}));
    assert(loadLevelCarriers(plain.outfile).empty());             // no seventh column

    // Signature loci, as the query counts them for the exclusion test: the level
    // sketch of the cap-1 build has A's signatures in position order, node-unique
    // ones over the genome and the filled ones in [250, 330), so with a gap of 500
    // all of A's signatures on the one contig are one locus; with a gap of 0 each
    // signature is its own locus. Within-group positions are merged in: one on
    // another contig adds a locus, one next to an existing signature does not. A
    // reference with no signatures has none.
    {
        std::vector<std::tuple<unsigned, std::string, unsigned>> none;
        std::vector<unsigned> l = countSignatureLoci(cap1.outfile, 3, none, 500);
        assert(l.size() == 3 && l[0] == 1 && l[1] == 1 && l[2] == 1);
        l = countSignatureLoci(cap1.outfile, 3, none, 0);
        std::set<unsigned> pA = positions(cap1.outfile, A, false);
        assert(l[0] == pA.size());
        std::vector<std::tuple<unsigned, std::string, unsigned>> extra = {
            {0, "seq", *pA.begin() + 1},          // beside a signature of A: same locus
            {0, "plasmid", 10},                   // another contig: a locus of its own
            {1, "seq", 5}, {1, "seq", 2000},      // B: 2 000 bases apart with a 500 gap
            {7, "seq", 1}};                       // out of range: ignored
        l = countSignatureLoci(cap1.outfile, 3, extra, 500);
        assert(l[0] == 2 && l[1] == 2 && l[2] == 1);
        l = countSignatureLoci("nc_cap1/no_such_sketch.tsv", 2, none, 500);
        assert(l.size() == 2 && l[0] == 0 && l[1] == 0);
    }

    for (const std::string& r : refs) std::remove(r.c_str());
    system("rm -rf nc_plain nc_cap0 nc_cap1");
    std::cerr << "PASSED: Test node signatures with carriers (" << rows << " level signatures, " << carried << " carried)\n";
}

int main() {
    testLoadTree();
    testBuildParamsRoundTrip();
    testTreeLeavesRoundTrip();
    testGetTreeLevelRefs();
    testTreeLevelPartition();
    testGroupSignatures();
    testBuildGroupSketch();
    testNodeCarriers();
    testBuildHierarchicalTree();
    std::cerr << "TreeSketch: All tests PASSED!\n\n";
    return 0;
}
