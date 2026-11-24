// The suites check with assert(); with NDEBUG defined (CMake's RelWithDebInfo and
// MinSizeRel, a packager's CPPFLAGS) every check compiled away, calls made inside
// an assert were not made at all, and the suite still printed PASSED and exited 0.
#undef NDEBUG
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <set>
#include <unordered_map>
#include <vector>
#include <cassert>

#include "SequenceUtil.hpp"
#include "TestUtil.hpp"
#include "QuerySketchUtil.hpp"
#include "GroupSketchUtil.hpp"

void testQuerySampleBatch(QuerySketchConfig cfg, const std::vector<std::string>& sampleFiles) {
    std::cerr << "START: Test querySample... \t";

    /*
    Sketch has 4 unique kmers from 2 references (ids 0, 1):
      ref 0 (expected_ref_test_1): AGGCTACAAA, AGACTGAACT
      ref 1 (expected_ref_test_2): CTCAGCTAAG, AGGCTACAAC

    r.fq has 3 reads of 11bp: AGACTGAACTC, AGACTGAACTG, AGACTGAACTT
    Querying 10bp kmers: AGACTGAACT hits ref 0 three times.
    Expected sketchCount: {3, 0}
    */
    SketchHash sketchHash{
        {"CTCAGCTAAG", 1}, {"AGGCTACAAC", 1},
        {"AGGCTACAAA", 0}, {"AGACTGAACT", 0}
    };
    std::vector<SketchHash> refSigCount{
        {{"AGGCTACAAA", 0}, {"AGACTGAACT", 0}},
        {{"CTCAGCTAAG", 0}, {"AGGCTACAAC", 0}}
    };
    std::vector<std::string> sketchRef{"expected_ref_test_1", "expected_ref_test_2"};
    std::vector<unsigned> sketchCount(sketchRef.size());
    std::vector<std::unordered_set<std::string>> refRead(sketchRef.size());

    loadSketch(cfg, cfg.ref, sketchHash, sketchRef, refSigCount);
    querySampleBatch(cfg, sampleFiles, sketchHash, sketchCount, refSigCount, refRead);

    std::vector<unsigned> expected{3, 0};
    assert(std::equal(sketchCount.begin(), sketchCount.end(), expected.begin()) == true);

    std::cerr << "PASSED: Test querySample\n";
}

void testGenerateQueryResult(const QuerySketchConfig& cfg) {
    std::cerr << "START: Test generateQueryResult... \t";

    /*
    After querying r.fq:
      expected_ref_test_1: count=3 (2 of its 2 signatures hit, 1 supporting read)
      expected_ref_test_2: count=0 -> not reported
    Output row: expected_ref_test_1  abundance=1  count=3  reads=1
                sigs_hit=2  sigs_total=2  resolution=strain
    */
    std::vector<unsigned> sketchCount = {3, 0};
    std::vector<std::string> sketchRef = {"expected_ref_test_1", "expected_ref_test_2"};
    std::vector<SketchHash> refSigCount{
        {{"AGGCTACAAA", 2}, {"AGACTGAACT", 1}},
        {{"CTCAGCTAAG", 0}, {"AGGCTACAAC", 0}}
    };
    // One entry per reference (the second has no supporting reads). Previously
    // sized 1 for 2 references -- an out-of-bounds access that happened not to
    // crash.
    std::vector<std::unordered_set<std::string>> refRead{{"1"}, {}};

    generateQueryResult(cfg, sketchCount, sketchRef, refSigCount, refRead, cfg.out);
    assert(compareFiles(cfg.out, "expected_out_querysketch.tsv") == true);

    std::cerr << "PASSED: Test generateQueryResult\n";
}

// Parse the primary output into ref -> abundance.
static std::unordered_map<std::string, double> readAbundances(const std::string& path) {
    std::unordered_map<std::string, double> ab;
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);  // header
    while (std::getline(in, line)) {
        std::istringstream ls(line);
        std::string ref;
        double abundance;
        ls >> ref >> abundance;
        ab[ref] = abundance;
    }
    return ab;
}

void testAbundanceIsHitsPerSignature(QuerySketchConfig cfg) {
    std::cerr << "START: Test abundance is hits per signature... \t";

    /*
    Two references at equal coverage, one carrying 4 signatures and the other 1
    (its whole candidate pool was smaller than -c). Equal coverage means equal hits
    per signature: 2 each. Raw counts (8 vs 2) would report an 80/20 split; hits
    per signature must report 50/50.
    */
    cfg.out = "out_querysketch_norm.tsv";
    std::vector<unsigned> sketchCount = {8, 2};
    std::vector<std::string> sketchRef = {"big_pool", "small_pool"};
    std::vector<SketchHash> refSigCount{
        {{"A", 2}, {"B", 2}, {"C", 2}, {"D", 2}},
        {{"E", 2}}
    };
    std::vector<std::unordered_set<std::string>> refRead{{"r1", "r2", "r3"}, {"r4", "r5", "r6"}};

    generateQueryResult(cfg, sketchCount, sketchRef, refSigCount, refRead, cfg.out);
    auto ab = readAbundances(cfg.out);
    assert(ab.size() == 2);
    assert(std::abs(ab["big_pool"] - 0.5) < 1e-9);
    assert(std::abs(ab["small_pool"] - 0.5) < 1e-9);

    // The tree writer applies the same normalisation, families included: a family
    // row with 10 signatures and 20 hits sits at the same coverage as the leaves.
    cfg.out = "out_querysketch_norm_tree.tsv";
    FamilyCall fam;
    fam.name = "family_x"; fam.level = 1; fam.count = 20; fam.reads = 20;
    fam.sigsHit = 10; fam.sigsTotal = 10; fam.coveredLeaves = {"absent_leaf"};
    std::vector<FamilyCall> families{fam};
    generateTreeQueryResult(cfg, sketchCount, sketchRef, refSigCount, refRead, families, cfg.out);
    ab = readAbundances(cfg.out);
    assert(ab.size() == 3);
    // The file carries six significant digits, hence the tolerance.
    assert(std::abs(ab["big_pool"] - 1.0 / 3) < 1e-5);
    assert(std::abs(ab["small_pool"] - 1.0 / 3) < 1e-5);
    assert(std::abs(ab["family_x"] - 1.0 / 3) < 1e-5);

    std::cerr << "PASSED: Test abundance is hits per signature\n";
}

void testFamilyBreadth(QuerySketchConfig cfg) {
    std::cerr << "START: Test family breadth floor... \t";

    /*
    A node of 2 000 signatures hit 2 000 times on 20 distinct ones: one locus. Against
    the node's total the ratio is 20 / (2000 (1 - e^-1)) = 0.016. Its members carry
    1 800 and 1 500 of the signatures; against the fewest, 1 500, the ratio is
    20 / (1500 (1 - e^-1.33)) = 0.018. A present member at 2 000 hits on 1 000 distinct
    signatures of the 1 500 it carries scores 1 000 / 1 105 = 0.90, where against the
    node's total it would be 0.79. The floor withholds the first and keeps the second;
    with familyBreadth off (an index built with the largest member) both are reported;
    --tree-breadth=0 reports both as well.
    */
    FamilyCall locus;
    locus.name = "family_locus"; locus.level = 1; locus.count = 2000; locus.reads = 300;
    locus.sigsHit = 20; locus.sigsTotal = 2000; locus.coveredLeaves = {"m1", "m2"};
    locus.memberCarried = {{"m1", 1800}, {"m2", 1500}};
    FamilyCall present = locus;
    present.name = "family_present"; present.sigsHit = 1000; present.coveredLeaves = {"m3", "m4"};
    present.memberCarried = {{"m3", 1800}, {"m4", 1500}};
    FamilyBreadth fl = familyBreadthRatio(locus), fp = familyBreadthRatio(present);
    assert(fl.carried == 1500 && fl.ratio > 0.017 && fl.ratio < 0.019);
    assert(fp.carried == 1500 && fp.ratio > 0.89 && fp.ratio < 0.91);
    assert(breadthRatio(present.count, present.sigsHit, present.sigsTotal) < 0.80);
    // A member's count is capped at the signatures in use; with no member carrying
    // anything the node's own total is the reference.
    FamilyCall bare = locus;
    bare.memberCarried = {{"m1", 0}, {"m2", 5000}};
    assert(familyBreadthRatio(bare).carried == 2000);
    bare.memberCarried.clear();
    assert(familyBreadthRatio(bare).carried == 2000);

    std::vector<unsigned> sketchCount = {0};
    std::vector<std::string> sketchRef = {"unrelated"};
    std::vector<SketchHash> refSigCount{{{"A", 0}}};
    std::vector<std::unordered_set<std::string>> refRead{{}};
    std::vector<FamilyCall> families{locus, present};
    auto names = [&](const std::string& path) {
        std::set<std::string> n;
        std::ifstream in(path);
        std::string line;
        std::getline(in, line);
        while (std::getline(in, line)) n.insert(line.substr(0, line.find('\t')));
        return n;
    };
    cfg.out = "out_querysketch_fam_breadth.tsv";
    cfg.familyBreadth = true;
    generateTreeQueryResult(cfg, sketchCount, sketchRef, refSigCount, refRead, families, cfg.out);
    assert(names(cfg.out) == std::set<std::string>{"family_present"});
    cfg.familyBreadth = false;
    generateTreeQueryResult(cfg, sketchCount, sketchRef, refSigCount, refRead, families, cfg.out);
    assert(names(cfg.out) == std::set<std::string>({"family_locus", "family_present"}));
    cfg.familyBreadth = true;
    cfg.treeBreadth = 0.0;
    generateTreeQueryResult(cfg, sketchCount, sketchRef, refSigCount, refRead, families, cfg.out);
    assert(names(cfg.out) == std::set<std::string>({"family_locus", "family_present"}));
    std::remove(cfg.out.c_str());
    std::remove("log_out_querysketch_fam_breadth.tsv");
    std::remove("logread_out_querysketch_fam_breadth.tsv");

    std::cerr << "PASSED: Test family breadth floor\n";
}

static std::string randDNA(size_t n, uint64_t seed) {
    std::string s(n, 'A');
    uint64_t x = seed * 6364136223846793005ull + 1442695040888963407ull;
    for (size_t i = 0; i < n; i++) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        s[i] = "ACGT"[x & 3];
    }
    return s;
}

void testChunkedLoad(QuerySketchConfig cfg) {
    std::cerr << "START: Test chunked sketch load... \t";

    /*
    A sketch of 7 references and 3 000 lines, written one reference after another
    as uniqsketch writes it, with a carriers column on every third line and two
    sequences that occur twice (the second line under another reference, as a
    carried k-mer that filled two sibling nodes does), loaded in one pass and in
    seven byte ranges on four threads: the maps, names, carriers and positions
    must come out the same, the duplicates resolved to their first line. Then the
    same for the within-group sketch loader.
    */
    const unsigned k = 21;
    std::vector<std::string> seqs;
    for (unsigned i = 0; i < 3000; i++) seqs.push_back(randDNA(k, 1000 + i));
    seqs[2500] = seqs[10];                 // duplicates: ref 5 repeats ref 0's k-mer, ref 6 ref 1's
    seqs[2900] = seqs[500];
    {
        std::ofstream out("chunked_sketch.tsv");
        for (unsigned i = 0; i < 3000; i++) {
            const unsigned ref = i / 430;      // refs 0..6
            out << seqs[i] << "\t" << ref << "\tref" << ref << "\t7\tcontig" << (i % 3) << "\t" << (i * 37 % 100000);
            if (i % 3 == 0) out << "\t" << (ref + 1) % 7 << "," << (ref + 2) % 7; else out << "\t-";
            out << "\n";
        }
    }
    auto load = [&](unsigned threads, uint64_t minChunk, SketchHash& hash, std::vector<std::string>& ref,
                    std::vector<SketchHash>& sig, std::unordered_map<std::string, std::vector<unsigned>>& car,
                    SignaturePositions& pos) {
        QuerySketchConfig c = cfg;
        loadSketch(c, "chunked_sketch.tsv", hash, ref, sig, &car, threads, &pos, minChunk);
        assert(c.kmerLen == k);
    };
    SketchHash h1, h7;
    std::vector<std::string> r1, r7;
    std::vector<SketchHash> s1, s7;
    std::unordered_map<std::string, std::vector<unsigned>> c1, c7;
    SignaturePositions p1, p7;
    load(1, 32u << 20, h1, r1, s1, c1, p1);
    load(4, 1, h7, r7, s7, c7, p7);          // every range at least one byte: as many as threads
    assert(h1.size() == 2998 && h1 == h7);
    assert(h1.at(seqs[10]) == 0 && h1.at(seqs[500]) == 1);   // first line wins
    assert(r1.size() == 7 && r1 == r7 && r1[6] == "ref6");
    assert(s1.size() == 7 && s1 == s7);
    size_t mapped = 0;
    for (const SketchHash& m : s1) mapped += m.size();
    assert(mapped == 3000);                  // each reference's map holds its own lines
    assert(c1.size() == 1000 && c1 == c7 && c1.at(seqs[3]) == std::vector<unsigned>({1, 2}));
    assert(p1.size() == 7 && p7.size() == 7);
    for (size_t r = 0; r < 7; r++) {
        auto a = p1[r], b = p7[r];
        std::sort(a.begin(), a.end()); std::sort(b.begin(), b.end());
        assert(a.size() == (r < 6 ? 430 : 420) && a == b);
    }
    std::remove("chunked_sketch.tsv");

    {
        std::ofstream out("chunked_group.tsv");
        for (unsigned i = 0; i < 2000; i++)
            out << seqs[i] << "\t" << i / 300 << "\tref" << i / 300 << "\t7\tcontig" << (i % 2) << "\t" << i * 11
                << "\t" << (i % 4 == 0 ? "-" : std::to_string(i % 7)) << "\tnode" << i / 600 << "\n";
    }
    std::vector<GroupSignature> g1 = loadGroupSketch("chunked_group.tsv"), g5 = loadGroupSketch("chunked_group.tsv", 5, 1);
    assert(g1.size() == 2000 && g5.size() == 2000);
    for (size_t i = 0; i < g1.size(); i++) {
        assert(g1[i].seq == g5[i].seq && g1[i].refId == g5[i].refId && g1[i].refName == g5[i].refName &&
               g1[i].contig == g5[i].contig && g1[i].pos == g5[i].pos && g1[i].carriers == g5[i].carriers &&
               g1[i].node == g5[i].node);
        assert(g1[i].seq == seqs[i] && g1[i].pos == i * 11);
    }
    std::remove("chunked_group.tsv");

    // The index's sort: distinct (hash, position) keys come out as std::sort orders
    // them, buckets or not.
    std::vector<std::pair<uint64_t, uint32_t>> keys(400000);
    uint64_t x = 88172645463325252ull;
    for (uint32_t i = 0; i < keys.size(); i++) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        keys[i] = {x >> (i % 5 == 0 ? 40 : 0), i};   // a few hashes repeat (small values), positions break ties
    }
    std::vector<std::pair<uint64_t, uint32_t>> expected = keys;
    std::sort(expected.begin(), expected.end());
    parallelSort(keys);
    assert(keys == expected);

    std::cerr << "PASSED: Test chunked sketch load\n";
}

int main() {
    QuerySketchConfig cfg;
    cfg.ref = "expected_out_uniqsketch.txt";
    cfg.bits = 64;
    cfg.dbfSize = 32;
    cfg.sbfSize = 16;
    cfg.nhash = 3;
    cfg.kmerLen = 10;
    cfg.sketchHit = 1;
    cfg.readCutoff = 0;
    cfg.out = "out_querysketch.tsv";

    std::vector<std::string> sampleFiles = {"r.fq"};

    // testQuerySampleBatch takes cfg by value: loadSketch sets kmerLen from the
    // index, and that must not leak into the second test.
    testQuerySampleBatch(cfg, sampleFiles);
    testGenerateQueryResult(cfg);
    testAbundanceIsHitsPerSignature(cfg);
    testFamilyBreadth(cfg);
    testChunkedLoad(cfg);

    std::cerr << "QuerySketch: All tests PASSED!\n\n";
    return 0;
}