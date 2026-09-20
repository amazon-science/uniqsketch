#include <fstream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <cassert>

#include "SequenceUtil.hpp"
#include "TestUtil.hpp"
#include "QuerySketchUtil.hpp"

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

    std::cerr << "QuerySketch: All tests PASSED!\n\n";
    return 0;
}