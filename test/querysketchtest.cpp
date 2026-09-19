#include <fstream>
#include <iterator>
#include <algorithm>
#include <string>
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

    std::cerr << "QuerySketch: All tests PASSED!\n\n";
    return 0;
}