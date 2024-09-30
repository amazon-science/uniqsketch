#include <fstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <vector>
#include <cassert>

#include "SequenceUtil.hpp"
#include "TestUtil.hpp"
#include "QuerySketchUtil.hpp"

// check if querySampleBatch generates the correct counts for the unique sketch in r.fq
void testQuerySampleBatch(const std::vector<std::string> &sampleFiles) {
    std::cerr << "START: Test querySample... \t";

    /*
    uniqsketch has four uniqe kmers from 2 references: expected_ref_test_1, expected_ref_test_2 with ref integer ids 0, 1.
    AGGCTACAAA, AGACTGAACT: from ref expected_ref_test_1 with id 0
    CTCAGCTAAG, AGGCTACAAC: from ref expected_ref_test_2 with id 1
    */
    SketchHash sketchHash {{"CTCAGCTAAG",1},{"AGGCTACAAC",1},{"AGGCTACAAA",0},{"AGACTGAACT",0}};
    std::vector<SketchHash> refSigCount {{{"AGGCTACAAA",0},{"AGACTGAACT",0}}, {{"CTCAGCTAAG",0},{"AGGCTACAAC",0}}};
    std::vector<std::string> sketchRef {"expected_ref_test_1", "expected_ref_test_2"};
    /*
    samplefile r.fq has three reads each 11bp. 
    AGACTGAACTC
    AGACTGAACTG
    AGACTGAACTT

    When querying the 10bp kmers in reads against uniqsketch, 
    AGACTGAACT observed three times in expected_ref_test_1 with id 0 and no kmers from reads are in expected_ref_test_2
    therefore we should have the sketchCount vector as the result: {3, 0}
    */
    std::vector<unsigned> sketchCount(sketchRef.size());    
    std::vector<std::unordered_set<std::string> > refRead(sketchRef.size());

    loadSketch(opt::ref, sketchHash, sketchRef, refSigCount);

    querySampleBatch(sampleFiles, sketchHash, sketchCount, refSigCount, refRead);

    std::vector<unsigned> expectedSketchCount {3, 0}; 
    
    assert(std::equal(sketchCount.begin(), sketchCount.end(), expectedSketchCount.begin())==true);

    std::cerr << "PASSED: Test querySample" << std::endl;
}

// check if generateQueryResult generated the correct tsv output as expected
void testGenerateQueryResult() {
    std::cerr << "START: Test generateQueryResult... \t";

    /*
    Querying kmers in r.fq resulted in the reference--count stat below:
    expected_ref_test_1: 3
    expected_ref_test_2: 0

    So the final tsv output of querysketch should be:
    ref	abundance	count
    expected_ref_test_1	1	3
    */
    std::vector<unsigned> sketchCount = {3, 0};
    std::vector<std::string> sketchRef = {"expected_ref_test_1", "expected_ref_test_2"};
    std::vector<SketchHash> refSigCount {{{"AGGCTACAAA",2},{"AGACTGAACT",1}}, {{"CTCAGCTAAG",0},{"AGGCTACAAC",0}}};
    std::vector<std::unordered_set<std::string> > refRead {{"1"}};

    generateQueryResult(sketchCount, sketchRef, refSigCount, refRead, opt::out);

    assert(compareFiles(opt::out, "expected_out_querysketch.tsv")==true);
    
    std::cerr << "PASSED: Test generateQueryResult" << std::endl;
}

int main() {
    opt::ref = "expected_out_uniqsketch.txt";
    opt::bits = 64;
    opt::dbfSize = 32;
    opt::sbfSize = 16;
    opt::nhash = 3;
    opt::kmerLen = 10;
    opt::sketchHit = 1;
    opt::readCutoff = 0;
    opt::out = "out_querysketch.tsv";

    std::vector<std::string> sampleFiles;
    sampleFiles.push_back("r.fq");

    testQuerySampleBatch(sampleFiles);
    testGenerateQueryResult();

    std::cerr << "QuerySketch: All tests PASSED!" << std::endl << std::endl;
    return 0;
}
