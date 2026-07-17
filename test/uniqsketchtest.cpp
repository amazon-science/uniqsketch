#include <fstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <vector>
#include <cassert>

#include "SequenceUtil.hpp"
#include "TestUtil.hpp"
#include "UniqSketchUtil.hpp"

void testIdentifyUniqKmers(const std::vector<std::string>& refFiles) {
    std::cerr << "START: Test identifyUniqKmers... \t";

    /*
    ref_test_1: {CAGGCTACAAA, GAGACTGAACT} -> 4 kmers of 10bp
    ref_test_2: {CTTAGCTGAGG, CAGGCTACAAC} -> 4 kmers of 10bp

    CAGGCTACAA is shared, so canonical unique sets are:
    ref_test_1: {AGGCTACAAA, GAGACTGAAC, AGACTGAACT}
    ref_test_2: {CTCAGCTAAG, CCTCAGCTAA, AGGCTACAAC}
    */
    BloomFilter dbFilter(opt::m1, opt::nhash1, opt::kmerLen);
    identifyUniqKmers(refFiles, dbFilter);

    assert(compareFiles("ref_test_1.tsv", "expected_ref_test_1.tsv") == true);
    assert(compareFiles("ref_test_2.tsv", "expected_ref_test_2.tsv") == true);

    std::cerr << "PASSED: Test identifyUniqKmers\n";
}

void testBuildSketch() {
    std::cerr << "START: Test buildSketch... \t";

    /*
    With opt::sketchnum = 2, each ref gets 2 unique kmers uniformly selected:
    ref_test_1: {AGACTGAACT, AGGCTACAAA}
    ref_test_2: {AGGCTACAAC, CTCAGCTAAG}
    Total sketch: 4 unique kmers
    */
    BloomFilter dbFilter(opt::m1, opt::nhash1, opt::kmerLen);
    buildSketch("expected_db_uniq_count.tsv", dbFilter);
    assert(compareFiles(opt::outfile, "expected_out_uniqsketch.txt") == true);

    std::cerr << "PASSED: Test buildSketch\n";
}

void testLowComplexity() {
    opt::maxEntropy = 24.56;
    std::cerr << "START: Test lowComplexity... \t";

    std::string sig1("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCT");
    std::string sig2("CGTTGCATGGCATGCGGTTCGTCTGCGCGTGCGTGCTGATGCAGTAAAAAAGCTGCATAAAATCGATCGAAAGTACATGCG");

    assert(lowComplexity(sig1) == true);
    assert(lowComplexity(sig2) == false);

    std::cerr << "PASSED: Test lowComplexity\n";
}

void testGetCanonical() {
    std::cerr << "START: Test get_canonical... \t";

    assert(get_canonical("GTTTATTTACC") == "GGTAAATAAAC");
    assert(get_canonical("AGGTTAGCTTTC") == "AGGTTAGCTTTC");

    std::cerr << "PASSED: Test get_canonical\n";
}

void testGetBaseId() {
    std::cerr << "START: Test getBaseId... \t";

    assert(getBaseId("/home/user/workspace/test/ref1.fasta") == "ref1");
    assert(getBaseId("ref2.fa") == "ref2");
    // Multiple dots — only strip last extension
    assert(getBaseId("/path/to/ref.test.fa") == "ref.test");
    // No directory separator
    assert(getBaseId("genome.fna") == "genome");
    // No extension
    assert(getBaseId("refonly") == "refonly");
    assert(getBaseId("/path/to/refonly") == "refonly");

    std::cerr << "PASSED: Test getBaseId\n";
}

int main() {
    opt::kmerLen = 10;
    opt::bits = 128;
    opt::dbfSize = 7;
    opt::sbfSize = 1;
    opt::m1 = opt::bits * opt::dbfSize;
    opt::m2 = opt::bits * opt::sbfSize;
    opt::numRef = 2;
    opt::sketchnum = 2;
    opt::outfile = "test_out_uniqsketch.txt";
    opt::outdir = ".";

    std::vector<std::string> refFiles = {"ref_test_1.fa", "ref_test_2.fa"};

    testIdentifyUniqKmers(refFiles);
    testBuildSketch();
    testLowComplexity();
    testGetCanonical();
    testGetBaseId();

    std::cerr << "UniqSketch: All tests PASSED!\n\n";
    return 0;
}