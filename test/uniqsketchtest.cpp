#include <fstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <vector>
#include <cassert>

#include "SequenceUtil.hpp"
#include "TestUtil.hpp"
#include "UniqSketchUtil.hpp"

void testIdentifyUniqKmers(const UniqSketchConfig& cfg,
                           const std::vector<std::string>& refFiles) {
    std::cerr << "START: Test identifyUniqKmers... \t";

    /*
    ref_test_1: {CAGGCTACAAA, GAGACTGAACT} -> 4 kmers of 10bp
    ref_test_2: {CTTAGCTGAGG, CAGGCTACAAC} -> 4 kmers of 10bp

    CAGGCTACAA is shared, so canonical unique sets are:
    ref_test_1: {AGGCTACAAA, GAGACTGAAC, AGACTGAACT}
    ref_test_2: {CTCAGCTAAG, CCTCAGCTAA, AGGCTACAAC}
    */
    BloomFilter dbFilter(cfg.m1(), cfg.nhash1, cfg.kmerLen);
    identifyUniqKmers(cfg, refFiles, dbFilter);

    assert(compareFiles("ref_test_1.tsv", "expected_ref_test_1.tsv") == true);
    assert(compareFiles("ref_test_2.tsv", "expected_ref_test_2.tsv") == true);

    std::cerr << "PASSED: Test identifyUniqKmers\n";
}

void testBuildSketch(const UniqSketchConfig& cfg) {
    std::cerr << "START: Test buildSketch... \t";

    /*
    With sketchnum = 2, each ref gets 2 unique kmers uniformly selected:
    ref_test_1: {AGACTGAACT, AGGCTACAAA}
    ref_test_2: {AGGCTACAAC, CTCAGCTAAG}
    Total sketch: 4 unique kmers
    */
    BloomFilter dbFilter(cfg.m1(), cfg.nhash1, cfg.kmerLen);
    buildSketch(cfg, "expected_db_uniq_count.tsv", dbFilter);
    assert(compareFiles(cfg.outfile, "expected_out_uniqsketch.txt") == true);

    std::cerr << "PASSED: Test buildSketch\n";
}

void testLowComplexity() {
    std::cerr << "START: Test lowComplexity... \t";

    // These are 81-mers, so the entropy ceiling is that of the default k=81
    // configuration. Previously the test had to hard-code the ceiling (24.56)
    // because buildSketch computed it as a side effect on a global; now it
    // follows from the config.
    UniqSketchConfig cfg;   // kmerLen = 81
    assert(std::abs(cfg.maxEntropy() - 24.56) < 0.01);

    std::string sig1("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCT");
    std::string sig2("CGTTGCATGGCATGCGGTTCGTCTGCGCGTGCGTGCTGATGCAGTAAAAAAGCTGCATAAAATCGATCGAAAGTACATGCG");

    assert(lowComplexity(cfg, sig1) == true);
    assert(lowComplexity(cfg, sig2) == false);

    std::cerr << "PASSED: Test lowComplexity\n";
}

void testHasLongHomopolymer() {
    std::cerr << "START: Test hasLongHomopolymer... \t";

    // Filter disabled (0 = off): even a pure homopolymer must pass.
    UniqSketchConfig off;
    off.maxHomopolymer = 0;
    assert(hasLongHomopolymer(off, "AAAAAAAAAAAA") == false);

    // Two configs side by side: no shared state, no restore step needed.
    UniqSketchConfig cfg;
    cfg.maxHomopolymer = 5;

    // Boundary: a run of exactly N is allowed, N+1 is rejected.
    assert(hasLongHomopolymer(cfg, "AAAAA") == false);
    assert(hasLongHomopolymer(cfg, "AAAAAA") == true);

    // Run position must not matter (start, middle, end).
    assert(hasLongHomopolymer(cfg, "GGGGGGCTACTA") == true);
    assert(hasLongHomopolymer(cfg, "CTAGGGGGGCTA") == true);
    assert(hasLongHomopolymer(cfg, "CTACTAGGGGGG") == true);

    // Run length resets on a base change: two runs of 5 are still acceptable.
    assert(hasLongHomopolymer(cfg, "AAAAACAAAAA") == false);

    // Mixed sequence with only short runs.
    assert(hasLongHomopolymer(cfg, "AACCGGTTAACCGGTT") == false);

    // Degenerate inputs.
    assert(hasLongHomopolymer(cfg, "") == false);
    assert(hasLongHomopolymer(cfg, "A") == false);

    std::cerr << "PASSED: Test hasLongHomopolymer\n";
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
    // Shared fixture for the two-reference, k=10 index tests. m1/m2 are now
    // derived (bits * size) rather than set by hand.
    UniqSketchConfig cfg;
    cfg.kmerLen = 10;
    cfg.bits = 128;
    cfg.dbfSize = 7;
    cfg.sbfSize = 1;
    cfg.numRef = 2;
    cfg.sketchnum = 2;
    cfg.outfile = "test_out_uniqsketch.txt";
    cfg.outdir = ".";

    std::vector<std::string> refFiles = {"ref_test_1.fa", "ref_test_2.fa"};

    testIdentifyUniqKmers(cfg, refFiles);
    testBuildSketch(cfg);
    testLowComplexity();
    testHasLongHomopolymer();
    testGetCanonical();
    testGetBaseId();

    std::cerr << "UniqSketch: All tests PASSED!\n\n";
    return 0;
}