#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "CompareSketchUtil.hpp"

// Test that loadBloomFilter populates the filter so known k-mers are found
void testLoadBloomFilter(const CompareSketchConfig& cfg) {
    std::cerr << "START: Test loadBloomFilter... \t";

    /*
    ref_test_1.fa has sequences CAGGCTACAAA and GAGACTGAACT (11bp each).
    With k=10, the 4 k-mers are: CAGGCTACAA, AGGCTACAAA, GAGACTGAAC, AGACTGAACT.
    After loading into a BF, querying these k-mers by hash should return true.
    */
    BloomFilter dbFilter(cfg.bits * cfg.dbfSize, cfg.nhash, cfg.kmerLen);
    loadBloomFilter(cfg, "ref_test_1.fa", dbFilter);

    // Verify a known k-mer is present by hashing it
    std::string kmer("CAGGCTACAA");
    assert(dbFilter.contains(kmer.c_str()) == true);

    // Verify a k-mer NOT in ref_test_1 is absent
    std::string absent("TTTTTTTTTT");
    assert(dbFilter.contains(absent.c_str()) == false);

    std::cerr << "PASSED: Test loadBloomFilter\n";
}

// Test checkBloomFilter computes correct similarity between identical references
void testCheckBloomFilterIdentical(const CompareSketchConfig& cfg) {
    std::cerr << "START: Test checkBloomFilter (identical)... \t";

    /*
    Load ref_test_1 into BF, then check ref_test_1 against itself.
    All k-mers should be found -> similarity = 1.0, countUnique = 0.
    */
    BloomFilter dbFilter(cfg.bits * cfg.dbfSize, cfg.nhash, cfg.kmerLen);
    loadBloomFilter(cfg, "ref_test_1.fa", dbFilter);

    size_t countAll = 0, countUnique = 0;
    double similarity = checkBloomFilter(cfg, "ref_test_1.fa", dbFilter, countAll, countUnique);

    assert(countAll == 4);      // 4 k-mers of length 10 from two 11bp sequences
    assert(countUnique == 0);   // all found in BF
    assert(std::abs(similarity - 1.0) < 1e-9);

    std::cerr << "PASSED: Test checkBloomFilter (identical)\n";
}

// Test checkBloomFilter between different references
void testCheckBloomFilterDifferent(const CompareSketchConfig& cfg) {
    std::cerr << "START: Test checkBloomFilter (different)... \t";

    /*
    Load ref_test_1 into BF, check ref_test_2 against it.
    ref_test_1 k-mers: {CAGGCTACAA, AGGCTACAAA, GAGACTGAAC, AGACTGAACT}
    ref_test_2 k-mers: {CTTAGCTGAG, TTAGCTGAGG, CAGGCTACAA, AGGCTACAAC}

    Only CAGGCTACAA is shared (by hash). So 1 out of 4 k-mers found.
    Note: BF uses ntHash (canonical), so we check by hash not exact string.
    countAll = 4, similarity should be ~0.25 (1 shared / 4 total).
    */
    BloomFilter dbFilter(cfg.bits * cfg.dbfSize, cfg.nhash, cfg.kmerLen);
    loadBloomFilter(cfg, "ref_test_1.fa", dbFilter);

    size_t countAll = 0, countUnique = 0;
    double similarity = checkBloomFilter(cfg, "ref_test_2.fa", dbFilter, countAll, countUnique);

    assert(countAll == 4);
    // At least some k-mers should be unique to ref_test_2
    assert(countUnique > 0);
    assert(similarity > 0.0 && similarity < 1.0);

    std::cerr << "PASSED: Test checkBloomFilter (different)\n";
}

int main() {
    CompareSketchConfig cfg;
    cfg.kmerLen = 10;
    cfg.bits = 128;
    cfg.nhash = 3;
    cfg.dbfSize = 100;

    testLoadBloomFilter(cfg);
    testCheckBloomFilterIdentical(cfg);
    testCheckBloomFilterDifferent(cfg);

    std::cerr << "CompareSketch: All tests PASSED!\n\n";
    return 0;
}