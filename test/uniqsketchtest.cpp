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

// NTMC64sub must reproduce NTMC64 bit for bit on every single-base substitution.
// isTwoSafe relies on that to test a candidate's 3k neighbours without
// rehashing, so an error here would change which signatures are selected.
// Exercised across k values with different remainders mod 31 and 33 (the
// rotation tables are indexed by those), at every position and for every
// base, including a no-op substitution, lowercase bases and N.
void testSubstitutionHash() {
    std::cerr << "START: Test NTMC64sub matches NTMC64... \t";
    const unsigned m = 5;
    const char* alphabet = "ACGTACGTACGTACGTacgtN";   // mostly upper case
    uint64_t x = 20260923;
    auto next = [&x]() { x = x * 6364136223846793005ULL + 1442695040888963407ULL; return x >> 33; };
    size_t checked = 0;
    for (unsigned k : {21u, 31u, 33u, 62u, 64u, 81u, 100u}) {
        for (int rep = 0; rep < 40; rep++) {
            std::string s(k, 'A');
            for (unsigned i = 0; i < k; i++) s[i] = alphabet[next() % 21];
            uint64_t fh = 0, rh = 0;
            NTC64(s.c_str(), k, fh, rh);
            uint64_t expect[m], got[m];
            for (unsigned i = 0; i < k; i++) {
                for (char b : {'A', 'C', 'G', 'T', 'N', 'a'}) {
                    std::string t = s;
                    t[i] = b;
                    NTMC64(t.c_str(), k, m, expect);
                    NTMC64sub(fh, rh, i, s[i], b, k, m, got);
                    for (unsigned j = 0; j < m; j++) assert(got[j] == expect[j]);
                    ++checked;
                }
            }
        }
    }
    std::cerr << "PASSED: Test NTMC64sub matches NTMC64 (" << checked << " substitutions)\n";
}

// isTwoSafe against a brute-force definition on a small universe: a candidate
// is 2-safe iff none of its Hamming-1 neighbours is in the universe. The
// filter is sized so that false positives cannot occur at this scale.
void testIsTwoSafe() {
    std::cerr << "START: Test isTwoSafe... \t";
    UniqSketchConfig cfg;
    cfg.kmerLen = 31;
    cfg.dbfSize = 100000;
    BloomFilter dbf(cfg.m1(), cfg.nhash1, cfg.kmerLen);

    uint64_t x = 777;
    auto next = [&x]() { x = x * 6364136223846793005ULL + 1442695040888963407ULL; return x >> 33; };
    auto randKmer = [&]() { std::string s(cfg.kmerLen, 'A'); for (auto& c : s) c = "ACGT"[next() % 4]; return s; };

    // Universe: 200 random k-mers plus, for 50 further candidates, exactly one
    // Hamming-1 neighbour each (at a random position, random other base).
    std::vector<std::string> safe, unsafe;
    for (int i = 0; i < 200; i++) { std::string u = randKmer(); dbf.insert(u.c_str()); safe.push_back(u); }
    for (int i = 0; i < 50; i++) {
        std::string c = randKmer(), n = c;
        unsigned pos = next() % cfg.kmerLen;
        char b; do { b = "ACGT"[next() % 4]; } while (b == c[pos]);
        n[pos] = b;
        dbf.insert(n.c_str());
        unsafe.push_back(c);
    }
    // A universe k-mer's own presence is not a neighbour; a random 31-mer has
    // no neighbour in a 250-k-mer universe, so the first set is 2-safe ...
    for (const std::string& s : safe) assert(isTwoSafe(cfg, s, dbf));
    // ... and every candidate with a planted neighbour is not.
    for (const std::string& s : unsafe) assert(!isTwoSafe(cfg, s, dbf));
    // Reverse complement of a planted neighbour is the same k-mer to ntHash.
    for (const std::string& s : unsafe) {
        std::string rc(s.rbegin(), s.rend());
        for (auto& c : rc) c = (c == 'A') ? 'T' : (c == 'T') ? 'A' : (c == 'C') ? 'G' : 'C';
        assert(!isTwoSafe(cfg, rc, dbf));
    }
    std::cerr << "PASSED: Test isTwoSafe\n";
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
    testSubstitutionHash();
    testIsTwoSafe();

    std::cerr << "UniqSketch: All tests PASSED!\n\n";
    return 0;
}