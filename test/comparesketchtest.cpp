// The suites check with assert(); with NDEBUG defined (CMake's RelWithDebInfo and
// MinSizeRel, a packager's CPPFLAGS) every check compiled away, calls made inside
// an assert were not made at all, and the suite still printed PASSED and exited 0.
#undef NDEBUG
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "CompareSketchUtil.hpp"
#include "BlockedBloomFilter.hpp"   // explicit, so both layouts are tested under -DUNIQSKETCH_STANDARD_BF too

#ifdef _OPENMP
#include <omp.h>
#endif

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

// unusableOptions names a k-mer length, a number of bits or a hash number of 0,
// which main() then refuses with exit status 1, and nothing else: the defaults
// and the smallest values a comparison can run with are accepted.
void testUnusableOptions() {
    std::cerr << "START: Test unusable option values are named... \t";
    CompareSketchConfig cfg;
    assert(unusableOptions(cfg).empty());
    cfg.kmerLen = 1;
    cfg.bits = 1;
    cfg.nhash = 1;
    assert(unusableOptions(cfg).empty());
    for (unsigned which = 0; which < 3; which++) {
        CompareSketchConfig zero;
        (which == 0 ? zero.kmerLen : which == 1 ? zero.bits : zero.nhash) = 0;
        assert(unusableOptions(zero) == "-k, -b and -d must be at least 1");
    }
    std::cerr << "PASSED: Test unusable option values are named\n";
}

// Deterministic sequence generator for the tests below.
static uint64_t rngState = 0x9E3779B97F4A7C15ull;
static uint64_t rnd() {
    rngState ^= rngState << 13; rngState ^= rngState >> 7; rngState ^= rngState << 17;
    return rngState;
}
static std::string dna(size_t n) {
    std::string d(n, 'A');
    for (char& c : d) c = "ACGT"[rnd() & 3];
    return d;
}
// A copy of `d` with one base redrawn every `every` positions.
static std::string mutate(std::string d, size_t every) {
    for (size_t i = every / 2; i < d.size(); i += every) d[i] = "ACGT"[rnd() & 3];
    return d;
}

// The comparison keeps one hash per k-mer position, hVal[0] of ntHashIterator, and
// the filters derive the other values from it with NTE64 (BloomFilter::countAbsent).
// The counts are today's only if those derived values are, bit for bit, the ones
// the iterator hands to insert() and contains(). Checked at every position of
// random sequences with non-ACGT characters and lower case, for k on both sides of
// 32 and 64 and for one to eight hash functions; the positions themselves do not
// depend on the number of hash functions either.
void testDerivedHashes() {
    std::cerr << "START: Test hashes derived from the leading hash equal ntHashIterator's... \t";
    const char* alphabet = "ACGTACGTACGTACGTACGTacgt";
    size_t positions = 0, values = 0;
    for (unsigned k : {1u, 11u, 31u, 32u, 33u, 63u, 64u, 65u, 81u, 127u}) {
        for (int rep = 0; rep < 6; rep++) {
            std::string seq(k + rnd() % 700, 'A');
            for (char& c : seq) c = alphabet[rnd() % 24];
            for (unsigned n = rep % 4; n > 0; n--) seq[rnd() % seq.size()] = "NnX-"[rnd() & 3];
            for (unsigned nhash : {1u, 2u, 3u, 4u, 5u, 8u}) {
                ntHashIterator one(seq, 1, k);
                ntHashIterator all(seq, nhash, k);
                for (; all != all.end(); ++all, ++one) {
                    assert(one != one.end() && one.get_pos() == all.get_pos());
                    const uint64_t h0 = (*one)[0];
                    assert(h0 == (*all)[0]);
                    for (unsigned h = 1; h < nhash; h++) {
                        assert(NTE64(h0, k, h) == (*all)[h]);
                        ++values;
                    }
                    ++positions;
                }
                assert(!(one != one.end()));
            }
        }
    }
    assert(positions > 60000 && values > 150000);
    std::cerr << "PASSED: Test hashes derived from the leading hash equal ntHashIterator's ("
              << positions << " positions, " << values << " derived values)\n";
}

// sortHashes against std::sort: empty and one-element vectors, values that differ
// in one byte only (every pass has to be stable), repeats, the extremes.
void testSortHashes() {
    std::cerr << "START: Test sortHashes... \t";
    std::vector<uint64_t> scratch;
    std::vector<std::vector<uint64_t>> inputs = {{}, {42}, {2, 1}, {0, ~uint64_t(0), 0, 1, ~uint64_t(0)}};
    for (unsigned byte = 0; byte < 8; byte++) {
        std::vector<uint64_t> v;
        for (unsigned i = 0; i < 600; i++) v.push_back((rnd() & 255) << (8 * byte) | (byte ? 0 : 0x5A00));
        inputs.push_back(v);
    }
    for (size_t n : {3u, 255u, 256u, 257u, 5000u, 70000u}) {
        std::vector<uint64_t> v(n);
        for (uint64_t& x : v) x = rnd();
        for (size_t i = 0; i + 7 < n; i += 7) v[i + 3] = v[i];       // repeats
        inputs.push_back(v);
    }
    for (std::vector<uint64_t> v : inputs) {
        std::vector<uint64_t> expect = v;
        std::sort(expect.begin(), expect.end());
        sortHashes(v, scratch);
        assert(v == expect);
    }
    std::cerr << "PASSED: Test sortHashes\n";
}

// Test references for the comparison tests: `a` has a repeat inside one record, a
// run of N, a lower-case stretch and a record shorter than k; `b` is `a`'s first
// record with one base in 300 redrawn, plus sequence of its own; `c` is unrelated
// and gzipped; `d` is a part of `a`; `e` holds no k-mer.
static const std::vector<std::string> cmpRefs = {"cmp_a.fa", "cmp_b.fa", "cmp_c.fa.gz", "cmp_d.fa", "cmp_e.fa"};
static void writeCmpRefs() {
    const std::string core = dna(40000);
    const std::string repeat = dna(600);
    std::string lower = dna(1500);
    for (char& c : lower) c = static_cast<char>(std::tolower(c));
    { std::ofstream o(cmpRefs[0]);
      o << ">a_core\n" << core << "\n>a_repeat\n" << dna(2500) << repeat << dna(1500) << repeat << "\n"
        << ">a_n\n" << dna(1200) << std::string(40, 'N') << lower << "\n>a_short\nACGTACGT\n"; }
    { std::ofstream o(cmpRefs[1]);
      o << ">b_core\n" << mutate(core, 300) << "\n>b_own\n" << dna(8000) << "\n"; }
    { std::ostringstream c;
      c << ">c\n" << dna(20000) << "\n";
      const std::string cs = c.str();
      gzFile g = gzopen(cmpRefs[2].c_str(), "wb");
      assert(g);
      assert(gzwrite(g, cs.data(), static_cast<unsigned>(cs.size())) == static_cast<int>(cs.size()));
      gzclose(g); }
    { std::ofstream o(cmpRefs[3]); o << ">d\n" << core.substr(12000, 9000) << "\n"; }
    { std::ofstream o(cmpRefs[4]); o << ">e\nACGTAC\n"; }
}
static void removeCmpRefs() {
    for (const std::string& f : cmpRefs) std::remove(f.c_str());
}

// countAbsent must return, for either filter layout and any number of hash
// functions, the number of k-mers for which contains() is false: the count the
// comparison reported before. The reference is built as before, from all nhash
// values of every position through insert() and contains(); the thin filters make
// false positives common, and the first setting is checked to have some, since a
// test in which every absent k-mer is reported absent could not tell a wrong bit
// from a right one. Three hash functions take the constant path of either
// filter, the other numbers its general one; the order of the leading hashes
// must not matter to the count.
template <class Filter>
static size_t checkCountAbsent(unsigned bits, size_t elements, unsigned nhash, unsigned k) {
    std::vector<uint64_t> scratch;
    std::vector<std::vector<uint64_t>> all(cmpRefs.size()), leading(cmpRefs.size());
    for (size_t r = 0; r < cmpRefs.size(); r++) {
        all[r] = loadKmerHashes(cmpRefs[r], nhash, k);
        leading[r] = loadLeadingHashes(cmpRefs[r], k);
        // The leading hashes are hVal[0] of every position, in file order.
        std::vector<uint64_t> first;
        for (size_t x = 0; x < all[r].size(); x += nhash) first.push_back(all[r][x]);
        assert(leading[r] == first);
    }
    size_t falsePositives = 0;
    for (size_t i = 0; i < cmpRefs.size(); i++) {
        Filter filter(bits * elements, nhash, k);
        for (size_t x = 0; x < all[i].size(); x += nhash) filter.insert(&all[i][x]);
        const std::unordered_set<uint64_t> held(leading[i].begin(), leading[i].end());
        for (size_t b = 0; b < cmpRefs.size(); b++) {
            size_t expect = 0, absent = 0;
            for (size_t x = 0; x < all[b].size(); x += nhash) {
                if (!filter.contains(&all[b][x])) ++expect;
                if (!held.count(all[b][x])) ++absent;
            }
            assert(filter.countAbsent(leading[b].data(), leading[b].size()) == expect);
            std::vector<uint64_t> sorted = leading[b];
            sortHashes(sorted, scratch);
            assert(filter.countAbsent(sorted.data(), sorted.size()) == expect);
            assert(expect <= absent);
            falsePositives += absent - expect;
        }
    }
    return falsePositives;
}
void testCountAbsent() {
    std::cerr << "START: Test countAbsent equals the count of contains() == false, both layouts... \t";
    struct Setting { unsigned bits; size_t elements; unsigned nhash, k; };
    const Setting settings[] = {
        {8, 50000, 3, 31}, {16, 50000, 3, 31}, {16, 5000000, 3, 81}, {8, 50000, 1, 31},
        {6, 50000, 2, 21}, {4, 50000, 4, 31}, {12, 50000, 5, 81}, {16, 20000, 8, 64},
    };
    for (const Setting& s : settings) {
        const size_t fpBlocked = checkCountAbsent<BlockedBloomFilter>(s.bits, s.elements, s.nhash, s.k);
        const size_t fpStandard = checkCountAbsent<StandardBloomFilter>(s.bits, s.elements, s.nhash, s.k);
        if (&s == settings) assert(fpBlocked > 100 && fpStandard > 100);
    }
    // No hash function: contains() holds for every k-mer.
    const uint64_t any[] = {0, 1, ~uint64_t(0)};
    assert(BlockedBloomFilter(1000, 0, 31).countAbsent(any, 3) == 0);
    assert(StandardBloomFilter(1000, 0, 31).countAbsent(any, 3) == 0);
    assert(BlockedBloomFilter(1000, 3, 31).countAbsent(any, 3) == 3);
    assert(StandardBloomFilter(1000, 3, 31).countAbsent(any, 3) == 3);
    assert(BlockedBloomFilter(1000, 3, 31).countAbsent(any, 0) == 0);
    std::cerr << "PASSED: Test countAbsent equals the count of contains() == false, both layouts\n";
}

// The -m budget counts a refSet2 genome as 8 bytes per byte of the largest refSet2
// file, whatever the hash number: the figures of the README (1 717 genomes of
// 5 000 000 bytes at the default 64 GB) and of the 120 E. coli at -m 1 (largest
// file 5 797 390 bytes: 23 a batch, 6 batches). The comparison tests below run at
// budgets 64 and 0, all of refSet2 or one genome, and cannot tell one accounting
// from another: with 8 bytes for each hash function, as a genome counted while
// all hash values were kept, or with one byte, they pass.
void testBatchSize() {
    std::cerr << "START: Test the budget counts 8 bytes per byte of the largest refSet2 file... \t";
    assert(refSet2BatchSize(64, 5000000, 100000) == 1717);
    assert(refSet2BatchSize(64, 5000000, 1717) == 1717);
    assert(refSet2BatchSize(64, 5000000, 1000) == 1000);       // all of refSet2
    assert(refSet2BatchSize(64, 6328931, 5493) == 1357);       // 5 493 E. coli, largest file 6 328 931 bytes
    assert(refSet2BatchSize(1, 5797390, 120) == 23);           // 120 E. coli at -m 1
    assert(refSet2BatchSize(1, 134217728, 120) == 1);          // one genome fills the budget
    assert(refSet2BatchSize(1, 134217729, 120) == 1);          // a genome above the budget: one a batch
    assert(refSet2BatchSize(1, 67108864, 120) == 2);
    assert(refSet2BatchSize(0, 5000000, 120) == 1);            // -m 0, --low-mem
    assert(refSet2BatchSize(64, 0, 7) == 7);                   // no file size: all of refSet2
    std::cerr << "PASSED: Test the budget counts 8 bytes per byte of the largest refSet2 file\n";
}

// One row of the similarity table as identifyDifference writes it.
static std::string similarityRow(const std::string& name2, const std::string& name1,
                                 size_t countUnique, size_t countAll, double score) {
    std::ostringstream row;
    row << name2 << "\t" << name1 << "\t" << countUnique << "\t" << countAll << "\t" << score;
    return row.str();
}

// The table's rows as the functions the comparison used before give them
// (loadBloomFilter, checkBloomFilter: all nhash values of every position, insert()
// and contains()), refSet1 outer, refSet2 inner.
static std::vector<std::string> expectedRows(const CompareSketchConfig& cfg,
                                             const std::vector<std::string>& refSet1,
                                             const std::vector<std::string>& refSet2) {
    std::vector<std::string> rows;
    for (const std::string& ref1 : refSet1) {
        BloomFilter dbFilter(cfg.bits * cfg.dbfSize, cfg.nhash, cfg.kmerLen);
        loadBloomFilter(cfg, ref1, dbFilter);
        for (const std::string& ref2 : refSet2) {
            size_t countAll = 0, countUnique = 0;
            const double score = checkBloomFilter(cfg, ref2, dbFilter, countAll, countUnique);
            rows.push_back(similarityRow(getBaseId(ref2), getBaseId(ref1), countUnique, countAll, score));
        }
    }
    return rows;
}

// Run identifyDifference and return the lines of its table; `log` gets what it
// wrote to stderr.
static std::vector<std::string> runIdentifyDifference(const CompareSketchConfig& cfg,
                                                      const std::vector<std::string>& refSet1,
                                                      const std::vector<std::string>& refSet2,
                                                      std::string& log) {
    std::ostringstream captured;
    std::streambuf* keep = std::cerr.rdbuf(captured.rdbuf());
    identifyDifference(cfg, refSet1, refSet2);
    std::cerr.rdbuf(keep);
    log = captured.str();
    std::vector<std::string> lines;
    std::ifstream in(cfg.outfile);
    for (std::string line; std::getline(in, line);) lines.push_back(line);
    return lines;
}

// The comparison's path (insertReference, checkBloomFilterFromLeading,
// identifyDifference) against the functions it used before, on the filter layout
// this suite is built with: the filter of a reference bit for bit, the three
// numbers of every pair, and the table, at the default sizing, with thin filters
// (false positives in most rows) and with other numbers of hash functions and
// other k. A budget of 0 takes refSet2 one genome a batch; the table is the same
// file. refSet2 names a reference twice and one that holds no k-mer. With one
// refSet1 genome identifyDifference leaves the leading hashes in file order, with
// four and with nine it sorts them (the blocked filter); the settings take the
// three in turn.
void testComparisonPath() {
    std::cerr << "START: Test the comparison gives the numbers of loadBloomFilter and checkBloomFilter... \t";
    const std::vector<std::string> four = {cmpRefs[0], cmpRefs[1], cmpRefs[2], cmpRefs[4]};
    const std::vector<std::string> nine = {cmpRefs[0], cmpRefs[1], cmpRefs[2], cmpRefs[4], cmpRefs[3],
                                           cmpRefs[3], cmpRefs[2], cmpRefs[4], cmpRefs[3]};
    const std::vector<std::string> one = {cmpRefs[1]};
    const std::vector<std::string> refSet2 = {cmpRefs[1], cmpRefs[0], cmpRefs[3], cmpRefs[4], cmpRefs[2], cmpRefs[3]};
    struct Setting { unsigned bits; size_t elements; unsigned nhash, k; };
    const Setting settings[] = {{16, 5000000, 3, 81}, {8, 50000, 3, 31}, {6, 40000, 2, 21}, {10, 50000, 5, 64}, {16, 50000, 1, 31}};
    std::vector<uint64_t> scratch;
    size_t lowRows = 0;
    for (const Setting& s : settings) {
        const size_t turn = (&s - settings) % 3;
        const std::vector<std::string>& refSet1 = turn == 0 ? four : turn == 1 ? nine : one;
        CompareSketchConfig cfg;
        cfg.bits = s.bits;
        cfg.dbfSize = s.elements;
        cfg.nhash = s.nhash;
        cfg.kmerLen = s.k;
        cfg.outfile = "cmp_out.tsv";
        for (const std::string& ref1 : four) {
            BloomFilter before(cfg.bits * cfg.dbfSize, cfg.nhash, cfg.kmerLen);
            BloomFilter now(cfg.bits * cfg.dbfSize, cfg.nhash, cfg.kmerLen);
            loadBloomFilter(cfg, ref1, before);
            insertReference(ref1, now);
            assert(before.sizeBytes() == now.sizeBytes());
            assert(std::equal(before.getFilter(), before.getFilter() + before.sizeBytes(), now.getFilter()));
            for (const std::string& ref2 : refSet2) {
                size_t allBefore = 0, uniqueBefore = 0, allNow = 0, uniqueNow = 0;
                const double scoreBefore = checkBloomFilter(cfg, ref2, before, allBefore, uniqueBefore);
                std::vector<uint64_t> leading = loadLeadingHashes(ref2, cfg.kmerLen);
                const double scoreNow = checkBloomFilterFromLeading(leading, now, allNow, uniqueNow);
                assert(allNow == allBefore && uniqueNow == uniqueBefore && scoreNow == scoreBefore);
                sortHashes(leading, scratch);
                assert(checkBloomFilterFromLeading(leading, now, allNow, uniqueNow) == scoreBefore);
                assert(allNow == allBefore && uniqueNow == uniqueBefore);
            }
        }
        const std::vector<std::string> expect = expectedRows(cfg, refSet1, refSet2);
        // Thin filters: a row that is neither "nothing found" nor "all found".
        for (const std::string& row : expect) {
            std::istringstream fields(row);
            std::string name2, name1;
            size_t countUnique = 0, countAll = 0;
            fields >> name2 >> name1 >> countUnique >> countAll;
            if (countUnique > 0 && countUnique < countAll) ++lowRows;
        }
        for (size_t budget : {64u, 0u}) {
            cfg.memBudgetGB = budget;
            std::string log;
            const std::vector<std::string> lines = runIdentifyDifference(cfg, refSet1, refSet2, log);
            assert(log.find(budget ? "refSet2 batch=6/6" : "refSet2 batch=1/6") != std::string::npos);
            assert(lines.size() == 1 + expect.size());
            assert(lines[0] == "ref1\tref2\tuniq_ref1\tall_ref1\tsimilarity");
            assert(std::equal(expect.begin(), expect.end(), lines.begin() + 1));
        }
    }
    assert(lowRows > 20);
    std::remove("cmp_out.tsv");
    std::cerr << "PASSED: Test the comparison gives the numbers of loadBloomFilter and checkBloomFilter\n";
}

// The table is in list order, refSet1 outer and refSet2 inner, whatever the number
// of threads and the budget. Before, the rows were written as they were computed:
// with one refSet2 genome a batch (a budget of 0) refSet2 became the outer loop,
// and with several threads a refSet1 genome's rows came when its thread finished.
// The first refSet1 genome here has 300 000 bases and the ones after it 6 to
// 48 000, so that with several threads they are complete, and have to wait, while
// it is still being read. Either of the old orders fails this test, the first of
// them on one thread.
void testRowOrder() {
    std::cerr << "START: Test rows are in list order on any thread count and budget... \t";
    const std::string big = "cmp_big.fa";
    { std::ofstream o(big); o << ">big\n" << dna(300000) << "\n"; }
    const std::vector<std::string> refSet1 = {big, cmpRefs[3], cmpRefs[2], cmpRefs[4], cmpRefs[3], cmpRefs[2], cmpRefs[3], cmpRefs[1]};
    const std::vector<std::string> refSet2 = {cmpRefs[1], cmpRefs[3], cmpRefs[4], cmpRefs[2], cmpRefs[0]};
    CompareSketchConfig cfg;
    cfg.bits = 8;
    cfg.dbfSize = 300000;
    cfg.kmerLen = 31;
    cfg.outfile = "cmp_order.tsv";
    const std::vector<std::string> expect = expectedRows(cfg, refSet1, refSet2);
    assert(expect.size() == refSet1.size() * refSet2.size());
#ifdef _OPENMP
    const int savedThreads = omp_get_max_threads();
    const std::vector<int> threadCounts = {1, 4, 8};
#else
    const std::vector<int> threadCounts = {1};      // built without OpenMP: the budget's order only
#endif
    for (int threads : threadCounts) {
#ifdef _OPENMP
        omp_set_num_threads(threads);
#endif
        for (size_t budget : {64u, 0u}) {
            cfg.memBudgetGB = budget;
            for (int run = 0; run < (threads > 1 ? 3 : 1); run++) {
                std::string log;
                const std::vector<std::string> lines = runIdentifyDifference(cfg, refSet1, refSet2, log);
                assert(log.find(budget ? "refSet2 batch=5/5" : "refSet2 batch=1/5") != std::string::npos);
                assert(lines.size() == 1 + expect.size());
                assert(std::equal(expect.begin(), expect.end(), lines.begin() + 1));
            }
        }
    }
#ifdef _OPENMP
    omp_set_num_threads(savedThreads);
#endif
    std::remove(big.c_str());
    std::remove(cfg.outfile.c_str());
    std::cerr << "PASSED: Test rows are in list order on any thread count and budget\n";
}

int main() {
#ifdef _OPENMP
    // One thread, unless a test asks for more: identifyDifference takes the number
    // OpenMP has, every core of the machine by default, and a filter for each.
    omp_set_num_threads(1);
#endif
    CompareSketchConfig cfg;
    cfg.kmerLen = 10;
    cfg.bits = 128;
    cfg.nhash = 3;
    cfg.dbfSize = 100;

    testLoadBloomFilter(cfg);
    testCheckBloomFilterIdentical(cfg);
    testCheckBloomFilterDifferent(cfg);
    testUnusableOptions();
    testDerivedHashes();
    testSortHashes();
    testBatchSize();
    writeCmpRefs();
    testCountAbsent();
    testComparisonPath();
    testRowOrder();
    removeCmpRefs();

    std::cerr << "CompareSketch: All tests PASSED!\n\n";
    return 0;
}