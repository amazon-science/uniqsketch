#include <iostream>
#include <cstdio>
#include <string>
#include <cassert>
#include <atomic>
#include <thread>
#include <vector>

#include "BloomFilter.hpp"

namespace opt {
size_t size;
unsigned hash;
unsigned k;
}

void initialize() {
    opt::size = 65536;
    opt::hash = 3;
    opt::k = 32;
}

void testInsertContain() {
    std::cerr << "START: Test insert, contain... \t";

    BloomFilter bFilter(opt::size, opt::hash, opt::k);
    uint64_t hVec1[] = {1, 3000, 50000000};
    uint64_t hVec2[] = {2000000, 700000000, 4000000000};

    bFilter.insert(hVec1);
    assert(bFilter.contains(hVec1) == true);
    assert(bFilter.contains(hVec2) == false);

    std::string kmer1("GCCTAGCTAGCTTTTAGCTGGGATTTTTT");
    std::string kmer2("ATTTAGCTAGCGCTAAGCTGGGACCCTGC");

    bFilter.insert(kmer1.c_str());
    assert(bFilter.contains(kmer1.c_str()) == true);
    assert(bFilter.contains(kmer2.c_str()) == false);

    std::cerr << "PASSED: Test insert, contain\n";
}

void testInsertMakeChange() {
    std::cerr << "START: Test insert_make_change... \t";

    BloomFilter bFilter(opt::size, opt::hash, opt::k);
    uint64_t hVec1[] = {34987645632, 89, 37692833};
    uint64_t hVec2[] = {34359738368, 4398046511238, 678};

    bFilter.insert(hVec1);
    assert(bFilter.insert_make_change(hVec1) == false);
    assert(bFilter.insert_make_change(hVec2) == true);

    std::cerr << "PASSED: Test insert_make_change\n";
}

// insert_make_change must report true to exactly one caller for a given k-mer,
// whatever the thread scheduling. uniqsketch's cascading distinct/solid filter
// pair depends on this: a k-mer shared by two references has to be seen as
// already present by the second caller, so it is recorded as a repeat and
// excluded from both references' signatures. Setting the individual bits
// atomically is not sufficient on its own, because two threads inserting the
// same k-mer can each win a subset of the bits and both conclude they inserted
// it -- which silently turns a shared k-mer into a signature for both
// references, and makes the index depend on thread timing.
void testInsertMakeChangeConcurrent() {
    std::cerr << "START: Test insert_make_change concurrency... \t";

    const unsigned NTHREADS = 16;
    const unsigned NTRIALS = 300;
    // Many hash functions widen the window between setting the first and last
    // bit, which is the interval in which the unfixed implementation lets two
    // threads each claim a subset of the bits.
    const unsigned NHASH = 64;

    for (unsigned trial = 0; trial < NTRIALS; trial++) {
        BloomFilter bFilter(opt::size, NHASH, opt::k);
        std::vector<uint64_t> hVec(NHASH);
        for (unsigned i = 0; i < NHASH; i++) {
            hVec[i] = (trial + 1) * 1299709ULL + i * 7919ULL;
        }

        // Release all threads at once so they collide inside the same call,
        // rather than each finishing before the next is spawned.
        std::atomic<unsigned> arrived(0);
        std::atomic<bool> go(false);
        std::vector<unsigned> won(NTHREADS, 0);
        std::vector<std::thread> pool;
        pool.reserve(NTHREADS);

        for (unsigned t = 0; t < NTHREADS; t++) {
            pool.emplace_back([&bFilter, &won, &hVec, &arrived, &go, t]() {
                arrived.fetch_add(1);
                while (!go.load()) {
                    std::this_thread::yield();
                }
                if (bFilter.insert_make_change(hVec.data())) {
                    won[t] = 1;
                }
            });
        }
        while (arrived.load() < NTHREADS) {
            std::this_thread::yield();
        }
        go.store(true);
        for (auto& th : pool) {
            th.join();
        }

        unsigned winners = 0;
        for (unsigned t = 0; t < NTHREADS; t++) {
            winners += won[t];
        }
        assert(winners == 1);
        assert(bFilter.contains(hVec.data()) == true);
    }

    std::cerr << "PASSED: Test insert_make_change concurrency\n";
}

// Test that an empty Bloom filter contains nothing
void testEmptyFilter() {
    std::cerr << "START: Test empty filter... \t";

    BloomFilter bFilter(opt::size, opt::hash, opt::k);
    uint64_t hVec[] = {42, 100, 999};
    assert(bFilter.contains(hVec) == false);
    assert(bFilter.get_pop() == 0);

    std::cerr << "PASSED: Test empty filter\n";
}

// Test get_pop returns correct population count
void testGetPop() {
    std::cerr << "START: Test get_pop... \t";

    BloomFilter bFilter(opt::size, opt::hash, opt::k);
    assert(bFilter.get_pop() == 0);

    uint64_t hVec1[] = {1, 3000, 50000000};
    bFilter.insert(hVec1);
    size_t pop1 = bFilter.get_pop();
    assert(pop1 > 0);
    // With 3 hash functions, at most 3 bits set (could be fewer if collisions)
    assert(pop1 <= opt::hash);

    // Inserting the same element again should not change population
    bFilter.insert(hVec1);
    assert(bFilter.get_pop() == pop1);

    // Inserting a new element should increase population
    uint64_t hVec2[] = {2000000, 700000000, 4000000000};
    bFilter.insert(hVec2);
    assert(bFilter.get_pop() >= pop1);

    std::cerr << "PASSED: Test get_pop\n";
}

// Test store and load roundtrip
void testStoreLoad() {
    std::cerr << "START: Test store/load... \t";

    BloomFilter bFilter(opt::size, opt::hash, opt::k);
    uint64_t hVec1[] = {1, 3000, 50000000};
    uint64_t hVec2[] = {2000000, 700000000, 4000000000};

    bFilter.insert(hVec1);
    bFilter.insert(hVec2);

    const char* tmpPath = "test_bf_roundtrip.bin";
    bFilter.storeFilter(tmpPath);

    // Load into a new filter and verify contents match
    BloomFilter loaded(opt::size, opt::hash, opt::k, tmpPath);
    assert(loaded.contains(hVec1) == true);
    assert(loaded.contains(hVec2) == true);

    uint64_t hVec3[] = {999999, 888888, 777777};
    assert(loaded.contains(hVec3) == false);

    assert(loaded.get_pop() == bFilter.get_pop());

    std::remove(tmpPath);

    std::cerr << "PASSED: Test store/load\n";
}

int main() {
    initialize();
    testInsertContain();
    testInsertMakeChange();
    testInsertMakeChangeConcurrent();
    testEmptyFilter();
    testGetPop();
    testStoreLoad();

    std::cerr << "BloomFilter: All tests PASSED!\n\n";
    return 0;
}
