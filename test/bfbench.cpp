// Microbenchmark: standard vs blocked Bloom filter at matched bits/element.
// Measures insert and probe throughput on random 81-mer hashes at a size large
// enough to defeat L3, and the empirical false-positive rate on k-mers that were
// never inserted. Not part of the test suite; built by the experiment script.
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

#include "BloomFilter.hpp"          // StandardBloomFilter, and BlockedBloomFilter unless opted out
#include "BlockedBloomFilter.hpp"   // explicit, so the bench builds under -DUNIQSKETCH_STANDARD_BF too

static double secs(std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
}

// Synthetic multi-hash: independent 64-bit values per "hash function", drawn once
// so both filters see identical input. Mirrors what ntHash hands the filter.
struct HashSet {
    std::vector<uint64_t> v;   // n * m, row-major
    unsigned m;
    const uint64_t* row(size_t i) const { return v.data() + i * m; }
};
static HashSet makeHashes(size_t n, unsigned m, uint64_t seed) {
    HashSet h; h.m = m; h.v.resize(n * m);
    std::mt19937_64 rng(seed);
    for (auto& x : h.v) x = rng();
    return h;
}

// Both filters get the SAME number of bits. The blocked filter rounds its block
// count up to a power of two, so pick nElem such that nElem*b is already one.
template <class BF>
static void bench(const char* name, size_t nElem, unsigned bitsPerElem, unsigned m,
                  const HashSet& ins, const HashSet& abs) {
    BF bf(static_cast<size_t>(nElem) * bitsPerElem, m, 81);
    auto t0 = std::chrono::steady_clock::now();
    for (size_t i = 0; i < nElem; i++) bf.insert(ins.row(i));
    auto t1 = std::chrono::steady_clock::now();
    size_t hitsIn = 0;
    for (size_t i = 0; i < nElem; i++) hitsIn += bf.contains(ins.row(i));
    auto t2 = std::chrono::steady_clock::now();
    size_t fp = 0;
    for (size_t i = 0; i < nElem; i++) fp += bf.contains(abs.row(i));
    auto t3 = std::chrono::steady_clock::now();
    double ins_ns = 1e9 * secs(t0, t1) / nElem, hit_ns = 1e9 * secs(t1, t2) / nElem,
           miss_ns = 1e9 * secs(t2, t3) / nElem;
    std::printf("%-9s b=%3u m=%u  insert %6.1f ns  probe-hit %6.1f ns  probe-miss %6.1f ns  "
                "recall %zu/%zu  FPR %.3e  (fill %.4f)\n",
                name, bitsPerElem, m, ins_ns, hit_ns, miss_ns, hitsIn, nElem,
                static_cast<double>(fp) / nElem,
                static_cast<double>(bf.get_pop()) / (static_cast<double>(nElem) * bitsPerElem));
}

int main(int argc, char** argv) {
    // 2^24 = 16.8M elements: at b=32/64/128 the bit count is 2^29/2^30/2^31,
    // a power of two, so the blocked filter's rounding is a no-op and both
    // filters occupy exactly 64 / 128 / 256 MB. All exceed the 32 MB L3.
    size_t n = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : (size_t(1) << 24);
    unsigned m = 5;
    std::printf("# %zu elements, %u hashes; equal memory for both filters: %.0f/%.0f/%.0f MB\n",
                n, m, n * 32.0 / 8 / 1e6, n * 64.0 / 8 / 1e6, n * 128.0 / 8 / 1e6);
    HashSet ins = makeHashes(n, m, 1), abs = makeHashes(n, m, 2);
    for (unsigned b : {32u, 64u, 128u}) {
        bench<StandardBloomFilter>("standard", n, b, m, ins, abs);
        bench<BlockedBloomFilter>("blocked", n, b, m, ins, abs);
    }
    return 0;
}
