/*
 * Copyright Hamid Mohamadi.
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License. See the LICENSE accompanying this file
 * for the specific language governing permissions and limitations under
 * the License.
 */

// Cache-line-blocked Bloom filter, API-compatible with StandardBloomFilter.
//
// The standard filter touches one random cache line per hash function and does
// a 64-bit modulo for each. Index building spends most of its time in that
// traffic: the filters are far larger than cache and 89% of k-mer sightings on
// a bacterial set are repeats, which probe all m lines and then insert m more
// into the solid filter. Here the leading hash selects one 512-bit block and the
// remaining hashes select bits inside it, so an insert or a hit is one cache
// line. A miss was already ~1 line on the standard filter (early exit on the
// first unset bit at low fill), so the miss-dominated phases, checkRef and the
// 243 margin probes per candidate, are unchanged.
//
// Block selection is a multiply-high, (h0 * nBlocks) >> 64, which is unbiased
// for a uniform 64-bit hash and needs no division. The filter therefore has
// exactly the requested size: no power-of-two rounding, no memory overhead.
//
// The trade-off is false-positive rate at equal memory: a k-mer's bits are
// confined to one block, so blocks fill unevenly. Measured at b=128, m=5:
// 4.8e-7 against 1.8e-7 for the standard filter. A false positive in the
// distinct or solid filter drops a candidate; on a 120-genome E. coli build
// it dropped 11 of 25.9 M, and the 280-mixture sweep was identical to the
// decimal (bbf_experiment/RESULTS.md).
//
// Partitions. Block locality also makes a filter divisible: the block count is
// a multiple of PARTITIONS, so the k-mers whose leading hash falls in the p-th
// 1/PARTITIONS of the hash space touch the p-th 1/PARTITIONS of the blocks and
// nothing else, in this filter and in every other filter (whatever its size).
// A scan that hands partition p to one thread therefore gives each block a
// single writer, and the thread can replay the serial order of operations on
// it exactly; see loadBFfaGroupParallel. The rounding costs at most
// PARTITIONS - 1 blocks (16 KB) per filter.

#ifndef BLOCKEDBLOOMFILTER_HPP_
#define BLOCKEDBLOOMFILTER_HPP_

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <new>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "nthash.hpp"

class BlockedBloomFilter {
public:
    static constexpr size_t BLOCK_BITS  = 512;          // one cache line
    static constexpr size_t BLOCK_WORDS = BLOCK_BITS / 64;
    static constexpr size_t BLOCK_BYTES = BLOCK_BITS / 8;
    static constexpr unsigned PARTITION_BITS = 8;       // leading bits of h0 that name a partition
    static constexpr unsigned PARTITIONS = 1u << PARTITION_BITS;
    static constexpr bool kPartitioned = true;          // partition() isolates blocks (see header)

    /**
     * Construct an empty filter of (at least) filterSize bits, rounded up to a
     * whole number of 512-bit blocks and then to a multiple of PARTITIONS blocks.
     */
    BlockedBloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize)
        : m_filter(nullptr), m_hashNum(hashNum), m_kmerSize(kmerSize) {
        m_blocks = (filterSize + BLOCK_BITS - 1) / BLOCK_BITS;
        m_blocks = (m_blocks + PARTITIONS - 1) / PARTITIONS * PARTITIONS;
        if (m_blocks == 0) m_blocks = PARTITIONS;
        m_size = m_blocks * BLOCK_BITS;
        const size_t bytes = m_blocks * BLOCK_BYTES;
        void* p = nullptr;
        if (posix_memalign(&p, BLOCK_BYTES, bytes) != 0 || !p)
            throw std::bad_alloc();
        m_filter = static_cast<uint64_t*>(p);
        // Zeroing a gigabyte-scale filter (a tree level's distinct filter, a
        // group's "seen" filter) on one thread held every other thread at the
        // start of each phase. Clear it in parallel when no parallel region is
        // already running; leaving the pages to the kernel's first-touch zeroing
        // instead was slower, since the faults then land inside the hot loop.
        // (The first-touch work is the same either way; MAP_POPULATE would just
        // move it back onto one thread.)
        constexpr size_t CHUNK = size_t(1) << 22;         // 4 MB
        const size_t chunks = (bytes + CHUNK - 1) / CHUNK;
#ifdef _OPENMP
        const bool par = chunks >= 8 && !omp_in_parallel();
#else
        const bool par = false;
#endif
        #pragma omp parallel for schedule(static) if(par)
        for (size_t c = 0; c < chunks; c++) {
            std::memset(reinterpret_cast<char*>(p) + c * CHUNK, 0, std::min(CHUNK, bytes - c * CHUNK));
        }
    }

    /** Load a filter from a file written by storeFilter. */
    BlockedBloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize, const char* fPath)
        : BlockedBloomFilter(filterSize, hashNum, kmerSize) {
        std::ifstream file(fPath, std::ios::in | std::ios::binary);
        file.read(reinterpret_cast<char*>(m_filter), m_blocks * BLOCK_BYTES);
    }

    BlockedBloomFilter(const BlockedBloomFilter&) = delete;
    BlockedBloomFilter& operator=(const BlockedBloomFilter&) = delete;
    ~BlockedBloomFilter() { std::free(m_filter); }

    /**
     * Insert and report whether the k-mer was new. Same contract as the standard
     * filter: exactly one caller observes true for a given k-mer, enforced by the
     * same striped lock around the test-then-insert sequence.
     */
    bool insert_make_change(const uint64_t* hVal) {
        if (contains(hVal)) return false;
        std::lock_guard<std::mutex> lock(stripe_lock(hVal[0]));
        if (contains(hVal)) return false;
        insert(hVal);
        return true;
    }

    void insert(const uint64_t* hVal) {
        uint64_t* blk = block(hVal[0]);
        for (unsigned i = 0; i < m_hashNum; i++) {
            unsigned b = bitIndex(hVal, i);
            __sync_or_and_fetch(&blk[b >> 6], uint64_t(1) << (b & 63));
        }
    }

    /**
     * insert_make_change / insert for a caller that is the only thread touching
     * this k-mer's block: no stripe lock, plain stores. The partitioned scan
     * (loadBFfaGroupParallel) gives each block one writer by construction, and
     * the stripe locks, keyed by the low bits of h0 while partitions are its
     * high bits, would otherwise bounce between every thread.
     */
    bool insert_make_change_exclusive(const uint64_t* hVal) {
        if (contains(hVal)) return false;
        insert_exclusive(hVal);
        return true;
    }

    void insert_exclusive(const uint64_t* hVal) {
        uint64_t* blk = block(hVal[0]);
        for (unsigned i = 0; i < m_hashNum; i++) {
            unsigned b = bitIndex(hVal, i);
            blk[b >> 6] |= uint64_t(1) << (b & 63);
        }
    }

    void insert(const char* kmer) {
        uint64_t hVal[m_hashNum];
        NTMC64(kmer, m_kmerSize, m_hashNum, hVal);
        insert(hVal);
    }

    bool contains(const uint64_t* hVal) const {
        const uint64_t* blk = block(hVal[0]);
        for (unsigned i = 0; i < m_hashNum; i++) {
            unsigned b = bitIndex(hVal, i);
            if ((blk[b >> 6] & (uint64_t(1) << (b & 63))) == 0) return false;
        }
        return true;
    }

    bool contains(const char* kmer) const {
        uint64_t hVal[m_hashNum];
        NTMC64(kmer, m_kmerSize, m_hashNum, hVal);
        return contains(hVal);
    }

    void storeFilter(const char* fPath) const {
        std::ofstream file(fPath, std::ios::out | std::ios::binary);
        file.write(reinterpret_cast<const char*>(m_filter), m_blocks * BLOCK_BYTES);
    }

    /** Number of set bits. */
    size_t get_pop() const {
        size_t pop = 0;
        const size_t words = m_blocks * BLOCK_WORDS;
        #pragma omp parallel for reduction(+:pop)
        for (size_t i = 0; i < words; i++) pop += __builtin_popcountll(m_filter[i]);
        return pop;
    }

    unsigned getHashNum() const { return m_hashNum; }
    unsigned getKmerSize() const { return m_kmerSize; }
    unsigned char* getFilter() { return reinterpret_cast<unsigned char*>(m_filter); }
    size_t sizeBits() const { return m_size; }
    size_t sizeBytes() const { return m_blocks * BLOCK_BYTES; }   // the allocation, not the request
    size_t blockCount() const { return m_blocks; }

    /**
     * Partition of a k-mer: the leading PARTITION_BITS of h0, i.e. the high 64
     * bits of h0 * PARTITIONS. Block selection is the high 64 bits of
     * h0 * nBlocks with nBlocks = c * PARTITIONS, so block b lies entirely in
     * partition floor(b / c): two k-mers that share a block share a partition,
     * in this filter and in any other.
     */
    static unsigned partition(uint64_t h0) {
        return static_cast<unsigned>(h0 >> (64 - PARTITION_BITS));
    }

private:
    /**
     * Block for a k-mer: the high 64 bits of h0 * nBlocks, which maps a uniform
     * 64-bit hash onto [0, nBlocks) without bias and without a division.
     */
    uint64_t* block(uint64_t h0) const {
        unsigned __int128 p = static_cast<unsigned __int128>(h0) * m_blocks;
        return m_filter + static_cast<size_t>(p >> 64) * BLOCK_WORDS;
    }
    /**
     * Bit i within the block. hVal[1..m-1] are independent ntHash multi-hash
     * values. hVal[0] chose the block from its high bits, so its low 9 bits are
     * still free to pick the first in-block position.
     */
    static unsigned bitIndex(const uint64_t* hVal, unsigned i) {
        return static_cast<unsigned>(hVal[i] & (BLOCK_BITS - 1));
    }
    static std::mutex& stripe_lock(uint64_t h) {
        static constexpr size_t STRIPES = 4096;
        static std::mutex locks[STRIPES];
        return locks[h & (STRIPES - 1)];
    }

    uint64_t* m_filter;
    size_t m_size;
    size_t m_blocks;
    unsigned m_hashNum;
    unsigned m_kmerSize;
};

#endif // BLOCKEDBLOOMFILTER_HPP_
