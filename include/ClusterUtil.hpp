#ifndef CLUSTERUTIL_HPP_
#define CLUSTERUTIL_HPP_

#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <climits>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <fstream>
#include <numeric>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "SequenceUtil.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"
#include "kseq_gz.h"

/**
 * Load all k-mer hashes from a FASTA/FASTA.gz reference into a flat vector.
 *
 * @param fPath  Path to a reference file (plain or gzipped).
 * @param nhash  Number of hash functions.
 * @param kmerLen  K-mer length.
 * @return Flattened vector of hash values (nhash values per k-mer).
 */
inline std::vector<uint64_t> clusterLoadKmerHashes(const std::string& fPath,
                                                   unsigned nhash, unsigned kmerLen) {
    std::vector<uint64_t> hashes;
    gzFile fp = gzopen(fPath.c_str(), "r");
    kseq_t* ks = kseq_init(fp);
    while (kseq_read(ks) >= 0) {
        std::string seq(ks->seq.s, ks->seq.l);
        ntHashIterator itr(seq, nhash, kmerLen);
        while (itr != itr.end()) {
            const uint64_t* hVec = *itr;
            for (unsigned h = 0; h < nhash; h++) {
                hashes.push_back(hVec[h]);
            }
            ++itr;
        }
    }
    kseq_destroy(ks);
    gzclose(fp);
    return hashes;
}

/**
 * Exact O(n^2) clustering: pairwise unique k-mer counts via a Bloom filter, merged
 * with single-linkage union-find. References with fewer than `threshold` unique
 * k-mers between them are merged. Precise, but quadratic in the number of
 * references (and it caches every reference's k-mer hashes).
 *
 * @param refFiles  All reference file paths.
 * @param threshold  Minimum unique k-mers to consider two references as distinct.
 * @param kmerLen  K-mer length.
 * @param nhash  Number of hash functions for Bloom filter.
 * @param bits  Bits per element in Bloom filter.
 * @param dbfSize  Approximate genome size for Bloom filter sizing.
 * @param clusterFile  Output TSV path mapping each reference to its cluster representative.
 * @return Vector of representative file paths (one per cluster).
 */
inline std::vector<std::string> clusterReferencesExact(
        const std::vector<std::string>& refFiles,
        unsigned threshold,
        unsigned kmerLen,
        unsigned nhash,
        unsigned bits,
        size_t dbfSize,
        const std::string& clusterFile) {
    size_t n = refFiles.size();
    if (n <= 1) return refFiles;

    std::cerr << "Clustering " << n << " references (threshold=" << threshold
              << " unique k-mers, k=" << kmerLen << ")...\n";

    // Pre-compute k-mer hashes for all references
    std::vector<std::vector<uint64_t>> hashSets(n);
    std::vector<std::string> names(n);

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < n; i++) {
        hashSets[i] = clusterLoadKmerHashes(refFiles[i], nhash, kmerLen);
        names[i] = getBaseId(refFiles[i]);
    }

    // Union-Find for single-linkage clustering
    std::vector<unsigned> parent(n);
    std::iota(parent.begin(), parent.end(), 0);

    std::function<unsigned(unsigned)> find = [&](unsigned x) -> unsigned {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };

    size_t bfBits = bits * dbfSize;

    #pragma omp parallel
    {
        BloomFilter dbFilter(bfBits, nhash, kmerLen);
        unsigned char* filterPtr = dbFilter.getFilter();
        size_t filterBytes = dbFilter.sizeBytes();   // the allocation: the blocked filter rounds up

        #pragma omp for schedule(dynamic)
        for (unsigned i = 0; i < n; i++) {
            std::memset(filterPtr, 0, filterBytes);
            for (size_t h = 0; h < hashSets[i].size(); h += nhash) {
                dbFilter.insert(&hashSets[i][h]);
            }

            for (unsigned j = i + 1; j < n; j++) {
                // Count how many k-mers in j are NOT in i's BF
                size_t countUnique = 0;
                size_t countAll = hashSets[j].size() / nhash;
                for (size_t h = 0; h < hashSets[j].size(); h += nhash) {
                    if (!dbFilter.contains(&hashSets[j][h])) {
                        ++countUnique;
                    }
                }

                if (countUnique < threshold) {
                    #pragma omp critical(union_op)
                    {
                        unsigned ri = find(i), rj = find(j);
                        if (ri != rj) parent[ri] = rj;
                    }
                }
            }
        }
    }

    // Build clusters
    std::unordered_map<unsigned, std::vector<unsigned>> clusters;
    for (unsigned i = 0; i < n; i++) {
        clusters[find(i)].push_back(i);
    }

    // Pick representative per cluster: member with the most k-mers
    std::vector<std::string> representatives;
    std::ofstream clusterOut(clusterFile);
    clusterOut << "cluster_id\trepresentative\tmember\tnum_kmers\n";

    unsigned clusterId = 0;
    for (auto& [root, members] : clusters) {
        unsigned repIdx = members[0];
        size_t maxKmers = hashSets[repIdx].size() / nhash;
        for (unsigned m : members) {
            size_t nk = hashSets[m].size() / nhash;
            if (nk > maxKmers) {
                maxKmers = nk;
                repIdx = m;
            }
        }
        representatives.push_back(refFiles[repIdx]);

        for (unsigned m : members) {
            clusterOut << clusterId << "\t" << names[repIdx] << "\t" << names[m]
                       << "\t" << hashSets[m].size() / nhash << "\n";
        }
        ++clusterId;
    }
    clusterOut.close();

    std::cerr << "Clustering complete: " << n << " references -> "
              << representatives.size() << " clusters\n";

    return representatives;
}

/**
 * Default size of the bottom-s MinHash sketch from which --cluster and --tree
 * estimate the distance between two references: the s smallest distinct k-mer
 * hashes of each. `uniqsketch --minhash-size=N` hands another size to the two
 * functions that sketch, clusterReferencesMinHash and buildHierarchicalTree, which
 * take this one when none is given; uniqsketch reads it as the option's default.
 *
 * The estimate of a pair moves in steps of one sketch hash, and one hash stands
 * for the smaller genome's k-mers divided by s: on the 120 E. coli genomes of the
 * validation index 2 507 k-mers at 2 000 (2 264 to 2 840 over the pairs) and 251
 * at 20 000, where the finest threshold of that index's tree is 10 000. A larger
 * sketch costs time in the tree step: on those genomes the sketch pass takes 35
 * to 37 ms per genome at 2 000 and 56 to 62 at 20 000 on one thread, the pair
 * counting 0.02 s and 0.21 to 0.24 s; a sketch is 8 bytes a hash.
 *
 * A reference with fewer distinct k-mers than s is sketched whole
 * (clusterMinhashSketch returns every hash it has), and its estimates then mean
 * something else. Between two full sketches the estimate follows the k-mers the
 * larger genome has and the smaller lacks; with both sketched whole it is the
 * k-mers the smaller has and the larger lacks, and a value in between with one of
 * them whole. A 3 000-base genome that is part of a 5 000-base one (k = 21,
 * testMinhashSketchSize) is estimated at about 1 200 k-mers from it at 2 000, at
 * 600 at 4 000 and at 0 at 10 000, where the two merge under any threshold. So
 * the size should stay below the distinct k-mers of the smallest reference; at
 * 2 000 that concerns references of less than 2 kb. The two functions say how
 * many references a size other than this one sketches whole, and nothing at
 * this one.
 */
static constexpr unsigned MINHASH_SKETCH_SIZE = 2000;

/**
 * Streaming bottom-s MinHash sketch of a reference (O(s) memory), plus its exact
 * k-mer count. Uses the same ntHash values as the rest of the pipeline.
 *
 * @param fPath  Path to a reference file (plain or gzipped).
 * @param kmerLen  K-mer length.
 * @param s  Sketch size (number of smallest distinct hashes to keep).
 * @param kmerCount  Output: total number of k-mer positions in the reference.
 * @return Sorted vector of the s smallest distinct k-mer hashes (every distinct
 *         hash of a reference that has fewer than s).
 */
inline std::vector<uint64_t> clusterMinhashSketch(const std::string& fPath,
                                                  unsigned kmerLen, unsigned s,
                                                  size_t& kmerCount) {
    std::priority_queue<uint64_t> heap;    // max-heap of the s smallest hashes seen
    std::unordered_set<uint64_t> inHeap;   // current heap contents (<= s), for dedup
    // Counted in a local and stored once at the end. `kmerCount` is the caller's
    // kmerCount[i], eight to a cache line, and the callers run one reference per
    // thread: incremented in place, once per k-mer, every increment took the line
    // away from the threads working on the neighbouring references.
    size_t counted = 0;
    gzFile fp = gzopen(fPath.c_str(), "r");
    kseq_t* ks = kseq_init(fp);
    while (kseq_read(ks) >= 0) {
        std::string seq(ks->seq.s, ks->seq.l);
        ntHashIterator itr(seq, 1, kmerLen);
        while (itr != itr.end()) {
            uint64_t h = (*itr)[0];
            ++counted;
            if (heap.size() < s) {
                if (inHeap.insert(h).second) heap.push(h);
            } else if (h < heap.top() && inHeap.find(h) == inHeap.end()) {
                inHeap.erase(heap.top());
                heap.pop();
                heap.push(h);
                inHeap.insert(h);
            }
            ++itr;
        }
    }
    kseq_destroy(ks);
    gzclose(fp);
    kmerCount = counted;
    std::vector<uint64_t> sketch(inHeap.begin(), inHeap.end());
    std::sort(sketch.begin(), sketch.end());
    return sketch;
}

/**
 * Shared sketch hashes per candidate pair: for every pair of references that share
 * at least one MinHash sketch hash, how many they share, keyed (i << 32 | j) with
 * i < j. The merges read the map in its iteration order and name nodes in the order
 * the merges produce, so the map must come out as the original loop built it: one
 * operator[] per (hash, pair) visit, over the inverted index in its order. libstdc++
 * lays an unordered_map out by its sequence of distinct-key insertions (the rehash
 * points and the within-bucket order follow from it), so the same sequence of first
 * insertions with the same final counts gives the same iteration order.
 *
 * On 5 493 E. coli the sketches share 136 k distinct hashes, 117 of them in nearly
 * every genome, and the visits number 10.4 billion: as map increments, 458 s of a
 * 500 s tree construction. Here the counts go into a dense n x n array, rows
 * partitioned over threads so no two threads touch a cell (each thread walks the
 * whole index and takes the pairs whose smaller member is one of its rows; a hash's
 * reference list is ascending, so the inner loop walks a row in order), and the
 * walk position at which a pair is first seen is kept beside the count. The pairs
 * sorted by (first position, i, j) are the original insertion sequence: within one
 * hash the original loop ran a then b, i.e. (i, j) ascending. The dense arrays are
 * 8 bytes per ordered pair (240 MB at 5 493 genomes); above `denseMaxRefs`
 * references the original loop runs instead.
 *
 * That loop does not visit a hash held by more than `cap` references pair by pair
 * (C(f, 2) map increments for a hash in f references). Such a hash used to be
 * skipped, which took it out of every pair's count. The hashes most genomes of a
 * species share are 6 % of an E. coli sketch at nine references in ten and 40 % at
 * one in three, so above 20 000 references every distance estimate rose by that
 * share of a genome and the tree lost its levels without a word: with the cap at
 * 950 of 1 000 E. coli the two finest levels are empty, at 500 the tree has five
 * nodes of two or three genomes, all at the coarsest level, in place of its 195.
 * They are counted now, through one bit per reference and popular hash: a pair's
 * AND gives how many of them it shares and, the bits being in walk order, the
 * first. The result is the uncapped loop's map in the uncapped loop's order
 * whatever `cap` and `denseMaxRefs` are, so both are costs only. No hash can
 * exceed the cap while the references do not, and such a set takes the path it
 * took before.
 */
inline std::unordered_map<uint64_t, uint32_t> countSharedSketchHashes(
        const std::vector<std::vector<uint64_t>>& sketches, size_t cap, unsigned denseMaxRefs = 20000) {
    const size_t n = sketches.size();
    std::unordered_map<uint64_t, std::vector<unsigned>> index;
    for (unsigned i = 0; i < n; i++)
        for (uint64_t h : sketches[i]) index[h].push_back(i);
    std::unordered_map<uint64_t, uint32_t> shared;
    if (n > denseMaxRefs) {
        // The hashes over the cap, with their positions in the walk (the hashes of
        // two or more references in the index's order, as the dense path numbers them).
        std::vector<const std::vector<unsigned>*> popular;
        std::vector<uint32_t> popularPos;
        uint32_t walkLen = 0;
        for (auto& kv : index) {
            if (kv.second.size() < 2) continue;
            if (kv.second.size() > cap) { popular.push_back(&kv.second); popularPos.push_back(walkLen); }
            ++walkLen;
        }
        if (popular.empty()) {
            for (auto& kv : index) {
                auto& refs = kv.second;
                if (refs.size() < 2) continue;
                for (size_t a = 0; a < refs.size(); a++)
                    for (size_t b = a + 1; b < refs.size(); b++) {
                        unsigned i = refs[a], j = refs[b];
                        if (i > j) std::swap(i, j);
                        shared[(static_cast<uint64_t>(i) << 32) | j]++;
                    }
            }
            return shared;
        }
        // The other hashes pair by pair as above, keeping where each pair was first
        // seen: the popular hashes can put a pair earlier in the walk, and the map
        // must receive the pairs in the order the uncapped loop would have.
        struct Cell { uint32_t count = 0, first = 0; };
        std::unordered_map<uint64_t, Cell> cells;
        uint32_t pos = 0;
        for (auto& kv : index) {
            auto& refs = kv.second;
            if (refs.size() < 2) continue;
            if (refs.size() <= cap)
                for (size_t a = 0; a < refs.size(); a++)
                    for (size_t b = a + 1; b < refs.size(); b++) {
                        Cell& c = cells[(static_cast<uint64_t>(refs[a]) << 32) | refs[b]];
                        if (c.count++ == 0) c.first = pos;
                    }
            ++pos;
        }
        // One bit per (reference, popular hash), bits in walk order.
        const size_t words = (popular.size() + 63) / 64;
        std::vector<uint64_t> bits(n * words, 0);
        for (size_t p = 0; p < popular.size(); p++)
            for (unsigned r : *popular[p]) bits[r * words + p / 64] |= uint64_t(1) << (p % 64);
        std::vector<unsigned> holders;   // references holding a popular hash, ascending
        for (unsigned i = 0; i < n; i++)
            for (size_t w = 0; w < words; w++)
                if (bits[i * words + w]) { holders.push_back(i); break; }
        for (size_t x = 0; x < holders.size(); x++) {
            const uint64_t* bx = &bits[holders[x] * words];
            for (size_t y = x + 1; y < holders.size(); y++) {
                const uint64_t* by = &bits[holders[y] * words];
                uint32_t common = 0, firstPos = 0;
                for (size_t w = 0; w < words; w++) {
                    const uint64_t both = bx[w] & by[w];
                    if (!both) continue;
                    if (common == 0) firstPos = popularPos[w * 64 + __builtin_ctzll(both)];
                    common += static_cast<uint32_t>(__builtin_popcountll(both));
                }
                if (common == 0) continue;
                Cell& c = cells[(static_cast<uint64_t>(holders[x]) << 32) | holders[y]];
                if (c.count == 0 || firstPos < c.first) c.first = firstPos;
                c.count += common;
            }
        }
        // The pairs in the uncapped loop's insertion order (first position, i, j).
        struct Entry { uint32_t pos, i, j, count; };
        std::vector<Entry> entries;
        entries.reserve(cells.size());
        for (auto& kv : cells)
            entries.push_back({kv.second.first, static_cast<uint32_t>(kv.first >> 32),
                               static_cast<uint32_t>(kv.first & 0xffffffffu), kv.second.count});
        cells = {};
        std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
            return a.pos != b.pos ? a.pos < b.pos : a.i != b.i ? a.i < b.i : a.j < b.j;
        });
        for (const Entry& e : entries)
            shared.emplace((static_cast<uint64_t>(e.i) << 32) | e.j, e.count);
        return shared;
    }
    // The hashes in the index's order, as the original loop visited them. Every
    // hash is counted here whatever the cap: a visit is an array increment.
    std::vector<const std::vector<unsigned>*> walk;
    walk.reserve(index.size());
    for (auto& kv : index)
        if (kv.second.size() >= 2) walk.push_back(&kv.second);
    std::vector<uint32_t> count(n * n, 0), first(n * n, UINT32_MAX);
    unsigned threads = 1;
#ifdef _OPENMP
    threads = static_cast<unsigned>(std::max(1, omp_get_max_threads()));
#endif
    #pragma omp parallel num_threads(static_cast<int>(threads))
    {
        unsigned t = 0, T = 1;
#ifdef _OPENMP
        t = static_cast<unsigned>(omp_get_thread_num()); T = static_cast<unsigned>(omp_get_num_threads());
#endif
        for (size_t pos = 0; pos < walk.size(); pos++) {
            const std::vector<unsigned>& refs = *walk[pos];
            for (size_t a = 0; a + 1 < refs.size(); a++) {
                const unsigned i = refs[a];
                if (i % T != t) continue;
                uint32_t* c = &count[static_cast<size_t>(i) * n];
                uint32_t* f = &first[static_cast<size_t>(i) * n];
                for (size_t b = a + 1; b < refs.size(); b++) {
                    const unsigned j = refs[b];
                    if (c[j]++ == 0) f[j] = static_cast<uint32_t>(pos);
                }
            }
        }
    }
    // The pairs in insertion order.
    struct Seen { uint32_t pos, i, j; };
    std::vector<Seen> seen;
    for (size_t i = 0; i < n; i++) {
        const uint32_t* f = &first[i * n];
        for (size_t j = i + 1; j < n; j++)
            if (f[j] != UINT32_MAX) seen.push_back({f[j], static_cast<uint32_t>(i), static_cast<uint32_t>(j)});
    }
    first = {};
    std::sort(seen.begin(), seen.end(), [](const Seen& a, const Seen& b) {
        return a.pos != b.pos ? a.pos < b.pos : a.i != b.i ? a.i < b.i : a.j < b.j;
    });
    for (const Seen& s : seen)
        shared.emplace((static_cast<uint64_t>(s.i) << 32) | s.j, count[static_cast<size_t>(s.i) * n + s.j]);
    return shared;
}

/**
 * Near-linear clustering via MinHash + inverted index + single-linkage union-find.
 * Instead of scoring all O(n^2) pairs, references are sketched (bottom-s MinHash)
 * and only pairs that share sketch hashes are scored, using an inverted index. The
 * estimated number of distinguishing k-mers is min(|i|,|j|)*(1 - Jaccard); pairs
 * below `threshold` are merged, matching the semantics of --cluster=N. Approximate
 * (sketch resolution ~ |genome|/s), but on well-separated near-identical strains it
 * reproduces the exact clustering at a fraction of the cost. Use clusterReferencesExact
 * for a precise result.
 *
 * @param refFiles  All reference file paths.
 * @param threshold  Distinguishing-k-mer threshold below which references merge.
 * @param kmerLen  K-mer length.
 * @param clusterFile  Output TSV path mapping each reference to its representative.
 * @param sketchSize  Hashes in a reference's sketch (--minhash-size; see MINHASH_SKETCH_SIZE).
 * @return Vector of representative file paths (one per cluster).
 */
inline std::vector<std::string> clusterReferencesMinHash(
        const std::vector<std::string>& refFiles,
        unsigned threshold,
        unsigned kmerLen,
        const std::string& clusterFile,
        unsigned sketchSize = MINHASH_SKETCH_SIZE) {
    size_t n = refFiles.size();
    if (n <= 1) return refFiles;

    const unsigned s = sketchSize;
    std::cerr << "Clustering " << n << " references (MinHash, threshold=" << threshold
              << " unique k-mers, k=" << kmerLen << ", sketch=" << s << ")...\n";

    std::vector<std::vector<uint64_t>> sketches(n);
    std::vector<size_t> kmerCount(n, 0);
    std::vector<std::string> names(n);

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < n; i++) {
        sketches[i] = clusterMinhashSketch(refFiles[i], kmerLen, s, kmerCount[i]);
        names[i] = getBaseId(refFiles[i]);
    }
    // A reference with fewer distinct k-mers than the sketch is sketched whole and its
    // estimates mean something else (see MINHASH_SKETCH_SIZE). Said for another size
    // than the default only, so that the stderr of a default build is what it was.
    if (s != MINHASH_SKETCH_SIZE) {
        size_t whole = 0;
        for (const std::vector<uint64_t>& sketch : sketches) whole += sketch.size() < s;
        if (whole)
            std::cerr << "Warning: " << whole << " of " << n << " references have fewer than " << s
                      << " distinct k-mers and are sketched whole (--minhash-size); the distance "
                         "estimates of their pairs mean something else (see README)\n";
    }

    // Shared-hash counts per candidate pair (countSharedSketchHashes). Above 20 000
    // references, hashes shared by more than CAP of them are counted through bit sets
    // rather than pair by pair; the counts are the same.
    const size_t CAP = 20000;
    std::unordered_map<uint64_t, uint32_t> shared = countSharedSketchHashes(sketches, CAP);

    // Single-linkage union-find.
    std::vector<unsigned> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](unsigned x) -> unsigned {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };

    for (auto& kv : shared) {
        unsigned i = static_cast<unsigned>(kv.first >> 32);
        unsigned j = static_cast<unsigned>(kv.first & 0xffffffffu);
        size_t denom = std::min(sketches[i].size(), sketches[j].size());
        if (denom == 0) continue;
        double jaccard = static_cast<double>(kv.second) / static_cast<double>(denom);
        double uniqueEst = static_cast<double>(std::min(kmerCount[i], kmerCount[j])) * (1.0 - jaccard);
        if (uniqueEst < static_cast<double>(threshold)) {
            unsigned ri = find(i), rj = find(j);
            if (ri != rj) parent[ri] = rj;
        }
    }

    // Build clusters.
    std::unordered_map<unsigned, std::vector<unsigned>> clusters;
    for (unsigned i = 0; i < n; i++) clusters[find(i)].push_back(i);

    // Representative per cluster: member with the most k-mers.
    std::vector<std::string> representatives;
    std::ofstream clusterOut(clusterFile);
    clusterOut << "cluster_id\trepresentative\tmember\tnum_kmers\n";
    unsigned clusterId = 0;
    for (auto& [root, members] : clusters) {
        unsigned repIdx = members[0];
        for (unsigned m : members)
            if (kmerCount[m] > kmerCount[repIdx]) repIdx = m;
        representatives.push_back(refFiles[repIdx]);
        for (unsigned m : members)
            clusterOut << clusterId << "\t" << names[repIdx] << "\t" << names[m]
                       << "\t" << kmerCount[m] << "\n";
        ++clusterId;
    }
    clusterOut.close();

    std::cerr << "Clustering complete: " << n << " references -> "
              << representatives.size() << " clusters\n";
    return representatives;
}

/**
 * Cluster references, dispatching to the near-linear MinHash method by default or
 * the exact O(n^2) method when useExact is set.
 */
inline std::vector<std::string> clusterReferences(
        const std::vector<std::string>& refFiles,
        unsigned threshold,
        unsigned kmerLen,
        unsigned nhash,
        unsigned bits,
        size_t dbfSize,
        const std::string& clusterFile,
        bool useExact = false,
        unsigned sketchSize = MINHASH_SKETCH_SIZE) {
    if (useExact)
        return clusterReferencesExact(refFiles, threshold, kmerLen, nhash, bits, dbfSize, clusterFile);
    return clusterReferencesMinHash(refFiles, threshold, kmerLen, clusterFile, sketchSize);
}

#endif // CLUSTERUTIL_HPP_