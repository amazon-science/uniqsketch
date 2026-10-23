#ifndef CLUSTERUTIL_HPP_
#define CLUSTERUTIL_HPP_

#include <algorithm>
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

// Default MinHash sketch size (bottom-s). Larger = finer similarity resolution.
static constexpr unsigned MINHASH_SKETCH_SIZE = 2000;

/**
 * Streaming bottom-s MinHash sketch of a reference (O(s) memory), plus its exact
 * k-mer count. Uses the same ntHash values as the rest of the pipeline.
 *
 * @param fPath  Path to a reference file (plain or gzipped).
 * @param kmerLen  K-mer length.
 * @param s  Sketch size (number of smallest distinct hashes to keep).
 * @param kmerCount  Output: total number of k-mer positions in the reference.
 * @return Sorted vector of the s smallest distinct k-mer hashes.
 */
inline std::vector<uint64_t> clusterMinhashSketch(const std::string& fPath,
                                                  unsigned kmerLen, unsigned s,
                                                  size_t& kmerCount) {
    std::priority_queue<uint64_t> heap;    // max-heap of the s smallest hashes seen
    std::unordered_set<uint64_t> inHeap;   // current heap contents (<= s), for dedup
    kmerCount = 0;
    gzFile fp = gzopen(fPath.c_str(), "r");
    kseq_t* ks = kseq_init(fp);
    while (kseq_read(ks) >= 0) {
        std::string seq(ks->seq.s, ks->seq.l);
        ntHashIterator itr(seq, 1, kmerLen);
        while (itr != itr.end()) {
            uint64_t h = (*itr)[0];
            ++kmerCount;
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
    std::vector<uint64_t> sketch(inHeap.begin(), inHeap.end());
    std::sort(sketch.begin(), sketch.end());
    return sketch;
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
 * @return Vector of representative file paths (one per cluster).
 */
inline std::vector<std::string> clusterReferencesMinHash(
        const std::vector<std::string>& refFiles,
        unsigned threshold,
        unsigned kmerLen,
        const std::string& clusterFile) {
    size_t n = refFiles.size();
    if (n <= 1) return refFiles;

    const unsigned s = MINHASH_SKETCH_SIZE;
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

    // Inverted index: sketch hash -> references that contain it.
    std::unordered_map<uint64_t, std::vector<unsigned>> idx;
    for (unsigned i = 0; i < n; i++)
        for (uint64_t h : sketches[i]) idx[h].push_back(i);

    // Accumulate shared-hash counts per candidate pair. Ultra-popular hashes
    // (shared by more than CAP references) are skipped to keep this near-linear.
    const size_t CAP = 20000;
    std::unordered_map<uint64_t, uint32_t> shared;
    for (auto& kv : idx) {
        auto& refs = kv.second;
        if (refs.size() < 2 || refs.size() > CAP) continue;
        for (size_t a = 0; a < refs.size(); a++)
            for (size_t b = a + 1; b < refs.size(); b++) {
                unsigned i = refs[a], j = refs[b];
                if (i > j) std::swap(i, j);
                shared[(static_cast<uint64_t>(i) << 32) | j]++;
            }
    }

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
        bool useExact = false) {
    if (useExact)
        return clusterReferencesExact(refFiles, threshold, kmerLen, nhash, bits, dbfSize, clusterFile);
    return clusterReferencesMinHash(refFiles, threshold, kmerLen, clusterFile);
}

#endif // CLUSTERUTIL_HPP_