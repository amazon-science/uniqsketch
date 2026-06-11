#ifndef COMPARESKETCHUTIL_HPP_
#define COMPARESKETCHUTIL_HPP_

#include <climits>
#include <cstring>

#include "SequenceUtil.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"
#include "kseq_gz.h"

namespace opt {
unsigned threads(1);                                    // number of threads
unsigned bits(16);                                      // bits per element in Bloom filter
unsigned nhash(3);                                      // number of hashes for Bloom filter
unsigned kmerLen(81);                                   // length of k-mer
size_t dbfSize(5000000);                                // distinct BF size (approx genome size)
std::string outfile("reference_similarity.tsv");        // output file name
bool autoSize(false);                                   // size BF from ntCard cardinality estimate
bool autoGsize(false);                                  // size BF from largest input genome
double targetFpr(0.0);                                  // target Bloom-filter FPR (0 = use --bit)
}

/**
 * Load all k-mer hashes from a FASTA/FASTA.gz reference into a vector.
 * Each entry stores nhash hash values for one k-mer.
 *
 * @param fPath  Path to a reference file (plain or gzipped).
 * @param nhash  Number of hash functions.
 * @param kmerLen  K-mer length.
 * @return Vector of hash value arrays (flattened: nhash values per k-mer).
 */
std::vector<uint64_t> loadKmerHashes(const std::string& fPath,
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
 * Insert precomputed k-mer hashes into a Bloom filter.
 *
 * @param hashes  Flattened hash array (nhash values per k-mer).
 * @param nhash  Number of hash functions.
 * @param dbFilter  Bloom filter to populate.
 */
void loadBloomFilterFromHashes(const std::vector<uint64_t>& hashes,
                               unsigned nhash, BloomFilter& dbFilter) {
    for (size_t i = 0; i < hashes.size(); i += nhash) {
        dbFilter.insert(&hashes[i]);
    }
}

/**
 * Check precomputed k-mer hashes against a Bloom filter.
 *
 * @param hashes  Flattened hash array (nhash values per k-mer).
 * @param nhash  Number of hash functions.
 * @param dbFilter  Bloom filter to query against.
 * @param countAll  Total number of k-mers (output).
 * @param countUnique  Number of k-mers not found in BF (output).
 * @return Similarity rate (fraction of k-mers found in BF).
 */
double checkBloomFilterFromHashes(const std::vector<uint64_t>& hashes,
                                  unsigned nhash, const BloomFilter& dbFilter,
                                  size_t& countAll, size_t& countUnique) {
    countAll = hashes.size() / nhash;
    countUnique = 0;
    for (size_t i = 0; i < hashes.size(); i += nhash) {
        if (!dbFilter.contains(&hashes[i])) {
            ++countUnique;
        }
    }
    if (countAll == 0) return 0.0;
    return static_cast<double>(countAll - countUnique) / countAll;
}

/**
 * Legacy interface: load k-mers from FASTA into a Bloom filter (for tests).
 */
void loadBloomFilter(const std::string& fPath, BloomFilter& dbFilter) {
    auto hashes = loadKmerHashes(fPath, opt::nhash, opt::kmerLen);
    loadBloomFilterFromHashes(hashes, opt::nhash, dbFilter);
}

/**
 * Legacy interface: check FASTA k-mers against a Bloom filter (for tests).
 */
double checkBloomFilter(const std::string& fPath, const BloomFilter& dbFilter,
                        size_t& countAll, size_t& countUnique) {
    auto hashes = loadKmerHashes(fPath, opt::nhash, opt::kmerLen);
    return checkBloomFilterFromHashes(hashes, opt::nhash, dbFilter, countAll, countUnique);
}

/**
 * Calculate pairwise similarity between two sets of references (M x N comparison).
 * Pre-caches all k-mer hashes to avoid redundant file I/O.
 *
 * @param refSet1  List of references in the first set.
 * @param refSet2  List of references in the second set.
 */
void identifyDifference(const std::vector<std::string>& refSet1,
                        const std::vector<std::string>& refSet2) {
    // Pre-compute hashes for refSet1 and refSet2 in parallel
    std::vector<std::vector<uint64_t>> hashSet1(refSet1.size());
    std::vector<std::vector<uint64_t>> hashSet2(refSet2.size());
    std::vector<std::string> names1(refSet1.size());
    std::vector<std::string> names2(refSet2.size());

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < refSet1.size(); i++) {
        hashSet1[i] = loadKmerHashes(refSet1[i], opt::nhash, opt::kmerLen);
        names1[i] = getBaseId(refSet1[i]);
    }

    #pragma omp parallel for schedule(dynamic)
    for (unsigned j = 0; j < refSet2.size(); j++) {
        hashSet2[j] = loadKmerHashes(refSet2[j], opt::nhash, opt::kmerLen);
        names2[j] = getBaseId(refSet2[j]);
    }

    std::ofstream simOut(opt::outfile);
    simOut << "ref1\tref2\tuniq_ref1\tall_ref1\tsimilarity\n";

    size_t bfBytes = opt::bits * opt::dbfSize;

    // Pairwise comparison: parallelize over refSet1
    #pragma omp parallel
    {
        // Each thread gets its own BF — allocated once, cleared per refSet1 entry
        BloomFilter dbFilter(bfBytes, opt::nhash, opt::kmerLen);
        unsigned char* filterPtr = dbFilter.getFilter();
        size_t filterBytes = (bfBytes + CHAR_BIT - 1) / CHAR_BIT;

        #pragma omp for schedule(dynamic)
        for (unsigned i = 0; i < refSet1.size(); i++) {
            // Clear and reload BF for this reference
            std::memset(filterPtr, 0, filterBytes);
            loadBloomFilterFromHashes(hashSet1[i], opt::nhash, dbFilter);

            for (unsigned j = 0; j < refSet2.size(); j++) {
                size_t countAll = 0, countUnique = 0;
                double score = checkBloomFilterFromHashes(
                    hashSet2[j], opt::nhash, dbFilter, countAll, countUnique);
                #pragma omp critical(simout)
                simOut << names2[j] << "\t" << names1[i]
                       << "\t" << countUnique << "\t" << countAll
                       << "\t" << score << "\n";
            }
        }
    }
    simOut.close();
}

#endif // COMPARESKETCHUTIL_HPP_
