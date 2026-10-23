#ifndef COMPARESKETCHUTIL_HPP_
#define COMPARESKETCHUTIL_HPP_

#include <algorithm>
#include <climits>
#include <cstring>

#include "SequenceUtil.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"
#include "FileUtil.hpp"
#include "kseq_gz.h"

/**
 * Settings for one pairwise comparison run. Passed explicitly; see
 * UniqSketchConfig for why these are no longer file-scope globals.
 */
struct CompareSketchConfig {
    std::string outfile = "reference_similarity.tsv";   // output similarity filename
    unsigned threads = 1;                               // number of threads
    unsigned kmerLen = 81;                              // k-mer length

    // Bloom filter sizing
    unsigned bits = 16;                                 // bits per element
    unsigned nhash = 3;                                 // hash functions
    size_t dbfSize = 5000000;                           // element count (approx. largest genome size)
    bool autoSize = false;                              // --auto: size from an ntCard cardinality estimate
    bool autoGsize = false;                             // --auto-gsize: size from the largest input genome
    double targetFpr = 0.0;                             // --fpr: target false-positive rate (0 = use bits)

    // Memory
    size_t memBudgetGB = 64;                            // budget for caching refSet2 hashes; 0 = stream
};

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
 * Load k-mers from a FASTA reference directly into a Bloom filter, reading the
 * file once and not retaining the hashes. Used by the unit tests.
 */
void loadBloomFilter(const CompareSketchConfig& cfg, const std::string& fPath, BloomFilter& dbFilter) {
    auto hashes = loadKmerHashes(fPath, cfg.nhash, cfg.kmerLen);
    loadBloomFilterFromHashes(hashes, cfg.nhash, dbFilter);
}

/**
 * Check a FASTA reference's k-mers against a Bloom filter, reading the file once
 * and not retaining the hashes. Used by the unit tests.
 */
double checkBloomFilter(const CompareSketchConfig& cfg, const std::string& fPath, const BloomFilter& dbFilter,
                        size_t& countAll, size_t& countUnique) {
    auto hashes = loadKmerHashes(fPath, cfg.nhash, cfg.kmerLen);
    return checkBloomFilterFromHashes(hashes, cfg.nhash, dbFilter, countAll, countUnique);
}

/**
 * Pairwise similarity between two reference sets (M x N), under a memory budget.
 *
 * refSet2 is cached in batches sized to fit cfg.memBudgetGB, and refSet1 is
 * streamed (loaded per entry, per batch). This one function spans the whole
 * memory/speed range:
 *   - a budget large enough to hold all of refSet2 -> a single batch, so every
 *     file is read once (M+N reads: fastest, highest memory);
 *   - cfg.memBudgetGB == 0 (or a genome bigger than the budget) -> batch size 1,
 *     i.e. one refSet2 genome cached at a time -> M x N reads (lowest memory).
 * The --low-mem flag is simply an alias for memBudgetGB == 0. Similarity values
 * are identical for any budget; only the memory/I/O profile changes.
 *
 * The batch size is derived from a per-genome hash-memory estimate
 * (largest input genome x nhash x 8 bytes), independent of the Bloom-filter
 * sizing, so --auto (which sets dbfSize to the whole-set cardinality) does not
 * distort it.
 *
 * @param refSet1  List of references in the first set.
 * @param refSet2  List of references in the second set.
 */
void identifyDifference(const CompareSketchConfig& cfg, const std::vector<std::string>& refSet1,
                        const std::vector<std::string>& refSet2) {
    std::vector<std::string> names1(refSet1.size());
    std::vector<std::string> names2(refSet2.size());
    for (unsigned i = 0; i < refSet1.size(); i++) names1[i] = getBaseId(refSet1[i]);
    for (unsigned j = 0; j < refSet2.size(); j++) names2[j] = getBaseId(refSet2[j]);

    std::ofstream simOut(cfg.outfile);
    simOut << "ref1\tref2\tuniq_ref1\tall_ref1\tsimilarity\n";

    size_t bfBits = cfg.bits * cfg.dbfSize;

    // Per-genome hash-memory estimate (bases ~ k-mers), independent of BF sizing.
    size_t perGenome = std::max<size_t>(1, estimateMaxGenomeSize(refSet2));
    size_t hashPerGenome = perGenome * cfg.nhash * sizeof(uint64_t);
    size_t memBudget = cfg.memBudgetGB * 1024ULL * 1024 * 1024;
    size_t batchSize = std::max<size_t>(1, memBudget / hashPerGenome);
    batchSize = std::min(batchSize, refSet2.size());

    std::cerr << "comparesketch: " << refSet1.size() << " x " << refSet2.size()
              << " pairwise, refSet2 batch=" << batchSize << "/" << refSet2.size()
              << " (budget " << cfg.memBudgetGB << " GB)\n";

    // Process refSet2 in batches that fit the budget; stream refSet1 within each.
    for (size_t bStart = 0; bStart < refSet2.size(); bStart += batchSize) {
        size_t bEnd = std::min(bStart + batchSize, refSet2.size());
        size_t bCount = bEnd - bStart;

        std::vector<std::vector<uint64_t>> batchHashes(bCount);
        #pragma omp parallel for schedule(dynamic)
        for (size_t b = 0; b < bCount; b++)
            batchHashes[b] = loadKmerHashes(refSet2[bStart + b], cfg.nhash, cfg.kmerLen);

        #pragma omp parallel
        {
            BloomFilter dbFilter(bfBits, cfg.nhash, cfg.kmerLen);
            unsigned char* filterPtr = dbFilter.getFilter();
            size_t filterBytes = dbFilter.sizeBytes();   // the allocation: the blocked filter rounds up

            #pragma omp for schedule(dynamic)
            for (unsigned i = 0; i < refSet1.size(); i++) {
                std::memset(filterPtr, 0, filterBytes);
                auto hashes1 = loadKmerHashes(refSet1[i], cfg.nhash, cfg.kmerLen);
                loadBloomFilterFromHashes(hashes1, cfg.nhash, dbFilter);
                hashes1.clear();
                hashes1.shrink_to_fit();

                for (size_t b = 0; b < bCount; b++) {
                    size_t countAll = 0, countUnique = 0;
                    double score = checkBloomFilterFromHashes(
                        batchHashes[b], cfg.nhash, dbFilter, countAll, countUnique);
                    #pragma omp critical(simout)
                    simOut << names2[bStart + b] << "\t" << names1[i]
                           << "\t" << countUnique << "\t" << countAll
                           << "\t" << score << "\n";
                }
            }
        }
    }
    simOut.close();
}

#endif // COMPARESKETCHUTIL_HPP_
