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
 * identifyDifference no longer goes through this function and the two after it
 * (see loadLeadingHashes); with loadBloomFilter and checkBloomFilter they are
 * what the unit tests hold its counts to.
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
 * Sort hash values ascending: an LSD radix sort, one byte a pass. The 5 M k-mer
 * positions of an E. coli genome take it 0.12 s (0.15 s when eight threads sort
 * at once), std::sort 0.33 s (0.34 s).
 *
 * @param hashes  Values to sort, in place.
 * @param scratch  The sort's second array, resized as needed; a caller that sorts
 *                 many vectors passes the same one so that it is allocated once.
 */
void sortHashes(std::vector<uint64_t>& hashes, std::vector<uint64_t>& scratch) {
    const size_t n = hashes.size();
    if (n < 2) return;
    scratch.resize(n);
    size_t count[8][256] = {};
    for (size_t i = 0; i < n; i++)
        for (unsigned d = 0; d < 8; d++) count[d][(hashes[i] >> (8 * d)) & 255]++;
    uint64_t* from = hashes.data();
    uint64_t* to = scratch.data();
    for (unsigned d = 0; d < 8; d++) {
        size_t next = 0;
        for (unsigned v = 0; v < 256; v++) {
            const size_t c = count[d][v];
            count[d][v] = next;
            next += c;
        }
        for (size_t i = 0; i < n; i++) to[count[d][(from[i] >> (8 * d)) & 255]++] = from[i];
        std::swap(from, to);
    }
    // Eight passes: the sorted values are back in `hashes`.
}

/**
 * The leading hash of every k-mer position of a reference, in file order: hVal[0]
 * of ntHashIterator, the canonical ntHash, one value per position, so a k-mer that
 * occurs m times in the reference is m values. This is what the comparison caches
 * of a refSet2 genome. The other nhash - 1 values of a k-mer are functions of the
 * leading one (NTE64) and are derived when the filter is probed
 * (BloomFilter::countAbsent), so a position takes 8 bytes whatever nhash is.
 *
 * @param fPath  Path to a reference file (plain or gzipped).
 * @param kmerLen  K-mer length.
 * @return Leading hash values, one per k-mer position.
 */
std::vector<uint64_t> loadLeadingHashes(const std::string& fPath, unsigned kmerLen) {
    std::vector<uint64_t> hashes;
    // The file's bytes, an upper bound of its k-mer positions. Grown by doubling
    // instead, the vector of a 5 Mb genome took 0.18 s to fill where it takes
    // 0.06 s, with eight threads loading at once (0.08 and 0.05 s with one).
    hashes.reserve(uncompressedSize(fPath));
    gzFile fp = gzopen(fPath.c_str(), "r");
    kseq_t* ks = kseq_init(fp);
    while (kseq_read(ks) >= 0) {
        std::string seq(ks->seq.s, ks->seq.l);
        for (ntHashIterator itr(seq, 1, kmerLen); itr != itr.end(); ++itr) {
            hashes.push_back((*itr)[0]);
        }
    }
    kseq_destroy(ks);
    gzclose(fp);
    return hashes;
}

/**
 * Insert every k-mer of a reference into a Bloom filter as the file is read; no
 * hash is kept. The filter gets the bits loadBloomFilter sets. It must be the
 * calling thread's own: the inserts are plain stores (insert_exclusive).
 *
 * @param fPath  Path to a reference file (plain or gzipped).
 * @param dbFilter  Bloom filter to populate; its hash number and k-mer size are used.
 */
void insertReference(const std::string& fPath, BloomFilter& dbFilter) {
    gzFile fp = gzopen(fPath.c_str(), "r");
    kseq_t* ks = kseq_init(fp);
    while (kseq_read(ks) >= 0) {
        std::string seq(ks->seq.s, ks->seq.l);
        for (ntHashIterator itr(seq, dbFilter.getHashNum(), dbFilter.getKmerSize()); itr != itr.end(); ++itr) {
            dbFilter.insert_exclusive(*itr);
        }
    }
    kseq_destroy(ks);
    gzclose(fp);
}

/**
 * Check a reference's k-mers, given by their leading hashes (loadLeadingHashes
 * with the filter's k-mer size, in any order), against a Bloom filter. The numbers
 * are those of checkBloomFilterFromHashes on the same reference.
 *
 * @param leading  Leading hash value of every k-mer position.
 * @param dbFilter  Bloom filter to query against.
 * @param countAll  Total number of k-mers (output).
 * @param countUnique  Number of k-mers not found in BF (output).
 * @return Similarity rate (fraction of k-mers found in BF).
 */
double checkBloomFilterFromLeading(const std::vector<uint64_t>& leading, const BloomFilter& dbFilter,
                                   size_t& countAll, size_t& countUnique) {
    countAll = leading.size();
    countUnique = dbFilter.countAbsent(leading.data(), leading.size());
    if (countAll == 0) return 0.0;
    return static_cast<double>(countAll - countUnique) / countAll;
}

/**
 * The refSet2 genomes identifyDifference caches at a time: as many as fit the
 * budget, at least one, at most all of them. A genome counts 8 bytes per byte
 * of the largest refSet2 file: one hash value per k-mer position, and a file
 * holds no more positions than bytes. The hash number does not enter: while all
 * hash values of a position were kept a genome counted 8 bytes for each hash
 * function, and the default 64 GB held 572 genomes of 5 000 000 bytes where it
 * holds 1 717.
 *
 * @param memBudgetGB  Budget for the cached hashes, in GB; 0 = one genome at a time.
 * @param largestGenome  Bytes of the largest refSet2 file (estimateMaxGenomeSize).
 * @param genomes  Number of refSet2 genomes.
 * @return Genomes per batch.
 */
size_t refSet2BatchSize(size_t memBudgetGB, size_t largestGenome, size_t genomes) {
    const size_t hashPerGenome = std::max<size_t>(1, largestGenome) * sizeof(uint64_t);
    const size_t memBudget = memBudgetGB * 1024ULL * 1024 * 1024;
    return std::min(std::max<size_t>(1, memBudget / hashPerGenome), genomes);
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
 * What is cached of a refSet2 genome is the leading hash of each k-mer position
 * (loadLeadingHashes), 8 bytes a position, and a genome is checked with one call
 * of BloomFilter::countAbsent; a refSet1 genome goes into its thread's filter as
 * its file is read (insertReference). Before, both were held as all nhash values
 * of every position, 24 bytes a position at the default, and each position was
 * checked with BloomFilter::contains.
 *
 * The batch size is derived from a per-genome hash-memory estimate
 * (largest input genome x 8 bytes), independent of the Bloom-filter
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
    const size_t batchSize = refSet2BatchSize(cfg.memBudgetGB, estimateMaxGenomeSize(refSet2), refSet2.size());

    std::cerr << "comparesketch: " << refSet1.size() << " x " << refSet2.size()
              << " pairwise, refSet2 batch=" << batchSize << "/" << refSet2.size()
              << " (budget " << cfg.memBudgetGB << " GB)\n";

    // Ascending leading hashes make a genome's probes walk the blocked filter front
    // to back, its block being chosen by the leading hash. Sorting takes 0.12 s a
    // genome and is made from two refSet1 genomes on. With the default filter
    // (10 MB) alone in a 32 MB cache, on one thread, it saves 15.6 ms a pair
    // (12.1 ms for the k-mers of a 5 Mb genome where the file's order takes
    // 27.7 ms, BlockedBloomFilter::countAbsent) and pays only from eight genomes
    // on: 2 x 16 E. coli take 3.4 s sorted and 2.0 s in file order, 4 x 16 4.0
    // and 3.1 s, 7 x 16 4.9 and 4.7 s. With one hash function, where one value a
    // position was stored before as well, the sorted run is slower than the run
    // before below seven genomes (-d 1: 2 x 16 take 3.2 s for 2.1 s, 4 x 16 3.5 s
    // for 2.9 s, 7 x 16 4.0 s either way). A filter that is not in the cache makes
    // the probes in the file's order cost more than those on the stored values
    // did before, and the sort pays from two genomes on. One thread with 80 MB
    // (-b 128): 2 x 16 take 3.9 s sorted and 5.4 s in file order (6.5 s before),
    // 7 x 16 6.5 and 16.4 s (15.1); with the default filter and eight threads on
    // the eight cores of one 32 MB cache 7 x 120 take 5.6 and 15.2 s (15.7 s
    // before). One refSet1 genome is probed in file order, the sort costing more
    // than one genome's probes save on one thread; with 80 MB (-b 128) and eight
    // threads that run is slower than before (1 x 120: 15.7 s for 14.7 s), with
    // 40 MB (-b 64) it is not (12.1 s for 12.8 s). The standard filter has no such
    // order to offer and keeps the file's.
    const bool sorted = BloomFilter::kPartitioned && refSet1.size() >= 2;

    // Process refSet2 in batches that fit the budget; stream refSet1 within each.
    for (size_t bStart = 0; bStart < refSet2.size(); bStart += batchSize) {
        size_t bEnd = std::min(bStart + batchSize, refSet2.size());
        size_t bCount = bEnd - bStart;

        std::vector<std::vector<uint64_t>> batchHashes(bCount);
        #pragma omp parallel
        {
            std::vector<uint64_t> scratch;                   // the sort's second array, one per thread
            #pragma omp for schedule(dynamic)
            for (size_t b = 0; b < bCount; b++) {
                batchHashes[b] = loadLeadingHashes(refSet2[bStart + b], cfg.kmerLen);
                if (sorted) sortHashes(batchHashes[b], scratch);
            }
        }

        #pragma omp parallel
        {
            BloomFilter dbFilter(bfBits, cfg.nhash, cfg.kmerLen);
            unsigned char* filterPtr = dbFilter.getFilter();
            size_t filterBytes = dbFilter.sizeBytes();   // the allocation: the blocked filter rounds up

            #pragma omp for schedule(dynamic)
            for (unsigned i = 0; i < refSet1.size(); i++) {
                std::memset(filterPtr, 0, filterBytes);
                insertReference(refSet1[i], dbFilter);

                for (size_t b = 0; b < bCount; b++) {
                    size_t countAll = 0, countUnique = 0;
                    double score = checkBloomFilterFromLeading(
                        batchHashes[b], dbFilter, countAll, countUnique);
                    #pragma omp critical(simout)
                    simOut << names2[bStart + b] << "\t" << names1[i]
                           << "\t" << countUnique << "\t" << countAll
                           << "\t" << score << "\n";
                }
            }
        }
    }
    simOut.close();
    if (!simOut) {
        std::cerr << "Error: could not write " << cfg.outfile << "\n";
        exit(EXIT_FAILURE);
    }
}

#endif // COMPARESKETCHUTIL_HPP_
