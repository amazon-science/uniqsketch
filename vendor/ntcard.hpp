/*
 * Copyright Hamid Mohamadi.
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License. See the LICENSE accompanying this file
 * for the specific language governing permissions and limitations under
 * the License.
 */

// Adapted from ntCard https://github.com/bcgsc/ntCard

#ifndef NTCARD_HPP_
#define NTCARD_HPP_

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "ntHashIterator.hpp"
#include "kseq_gz.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace opt {
size_t rBuck;
unsigned rBits(27);
unsigned sBits(11);
unsigned sMask;
unsigned covMax(65536);
size_t nSamp(2);
size_t nK(1);
bool samH(true);
}

/**
 * Classify a hash value into a sampling bucket and increment the counter.
 *
 * @param hVal  Hash value to classify.
 * @param t_Counter  Counter array for sampling buckets.
 */
inline void ntComp(uint64_t hVal, uint16_t* t_Counter) {
    uint64_t indBit = opt::nSamp;
    if (hVal >> (64 - opt::sBits) == 1) {
        indBit = 0;
    }
    if (hVal >> (64 - opt::sBits) == opt::sMask) {
        indBit = 1;
    }
    if (indBit < opt::nSamp) {
        size_t shVal = hVal & (opt::rBuck - 1);
        #pragma omp atomic
        ++t_Counter[indBit * opt::rBuck + shVal];
    }
}

/**
 * Hash all k-mers in a sequence and update sampling counters.
 *
 * @param seq  DNA sequence.
 * @param kList  List of k-mer sizes.
 * @param t_Counter  Counter array for sampling buckets.
 * @param totKmer  Per-k total k-mer counts (output).
 */
inline void ntRead(const std::string& seq, const std::vector<unsigned>& kList,
                   uint16_t* t_Counter, size_t totKmer[]) {
    for (unsigned k = 0; k < kList.size(); k++) {
        ntHashIterator itr(seq, 1, kList[k]);
        while (itr != itr.end()) {
            ntComp((*itr)[0], t_Counter + k * opt::nSamp * opt::rBuck);
            ++itr;
            ++totKmer[k];
        }
    }
}

/**
 * Compute cardinality estimates from sampling counters.
 *
 * Estimates F0 (distinct k-mers) and f1 (k-mers seen once), which is all
 * getCardinality reports: the Bloom filters are sized from F0 and F0 - f1.
 * Upstream ntCard continues from here to the whole frequency spectrum
 * f2..f_covMax with an O(covMax^2) recurrence. With covMax = 65536 that took
 * about 1.9 s per call, for output no caller read, on every index build and
 * every query, and once per group in the within-group signature build (60 s of
 * a 100 s phase on 120 E. coli genomes). F0 and f1 are computed exactly as
 * before, so filter sizes, and with them every index and result, are unchanged.
 *
 * @param t_Counter  Counter array from sampling.
 * @param F0Mean  Estimated number of distinct k-mers (output).
 * @param fMean  Frequency estimates (output); only fMean[1] is filled in.
 */
void compEst(const uint16_t* t_Counter, double& F0Mean, std::vector<double>& fMean) {
    std::vector<std::vector<unsigned>> p(opt::nSamp, std::vector<unsigned>(opt::covMax, 0));
    for (size_t i = 0; i < opt::nSamp; i++) {
        for (size_t j = 0; j < opt::rBuck; j++) {
            ++p[i][t_Counter[i * opt::rBuck + j]];
        }
    }

    // pMean[i]: mean over samples of the number of counters holding value i.
    double pMean0 = 0.0, pMean1 = 0.0;
    for (size_t j = 0; j < opt::nSamp; j++) {
        pMean0 += p[j][0];
        pMean1 += p[j][1];
    }
    pMean0 /= static_cast<double>(opt::nSamp);
    pMean1 /= static_cast<double>(opt::nSamp);

    F0Mean = static_cast<double>(
        (opt::rBits * std::log(2) - std::log(pMean0)) *
        (static_cast<size_t>(1) << (opt::sBits + opt::rBits)));

    fMean[1] = -1.0 * pMean1 / (pMean0 * (std::log(pMean0) - opt::rBits * std::log(2)));
    fMean[1] = std::abs(static_cast<ssize_t>(fMean[1] * F0Mean));
}

/**
 * Estimate k-mer cardinality from input sequence files.
 *
 * @param kmerNum  Total number of k-mers processed (output).
 * @param dbsize  Estimated number of distinct k-mers (output).
 * @param sbsize  Estimated number of solid k-mers (output).
 * @param kmer_len  K-mer length.
 * @param num_thread  Number of threads.
 * @param inFiles  Input sequence file paths.
 * @return Always returns false.
 */
bool getCardinality(size_t& kmerNum, size_t& dbsize, size_t& sbsize,
                    unsigned kmer_len, unsigned num_thread,
                    const std::vector<std::string>& inFiles) {
    std::vector<unsigned> kList = {kmer_len};
    std::vector<size_t> totalKmers(kList.size(), 0);
    opt::rBuck = static_cast<size_t>(1) << opt::rBits;
    opt::sMask = (static_cast<size_t>(1) << (opt::sBits - 1)) - 1;
    uint16_t* t_Counter = new uint16_t[opt::nK * opt::nSamp * opt::rBuck]();

#ifdef _OPENMP
    omp_set_num_threads(num_thread);
#endif

    #pragma omp parallel for schedule(dynamic)
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
        size_t totKmer[kList.size()];
        for (unsigned k = 0; k < kList.size(); k++) {
            totKmer[k] = 0;
        }
        gzFile fp = gzopen(inFiles[file_i].c_str(), "r");
        kseq_t* record = kseq_init(fp);
        while (kseq_read(record) >= 0) {
            std::string seq(record->seq.s, record->seq.l);
            ntRead(seq, kList, t_Counter, totKmer);
        }
        kseq_destroy(record);
        gzclose(fp);
        for (unsigned k = 0; k < kList.size(); k++) {
            #pragma omp atomic
            totalKmers[k] += totKmer[k];
        }
    }

    double F0Mean = 0.0;
    std::vector<double> fMean(opt::covMax, 0.0);
    compEst(t_Counter, F0Mean, fMean);

    kmerNum = totalKmers[0];
    dbsize = static_cast<size_t>(F0Mean);
    sbsize = dbsize - static_cast<size_t>(fMean[1]);

    delete[] t_Counter;
    return false;
}

/**
 * The sample ntCard draws from one file, kept instead of counted: the counter
 * slots (indBit * rBuck + shVal, as ntComp indexes them) its sampled k-mers fall
 * into, with multiplicity, sorted, and the file's k-mer count. The counter table is
 * additive over files, so the sample of a set of files is the concatenation of
 * their samples, and the estimate getCardinality would compute for the set
 * (estimateCardinality) comes from the merged slots alone. About one k-mer in a
 * thousand is sampled (two of 2^sBits top-bit patterns), 10 k slots for a 5 Mb
 * genome.
 */
struct CardinalitySample {
    std::vector<uint32_t> slots;
    size_t kmers = 0;
};

inline CardinalitySample ntSample(const std::string& file, unsigned kmerLen) {
    opt::rBuck = static_cast<size_t>(1) << opt::rBits;
    opt::sMask = (static_cast<size_t>(1) << (opt::sBits - 1)) - 1;
    CardinalitySample out;
    gzFile fp = gzopen(file.c_str(), "r");
    kseq_t* record = kseq_init(fp);
    while (kseq_read(record) >= 0) {
        std::string seq(record->seq.s, record->seq.l);
        ntHashIterator itr(seq, 1, kmerLen);
        while (itr != itr.end()) {
            const uint64_t hVal = (*itr)[0];
            uint64_t indBit = opt::nSamp;
            if (hVal >> (64 - opt::sBits) == 1) indBit = 0;
            if (hVal >> (64 - opt::sBits) == opt::sMask) indBit = 1;
            if (indBit < opt::nSamp)
                out.slots.push_back(static_cast<uint32_t>(indBit * opt::rBuck + (hVal & (opt::rBuck - 1))));
            ++itr;
            ++out.kmers;
        }
    }
    kseq_destroy(record);
    gzclose(fp);
    std::sort(out.slots.begin(), out.slots.end());
    return out;
}

/**
 * getCardinality's estimate for the union of the files behind `parts`, from their
 * samples: the same F0 and f1 as compEst reads off the counter table (a slot's
 * count taken as the uint16 the table would hold), so the same filter sizes.
 * Linear in the samples rather than in the 2^28-entry table, which getCardinality
 * allocates, zeroes and walks on every call: the within-group build on 5 493
 * genomes made 1 125 such calls, 419 s, now one pass over the genomes and a merge
 * per group.
 */
inline void estimateCardinality(const std::vector<const CardinalitySample*>& parts,
                                size_t& kmerNum, size_t& dbsize, size_t& sbsize) {
    opt::rBuck = static_cast<size_t>(1) << opt::rBits;
    std::vector<uint32_t> slots;
    size_t total = 0;
    kmerNum = 0;
    for (const CardinalitySample* p : parts) { total += p->slots.size(); kmerNum += p->kmers; }
    slots.reserve(total);
    for (const CardinalitySample* p : parts) slots.insert(slots.end(), p->slots.begin(), p->slots.end());
    std::sort(slots.begin(), slots.end());
    // Per sample index: counters at zero (every slot not seen, and any whose uint16
    // count wrapped to zero) and counters at one.
    std::vector<size_t> zeros(opt::nSamp, opt::rBuck), ones(opt::nSamp, 0);
    for (size_t i = 0; i < slots.size();) {
        size_t j = i + 1;
        while (j < slots.size() && slots[j] == slots[i]) ++j;
        const size_t s = slots[i] / opt::rBuck;
        const uint16_t c = static_cast<uint16_t>(j - i);
        if (c != 0) --zeros[s];
        if (c == 1) ++ones[s];
        i = j;
    }
    double pMean0 = 0.0, pMean1 = 0.0;
    for (size_t j = 0; j < opt::nSamp; j++) {
        pMean0 += static_cast<unsigned>(zeros[j]);
        pMean1 += static_cast<unsigned>(ones[j]);
    }
    pMean0 /= static_cast<double>(opt::nSamp);
    pMean1 /= static_cast<double>(opt::nSamp);
    const double F0Mean = static_cast<double>(
        (opt::rBits * std::log(2) - std::log(pMean0)) *
        (static_cast<size_t>(1) << (opt::sBits + opt::rBits)));
    double f1 = -1.0 * pMean1 / (pMean0 * (std::log(pMean0) - opt::rBits * std::log(2)));
    f1 = std::abs(static_cast<ssize_t>(f1 * F0Mean));
    dbsize = static_cast<size_t>(F0Mean);
    sbsize = dbsize - static_cast<size_t>(f1);
}
#endif // NTCARD_HPP_
