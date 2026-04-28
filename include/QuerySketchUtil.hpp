#ifndef QUERYSKETCHUTIL_HPP_
#define QUERYSKETCHUTIL_HPP_

#include <algorithm>
#include <fstream>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "SequenceUtil.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"
#include "FileUtil.hpp"

#include "kseq_gz.h"

namespace opt {
unsigned kmerLen;                                       // k-mer length
unsigned threads(1);                                    // number of threads
std::string out("out_sketch.tsv");                      // output file name
std::string r1;                                         // read 1 input
std::string r2;                                         // read 2 input
std::string ref;                                        // input uniqsketch reference
unsigned nhash(3);                                      // number of hashes for Bloom filter
size_t dbfSize;                                         // distinct Bloom filter size
size_t sbfSize;                                         // solid Bloom filter size
unsigned bits(64);                                      // bits per element in Bloom filter
int sketchHit(10);                                      // unique k-mer hits to call a reference
int solid(0);                                           // flag: only use solid k-mers from reads
double abundanceCutoff(0.0);                            // abundance cutoff to report
int readCutoff(2);                                      // read number cutoff to report
}

using SketchHash = std::unordered_map<std::string, unsigned>;

/**
 * Load the uniqsketch index into a hash table.
 *
 * @param fPath  Path to the uniqsketch file.
 * @param sketchHash  Hash table for the uniqsketch.
 * @param sketchRef  Vector mapping refId to refName.
 * @param refSigCount  Per-reference signature count map.
 */
void loadSketch(const std::string& fPath, SketchHash& sketchHash,
                std::vector<std::string>& sketchRef,
                std::vector<SketchHash>& refSigCount) {
    std::ifstream sketchFile(fPath);
    bool init = true;
    std::string seq;
    while (std::getline(sketchFile, seq)) {
        std::string uniqSeq, refName, contig;
        unsigned refId, refNum, pos;
        std::istringstream seqstm(seq);
        seqstm >> uniqSeq >> refId >> refName >> refNum >> contig >> pos;
        if (seqstm.fail()) {
            std::cerr << "Error in parsing uniqsketch entry:\n" << seq << "\n";
            exit(EXIT_FAILURE);
        }
        sketchHash.emplace(uniqSeq, refId);
        if (init) {
            opt::kmerLen = uniqSeq.length();
            sketchRef.resize(refNum);
            refSigCount.resize(refNum);
            init = false;
        }
        sketchRef[refId] = refName;
        refSigCount[refId].emplace(uniqSeq, 0);
    }
    sketchFile.close();
}

/**
 * Query an input read file against the uniqsketch index.
 *
 * @param in  Path to the input read file.
 * @param sketchHash  Hash table for the uniqsketch.
 * @param dbFilter  Bloom filter for tracking distinct k-mers.
 * @param sbFilter  Bloom filter for tracking solid k-mers.
 * @param sketchCount  Per-reference unique k-mer hit counts.
 * @param refSigCount  Per-reference per-signature counts.
 * @param refRead  Per-reference distinct read id sets.
 */
void querySample(const std::string& in, const SketchHash& sketchHash,
                 BloomFilter& dbFilter, BloomFilter& sbFilter,
                 std::vector<unsigned>& sketchCount,
                 std::vector<SketchHash>& refSigCount,
                 std::vector<std::unordered_set<std::string>>& refRead) {
    static const size_t BATCH_SIZE = 4096;
    gzFile fp = gzopen(in.c_str(), "r");
    kseq_t* ks = kseq_init(fp);

    std::vector<std::string> seqBatch;
    std::vector<std::string> nameBatch;
    seqBatch.reserve(BATCH_SIZE);
    nameBatch.reserve(BATCH_SIZE);

    bool eof = false;
    while (!eof) {
        seqBatch.clear();
        nameBatch.clear();
        for (size_t i = 0; i < BATCH_SIZE; ++i) {
            if (kseq_read(ks) < 0) {
                eof = true;
                break;
            }
            seqBatch.emplace_back(ks->seq.s, ks->seq.l);
            nameBatch.emplace_back(ks->name.s, ks->name.l);
        }

        #pragma omp parallel for schedule(dynamic, 64)
        for (size_t r = 0; r < seqBatch.size(); ++r) {
            const std::string& seq = seqBatch[r];
            if (seq.length() < opt::kmerLen) continue;

            ntHashIterator itr(seq, opt::nhash, opt::kmerLen);
            while (itr != itr.end()) {
                if (opt::solid) {
                    if (dbFilter.insert_make_change(*itr)) {
                        ++itr;
                        continue;
                    }
                }

                std::string kmer = get_canonical(seq.substr(itr.get_pos(), opt::kmerLen));
                auto query = sketchHash.find(kmer);
                if (query != sketchHash.end()) {
                    unsigned refId = query->second;

                    #pragma omp atomic
                    ++sketchCount[refId];

                    #pragma omp critical(refout)
                    refRead[refId].insert(nameBatch[r]);

                    #pragma omp atomic
                    refSigCount[refId][kmer]++;
                }
                ++itr;
            }
        }
    }
    kseq_destroy(ks);
    gzclose(fp);
}

/**
 * Query batch of input read files against the uniqsketch index.
 *
 * @param sampleFiles  All input read file paths.
 * @param sketchHash  Hash table for the uniqsketch.
 * @param sketchCount  Per-reference unique k-mer hit counts.
 * @param refSigCount  Per-reference per-signature counts.
 * @param refRead  Per-reference distinct read id sets.
 */
void querySampleBatch(const std::vector<std::string>& sampleFiles,
                      const SketchHash& sketchHash,
                      std::vector<unsigned>& sketchCount,
                      std::vector<SketchHash>& refSigCount,
                      std::vector<std::unordered_set<std::string>>& refRead) {
    BloomFilter dbFilter(opt::dbfSize * opt::bits, opt::nhash, opt::kmerLen);
    BloomFilter sbFilter(opt::sbfSize * opt::bits, opt::nhash, opt::kmerLen);

    for (unsigned i = 0; i < sampleFiles.size(); i++) {
        querySample(sampleFiles[i], sketchHash, dbFilter, sbFilter,
                    sketchCount, refSigCount, refRead);
    }

    // Adjust for distinct/solid BF masking
    if (opt::solid) {
        for (auto& count : sketchCount) {
            if (count) {
                ++count;
            }
        }
    }
}

/**
 * Generate the query result and output to TSV files.
 *
 * @param sketchCount  Per-reference unique k-mer hit counts.
 * @param sketchRef  Vector mapping refId to refName.
 * @param refSigCount  Per-reference per-signature counts.
 * @param refRead  Per-reference distinct read id sets.
 * @param fPath  Path to the output result TSV file.
 */
void generateQueryResult(std::vector<unsigned>& sketchCount,
                         const std::vector<std::string>& sketchRef,
                         const std::vector<SketchHash>& refSigCount,
                         const std::vector<std::unordered_set<std::string>>& refRead,
                         const std::string& fPath) {
    // Sort index by descending count
    std::vector<unsigned> sketchIndex(sketchCount.size());
    std::iota(sketchIndex.begin(), sketchIndex.end(), 0);
    std::sort(sketchIndex.begin(), sketchIndex.end(),
              [&](unsigned i, unsigned j) { return sketchCount[i] > sketchCount[j]; });

    size_t sketchSum = 0;
    for (unsigned i = 0; i < sketchCount.size(); i++) {
        if (refRead[i].size() > static_cast<size_t>(opt::readCutoff) &&
            sketchCount[i] >= static_cast<unsigned>(opt::sketchHit)) {
            sketchSum += sketchCount[i];
        }
    }

    // Primary output: reference abundance table
    std::ofstream outFile(opt::out);
    outFile << "ref\tabundance\tcount\n";
    if (sketchSum != 0) {
        for (unsigned i = 0; i < sketchCount.size(); i++) {
            unsigned idx = sketchIndex[i];
            if (refRead[idx].size() > static_cast<size_t>(opt::readCutoff) &&
                sketchCount[idx] >= static_cast<unsigned>(opt::sketchHit)) {
                double abundance = static_cast<double>(sketchCount[idx]) / sketchSum;
                if (abundance > opt::abundanceCutoff) {
                    outFile << sketchRef[idx] << "\t" << abundance
                            << "\t" << sketchCount[idx] << "\n";
                }
            }
        }
    }
    outFile.close();

    // Detailed signature log
    std::ofstream logSig(prefixBasename(opt::out, "log_"));
    logSig << "ref\ttotal_reads\tsignature\tcount\n";
    for (unsigned i = 0; i < refSigCount.size(); i++) {
        for (const auto& [sig, count] : refSigCount[i]) {
            if (count) {
                logSig << sketchRef[i] << "\t" << refRead[i].size()
                       << "\t" << sig << "\t" << count << "\n";
            }
        }
    }

    // Read-level log: matched read ids per reference
    std::ofstream logRead(prefixBasename(opt::out, "logread_"));
    logRead << "ref\tnum_reads\treads\n";
    for (unsigned i = 0; i < refRead.size(); i++) {
        if (refRead[i].size() > 1) {
            logRead << sketchRef[i] << "\t" << refRead[i].size() << "\t";
            for (const auto& rd : refRead[i]) {
                logRead << rd << "\t";
            }
            logRead << "\n";
        }
    }
}

#endif // QUERYSKETCHUTIL_HPP_