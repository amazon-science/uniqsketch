#ifndef QUERYSKETCHUTIL_HPP_
#define QUERYSKETCHUTIL_HPP_

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "SequenceUtil.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"
#include "FileUtil.hpp"

#include "kseq_gz.h"

/**
 * Settings for one query run. Passed explicitly; see UniqSketchConfig for why
 * these are no longer file-scope globals.
 *
 * kmerLen is not a user option: loadSketch sets it from the signature length in
 * the index, which is why loadSketch takes the config by non-const reference and
 * every other function takes it by const reference.
 */
struct QuerySketchConfig {
    // Inputs and outputs
    std::string out = "out_sketch.tsv";                 // output filename
    std::string r1;                                     // read 1 input
    std::string r2;                                     // read 2 input
    std::string ref;                                    // uniqsketch index (leaf sketch in tree mode)

    // Core parameters
    unsigned kmerLen = 0;                               // k-mer length, read from the index by loadSketch
    unsigned threads = 1;                               // number of threads

    // Read Bloom filters (sized per sample by ntCard)
    unsigned bits = 64;                                 // bits per element
    unsigned nhash = 3;                                 // hash functions
    size_t dbfSize = 0;                                 // distinct BF element count
    size_t sbfSize = 0;                                 // solid BF element count
    int solid = 0;                                      // --solid: query only k-mers seen at least twice

    // Reporting thresholds
    int sketchHit = 10;                                 // min signature hits to call a reference (-h)
    int readCutoff = 2;                                 // min distinct supporting reads, exclusive (-s)
    double abundanceCutoff = 0.0;                       // min abundance to report, exclusive (-a)

    // Tree query
    int treePruneHit = 1;                               // min coarse-level hits to keep a sub-tree in play
};

using SketchHash = std::unordered_map<std::string, unsigned>;

/**
 * Load the uniqsketch index into a hash table.
 *
 * @param fPath  Path to the uniqsketch file.
 * @param sketchHash  Hash table for the uniqsketch.
 * @param sketchRef  Vector mapping refId to refName.
 * @param refSigCount  Per-reference signature count map.
 */
void loadSketch(QuerySketchConfig& cfg, const std::string& fPath, SketchHash& sketchHash,
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
            cfg.kmerLen = uniqSeq.length();
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
/**
 * Does the k-mer at `kp` equal the canonical signature `sig`, in either
 * orientation? Compares in place -- no substring or reverse-complement
 * allocation -- which matters because this runs on every hash-probe hit.
 *
 * @param kp  Pointer to the first base of the read k-mer (sig.size() bases).
 * @param sig  Canonical signature sequence.
 * @return True if the read k-mer matches sig forward or as its reverse complement.
 */
inline bool kmerMatches(const char* kp, const std::string& sig) {
    const size_t k = sig.size();
    if (std::memcmp(kp, sig.data(), k) == 0) return true;
    for (size_t i = 0; i < k; i++) {
        char c;
        switch (kp[k - 1 - i]) {
            case 'A': case 'a': c = 'T'; break;
            case 'C': case 'c': c = 'G'; break;
            case 'G': case 'g': c = 'C'; break;
            case 'T': case 't': c = 'A'; break;
            default: return false;
        }
        if (c != sig[i]) return false;
    }
    return true;
}

/**
 * Hash-keyed view of a sketch for the query hot loop.
 *
 * Each signature is stored under its canonical ntHash value, the same value the
 * read iterator produces for every k-mer position, so a position can be tested
 * against the whole sketch with one 64-bit probe and no string work. On a probe
 * hit the stored sequence is compared to the read's k-mer before it counts, so a
 * hash collision can never produce a false hit: results are identical to a
 * string-keyed lookup, only the misses got cheap -- and in a typical query more
 * than 99.99% of positions are misses.
 */
struct SketchEntry {
    std::string seq;        // canonical signature sequence, for collision verification
    unsigned refId;         // reference the signature belongs to
    unsigned sigIdx;        // dense index of this signature within its reference
};
using SketchIndex = std::unordered_map<uint64_t, std::vector<SketchEntry>>;

/**
 * Build the hash-keyed view of a sketch and a dense per-reference signature
 * ordering, so per-thread hit tallies can be plain integer arrays rather than
 * string-keyed maps.
 *
 * @param sketchHash  Signature -> refId map, as loaded by loadSketch.
 * @param refSigCount  Per-reference signature -> count maps (defines the order).
 * @param sigNames  Output: for each reference, the signature at each dense index.
 * @return The hash-keyed index.
 */
inline SketchIndex buildSketchIndex(const SketchHash& sketchHash,
                                    const std::vector<SketchHash>& refSigCount,
                                    std::vector<std::vector<std::string>>& sigNames) {
    std::vector<std::unordered_map<std::string, unsigned>> sigIndexOf(refSigCount.size());
    sigNames.assign(refSigCount.size(), {});
    for (size_t ref = 0; ref < refSigCount.size(); ref++) {
        for (const auto& [sig, c] : refSigCount[ref]) {
            sigIndexOf[ref].emplace(sig, static_cast<unsigned>(sigNames[ref].size()));
            sigNames[ref].push_back(sig);
        }
    }

    SketchIndex index;
    index.reserve(sketchHash.size());
    for (const auto& [sig, refId] : sketchHash) {
        uint64_t h = NTC64(sig.c_str(), static_cast<unsigned>(sig.size()));
        auto it = sigIndexOf[refId].find(sig);
        // A signature present in the sketch but absent from refSigCount (e.g. a
        // filtered tree-level view) is still counted against its reference.
        unsigned sigIdx = (it != sigIndexOf[refId].end()) ? it->second : UINT32_MAX;
        index[h].push_back({sig, refId, sigIdx});
    }
    return index;
}

/**
 * Query an input read file against the uniqsketch index.
 *
 * Threads accumulate hits into private tallies and merge them once per batch,
 * so the hot loop takes no locks and performs no atomics. The previous version
 * took a global critical section on every hit, and applied `omp atomic` to a
 * map subscript -- which only covers the final increment, leaving the map
 * lookup itself unsynchronised.
 *
 * @param in  Path to the input read file.
 * @param index  Hash-keyed sketch view from buildSketchIndex.
 * @param sigNames  Dense signature ordering from buildSketchIndex.
 * @param dbFilter  Bloom filter for tracking distinct k-mers.
 * @param sbFilter  Bloom filter for tracking solid k-mers.
 * @param sketchCount  Per-reference unique k-mer hit counts.
 * @param refSigCount  Per-reference per-signature counts.
 * @param refRead  Per-reference distinct read id sets.
 */
void querySample(const QuerySketchConfig& cfg, const std::string& in, const SketchIndex& index,
                 const std::vector<std::vector<std::string>>& sigNames,
                 BloomFilter& dbFilter, BloomFilter& sbFilter,
                 std::vector<unsigned>& sketchCount,
                 std::vector<SketchHash>& refSigCount,
                 std::vector<std::unordered_set<std::string>>& refRead) {
    (void) sbFilter;
    static const size_t BATCH_SIZE = 4096;
    gzFile fp = gzopen(in.c_str(), "r");
    kseq_t* ks = kseq_init(fp);

    const size_t numRef = sketchCount.size();
    std::vector<std::string> seqBatch;
    std::vector<std::string> nameBatch;
    seqBatch.reserve(BATCH_SIZE);
    nameBatch.reserve(BATCH_SIZE);

    // Per-thread tallies, allocated once and reused across batches.
    struct Tally {
        std::vector<unsigned> count;                          // per reference
        std::vector<std::vector<unsigned>> sig;               // per reference, per dense sig index
        std::vector<std::vector<unsigned>> readIdx;           // per reference, batch-local read indices
    };
    int maxThreads = 1;
#ifdef _OPENMP
    maxThreads = omp_get_max_threads();
#endif
    std::vector<Tally> tallies(maxThreads);
    for (Tally& t : tallies) {
        t.count.assign(numRef, 0);
        t.sig.resize(numRef);
        for (size_t ref = 0; ref < numRef; ref++) t.sig[ref].assign(sigNames[ref].size(), 0);
        t.readIdx.resize(numRef);
    }

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

        #pragma omp parallel
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            Tally& t = tallies[tid];

            #pragma omp for schedule(dynamic, 64)
            for (size_t r = 0; r < seqBatch.size(); ++r) {
                const std::string& seq = seqBatch[r];
                if (seq.length() < cfg.kmerLen) continue;

                ntHashIterator itr(seq, cfg.nhash, cfg.kmerLen);
                while (itr != itr.end()) {
                    if (cfg.solid) {
                        if (dbFilter.insert_make_change(*itr)) {
                            ++itr;
                            continue;
                        }
                    }

                    // One 64-bit probe decides the (overwhelmingly common) miss.
                    auto bucket = index.find((*itr)[0]);
                    if (bucket != index.end()) {
                        // Probe hit: confirm the sequence before counting.
                        const char* kp = seq.c_str() + itr.get_pos();
                        for (const SketchEntry& e : bucket->second) {
                            if (kmerMatches(kp, e.seq)) {
                                ++t.count[e.refId];
                                if (e.sigIdx != UINT32_MAX) ++t.sig[e.refId][e.sigIdx];
                                t.readIdx[e.refId].push_back(static_cast<unsigned>(r));
                                break;
                            }
                        }
                    }
                    ++itr;
                }
            }
        }

        // Merge this batch's per-thread tallies into the shared results.
        for (Tally& t : tallies) {
            for (size_t ref = 0; ref < numRef; ref++) {
                if (t.count[ref] == 0) continue;
                sketchCount[ref] += t.count[ref];
                t.count[ref] = 0;
                for (size_t s = 0; s < t.sig[ref].size(); s++) {
                    if (t.sig[ref][s]) {
                        refSigCount[ref][sigNames[ref][s]] += t.sig[ref][s];
                        t.sig[ref][s] = 0;
                    }
                }
                for (unsigned r : t.readIdx[ref]) refRead[ref].insert(nameBatch[r]);
                t.readIdx[ref].clear();
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
void querySampleBatch(const QuerySketchConfig& cfg, const std::vector<std::string>& sampleFiles,
                      const SketchHash& sketchHash,
                      std::vector<unsigned>& sketchCount,
                      std::vector<SketchHash>& refSigCount,
                      std::vector<std::unordered_set<std::string>>& refRead) {
    BloomFilter dbFilter(cfg.dbfSize * cfg.bits, cfg.nhash, cfg.kmerLen);
    BloomFilter sbFilter(cfg.sbfSize * cfg.bits, cfg.nhash, cfg.kmerLen);

    // Build the hash-keyed view once; it is shared read-only across all files.
    std::vector<std::vector<std::string>> sigNames;
    SketchIndex index = buildSketchIndex(sketchHash, refSigCount, sigNames);

    for (unsigned i = 0; i < sampleFiles.size(); i++) {
        querySample(cfg, sampleFiles[i], index, sigNames, dbFilter, sbFilter,
                    sketchCount, refSigCount, refRead);
    }

    // Adjust for distinct/solid BF masking
    if (cfg.solid) {
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
/**
 * A detection made at an internal (coarse) tree level: a cluster representative
 * whose signatures cleared the reporting thresholds. Used to report a family-level
 * call when the evidence identifies a group of near-identical references but no
 * single leaf in that group can be resolved.
 */
struct FamilyCall {
    std::string name;                                   // internal node name, e.g. cluster_L1_AWS_Ec377
    unsigned level;                                     // tree level the detection was made at
    unsigned count;                                     // signature hits at that level
    size_t reads;                                       // distinct reads carrying those hits
    size_t sigsHit;                                     // distinct signatures observed at least once
    size_t sigsTotal;                                   // signatures the representative carries at that level
    std::unordered_set<std::string> coveredLeaves;      // leaves descended from this node
};

/**
 * Count how many of a reference's signatures were observed at least once.
 *
 * @param sigCount  Signature -> hit count map for one reference (all signatures
 *                  present, unobserved ones at zero).
 * @return Number of entries with a non-zero count.
 */
inline size_t countSigsHit(const SketchHash& sigCount) {
    size_t n = 0;
    for (const auto& [sig, c] : sigCount) {
        if (c) ++n;
    }
    return n;
}

/**
 * Signature hits per signature: the quantity abundance is computed from.
 *
 * A reference's total hit count grows with the number of signatures it carries
 * as well as with its coverage. When every reference has the same number of
 * signatures the two are interchangeable and a share of raw counts is a share of
 * coverage. They stop being interchangeable as soon as the counts differ: a
 * close relative whose whole candidate pool is smaller than -c gets fewer
 * signatures than a distant reference, and a tree index built with a large -c
 * makes that the normal case. On the 8-phage benchmark with -c500 (119-500
 * signatures per reference) a 10% contaminant was reported anywhere between
 * 2.6% and 32% from raw counts, and between 9.3% and 11.3% from hits per
 * signature.
 * Dividing by the signature count removes the dependence; with equal counts it
 * cancels out and the result is unchanged.
 */
inline double hitsPerSignature(size_t count, size_t sigsTotal) {
    return sigsTotal ? static_cast<double>(count) / static_cast<double>(sigsTotal) : 0.0;
}

/**
 * Write the per-signature hit log (`log_<out>`): one row per signature that was
 * seen at least once, with the number of distinct reads supporting its reference.
 * Shared by the flat and tree writers so both produce the same evidence trail.
 */
inline void writeSignatureLog(const QuerySketchConfig& cfg,
                              const std::vector<std::string>& sketchRef,
                              const std::vector<SketchHash>& refSigCount,
                              const std::vector<std::unordered_set<std::string>>& refRead) {
    std::ofstream logSig(prefixBasename(cfg.out, "log_"));
    logSig << "ref\ttotal_reads\tsignature\tcount\n";
    for (unsigned i = 0; i < refSigCount.size(); i++) {
        for (const auto& [sig, count] : refSigCount[i]) {
            if (count) {
                logSig << sketchRef[i] << "\t" << refRead[i].size()
                       << "\t" << sig << "\t" << count << "\n";
            }
        }
    }
}

/**
 * Write the read-level log (`logread_<out>`): for every reference with at least
 * one matching read, the ids of those reads. This is the only place the tool
 * exposes which reads carried a call, so a reference supported by a single read
 * is listed too -- at low abundance that one read is the whole case for or
 * against the call, and hiding it would leave the log unable to explain the
 * `reads` column of the result table.
 */
inline void writeReadLog(const QuerySketchConfig& cfg,
                         const std::vector<std::string>& sketchRef,
                         const std::vector<std::unordered_set<std::string>>& refRead) {
    std::ofstream logRead(prefixBasename(cfg.out, "logread_"));
    logRead << "ref\tnum_reads\treads\n";
    for (unsigned i = 0; i < refRead.size(); i++) {
        if (refRead[i].empty()) continue;
        logRead << sketchRef[i] << "\t" << refRead[i].size() << "\t";
        for (const auto& rd : refRead[i]) {
            logRead << rd << "\t";
        }
        logRead << "\n";
    }
}

/**
 * Write the query result for a tree index.
 *
 * Leaf references are reported exactly as in the flat case. In addition, a
 * family-level row is emitted for any coarse-level detection that cleared the
 * reporting thresholds but whose descendant leaves were all left unreported.
 *
 * The motivation is sensitivity. A cluster representative is chosen among distant
 * references, so its signatures are drawn from a far larger pool of candidates than
 * any of its near-identical members can offer, and it stays detectable at
 * abundances where no individual member is separable. Reporting only leaves throws
 * that evidence away and returns nothing; reporting the family instead degrades
 * gracefully, naming the group when the member cannot be pinned down. Family rows
 * carry the internal node name from the index's tree_index.tsv, which lists the
 * members they stand for.
 *
 * @param sketchCount  Per-leaf signature hit counts.
 * @param sketchRef  Leaf reference names.
 * @param refSigCount  Per-leaf signature:count maps, for the detailed log.
 * @param refRead  Per-leaf sets of matching read ids.
 * @param families  Coarse-level detections, ordered coarsest first.
 * @param fPath  Unused; retained to mirror generateQueryResult's signature.
 */
void generateTreeQueryResult(const QuerySketchConfig& cfg, std::vector<unsigned>& sketchCount,
                             const std::vector<std::string>& sketchRef,
                             const std::vector<SketchHash>& refSigCount,
                             const std::vector<std::unordered_set<std::string>>& refRead,
                             const std::vector<FamilyCall>& families,
                             const std::string& fPath) {
    (void) fPath;

    auto leafReported = [&](unsigned i) {
        return refRead[i].size() > static_cast<size_t>(cfg.readCutoff) &&
               sketchCount[i] >= static_cast<unsigned>(cfg.sketchHit);
    };

    // Which leaves clear the reporting thresholds.
    std::unordered_set<std::string> reportedLeaves;
    for (unsigned i = 0; i < sketchCount.size(); i++) {
        if (leafReported(i)) reportedLeaves.insert(sketchRef[i]);
    }

    // A family is only worth reporting when it is itself confident and none of the
    // leaves it covers was resolved; otherwise the leaf rows already say more.
    std::vector<const FamilyCall*> reportedFamilies;
    for (const FamilyCall& f : families) {
        if (f.count < static_cast<unsigned>(cfg.sketchHit) ||
            f.reads <= static_cast<size_t>(cfg.readCutoff)) {
            continue;
        }
        bool anyLeafResolved = false;
        for (const std::string& leaf : f.coveredLeaves) {
            if (reportedLeaves.count(leaf)) { anyLeafResolved = true; break; }
        }
        if (!anyLeafResolved) reportedFamilies.push_back(&f);
    }

    // Merge leaf and family rows, then order by descending count.
    struct Row {
        std::string name;
        unsigned count;
        size_t reads;
        size_t sigsHit;
        size_t sigsTotal;
        bool isFamily;
    };
    std::vector<Row> rows;
    for (unsigned i = 0; i < sketchCount.size(); i++) {
        if (leafReported(i)) {
            rows.push_back({sketchRef[i], sketchCount[i], refRead[i].size(),
                            countSigsHit(refSigCount[i]), refSigCount[i].size(), false});
        }
    }
    for (const FamilyCall* f : reportedFamilies) {
        rows.push_back({f->name, f->count, f->reads, f->sigsHit, f->sigsTotal, true});
    }
    std::sort(rows.begin(), rows.end(),
              [](const Row& a, const Row& b) { return a.count > b.count; });

    // Abundance is each row's hits per signature, relative to everything reported,
    // families included, so the column still sums to 1 over the rows in the file.
    // See hitsPerSignature for why the count is normalised.
    double densitySum = 0.0;
    for (const Row& r : rows) densitySum += hitsPerSignature(r.count, r.sigsTotal);

    // Same schema as the flat writer; `resolution` distinguishes individually
    // resolved references (`strain`) from group calls (`family`), so a family
    // row can never be mistaken for a strain call downstream.
    std::ofstream outFile(cfg.out);
    outFile << "ref\tabundance\tcount\treads\tsigs_hit\tsigs_total\tresolution\n";
    if (densitySum > 0.0) {
        for (const Row& r : rows) {
            double abundance = hitsPerSignature(r.count, r.sigsTotal) / densitySum;
            if (abundance > cfg.abundanceCutoff) {
                outFile << r.name << "\t" << abundance << "\t" << r.count
                        << "\t" << r.reads << "\t" << r.sigsHit
                        << "\t" << r.sigsTotal
                        << "\t" << (r.isFamily ? "family" : "strain") << "\n";
            }
        }
    }
    outFile.close();

    if (!reportedFamilies.empty()) {
        std::cerr << "Tree query: " << reportedFamilies.size()
                  << " family-level call(s) reported where no leaf could be resolved:";
        for (const FamilyCall* f : reportedFamilies) {
            std::cerr << " " << f->name << "(" << f->coveredLeaves.size() << " members)";
        }
        std::cerr << "\n";
    }

    // Signature and read logs over leaves, same files as the flat writer.
    writeSignatureLog(cfg, sketchRef, refSigCount, refRead);
    writeReadLog(cfg, sketchRef, refRead);
}

void generateQueryResult(const QuerySketchConfig& cfg, std::vector<unsigned>& sketchCount,
                         const std::vector<std::string>& sketchRef,
                         const std::vector<SketchHash>& refSigCount,
                         const std::vector<std::unordered_set<std::string>>& refRead,
                         const std::string& fPath) {
    // Sort index by descending count
    std::vector<unsigned> sketchIndex(sketchCount.size());
    std::iota(sketchIndex.begin(), sketchIndex.end(), 0);
    std::sort(sketchIndex.begin(), sketchIndex.end(),
              [&](unsigned i, unsigned j) { return sketchCount[i] > sketchCount[j]; });

    auto reported = [&](unsigned i) {
        return refRead[i].size() > static_cast<size_t>(cfg.readCutoff) &&
               sketchCount[i] >= static_cast<unsigned>(cfg.sketchHit);
    };

    // Abundance is each reference's hits per signature relative to all reported
    // references; see hitsPerSignature.
    double densitySum = 0.0;
    for (unsigned i = 0; i < sketchCount.size(); i++) {
        if (reported(i)) densitySum += hitsPerSignature(sketchCount[i], refSigCount[i].size());
    }

    // Primary output: reference abundance table. Columns beyond the third expose
    // the evidence behind each call: `reads` is the number of distinct reads
    // supporting it (two references can share an abundance while one rests on
    // thousands of reads and the other on three), and `sigs_hit`/`sigs_total`
    // show how many of the reference's signatures were seen at least once -- a
    // call carried by a handful of its signatures suggests a repeat or
    // cross-mapping, one supported across most of them does not. `resolution` is
    // `strain` for individually resolved references; the tree query may also
    // emit `family` rows (see generateTreeQueryResult). Columns are appended so
    // parsers of the original three are unaffected.
    std::ofstream outFile(cfg.out);
    outFile << "ref\tabundance\tcount\treads\tsigs_hit\tsigs_total\tresolution\n";
    if (densitySum > 0.0) {
        for (unsigned i = 0; i < sketchCount.size(); i++) {
            unsigned idx = sketchIndex[i];
            if (reported(idx)) {
                double abundance =
                    hitsPerSignature(sketchCount[idx], refSigCount[idx].size()) / densitySum;
                if (abundance > cfg.abundanceCutoff) {
                    outFile << sketchRef[idx] << "\t" << abundance
                            << "\t" << sketchCount[idx]
                            << "\t" << refRead[idx].size()
                            << "\t" << countSigsHit(refSigCount[idx])
                            << "\t" << refSigCount[idx].size()
                            << "\tstrain\n";
                }
            }
        }
    }
    outFile.close();

    writeSignatureLog(cfg, sketchRef, refSigCount, refRead);
    writeReadLog(cfg, sketchRef, refRead);
}

#endif // QUERYSKETCHUTIL_HPP_