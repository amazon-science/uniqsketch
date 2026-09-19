#ifndef UNIQSKETCHUTIL_HPP_
#define UNIQSKETCHUTIL_HPP_

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <numeric>
#include <random>
#include <sstream>
#include <unordered_map>

#include "SequenceUtil.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"

/**
 * Settings for one index build.
 *
 * Passed explicitly to every function that needs it. Previously these lived in a
 * mutable `namespace opt` defined in this header, which meant three things: no
 * two *Util headers could be included in one translation unit (each defined the
 * same names), no function could be tested under two settings in one process,
 * and the tree path had to save and restore nine globals around every level
 * because buildSketch mutated one of them as a side effect. A value type fixes
 * all three: copy it, change the copy, pass the copy.
 */
struct UniqSketchConfig {
    // Inputs and outputs
    std::string outdir = "outdir_uniqsketch";           // dir for per-reference unique k-mer candidates
    std::string refstat = "db_uniq_count.tsv";          // unique k-mer stat filename
    std::string outfile = "sketch_uniq.tsv";            // output sketch filename
    unsigned numRef = 0;                                // number of references in the universe

    // Core parameters
    unsigned threads = 1;                               // number of threads
    unsigned kmerLen = 81;                              // k-mer length
    int sketchnum = 100;                                // signatures per reference (-c)

    // Bloom filters
    unsigned bits = 128;                                // bits per element
    unsigned nhash1 = 5;                                // hashes, distinct BF
    unsigned nhash2 = 5;                                // hashes, solid BF
    size_t dbfSize = 0;                                 // distinct BF element count (from ntCard)
    size_t sbfSize = 0;                                 // solid BF element count (from ntCard)

    // Signature filters
    unsigned kRange = 5;                                // k'-mer spectrum range for the entropy filter
    double entropyThreshold = 0.65;                     // low-complexity threshold (fraction of max entropy)
    unsigned minMargin = 1;                             // min Hamming distance to other references (1 = off)
    unsigned maxHomopolymer = 0;                        // max homopolymer run in a signature (0 = off)
    bool strictMargin = false;                          // --min-margin is a hard requirement (no fallback)

    // Derived Bloom-filter bit sizes.
    size_t m1() const { return static_cast<size_t>(bits) * dbfSize; }
    size_t m2() const { return static_cast<size_t>(bits) * sbfSize; }

    /**
     * Maximum aggregate entropy A_K for a signature of length kmerLen: 2k' for
     * k' <= 3 (all 4^k' k'-mers can appear), log2(kmerLen - k' + 1) beyond that
     * (bounded by the number of k'-mer positions). Was previously accumulated into
     * a global by buildSketch, which is why the tree path had to snapshot it.
     */
    double maxEntropy() const {
        double e = 12.0;                                // 2 + 4 + 6 for k' = 1, 2, 3
        for (unsigned i = 4; i <= kRange; i++) e += std::log2(kmerLen - i + 1);
        return e;
    }
};

using SketchHash = std::unordered_map<std::string, unsigned>;

/**
 * Insert k-mers from a FASTA reference into distinct and solid Bloom filters.
 *
 * @param fPath  Path to a reference file.
 * @param dbFilter  Distinct Bloom filter for tracking distinct k-mers.
 * @param sbFilter  Solid Bloom filter for tracking repeat k-mers.
 */
void loadBFfa(const UniqSketchConfig& cfg, const std::string& fPath, BloomFilter& dbFilter, BloomFilter& sbFilter) {
    std::ifstream refFile(fPath);
    std::string seq;
    bool good = static_cast<bool>(std::getline(refFile, seq));
    while (good) {
        std::string faSeq;
        good = static_cast<bool>(std::getline(refFile, seq));
        while (good && seq[0] != '>') {
            faSeq += seq;
            good = static_cast<bool>(std::getline(refFile, seq));
        }
        ntHashIterator itr(faSeq, std::max(cfg.nhash1, cfg.nhash2), cfg.kmerLen);
        while (itr != itr.end()) {
            if (!dbFilter.insert_make_change(*itr)) {
                sbFilter.insert(*itr);
            }
            ++itr;
        }
    }
}

/**
 * Scan all k-mers in a reference and record unique k-mers.
 *
 * @param fPath  Path to a reference file.
 * @param refId  Integer id assigned to the reference.
 * @param sbFilter  Solid Bloom filter for tracking repeat k-mers.
 * @param refCount  Vector of unique k-mer counts per reference.
 */
void checkRef(const UniqSketchConfig& cfg, const std::string& fPath, int refId,
              const BloomFilter& sbFilter, std::vector<int>& refCount) {
    std::ifstream refFile(fPath);
    std::string sampleId = getBaseId(fPath);
    std::ofstream uniqFile(cfg.outdir + "/" + sampleId + ".tsv");
    std::string seq;
    size_t countUnique = 0;

    bool good = static_cast<bool>(std::getline(refFile, seq));
    std::string headerSeq = seq.substr(1);
    std::replace(headerSeq.begin(), headerSeq.end(), ' ', '_');

    while (good) {
        std::string faSeq;
        good = static_cast<bool>(std::getline(refFile, seq));
        while (good && seq[0] != '>') {
            faSeq += seq;
            good = static_cast<bool>(std::getline(refFile, seq));
        }
        ntHashIterator itr(faSeq, cfg.nhash2, cfg.kmerLen);
        while (itr != itr.end()) {
            if (!sbFilter.contains(*itr)) {
                ++countUnique;
                size_t seqLoc = itr.get_pos();
                uniqFile << get_canonical(faSeq.substr(seqLoc, cfg.kmerLen))
                         << "\t" << seqLoc << "\t" << headerSeq << "\n";
            }
            ++itr;
        }
        if (seq.size() > 1) {
            headerSeq = seq.substr(1);
            std::replace(headerSeq.begin(), headerSeq.end(), ' ', '_');
        }
    }
    uniqFile.close();
    refCount[refId] = countUnique;
}

/**
 * Identify all unique k-mers in the reference universe.
 *
 * @param refFiles  Vector of all reference file paths.
 * @param dbFilter  Distinct Bloom filter, owned by the caller so it can be
 *                  reused for the cross-reference margin check in buildSketch.
 */
void identifyUniqKmers(const UniqSketchConfig& cfg, const std::vector<std::string>& refFiles, BloomFilter& dbFilter) {
    std::cout << "k-mer length: " << cfg.kmerLen << "\n"
              << "Number of distinct k-mers: " << cfg.dbfSize << "\n"
              << "Number of repeat k-mers: " << cfg.sbfSize << "\n"
              << "Number of unique k-mers: " << cfg.dbfSize - cfg.sbfSize << "\n";

    BloomFilter sbFilter(cfg.m2(), cfg.nhash2, cfg.kmerLen);

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < refFiles.size(); ++i) {
        loadBFfa(cfg, refFiles[i], dbFilter, sbFilter);
    }

    std::vector<int> refCount(cfg.numRef);
    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < refFiles.size(); ++i) {
        checkRef(cfg, refFiles[i], i, sbFilter, refCount);
    }

    std::ofstream tsvOut(cfg.refstat);
    for (size_t i = 0; i < refCount.size(); i++) {
        tsvOut << cfg.outdir << "/" << getBaseId(refFiles[i]) << ".tsv\t"
               << refCount[i] << "\t" << i << "\n";
    }
    tsvOut.close();

    std::cout << "Distinct BF actual fpr: " << std::setprecision(4) << std::fixed
              << std::pow(static_cast<double>(dbFilter.get_pop()) / cfg.m1(), cfg.nhash1) << "\n"
              << "Solid BF actual fpr: " << std::setprecision(4) << std::fixed
              << std::pow(static_cast<double>(sbFilter.get_pop()) / cfg.m2(), cfg.nhash2) << "\n";
}

/**
 * Digitize a DNA sequence into an integer representation.
 *
 * @param seq  DNA sequence to digitize.
 * @return Integer representation of the DNA sequence.
 */
unsigned digitize(const UniqSketchConfig& cfg, const std::string& seq) {
    if (seq.length() > cfg.kRange) {
        return 0;
    }
    uint64_t hVal = 0;
    for (unsigned i = 0; i < seq.length(); i++) {
        hVal = (hVal << 2) | b2f[static_cast<unsigned char>(seq[i])];
    }
    return static_cast<unsigned>(hVal);
}

/**
 * Check whether a signature contains a homopolymer run longer than the allowed
 * maximum. Aggregate entropy dilutes a short local run across the full k-mer, so
 * this is a direct local-run-length gate for indel-prone tracts (e.g. poly-G).
 *
 * @param signature  Signature sequence to check.
 * @return True if the longest single-base run exceeds cfg.maxHomopolymer.
 */
bool hasLongHomopolymer(const UniqSketchConfig& cfg, const std::string& signature) {
    if (cfg.maxHomopolymer == 0) return false;   // filter disabled
    unsigned run = 1;
    for (size_t i = 1; i < signature.size(); i++) {
        if (signature[i] == signature[i - 1]) {
            if (++run > cfg.maxHomopolymer) return true;
        } else {
            run = 1;
        }
    }
    return false;
}

/**
 * Check a signature for low-complexity content using aggregate Shannon entropy.
 *
 * @param signature  Signature sequence to check.
 * @return True if the signature is low-complexity, false otherwise.
 */
bool lowComplexity(const UniqSketchConfig& cfg, const std::string& signature) {
    if (cfg.kmerLen <= cfg.kRange) {
        return true;
    }

    // Max spectrum size for kRange=5 is 4^5 = 1024
    static constexpr unsigned MAX_SPEC = 1024;
    unsigned kmerSpec[MAX_SPEC];

    double aggEntropy = 0.0;
    for (unsigned k = 1; k <= cfg.kRange; k++) {
        unsigned specSize = 1u << (2 * k);
        unsigned mask = specSize - 1;
        std::memset(kmerSpec, 0, specSize * sizeof(unsigned));

        // Seed the rolling hash with the first k-mer
        unsigned hVal = 0;
        for (unsigned j = 0; j < k; j++) {
            hVal = (hVal << 2) | b2f[static_cast<unsigned char>(signature[j])];
        }
        kmerSpec[hVal]++;

        // Roll through remaining positions: shift out old base, shift in new base
        for (unsigned i = 1; i <= signature.length() - k; i++) {
            hVal = ((hVal << 2) | b2f[static_cast<unsigned char>(signature[i + k - 1])]) & mask;
            kmerSpec[hVal]++;
        }

        double entropy = 0.0;
        double denom = signature.size() - k + 1;
        for (unsigned i = 0; i < specSize; i++) {
            if (kmerSpec[i]) {
                double p = static_cast<double>(kmerSpec[i]) / denom;
                entropy -= p * std::log2(p);
            }
        }
        aggEntropy += entropy;
    }
    return (aggEntropy / cfg.maxEntropy() <= cfg.entropyThreshold);
}

/**
 * Selection statistics returned by getUniqSet for a single reference.
 */
struct SketchStat {
    size_t selected;   // signatures written for this reference
    size_t safe;       // of those, how many satisfied the cross-reference margin
};

/**
 * Test whether a candidate signature is "2-safe": no k-mer in the reference
 * universe lies at Hamming distance 1 from it. Signatures are already unique, so
 * this makes the minimum distance to any foreign k-mer >= 2, meaning a single
 * substitution read error cannot turn a foreign k-mer into this signature.
 *
 * The distinct Bloom filter records only presence, not which reference a k-mer
 * belongs to, so a neighbor that occurs solely within this signature's own
 * reference also triggers a reject. This is deliberately conservative: it may
 * drop an otherwise-good signature but never keeps an unsafe one. ntHash is
 * canonical, so probing a mutated neighbor is reverse-complement-correct.
 *
 * @param sig  Candidate signature (length cfg.kmerLen).
 * @param dbFilter  Distinct Bloom filter over all universe k-mers.
 * @return True if no Hamming-1 neighbor is present in the universe.
 */
bool isTwoSafe(const UniqSketchConfig& cfg, const std::string& sig, const BloomFilter& dbFilter) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    std::string neighbor = sig;
    uint64_t hVal[cfg.nhash1];
    for (size_t i = 0; i < neighbor.size(); ++i) {
        const char orig = neighbor[i];
        for (char b : BASES) {
            if (b == orig) continue;
            neighbor[i] = b;
            NTMC64(neighbor.c_str(), cfg.kmerLen, cfg.nhash1, hVal);
            if (dbFilter.contains(hVal)) {
                neighbor[i] = orig;
                return false;
            }
        }
        neighbor[i] = orig;
    }
    return true;
}

/**
 * Extract the unique k-mer set for a given reference.
 *
 * @param fPath  Path to a reference's candidate unique k-mers.
 * @param refId  Integer id assigned to the reference.
 * @param sketchFilter  Bloom filter for the unique sketch.
 * @param out  Stream receiving this reference's signatures. buildSketch passes a
 *             per-reference buffer so output order is independent of thread
 *             scheduling.
 * @param dbFilter  Distinct Bloom filter over the whole universe (margin check).
 * @return Number of signatures selected and how many met the margin.
 */
SketchStat getUniqSet(const UniqSketchConfig& cfg, const std::string& fPath, unsigned refId,
                      BloomFilter& sketchFilter, std::ostream& out,
                      const BloomFilter& dbFilter) {
    // Read all candidate lines
    std::vector<std::string> lines;
    std::ifstream sketchFile(fPath);
    std::string line;
    while (std::getline(sketchFile, line)) {
        lines.push_back(std::move(line));
    }
    sketchFile.close();

    if (lines.empty()) return {0, 0};

    // Build a lightweight index: extract (contig, pos) for sorting without extra string copies
    size_t total = lines.size();
    std::vector<unsigned> indices(total);
    std::iota(indices.begin(), indices.end(), 0);

    // Parse position from each line for sorting (format: seq \t pos \t contig)
    // Use a fast scan: skip first field, read pos, read contig
    struct PosKey {
        unsigned pos;
        const char* contigStart;
        size_t contigLen;
    };
    std::vector<PosKey> keys(total);
    for (size_t i = 0; i < total; i++) {
        const char* s = lines[i].c_str();
        // skip seq field
        while (*s && *s != '\t') ++s;
        if (*s) ++s;
        // parse pos
        keys[i].pos = 0;
        while (*s >= '0' && *s <= '9') {
            keys[i].pos = keys[i].pos * 10 + (*s - '0');
            ++s;
        }
        if (*s) ++s;
        // contig starts here
        keys[i].contigStart = s;
        keys[i].contigLen = lines[i].c_str() + lines[i].size() - s;
    }

    std::sort(indices.begin(), indices.end(),
              [&keys](unsigned a, unsigned b) {
                  int cmp = std::strncmp(keys[a].contigStart, keys[b].contigStart,
                                         std::min(keys[a].contigLen, keys[b].contigLen));
                  if (cmp != 0) return cmp < 0;
                  if (keys[a].contigLen != keys[b].contigLen) return keys[a].contigLen < keys[b].contigLen;
                  return keys[a].pos < keys[b].pos;
              });

    // Select evenly spaced candidates, checking low-complexity (and, when
    // cfg.minMargin >= 2, the cross-reference Hamming margin) lazily.
    int count = 0;
    size_t safeCount = 0;
    size_t needed = std::min(static_cast<size_t>(cfg.sketchnum), total);
    double step = static_cast<double>(total) / needed;
    const bool useMargin = (cfg.minMargin >= 2);

    // Commit the first usable candidate in [a, b). When requireSafe is set, only
    // 2-safe candidates are accepted. Returns true if a signature was committed.
    auto trySlot = [&](size_t a, size_t b, bool requireSafe) -> bool {
        for (size_t i = a; i < b; i++) {
            unsigned idx = indices[i];
            std::istringstream seqstm(lines[idx]);
            std::string useq, contig;
            unsigned pos;
            seqstm >> useq >> pos >> contig;

            if (lowComplexity(cfg, useq)) continue;
            if (hasLongHomopolymer(cfg, useq)) continue;
            if (requireSafe && !isTwoSafe(cfg, useq, dbFilter)) continue;

            ntHashIterator itr(useq, cfg.nhash1, cfg.kmerLen);
            bool found = false;
            while (itr != itr.end()) {
                if (sketchFilter.insert_make_change(*itr)) {
                    ++count;
                    found = true;
                    // `out` is this reference's private buffer, so no lock is needed.
                    out << useq << "\t" << refId << "\t" << getBaseId(fPath)
                        << "\t" << cfg.numRef << "\t" << contig << "\t" << pos << "\n";
                }
                ++itr;
            }
            if (found) return true;
        }
        return false;
    };

    for (size_t s = 0; s < needed && count < cfg.sketchnum; s++) {
        size_t slotStart = static_cast<size_t>(s * step);
        size_t slotEnd = (s + 1 < needed) ? static_cast<size_t>((s + 1) * step) : total;

        if (useMargin) {
            // Prefer a 2-safe candidate. In soft mode (default) fall back to a
            // margin-1 candidate if the slot has no 2-safe option; in strict mode
            // leave the slot empty rather than admit a margin-1 signature.
            if (trySlot(slotStart, slotEnd, true)) ++safeCount;
            else if (!cfg.strictMargin) trySlot(slotStart, slotEnd, false);
        } else {
            trySlot(slotStart, slotEnd, false);
        }
    }

    return SketchStat{static_cast<size_t>(count), safeCount};
}

/**
 * Get the unique k-mer count from a stat line for sorting.
 *
 * @param sketchLine  Tab-separated line: path, count, refId.
 * @return Number of unique k-mers.
 */
inline int sketchRank(const std::string& sketchLine) {
    std::istringstream ss(sketchLine);
    std::string sketchPath;
    int sketchCount;
    ss >> sketchPath >> sketchCount;
    return sketchCount;
}

/**
 * Construct the unique sketch for the reference universe.
 *
 * @param statFile  Path to the unique k-mer stat file for all references.
 * @param dbFilter  Distinct Bloom filter over the whole universe (margin check).
 */
void buildSketch(const UniqSketchConfig& cfg, const std::string& statFile, const BloomFilter& dbFilter) {
    std::ifstream tsvIn(statFile);
    std::vector<std::string> uniqStats;
    std::string line;
    while (std::getline(tsvIn, line)) {
        uniqStats.push_back(line);
    }

    std::sort(uniqStats.begin(), uniqStats.end(),
              [](const std::string& s1, const std::string& s2) {
                  return sketchRank(s1) < sketchRank(s2);
              });

    BloomFilter sketchFilter(cfg.bits * cfg.sketchnum * uniqStats.size(),
                             cfg.nhash1, cfg.kmerLen);

    // Each reference writes its signatures into its own buffer, and the buffers
    // are emitted in reference order once the parallel loop is done. Writing
    // straight to the file from inside the loop (as this used to) made the block
    // order follow thread completion, so two builds of the same index produced
    // the same signatures in a different line order -- a real reproducibility
    // gap, invisible to any check that sorted the file first.
    std::vector<std::ostringstream> perRef(uniqStats.size());

    size_t totalSel = 0, totalSafe = 0;
    #pragma omp parallel for schedule(dynamic) reduction(+:totalSel, totalSafe)
    for (unsigned i = 0; i < uniqStats.size(); i++) {
        std::istringstream uniqstm(uniqStats[i]);
        std::string sketchPath;
        size_t sketchCount, sketchRefId;
        uniqstm >> sketchPath >> sketchCount >> sketchRefId;
        SketchStat st = getUniqSet(cfg, sketchPath, sketchRefId, sketchFilter, perRef[i], dbFilter);
        totalSel += st.selected;
        totalSafe += st.safe;
    }

    std::ofstream out(cfg.outfile);
    for (const auto& buf : perRef) out << buf.str();
    out.close();

    if (cfg.minMargin >= 2) {
        std::cout << "Signatures satisfying margin >= " << cfg.minMargin << ": "
                  << totalSafe << " / " << totalSel;
        if (totalSel) {
            std::cout << " (" << std::setprecision(1) << std::fixed
                      << 100.0 * static_cast<double>(totalSafe) / totalSel << "%)";
        }
        std::cout << "\n";
    }
}

#endif // UNIQSKETCHUTIL_HPP_