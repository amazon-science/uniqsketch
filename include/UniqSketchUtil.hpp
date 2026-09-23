#ifndef UNIQSKETCHUTIL_HPP_
#define UNIQSKETCHUTIL_HPP_

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <numeric>
#include <random>
#include <sstream>
#include <unordered_map>
#include <vector>

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
 * Walk a FASTA file contig by contig: fn(header, sequence). `header` is the
 * '>' line without the '>' and with spaces replaced by '_', the form the
 * sketch's contig column uses; `sequence` is the concatenation of the lines up
 * to the next '>' line, passed by non-const reference so a caller may move it.
 * Every scanner of a reference goes through this, so a position recorded by
 * one of them means the same base to all of them. The line handling is the one
 * the scanners had when each carried its own copy of the loop, kept exactly:
 * a header line of just ">" leaves the previous header in force, and an empty
 * line inside a record is skipped.
 */
template <typename Fn>
inline void scanFasta(const std::string& fPath, Fn fn) {
    std::ifstream refFile(fPath);
    std::string line;
    bool good = static_cast<bool>(std::getline(refFile, line));
    std::string header = line.empty() ? std::string() : line.substr(1);
    std::replace(header.begin(), header.end(), ' ', '_');
    while (good) {
        std::string faSeq;
        good = static_cast<bool>(std::getline(refFile, line));
        while (good && line[0] != '>') {
            faSeq += line;
            good = static_cast<bool>(std::getline(refFile, line));
        }
        fn(header, faSeq);
        if (line.size() > 1) {
            header = line.substr(1);
            std::replace(header.begin(), header.end(), ' ', '_');
        }
    }
}

/** Candidate file of a reference: <outdir>/<basename without extension>.cand. */
inline std::string candidatePath(const std::string& outdir, const std::string& refPath) {
    return outdir + "/" + getBaseId(refPath) + ".cand";
}

/**
 * Per-reference candidate file: the reference's unique k-mers, as positions in
 * its own FASTA file, one run per contig in file order.
 *
 *   "UQSKCAND"  u32 version  u32 k  u32 pathLen  path             header
 *   u32 nameLen  name  u64 seqLen  u64 count  u32 pos[count]       per contig
 *
 * `path` is the reference FASTA as given to uniqsketch; `name` is the contig
 * header as scanFasta reports it; `seqLen` lets a reader notice a genome that
 * changed under it. A candidate's sequence is the canonical k-mer at `pos` of
 * its contig, which the reader recomputes for the candidates it examines, a
 * few thousand of the millions in the file. Native byte order: the file is
 * scratch for the build that writes it.
 *
 * These used to be text, one line per k-mer with the sequence spelled out and
 * the contig header repeated: 162 bytes each, 22 GB written, read back, parsed
 * and sorted for a 120-genome tree index, about a quarter of the build and a
 * hundred megabytes of strings per thread during selection. This is 4 bytes
 * per k-mer.
 */
class CandidateWriter {
public:
    static constexpr const char* MAGIC = "UQSKCAND";
    static constexpr uint32_t FORMAT_VERSION = 1;

    CandidateWriter(const std::string& path, unsigned k, const std::string& refPath)
        : out_(path, std::ios::binary | std::ios::trunc) {
        out_.write(MAGIC, 8);
        put32(FORMAT_VERSION);
        put32(k);
        put32(static_cast<uint32_t>(refPath.size()));
        out_.write(refPath.data(), static_cast<std::streamsize>(refPath.size()));
    }

    void beginContig(const std::string& name, uint64_t seqLen) {
        put32(static_cast<uint32_t>(name.size()));
        out_.write(name.data(), static_cast<std::streamsize>(name.size()));
        put64(seqLen);
        countAt_ = out_.tellp();
        put64(0);                       // patched by endContig
        count_ = 0;
    }

    void add(uint32_t pos) {
        buf_.push_back(pos);
        ++count_;
        if (buf_.size() == BUF_ENTRIES) flush();
    }

    void endContig() {
        flush();
        std::streampos end = out_.tellp();
        out_.seekp(countAt_);
        put64(count_);
        out_.seekp(end);
    }

    bool good() const { return out_.good(); }

private:
    static constexpr size_t BUF_ENTRIES = size_t(1) << 14;

    void put32(uint32_t v) { out_.write(reinterpret_cast<const char*>(&v), sizeof v); }
    void put64(uint64_t v) { out_.write(reinterpret_cast<const char*>(&v), sizeof v); }
    void flush() {
        if (!buf_.empty())
            out_.write(reinterpret_cast<const char*>(buf_.data()),
                       static_cast<std::streamsize>(buf_.size() * sizeof(uint32_t)));
        buf_.clear();
    }

    std::ofstream out_;
    std::vector<uint32_t> buf_;
    std::streampos countAt_{};
    uint64_t count_ = 0;
};

struct CandidateRun {
    std::string name;                   // contig header, as scanFasta reports it
    uint64_t seqLen = 0;                // contig length when the candidates were written
    std::vector<uint32_t> pos;          // positions of the unique k-mers, ascending
};

struct CandidateSet {
    unsigned k = 0;
    std::string refPath;                // the reference FASTA the positions index into
    std::vector<CandidateRun> runs;     // one per contig, in file order
    size_t total = 0;
};

/** Read a candidate file. False, with a message on stderr, if it is not one. */
inline bool readCandidates(const std::string& path, CandidateSet& out) {
    std::ifstream in(path, std::ios::binary);
    auto get32 = [&](uint32_t& v) { return static_cast<bool>(in.read(reinterpret_cast<char*>(&v), sizeof v)); };
    auto get64 = [&](uint64_t& v) { return static_cast<bool>(in.read(reinterpret_cast<char*>(&v), sizeof v)); };
    auto getStr = [&](std::string& s) {
        uint32_t n = 0;
        if (!get32(n)) return false;
        s.resize(n);
        return n == 0 || static_cast<bool>(in.read(&s[0], n));
    };
    auto bad = [&](const char* why) {
        std::cerr << "Error: " << path << ": " << why << "\n";
        return false;
    };

    char magic[8];
    uint32_t version = 0, k = 0;
    if (!in.read(magic, 8) || std::memcmp(magic, CandidateWriter::MAGIC, 8) != 0)
        return bad("not a uniqsketch candidate file");
    if (!get32(version) || version != CandidateWriter::FORMAT_VERSION)
        return bad("unsupported candidate file version");
    if (!get32(k) || !getStr(out.refPath)) return bad("truncated header");
    out.k = k;
    out.runs.clear();
    out.total = 0;
    while (in.peek() != std::char_traits<char>::eof()) {
        CandidateRun run;
        uint64_t count = 0;
        if (!getStr(run.name) || !get64(run.seqLen) || !get64(count)) return bad("truncated contig record");
        run.pos.resize(count);
        if (count && !in.read(reinterpret_cast<char*>(run.pos.data()),
                              static_cast<std::streamsize>(count * sizeof(uint32_t))))
            return bad("truncated position list");
        out.total += count;
        out.runs.push_back(std::move(run));
    }
    return true;
}

/**
 * Insert k-mers from a FASTA reference into distinct and solid Bloom filters.
 *
 * @param fPath  Path to a reference file.
 * @param dbFilter  Distinct Bloom filter for tracking distinct k-mers.
 * @param sbFilter  Solid Bloom filter for tracking repeat k-mers.
 */
void loadBFfa(const UniqSketchConfig& cfg, const std::string& fPath, BloomFilter& dbFilter, BloomFilter& sbFilter) {
    scanFasta(fPath, [&](const std::string&, const std::string& faSeq) {
        ntHashIterator itr(faSeq, std::max(cfg.nhash1, cfg.nhash2), cfg.kmerLen);
        while (itr != itr.end()) {
            if (!dbFilter.insert_make_change(*itr)) {
                sbFilter.insert(*itr);
            }
            ++itr;
        }
    });
}

/**
 * Scan a reference and record its unique k-mers: every k-mer not in the solid
 * filter, as a position in the reference's candidate file (see CandidateWriter).
 *
 * @param fPath  Path to a reference file.
 * @param refId  Integer id assigned to the reference.
 * @param sbFilter  Solid Bloom filter for tracking repeat k-mers.
 * @param refCount  Vector of unique k-mer counts per reference.
 */
void checkRef(const UniqSketchConfig& cfg, const std::string& fPath, int refId,
              const BloomFilter& sbFilter, std::vector<int>& refCount) {
    const std::string candPath = candidatePath(cfg.outdir, fPath);
    CandidateWriter cand(candPath, cfg.kmerLen, fPath);
    size_t countUnique = 0;
    scanFasta(fPath, [&](const std::string& header, const std::string& faSeq) {
        cand.beginContig(header, faSeq.size());
        ntHashIterator itr(faSeq, cfg.nhash2, cfg.kmerLen);
        while (itr != itr.end()) {
            if (!sbFilter.contains(*itr)) {
                ++countUnique;
                cand.add(static_cast<uint32_t>(itr.get_pos()));
            }
            ++itr;
        }
        cand.endContig();
    });
    if (!cand.good()) std::cerr << "Error: could not write " << candPath << "\n";
    refCount[refId] = static_cast<int>(countUnique);
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
        tsvOut << candidatePath(cfg.outdir, refFiles[i]) << "\t"
               << refCount[i] << "\t" << i << "\n";
    }
    tsvOut.close();

    std::cout << "Distinct BF actual fpr: " << std::setprecision(4) << std::fixed
              << std::pow(static_cast<double>(dbFilter.get_pop()) / cfg.m1(), cfg.nhash1) << "\n"
              << "Solid BF actual fpr: " << std::setprecision(4) << std::fixed
              << std::pow(static_cast<double>(sbFilter.get_pop()) / cfg.m2(), cfg.nhash2) << "\n";
}

/**
 * Insert the k-mers of a group of genomes into the distinct/solid filters,
 * counting each k-mer at most once for the whole group.
 *
 * A tree node stands for several genomes. Its signatures should be k-mers that
 * occur in the node's members and nowhere else, and a k-mer shared by all the
 * members is the best kind. Fed through loadBFfa one file at a time, such a
 * k-mer is seen twice and lands in the solid filter as a repeat. Here a
 * group-private "seen" filter admits each k-mer to the global filters once, so
 * the solid filter ends up holding k-mers present in two or more *groups*, and
 * "distinct but not solid" means unique to one group.
 *
 * The seen filter is sized generously (64 bits per base of input, five hashes:
 * false-positive rate ~2e-6) because a false positive here is not conservative:
 * it drops the group's only insertion of a k-mer, so a k-mer this group shares
 * with one other group is counted once and handed to that group as a signature.
 */
void loadBFfaGroup(const UniqSketchConfig& cfg, const std::vector<std::string>& files,
                   BloomFilter& dbFilter, BloomFilter& sbFilter) {
    if (files.size() == 1) {
        loadBFfa(cfg, files[0], dbFilter, sbFilter);
        return;
    }
    size_t bases = 0;
    for (const auto& f : files) {
        std::ifstream in(f, std::ios::binary | std::ios::ate);
        if (in) bases += static_cast<size_t>(in.tellg());
    }
    if (bases == 0) bases = 1;
    BloomFilter seen(bases * 64, cfg.nhash1, cfg.kmerLen);
    for (const auto& fPath : files) {
        scanFasta(fPath, [&](const std::string&, const std::string& faSeq) {
            ntHashIterator itr(faSeq, std::max(cfg.nhash1, cfg.nhash2), cfg.kmerLen);
            while (itr != itr.end()) {
                if (seen.insert_make_change(*itr)) {
                    if (!dbFilter.insert_make_change(*itr)) {
                        sbFilter.insert(*itr);
                    }
                }
                ++itr;
            }
        });
    }
}

/**
 * identifyUniqKmers for a universe of genome groups.
 *
 * Each group is a tree-level entry: its first file is the representative whose
 * k-mers are the candidates, the rest are the other members of the node (none
 * for a leaf standing for itself). A candidate is a k-mer of the representative
 * that occurs in no other group, members included (see loadBFfaGroup). Built
 * over the complete partition of the leaves at a level, this makes every
 * coarse-level signature unique against every genome in the index, not just
 * against the other representatives.
 *
 * Output layout matches identifyUniqKmers: one candidate file per group named
 * after the representative, and a stat file buildSketch can consume.
 */
void identifyGroupUniqKmers(const UniqSketchConfig& cfg,
                            const std::vector<std::vector<std::string>>& groups,
                            BloomFilter& dbFilter) {
    size_t members = 0;
    for (const auto& g : groups) members += g.size();
    std::cout << "k-mer length: " << cfg.kmerLen << "\n"
              << "Groups: " << groups.size() << " over " << members << " genomes\n";

    BloomFilter sbFilter(cfg.m2(), cfg.nhash2, cfg.kmerLen);

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < groups.size(); ++i) {
        loadBFfaGroup(cfg, groups[i], dbFilter, sbFilter);
    }

    std::vector<int> refCount(groups.size());
    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < groups.size(); ++i) {
        checkRef(cfg, groups[i][0], i, sbFilter, refCount);
    }

    std::ofstream tsvOut(cfg.refstat);
    for (size_t i = 0; i < refCount.size(); i++) {
        tsvOut << candidatePath(cfg.outdir, groups[i][0]) << "\t"
               << refCount[i] << "\t" << i << "\n";
    }
    tsvOut.close();
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
 * The k-mer is hashed once; each of the 3k neighbors' hashes is then derived
 * from it in O(1) with NTMC64sub (ntHash is an XOR of per-base terms), so a
 * candidate costs one O(k) hash plus 3k Bloom probes instead of 3k O(k) hashes.
 *
 * @param sig  Candidate signature (length cfg.kmerLen).
 * @param dbFilter  Distinct Bloom filter over all universe k-mers.
 * @return True if no Hamming-1 neighbor is present in the universe.
 */
bool isTwoSafe(const UniqSketchConfig& cfg, const std::string& sig, const BloomFilter& dbFilter) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    uint64_t fhVal = 0, rhVal = 0;
    NTC64(sig.c_str(), cfg.kmerLen, fhVal, rhVal);
    uint64_t hVal[cfg.nhash1];
    for (unsigned i = 0; i < cfg.kmerLen; ++i) {
        const char orig = sig[i];
        for (char b : BASES) {
            if (b == orig) continue;
            NTMC64sub(fhVal, rhVal, i, orig, b, cfg.kmerLen, cfg.nhash1, hVal);
            if (dbFilter.contains(hVal)) return false;
        }
    }
    return true;
}

/**
 * Extract the unique k-mer set for a given reference.
 *
 * @param fPath  Path to a reference's candidate file (see CandidateWriter).
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
    CandidateSet cand;
    if (!readCandidates(fPath, cand)) return {0, 0};
    if (cand.total == 0) return {0, 0};
    if (cand.k != cfg.kmerLen) {
        std::cerr << "Error: " << fPath << " was written for k=" << cand.k
                  << ", not k=" << cfg.kmerLen << "\n";
        return {0, 0};
    }

    // The contigs, split as checkRef split them. A genome that changed since
    // the candidates were written is refused rather than mis-read.
    std::vector<std::string> contigs;
    std::vector<std::string> names;
    scanFasta(cand.refPath, [&](const std::string& header, std::string& seq) {
        names.push_back(header);
        contigs.push_back(std::move(seq));
    });
    bool same = contigs.size() == cand.runs.size();
    for (size_t r = 0; same && r < cand.runs.size(); r++)
        same = contigs[r].size() == cand.runs[r].seqLen && names[r] == cand.runs[r].name;
    if (!same) {
        std::cerr << "Error: " << cand.refPath << " does not match the candidates in " << fPath
                  << " (changed since the build started?)\n";
        return {0, 0};
    }
    names = std::vector<std::string>();

    // Candidates in the order selection walks them: by contig name, as the
    // sketch's contig column compares (byte order, shorter first on a tie),
    // then by position. With distinct names that is the runs in name order.
    // Two contigs with the same header interleave by position; that case
    // starts from file order and sorts, which is what always happened.
    const size_t nRuns = cand.runs.size();
    std::vector<uint32_t> runOrder(nRuns);
    std::iota(runOrder.begin(), runOrder.end(), 0u);
    auto nameLess = [&](uint32_t a, uint32_t b) {
        const std::string& x = cand.runs[a].name;
        const std::string& y = cand.runs[b].name;
        int c = std::strncmp(x.c_str(), y.c_str(), std::min(x.size(), y.size()));
        if (c != 0) return c < 0;
        return x.size() < y.size();
    };
    std::stable_sort(runOrder.begin(), runOrder.end(), nameLess);
    bool duplicateNames = false;
    for (size_t i = 1; i < nRuns; i++)
        if (!nameLess(runOrder[i - 1], runOrder[i])) duplicateNames = true;

    struct Cand { uint32_t run, pos; };
    std::vector<Cand> cands;
    cands.reserve(cand.total);
    if (!duplicateNames) {
        for (uint32_t r : runOrder)
            for (uint32_t p : cand.runs[r].pos) cands.push_back({r, p});
    } else {
        std::vector<uint32_t> rank(nRuns);
        for (size_t i = 0, k = 0; i < nRuns; i++) {
            if (i && nameLess(runOrder[i - 1], runOrder[i])) ++k;
            rank[runOrder[i]] = static_cast<uint32_t>(k);
        }
        for (uint32_t r = 0; r < nRuns; r++)
            for (uint32_t p : cand.runs[r].pos) cands.push_back({r, p});
        std::sort(cands.begin(), cands.end(), [&](const Cand& a, const Cand& b) {
            if (rank[a.run] != rank[b.run]) return rank[a.run] < rank[b.run];
            return a.pos < b.pos;
        });
    }
    for (CandidateRun& r : cand.runs) r.pos = std::vector<uint32_t>();

    // The contig column of the sketch: the header's first whitespace-delimited
    // token, which is what parsing the text candidate lines yielded.
    std::vector<std::string> contigToken(nRuns);
    for (size_t r = 0; r < nRuns; r++) {
        const std::string& n = cand.runs[r].name;
        size_t b = 0;
        while (b < n.size() && std::isspace(static_cast<unsigned char>(n[b]))) ++b;
        size_t e = b;
        while (e < n.size() && !std::isspace(static_cast<unsigned char>(n[e]))) ++e;
        contigToken[r] = n.substr(b, e - b);
    }
    const std::string baseId = getBaseId(fPath);

    // Select evenly spaced candidates, checking low-complexity (and, when
    // cfg.minMargin >= 2, the cross-reference Hamming margin) lazily. A
    // candidate's sequence is computed only here, when its slot is examined.
    const size_t total = cands.size();
    int count = 0;
    size_t safeCount = 0;
    size_t needed = std::min(static_cast<size_t>(cfg.sketchnum), total);
    double step = static_cast<double>(total) / needed;
    const bool useMargin = (cfg.minMargin >= 2);

    // Commit the first usable candidate in [a, b). When requireSafe is set, only
    // 2-safe candidates are accepted. Returns true if a signature was committed.
    auto trySlot = [&](size_t a, size_t b, bool requireSafe) -> bool {
        for (size_t i = a; i < b; i++) {
            const Cand& c = cands[i];
            std::string useq = get_canonical(contigs[c.run].substr(c.pos, cfg.kmerLen));

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
                    out << useq << "\t" << refId << "\t" << baseId
                        << "\t" << cfg.numRef << "\t" << contigToken[c.run] << "\t" << c.pos << "\n";
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