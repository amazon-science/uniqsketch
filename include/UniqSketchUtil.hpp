#ifndef UNIQSKETCHUTIL_HPP_
#define UNIQSKETCHUTIL_HPP_

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "SequenceUtil.hpp"
#include "FileUtil.hpp"
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
    bool circular = false;                              // --circular: references are circular molecules

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
inline void scanFasta(const std::string& fPath, Fn fn, unsigned padK = 0) {
    GzLineReader refFile(fPath);          // plain or gzip
    std::string line;
    bool good = refFile.getline(line);
    std::string header = line.empty() ? std::string() : line.substr(1);
    std::replace(header.begin(), header.end(), ' ', '_');
    while (good) {
        std::string faSeq;
        good = refFile.getline(line);
        while (good && line[0] != '>') {
            faSeq += line;
            good = refFile.getline(line);
        }
        // A circular molecule rendered linearly loses the k-1 k-mers that span
        // its end-to-start join: they are in the genome but not in the file.
        // Appending the first k-1 bases restores exactly those k-mers and no
        // others, so a contig of length L yields L k-mers, which is what a
        // circular sequence of length L has. Their start positions are all <= L,
        // so the candidate positions stay in the file's own coordinate frame and
        // nothing downstream needs to translate them. The true length is passed
        // separately so the padding cannot leak into a genome length.
        const size_t trueLen = faSeq.size();
        if (padK != 0 && trueLen != 0) faSeq += faSeq.substr(0, std::min<size_t>(padK, trueLen));
        fn(header, faSeq, trueLen);
        if (line.size() > 1) {
            header = line.substr(1);
            std::replace(header.begin(), header.end(), ' ', '_');
        }
    }
}

/** Wrap padding for a build: k-1 bases with --circular, none without. */
inline unsigned wrapPad(const UniqSketchConfig& cfg) {
    return cfg.circular && cfg.kmerLen > 1 ? cfg.kmerLen - 1 : 0;
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
    scanFasta(fPath, [&](const std::string&, const std::string& faSeq, size_t) {
        ntHashIterator itr(faSeq, std::max(cfg.nhash1, cfg.nhash2), cfg.kmerLen);
        while (itr != itr.end()) {
            if (!dbFilter.insert_make_change(*itr)) {
                sbFilter.insert(*itr);
            }
            ++itr;
        }
    }, wrapPad(cfg));
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
    size_t countUnique = 0, countAll = 0;
    scanFasta(fPath, [&](const std::string& header, const std::string& faSeq, size_t trueLen) {
        // The contig length recorded is the genome's, not the padded sequence's,
        // so a wrap pad can never be mistaken for sequence by a reader of the
        // candidate file or by the window layout in getUniqSet.
        cand.beginContig(header, trueLen);
        ntHashIterator itr(faSeq, cfg.nhash2, cfg.kmerLen);
        while (itr != itr.end()) {
            ++countAll;
            if (!sbFilter.contains(*itr)) {
                ++countUnique;
                cand.add(static_cast<uint32_t>(itr.get_pos()));
            }
            ++itr;
        }
        cand.endContig();
    }, wrapPad(cfg));
    // A reference that contributes nothing is almost always an input problem
    // (not FASTA, empty, sequences shorter than k) or a duplicate of another
    // reference. Either way it gets no signatures and cannot be detected, which
    // used to pass in silence. One message each, built first: this runs inside
    // a parallel loop.
    if (!cand.good()) {
        std::cerr << "Error: could not write " + candPath + "\n";
    } else if (countAll == 0) {
        std::cerr << "Warning: " + fPath + ": no k-mers of length " + std::to_string(cfg.kmerLen) +
                     " (empty, not FASTA, or every sequence shorter than k); it will get no signatures\n";
    } else if (countUnique == 0) {
        std::cerr << "Warning: " + fPath + ": no unique k-mers (identical to another reference?); "
                     "it will get no signatures\n";
    }
    refCount[refId] = static_cast<int>(countUnique);
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
 * The seen filter is sized generously (64 bits per element it can hold, five
 * hashes: false-positive rate ~2e-6; seenFilterBits) because a false positive
 * here is not conservative: it drops the group's only insertion of a k-mer, so
 * a k-mer this group shares with one other group is counted once and handed to
 * that group as a signature.
 *
 * Counting once per group must not also hide a k-mer that occurs twice inside a
 * single genome. loadBFfa marks such a k-mer solid on its second sighting and a
 * leaf therefore never gets one as a signature; left to the seen filter alone
 * the second sighting never reaches the solid filter and the repeat survives as
 * a group signature. Multi-copy k-mers sit in mobile elements, so they are the
 * sequence most likely to be present in organisms outside the index entirely: a
 * two-copy erm(B) in one representative took ~250 reads of a human gut
 * metagenome, from no Escherichia at all, and turned them into a family call.
 * Candidates are drawn from the representative only (see identifyGroupUniqKmers
 * and the checkRef call there), so the representative is the one genome whose
 * internal multiplicity can reach a signature, and repSeen restores loadBFfa's
 * rule for it. A false positive in repSeen costs a signature rather than
 * admitting a bad one, the opposite of the seen filter's exposure. Marking the
 * k-mer solid cannot harm another group: one present in two groups is already
 * solid by way of the seen/distinct path above.
 */

/**
 * Size of a group's seen filter, in bits: 64 per element it can be asked to
 * hold. The elements are the group's distinct k-mers, which are no more than
 * the group's bases and no more than the distinct k-mers of the whole index
 * (cfg.dbfSize, ntCard's estimate, which every tree level scans against), so
 * the smaller of the two sizes it; a configuration without a cardinality
 * (dbfSize 0, as in a unit test) keeps the base count. Sized by bases alone,
 * the filter of the largest node of a 5 493-genome E. coli index, 896 genomes
 * of 5.1 Mb, was 36.6 GB and the build's memory peak (42.6 GB), where the
 * index has 332 M distinct k-mers in all: 2.7 GB. The bound binds for a node
 * whose bases exceed the universe's distinct k-mers (on the 120-genome index
 * its 20-genome node alone; 25 of 1 143 nodes at 5 493 genomes), and there it
 * moves the filter's false positives, as any resizing does; the rate stays at
 * or below the 64-bit one.
 */
inline size_t seenFilterBits(const UniqSketchConfig& cfg, size_t bases) {
    size_t elems = std::max<size_t>(bases, 1);
    if (cfg.dbfSize > 0) elems = std::min(elems, cfg.dbfSize);
    return elems * 64;
}

void loadBFfaGroup(const UniqSketchConfig& cfg, const std::vector<std::string>& files,
                   BloomFilter& dbFilter, BloomFilter& sbFilter) {
    if (files.size() == 1) {
        loadBFfa(cfg, files[0], dbFilter, sbFilter);
        return;
    }
    size_t bases = 0;
    for (const auto& f : files) bases += uncompressedSize(f);
    BloomFilter seen(seenFilterBits(cfg, bases), cfg.nhash1, cfg.kmerLen);
    BloomFilter repSeen(std::max<size_t>(uncompressedSize(files[0]), 1) * 64,
                        cfg.nhash1, cfg.kmerLen);
    for (size_t fi = 0; fi < files.size(); ++fi) {
        const bool isRep = (fi == 0);
        scanFasta(files[fi], [&](const std::string&, const std::string& faSeq, size_t) {
            ntHashIterator itr(faSeq, std::max(cfg.nhash1, cfg.nhash2), cfg.kmerLen);
            while (itr != itr.end()) {
                if (isRep && !repSeen.insert_make_change(*itr)) {
                    sbFilter.insert(*itr);
                }
                if (seen.insert_make_change(*itr)) {
                    if (!dbFilter.insert_make_change(*itr)) {
                        sbFilter.insert(*itr);
                    }
                }
                ++itr;
            }
        }, wrapPad(cfg));
    }
}

/**
 * Scan groups of genomes into a distinct/solid filter pair in a fixed order,
 * with the result the serial scan would produce, on any number of threads.
 *
 * The serial scan is loadBFfaGroup on each group in turn (loadBFfa for a group
 * of one). Running the groups one per thread was simple but not equivalent:
 * which k-mers a Bloom filter mistakes for present depends on which bits were
 * set before them, so the filters' decisions, and through them the candidates
 * and the signatures, depend on the order of insertion, and groups scanned
 * concurrently have no fixed order. Measured on the 120-genome E. coli tree
 * build, one run in three differed from the other two in one reference's
 * candidate count (one k-mer), which selection turned into a thousand or two
 * redrawn signatures because the slot width is total / needed; -t1 and -t32
 * differed in two to five references per level. It also left a level waiting
 * on its largest group: the top level's 20-genome cluster took 11 s of a 12 s
 * phase on one thread while 31 threads finished the other 58 groups in 3 s.
 *
 * Here the k-mers are split by hash instead. Every filter is blocked
 * (BlockedBloomFilter): a k-mer touches one 512-bit block per filter, chosen by
 * its leading hash h0, and block counts are multiples of BloomFilter::PARTITIONS,
 * so the k-mers with h0 in the p-th slice of the hash space own the p-th slice
 * of the blocks of every filter, shared or group-private. Contigs are read in
 * file order and hashed in parallel pieces, each piece's h0 values sorted into
 * their partitions in position order; then one thread per partition replays its
 * k-mers piece by piece, deriving the remaining hashes from h0 as ntHashIterator
 * does (NTE64). Every block therefore sees exactly the sequence of operations
 * the serial scan performs on it, from a single writer, so the stripe locks are
 * skipped (insert_make_change_exclusive) and the filters come out bit for bit
 * the same whatever the thread count. A flush at the end of a group's
 * representative keeps the "second sighting in the representative" rule on
 * exactly the representative's k-mers. Each thread's slice of the filters is
 * 1/PARTITIONS of their size, a few megabytes on a bacterial set, so the replay
 * also runs near cache speed where the one-genome-per-thread scan paid a memory
 * round trip per probe.
 *
 * With the standard filter a k-mer's bits are spread over the whole array,
 * partitions isolate nothing, and the groups are scanned one per thread as
 * before: fast, but with the order dependence above.
 *
 * The scratch (one 8-byte h0 per k-mer of a flush, 256 MB at the default flush
 * size) belongs to the object so that repeated loads, as in the within-group
 * build, reuse it.
 */
class GroupScanner {
public:
    explicit GroupScanner(const UniqSketchConfig& cfg, size_t flushBases = size_t(32) << 20)
        : cfg_(cfg), flushBases_(flushBases) {}

    /**
     * Insert the k-mers of every group, in list order, into dbFilter/sbFilter.
     * Nothing else may touch the two filters during the call.
     */
    void load(const std::vector<std::vector<std::string>>& groups,
              BloomFilter& dbFilter, BloomFilter& sbFilter) {
        if (!BloomFilter::kPartitioned) {
            loadUnordered(groups, dbFilter, sbFilter);
            return;
        }
        for (const auto& files : groups) loadGroup(files, dbFilter, sbFilter);
    }

private:
    static constexpr unsigned P = BloomFilter::PARTITIONS;
    static constexpr size_t PIECE_KMERS = size_t(1) << 18;     // hashing work unit

    /** Standard filter: one group per thread, largest first. */
    void loadUnordered(const std::vector<std::vector<std::string>>& groups,
                       BloomFilter& dbFilter, BloomFilter& sbFilter) {
        std::vector<unsigned> order(groups.size());
        std::iota(order.begin(), order.end(), 0u);
        std::stable_sort(order.begin(), order.end(), [&](unsigned a, unsigned b) {
            return groups[a].size() > groups[b].size();
        });
        #pragma omp parallel for schedule(dynamic)
        for (size_t j = 0; j < order.size(); ++j) {
            loadBFfaGroup(cfg_, groups[order[j]], dbFilter, sbFilter);
        }
    }

    void loadGroup(const std::vector<std::string>& files, BloomFilter& dbFilter, BloomFilter& sbFilter) {
        const unsigned k = cfg_.kmerLen;
        // Group-private filters, sized as loadBFfaGroup sizes them; none for a
        // group of one, which loadBFfa scans without them.
        std::unique_ptr<BloomFilter> seen, repSeen;
        if (files.size() > 1) {
            size_t bases = 0;
            for (const auto& f : files) bases += uncompressedSize(f);
            seen.reset(new BloomFilter(seenFilterBits(cfg_, bases), cfg_.nhash1, k));
            repSeen.reset(new BloomFilter(std::max<size_t>(uncompressedSize(files[0]), 1) * 64,
                                          cfg_.nhash1, k));
        }
        for (size_t fi = 0; fi < files.size(); ++fi) {
            const bool isRep = (fi == 0);
            scanFasta(files[fi], [&](const std::string&, std::string& faSeq, size_t) {
                pendingBases_ += faSeq.size();
                contigs_.push_back(std::move(faSeq));
                if (pendingBases_ >= flushBases_) flush(isRep, seen.get(), repSeen.get(), dbFilter, sbFilter);
            }, wrapPad(cfg_));
            if (isRep || fi + 1 == files.size()) flush(isRep, seen.get(), repSeen.get(), dbFilter, sbFilter);
        }
    }

    /**
     * Hash the pending contigs and replay their k-mers. `seen` is null for a
     * group of one (loadBFfa's rule applies); otherwise `isRep` says whether the
     * pending contigs are the representative's, to which repSeen applies.
     */
    void flush(bool isRep, BloomFilter* seen, BloomFilter* repSeen,
               BloomFilter& dbFilter, BloomFilter& sbFilter) {
        const unsigned k = cfg_.kmerLen;
        const unsigned nh = std::max(cfg_.nhash1, cfg_.nhash2);
        // Pieces of at most PIECE_KMERS k-mer start positions, in file order. Each
        // owns a region of h0_ with one slot per position (fewer are used where
        // the sequence has N), holding its k-mers partition by partition.
        struct Piece { size_t contig, begin, end, base; };
        std::vector<Piece> pieces;
        size_t slots = 0;
        for (size_t c = 0; c < contigs_.size(); c++) {
            if (contigs_[c].size() < k) continue;
            const size_t n = contigs_[c].size() - k + 1;
            for (size_t b = 0; b < n; b += PIECE_KMERS) {
                const size_t e = std::min(n, b + PIECE_KMERS);
                pieces.push_back({c, b, e, slots});
                slots += e - b;
            }
        }
        if (h0_.size() < slots) h0_.resize(slots);
        off_.assign(pieces.size() * (P + 1), 0);   // off_[i*(P+1)+p]: partition p's start in piece i

        #pragma omp parallel
        {
            std::vector<uint64_t> inOrder;                 // one piece, position order
            inOrder.reserve(PIECE_KMERS);
            uint32_t count[P], cursor[P];
            #pragma omp for schedule(dynamic)
            for (size_t i = 0; i < pieces.size(); i++) {
                const Piece& pc = pieces[i];
                inOrder.clear();
                ntHashIterator itr(contigs_[pc.contig].substr(pc.begin, pc.end - pc.begin + k - 1), 1, k);
                while (itr != itr.end()) {
                    inOrder.push_back((*itr)[0]);
                    ++itr;
                }
                std::memset(count, 0, sizeof count);
                for (uint64_t h0 : inOrder) ++count[BloomFilter::partition(h0)];
                uint32_t* off = &off_[i * (P + 1)];
                off[0] = 0;
                for (unsigned p = 0; p < P; p++) off[p + 1] = off[p] + count[p];
                std::memcpy(cursor, off, sizeof cursor);
                uint64_t* region = h0_.data() + pc.base;
                for (uint64_t h0 : inOrder) region[cursor[BloomFilter::partition(h0)]++] = h0;
            }
        }
        // One thread per partition: its blocks of every filter, nobody else's.
        #pragma omp parallel for schedule(dynamic)
        for (unsigned p = 0; p < P; p++) {
            uint64_t hVal[nh];
            for (size_t i = 0; i < pieces.size(); i++) {
                const uint32_t* off = &off_[i * (P + 1)];
                const uint64_t* from = h0_.data() + pieces[i].base + off[p];
                const uint64_t* to = h0_.data() + pieces[i].base + off[p + 1];
                for (; from != to; ++from) {
                    const uint64_t h0 = *from;
                    hVal[0] = h0;
                    for (unsigned j = 1; j < nh; j++) hVal[j] = NTE64(h0, k, j);
                    if (seen) {
                        if (isRep && !repSeen->insert_make_change_exclusive(hVal)) {
                            sbFilter.insert_exclusive(hVal);
                        }
                        if (seen->insert_make_change_exclusive(hVal)) {
                            if (!dbFilter.insert_make_change_exclusive(hVal)) {
                                sbFilter.insert_exclusive(hVal);
                            }
                        }
                    } else if (!dbFilter.insert_make_change_exclusive(hVal)) {
                        sbFilter.insert_exclusive(hVal);
                    }
                }
            }
        }
        contigs_.clear();
        pendingBases_ = 0;
    }

    const UniqSketchConfig& cfg_;
    size_t flushBases_;
    std::vector<std::string> contigs_;      // pending, in file order
    size_t pendingBases_ = 0;
    std::vector<uint64_t> h0_;              // a flush's leading hashes, by piece then partition
    std::vector<uint32_t> off_;             // partition offsets within each piece's region
};

/** GroupScanner::load for a single call site. */
void loadBFfaGroups(const UniqSketchConfig& cfg, const std::vector<std::vector<std::string>>& groups,
                    BloomFilter& dbFilter, BloomFilter& sbFilter,
                    size_t flushBases = size_t(32) << 20) {
    GroupScanner(cfg, flushBases).load(groups, dbFilter, sbFilter);
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

    // Every reference is a group of one; GroupScanner says why they are not
    // simply scanned one per thread.
    std::vector<std::vector<std::string>> singles;
    singles.reserve(refFiles.size());
    for (const std::string& f : refFiles) singles.push_back({f});
    loadBFfaGroups(cfg, singles, dbFilter, sbFilter);

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
 * identifyUniqKmers for a universe of genome groups.
 *
 * Each group is a tree-level entry: its first file is the representative whose
 * k-mers are the candidates, the rest are the other members of the node (none
 * for a leaf standing for itself). A candidate is a k-mer that occurs once in
 * the representative and in no other group, members included (see
 * loadBFfaGroup). Built over the complete partition of the leaves at a level,
 * this makes every coarse-level signature unique against every genome in the
 * index, not just against the other representatives.
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

    loadBFfaGroups(cfg, groups, dbFilter, sbFilter);

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
 * @param out  Stream receiving this reference's signatures. buildSketch passes a
 *             per-reference buffer so output order is independent of thread
 *             scheduling.
 * @param dbFilter  Distinct Bloom filter over the whole universe (margin check).
 * @return Number of signatures selected and how many met the margin.
 */
SketchStat getUniqSet(const UniqSketchConfig& cfg, const std::string& fPath, unsigned refId,
                      std::ostream& out, const BloomFilter& dbFilter) {
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
    std::vector<size_t> trueLens;
    scanFasta(cand.refPath, [&](const std::string& header, std::string& seq, size_t trueLen) {
        names.push_back(header);
        trueLens.push_back(trueLen);
        contigs.push_back(std::move(seq));
    }, wrapPad(cfg));
    bool same = contigs.size() == cand.runs.size();
    for (size_t r = 0; same && r < cand.runs.size(); r++)
        same = trueLens[r] == cand.runs[r].seqLen && names[r] == cand.runs[r].name;
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
    if (cands.empty()) return {0, 0};

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

    // Selection: spread the signatures along the genome.
    //
    // The genome (its contigs laid end to end in the order above) is cut into
    // `quota` windows of equal length, and the first usable candidate of every
    // window is taken. Windows left empty, because they hold no candidate or no
    // candidate passes the margin, do not hand their place to anyone; the quota
    // that remains goes to the windows that did supply one, each in proportion
    // to its candidates and capped at a common level chosen so that the shares
    // add up to the quota (pass 2 below). The cap binds only when the supplied
    // windows' shares exceed the quota, that is for a reference that reaches
    // it: such a reference gets a flat density, one signature per window when
    // the windows outnumber the quota. Below the quota the cap is idle and a
    // window's share is its candidate share, which is the old allocation with
    // the slots the old grid left empty now filled: on a 120-genome E. coli
    // index at -c5000 no reference reaches the quota, the sketch grows from
    // 352 797 to 373 001 signatures, and the share of each reference's
    // signatures in its densest 2-kb window is unchanged (median 0.027).
    //
    // Selection used to space its slots over the candidate list instead, one per
    // total/quota candidates, and fill rates differ by locus: in the core of a
    // bacterial genome most candidates are a substitution away from a relative
    // and fail the strict margin, while every position of a foreign insertion
    // is a candidate and every one passes. On that 120-genome index a 785-bp
    // insertion holding 3 % of one genome's candidates took 21 % of its
    // chromosomal signatures (191 per kb against 0.15 elsewhere), a 65-kb
    // plasmid took 72 % of the whole set, and 717 reads of a human gut
    // metagenome with no Escherichia in it hit 155 of those signatures and
    // called the strain. Mobile sequence is both the most candidate-dense and
    // the most likely to turn up in organisms outside the index, which is
    // the same mechanism as the family-level erm(B) call fixed in the group
    // scan. Spacing by position does not cure that below the quota (the
    // insertion keeps its 150 signatures; that genome has 2-safe candidates in
    // 63 of 5 000 windows, so any rule that honours the quota gives each of
    // them dozens), and the two rules that flattened regardless of the quota
    // were rejected on the sweep (pass 2 below); the query's breadth floor
    // (QuerySketchUtil.hpp, breadthRatio) is what closes the gut call. What
    // selection by position guarantees is a floor, one signature in every
    // window that has a usable candidate, and a flat layout for a reference at
    // its quota.
    //
    // Candidates are examined lazily, and only as far as needed: a window is
    // walked in position order until one candidate is usable (low-complexity
    // and homopolymer filters, then the margin when cfg.minMargin >= 2), so a
    // dense window costs one margin check and a sparse one at most its size.
    // In soft-margin mode a window with no 2-safe candidate takes its first
    // margin-1 one; strict mode leaves it empty.
    const size_t total = cands.size();
    const bool useMargin = (cfg.minMargin >= 2);
    const size_t quota = static_cast<size_t>(std::max(cfg.sketchnum, 1));

    // Global coordinate: contigs end to end in sorted order. Same-named contigs
    // share one span (their positions interleave in `cands`), so rank, not run,
    // indexes the offsets; with distinct names rank is the sorted position.
    std::vector<uint32_t> rankOf(nRuns);
    std::vector<uint64_t> rankLen;
    for (size_t i = 0; i < nRuns; i++) {
        if (i == 0 || nameLess(runOrder[i - 1], runOrder[i])) rankLen.push_back(0);
        rankOf[runOrder[i]] = static_cast<uint32_t>(rankLen.size() - 1);
        rankLen.back() = std::max<uint64_t>(rankLen.back(), cand.runs[runOrder[i]].seqLen);
    }
    std::vector<uint64_t> rankOffset(rankLen.size(), 0);
    for (size_t r = 1; r < rankLen.size(); r++) rankOffset[r] = rankOffset[r - 1] + rankLen[r - 1];
    const uint64_t genomeLen = std::max<uint64_t>(rankOffset.back() + rankLen.back(), 1);
    const size_t windows = static_cast<size_t>(std::min<uint64_t>(quota, genomeLen));
    auto windowOf = [&](const Cand& c) -> size_t {
        const uint64_t g = rankOffset[rankOf[c.run]] + c.pos;
        return static_cast<size_t>(std::min<uint64_t>(windows - 1, g * windows / genomeLen));
    };
    // cands is sorted by (rank, pos), so windows are contiguous ranges of it.
    std::vector<size_t> wBegin(windows + 1, total);
    for (size_t i = 0, w = 0; i < total; i++) {
        const size_t wi = windowOf(cands[i]);
        while (w <= wi) wBegin[w++] = i;
    }

    struct Pick { size_t index; std::string seq; bool safe; };
    std::vector<Pick> picks;
    picks.reserve(std::min(quota, total));
    std::vector<char> noSafeLeft(windows, 0);       // a full walk of the window found no 2-safe candidate

    // 0: unusable; 1: usable, margin 1 (or margin not required); 2: usable, 2-safe.
    auto examine = [&](size_t i, bool checkMargin, std::string& useq) -> int {
        const Cand& c = cands[i];
        useq = get_canonical(contigs[c.run].substr(c.pos, cfg.kmerLen));
        if (lowComplexity(cfg, useq) || hasLongHomopolymer(cfg, useq)) return 0;
        if (!checkMargin) return 1;
        return isTwoSafe(cfg, useq, dbFilter) ? 2 : 1;
    };
    // Take the first usable candidate in [a, b), a range inside window w, or
    // report that it has none. On the window's first call, in soft-margin mode,
    // a window with no 2-safe candidate yields its first margin-1 one; strict
    // mode leaves it empty. Either way a window once found to have no 2-safe
    // candidate is not probed for the margin again.
    auto takeIn = [&](size_t w, size_t a, size_t b) -> bool {
        const bool wantSafe = useMargin && !noSafeLeft[w];
        size_t fallback = total;
        std::string fallbackSeq, useq;
        for (size_t i = a; i < b; i++) {
            const int v = examine(i, wantSafe, useq);
            if (v == 2 || (v == 1 && !wantSafe)) {
                picks.push_back({i, useq, v == 2});
                return true;
            }
            if (v == 1 && fallback == total) { fallback = i; fallbackSeq = useq; }
        }
        if (wantSafe && a == wBegin[w] && b == wBegin[w + 1]) {
            noSafeLeft[w] = 1;
            if (!cfg.strictMargin && fallback != total) {
                picks.push_back({fallback, fallbackSeq, false});
                return true;
            }
        }
        return false;
    };

    // Pass 1: the first usable candidate of every window.
    std::vector<size_t> active;             // windows that supplied
    std::vector<size_t> firstPick(windows, total);
    for (size_t w = 0; w < windows && picks.size() < quota; w++) {
        if (wBegin[w] == wBegin[w + 1]) continue;
        if (takeIn(w, wBegin[w], wBegin[w + 1])) { active.push_back(w); firstPick[w] = picks.back().index; }
    }
    // Pass 2: the unfilled quota goes to the windows that supplied. A window's
    // share is its share of the candidates, quota * (its candidates / total),
    // capped at a common level chosen so that the shares add up to the quota.
    // The share is laid out as equal runs of the window's candidate list, one
    // pick per run (the first usable), which inside a window is what the old
    // slots did over the whole list.
    //
    // The candidate share alone was the old rule, and it concentrated: a
    // 785-bp foreign insertion that held 3 % of a genome's candidates took 21 %
    // of its chromosomal signatures, because every one of its positions is a
    // candidate and every one passes the margin while the core's mostly do not.
    // An equal share for every supplied window (tried first) concentrated the
    // other way: a window whose supply is one 74-bp cluster of
    // substitution-overlapping k-mers got the same dozen picks as a window with
    // candidates throughout, 12 to 17 signatures within one read length, where
    // the candidate share had given it 2 to 4. A read from a close relative
    // carrying the two coincident errors that make this genome's k-mers look
    // like its own hits every signature whose k-mer covers both sites, so a
    // dozen in one cluster is a single-read call: four false rows on the
    // 120-genome sweep that had none. Taking a window's share as its next
    // usable candidates (tried second) put the share at consecutive positions
    // (18 such rows). Capping the candidate share at a common level is the old
    // allocation for a reference below its quota, and a flat one for a
    // reference at it. The picks move either way, since a run is laid over a
    // window's candidates where a slot was laid over the whole list: on the
    // 8-phage index at -c500 the four phages that take their whole pool are
    // unchanged and the other four keep 1 to 101 of 500 (1 834 of 2 835
    // signatures change); on the 120-genome index 7 of 119 references are
    // unchanged and 24 % of the old signatures are kept. No row of either
    // 280-mixture sweep changes, nor the abundance error.
    if (!active.empty() && picks.size() < quota) {
        // Each window's candidate share, then the water level: the largest cap L
        // with sum(min(candShare, L)) <= quota, and the remainder as one more to
        // the first windows (in position order) still above L. The sum reaches
        // the quota whenever the candidate shares do, which they always do (they
        // are rounded up), so a reference falls short of the quota only when its
        // windows run out of usable candidates, never through rounding.
        std::vector<size_t> candShare(windows, 0);
        size_t maxShare = 1;
        for (size_t w : active) {
            const size_t nw = wBegin[w + 1] - wBegin[w];
            candShare[w] = std::max<size_t>(1, static_cast<size_t>(
                (static_cast<unsigned __int128>(nw) * quota + total - 1) / total));
            maxShare = std::max(maxShare, candShare[w]);
        }
        auto sumAt = [&](size_t L) {
            size_t s = 0;
            for (size_t w : active) s += std::min(candShare[w], L);
            return s;
        };
        size_t lo = 1, hi = maxShare;                     // largest L with sumAt(L) <= quota
        while (lo < hi) {
            const size_t mid = lo + (hi - lo + 1) / 2;
            if (sumAt(mid) <= quota) lo = mid; else hi = mid - 1;
        }
        const size_t level = lo;
        size_t remainder = quota - std::min(quota, sumAt(level));
        std::vector<size_t> share(windows, 0);
        for (size_t w : active) {
            share[w] = std::min(candShare[w], level);
            if (remainder && candShare[w] > level) { ++share[w]; --remainder; }
        }
        for (size_t sub = 1, more = true; more && picks.size() < quota; sub++) {
            more = false;
            for (size_t w : active) {
                if (picks.size() >= quota) break;
                if (sub >= share[w]) continue;
                more = true;
                const size_t wb = wBegin[w], nw = wBegin[w + 1] - wb;
                const size_t a = wb + sub * nw / share[w];
                const size_t b = wb + (sub + 1) * nw / share[w];
                if (b <= firstPick[w] || a >= b) continue;          // pass 1's run, or an empty one
                takeIn(w, std::max(a, firstPick[w] + 1), b);
            }
        }
    }

    // Position order, as the candidates are, so two builds of the same index
    // are comparable line by line.
    std::sort(picks.begin(), picks.end(), [](const Pick& a, const Pick& b) { return a.index < b.index; });
    size_t safeCount = 0;
    for (const Pick& p : picks) {
        const Cand& c = cands[p.index];
        safeCount += p.safe;
        // `out` is this reference's private buffer, so no lock is needed. No
        // claim filter across references either: a candidate is a k-mer the
        // solid filter does not contain, and a Bloom filter has no false
        // negatives, so a k-mer in two references, or twice in one, is a
        // candidate for neither (see GroupScanner for why the filter went).
        out << p.seq << "\t" << refId << "\t" << baseId
            << "\t" << cfg.numRef << "\t" << contigToken[c.run] << "\t" << c.pos << "\n";
    }
    return SketchStat{picks.size(), safeCount};
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

    // Each reference writes its signatures into its own buffer, and the buffers
    // are emitted in reference order once the parallel loop is done. Writing
    // straight to the file from inside the loop (as this used to) made the block
    // order follow thread completion, so two builds of the same index produced
    // the same signatures in a different line order -- a real reproducibility
    // gap, invisible to any check that sorted the file first.
    std::vector<std::ostringstream> perRef(uniqStats.size());

    // A reference that ends the build with no signatures can never be reported by
    // any query, at any threshold. The candidate scan already warns when the pool
    // itself is empty (an identical reference, or a file yielding no k-mers), but a
    // reference whose pool is non-empty and is then consumed entirely by the
    // entropy, homopolymer and margin filters used to pass in silence. On an
    // internal 40-reference set 14 references held no signatures and only 7 said
    // so; the 7 that did not had pools of 56 to 560 candidates, and one of them was
    // the dominant vector of a sample the screen then reported nothing for.
    size_t totalSel = 0, totalSafe = 0;
    std::vector<size_t> refSelected(uniqStats.size(), 0);
    std::vector<size_t> refCandidates(uniqStats.size(), 0);
    std::vector<std::string> refPath(uniqStats.size());
    #pragma omp parallel for schedule(dynamic) reduction(+:totalSel, totalSafe)
    for (unsigned i = 0; i < uniqStats.size(); i++) {
        std::istringstream uniqstm(uniqStats[i]);
        std::string sketchPath;
        size_t sketchCount, sketchRefId;
        uniqstm >> sketchPath >> sketchCount >> sketchRefId;
        SketchStat st = getUniqSet(cfg, sketchPath, sketchRefId, perRef[i], dbFilter);
        refSelected[i] = st.selected;
        refCandidates[i] = sketchCount;
        refPath[i] = sketchPath;
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

    size_t noSignature = 0;
    for (unsigned i = 0; i < refSelected.size(); i++) {
        if (refSelected[i] != 0) continue;
        noSignature++;
        // Pools that were empty to begin with were already named by the candidate
        // scan; naming them twice would bury the ones nobody has heard about.
        if (refCandidates[i] == 0) continue;
        std::cerr << "Warning: " << getBaseId(refPath[i]) << ": no signatures selected from "
                  << refCandidates[i] << " unique k-mer(s) (entropy, homopolymer or margin "
                     "filter consumed the pool); it can never be reported\n";
    }
    if (noSignature != 0) {
        std::cerr << noSignature << " of " << refSelected.size()
                  << " references hold no signatures and can never be reported";
        if (cfg.strictMargin || cfg.minMargin >= 2) {
            std::cerr << "; a lower --min-margin, dropping --strict-margin or a --tree index"
                         " (which adds within-group signatures and family-level calls) would reach them";
        }
        std::cerr << "\n";
    }
}

#endif // UNIQSKETCHUTIL_HPP_