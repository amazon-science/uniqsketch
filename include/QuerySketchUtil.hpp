#ifndef QUERYSKETCHUTIL_HPP_
#define QUERYSKETCHUTIL_HPP_

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <memory>
#include <condition_variable>
#include <deque>
#include <mutex>
#include <numeric>
#include <thread>
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

    // Reporting thresholds, flat query
    double breadth = 0.0;                               // --breadth: min observed/expected signatures hit (0 = off)

    // Tree query
    int treePruneHit = 1;                               // min coarse-level hits to keep a sub-tree in play
    bool noGroupSignatures = false;                     // --no-group-signatures: leaf level uses global signatures only
    double treeBreadth = 0.5;                           // --tree-breadth: min observed/expected signatures hit (0 = off)
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
    std::string line;
    // Fields: seq, refId, refName, numRef, contig, pos. Split on whitespace by
    // hand: an istringstream per line cost 1.2 us, which on a 120-genome tree
    // index (1.5 M signatures over five files) was 1.8 s of the query.
    std::vector<std::pair<size_t, size_t>> field(6);
    while (std::getline(sketchFile, line)) {
        size_t n = 0, i = 0;
        while (n < 6) {
            while (i < line.size() && (line[i] == '\t' || line[i] == ' ')) ++i;
            if (i >= line.size()) break;
            size_t b = i;
            while (i < line.size() && line[i] != '\t' && line[i] != ' ') ++i;
            field[n++] = {b, i - b};
        }
        auto number = [&](size_t f, unsigned& out) {
            char* end = nullptr;
            unsigned long v = std::strtoul(line.c_str() + field[f].first, &end, 10);
            if (end != line.c_str() + field[f].first + field[f].second) return false;
            out = static_cast<unsigned>(v);
            return true;
        };
        unsigned refId = 0, refNum = 0, pos = 0;
        if (n < 6 || !number(1, refId) || !number(3, refNum) || !number(5, pos)) {
            std::cerr << "Error in parsing uniqsketch entry:\n" << line << "\n";
            exit(EXIT_FAILURE);
        }
        std::string uniqSeq = line.substr(field[0].first, field[0].second);
        if (init) {
            cfg.kmerLen = uniqSeq.length();
            sketchRef.resize(refNum);
            refSigCount.resize(refNum);
            init = false;
        }
        sketchRef[refId].assign(line, field[2].first, field[2].second);
        refSigCount[refId].emplace(uniqSeq, 0);
        sketchHash.emplace(std::move(uniqSeq), refId);
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
 * The query engine, shared by the flat and the tree query.
 *
 * Signatures are stored under their canonical ntHash value, the same value the
 * read iterator produces at every k-mer position, so a position is tested
 * against every signature in the index with one 64-bit probe and no string
 * work. On a probe hit the stored sequence is compared to the read's k-mer
 * before anything counts, so a hash collision can never produce a false hit:
 * results are identical to a string-keyed lookup, only the misses got cheap,
 * and in a typical query more than 99.99% of positions are misses.
 *
 * A hit counts toward a *target* and, optionally, a per-signature *slot*. The
 * flat query has one target per reference and one slot per signature. The tree
 * query puts every level's sketch and the within-group signatures into one index
 * (targets are (level, reference) pairs and individual group signatures), reads
 * the sample once, and applies the level-by-level narrowing to the tallies
 * afterwards. It used to run one pass over the reads per level: on a 120-genome
 * index that was six parses and six hashings of the same 3 M reads for one
 * result. The same k-mer may be a signature at several levels; every entry that
 * matches the read's k-mer counts, which is exactly what the separate passes
 * added up to.
 */
struct QueryEntry {
    const std::string* seq; // canonical signature, for collision verification; owned by the caller
    unsigned target;        // what a hit counts toward
    unsigned slot;          // per-signature counter, or NO_SLOT
};

static const unsigned NO_SLOT = UINT32_MAX;

/**
 * The index does not own its sequences: `seq` points into the caller's
 * signature maps (std::unordered_map keys never move) or containers with stable
 * addresses, which must outlive the index. On a 120-genome tree index that is
 * 1.8 M signatures of 81 bases; one copy instead of four is ~170 MB.
 *
 * Layout: the entries sorted by hash, and an open-addressing table from hash to
 * the run of entries carrying it. A probe touches one 16-byte slot, and for the
 * miss that almost every k-mer position is, that is the whole cost. The
 * std::unordered_map of vectors this replaced touched a bucket, a node and a
 * vector's storage, and took 1.7 M allocations and three seconds to build for
 * that index; this builds in a sort.
 */
struct QueryIndex {
    struct Slot {
        uint64_t hash;
        uint32_t first;         // index of the first entry with this hash
        uint32_t count;         // number of entries with this hash; 0 marks an empty slot
    };
    std::vector<QueryEntry> entries;    // sorted by hash after finalize()
    std::vector<Slot> slots;            // size is a power of two, at most half full
    uint64_t mask = 0;
    // One-hash Bloom filter over the keys, indexed by the hash's high bits (the
    // table uses the low ones). At 16+ bits per key it answers ~94% of the
    // misses from a table small enough to stay in cache (4 MB for a 120-genome
    // tree index), so the common case never touches the 64 MB slot table.
    std::vector<uint64_t> prefilter;
    unsigned prefilterShift = 64;
    unsigned numTargets = 0;
    unsigned numSlots = 0;

    /** Register a signature. Targets and slots are numbered by the caller. */
    void add(const std::string& seq, unsigned target, unsigned slot) {
        entries.push_back({&seq, target, slot});
    }

    /**
     * Sort the entries and build the table. Call once, after the last add().
     * Entries with equal hashes keep their add() order, so a duplicate line in
     * a sketch resolves the way the string-keyed maps resolved it.
     */
    void finalize() {
        std::vector<std::pair<uint64_t, uint32_t>> keyed(entries.size());
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < entries.size(); i++)
            keyed[i] = {NTC64(entries[i].seq->c_str(), static_cast<unsigned>(entries[i].seq->size())),
                        static_cast<uint32_t>(i)};
        std::sort(keyed.begin(), keyed.end());   // by hash, then by add() order

        std::vector<QueryEntry> sorted(entries.size());
        for (size_t i = 0; i < keyed.size(); i++) sorted[i] = entries[keyed[i].second];
        entries.swap(sorted);

        size_t n = 1;
        while (n < 2 * entries.size() + 1) n <<= 1;
        slots.assign(n, Slot{0, 0, 0});
        mask = n - 1;
        for (size_t i = 0; i < keyed.size();) {
            size_t j = i + 1;
            while (j < keyed.size() && keyed[j].first == keyed[i].first) ++j;
            size_t s = keyed[i].first & mask;
            while (slots[s].count) s = (s + 1) & mask;
            slots[s] = Slot{keyed[i].first, static_cast<uint32_t>(i), static_cast<uint32_t>(j - i)};
            i = j;
        }

        unsigned bits = 12;                                   // at least 4096 bits
        while ((size_t(1) << bits) < 16 * entries.size()) ++bits;
        prefilterShift = 64 - bits;
        prefilter.assign((size_t(1) << bits) / 64, 0);
        for (const auto& kv : keyed) {
            uint64_t b = kv.first >> prefilterShift;
            prefilter[b >> 6] |= uint64_t(1) << (b & 63);
        }
    }

    /**
     * The entries whose signature hashes to `h`: a pointer to the first and their
     * number, or nullptr. Linear probing; the table is at most half full.
     */
    const QueryEntry* find(uint64_t h, uint32_t& count) const {
        uint64_t b = h >> prefilterShift;
        if (!(prefilter[b >> 6] & (uint64_t(1) << (b & 63)))) return nullptr;
        size_t s = h & mask;
        while (slots[s].count) {
            if (slots[s].hash == h) {
                count = slots[s].count;
                return &entries[slots[s].first];
            }
            s = (s + 1) & mask;
        }
        return nullptr;
    }
};

/** What one pass over the reads produced. */
struct QueryTally {
    std::vector<unsigned> count;                            // per target: signature hits
    std::vector<unsigned> slotCount;                        // per slot: hits on that signature
    std::vector<std::unordered_set<std::string>> reads;     // per target: distinct read ids

    explicit QueryTally(const QueryIndex& index)
        : count(index.numTargets, 0), slotCount(index.numSlots, 0), reads(index.numTargets) {}
};

/**
 * One pass over the sample: every k-mer of every read against the index.
 *
 * Reads are parsed by a thread of their own per file into batches that the
 * OpenMP team consumes, so parsing and hashing overlap; before, the team sat
 * idle while the main thread parsed the next batch and then parsed nothing
 * while the team hashed, which capped a 3 M-read query at six times the
 * single-thread speed no matter how many threads it had.
 *
 * The team is one parallel region for the whole sample. Each thread takes a
 * batch from the queue, hashes it, and merges its hits under a lock. Earlier
 * each batch was its own parallel region: the main thread took the batch and
 * merged alone, and the other threads, with nothing to do, spun at the barrier
 * until the next batch arrived. With gzip input the readers deliver a batch
 * every few milliseconds and the hashing takes a fraction of that, so the spin
 * was most of the process: 436 s of CPU for a 25 s query at 32 threads (27 s
 * with OMP_WAIT_POLICY=passive). A thread that finds the queue empty now sleeps
 * on the condition variable like the readers do. The merge order follows
 * scheduling, as the per-thread merge order did before; the results it feeds
 * are counts and, per target, the *set* of read names, which do not depend on
 * it.
 *
 * The hot loop records a hit as (target, slot, read); the merge counts and
 * remembers the read's name once. The per-target sets of distinct read names
 * are built at the end, in parallel over targets, from those names: inserting
 * into the sets from the merge cost as much as the hashing did, and held the
 * dominant reference's names once per tree level.
 *
 * The read filter for --solid is allocated only when --solid is given; it is
 * sized for the sample's distinct k-mers (138 M on a 3 M-read E. coli sample,
 * 1.1 GB at 64 bits each) and used to be allocated, and zeroed, on every query.
 *
 * The --solid adjustment (a reference with any hit gets one more, for the first
 * occurrence the filter swallowed) is the caller's, because it applies per
 * reference and a tree target may be one signature of a reference.
 */
void queryReads(const QuerySketchConfig& cfg, const std::vector<std::string>& sampleFiles,
                const QueryIndex& index, QueryTally& tally) {
    static const size_t BATCH_SIZE = 4096;
    static const size_t QUEUE_DEPTH = 4;
    struct Hit { unsigned target, slot, read; };
    struct Batch {
        std::vector<std::string> seq, name;
        bool last = false;
    };

    std::unique_ptr<BloomFilter> solidFilter;
    if (cfg.solid) {
        solidFilter.reset(new BloomFilter(cfg.dbfSize * cfg.bits, cfg.nhash, cfg.kmerLen));
    }

    int maxThreads = 1;
#ifdef _OPENMP
    maxThreads = omp_get_max_threads();
#endif
    std::vector<std::vector<Hit>> hits(maxThreads);

    // Names of the reads that hit anything, each once, and per target the ids of
    // its hit reads (with runs of the same read collapsed).
    std::vector<std::string> hitNames;
    std::vector<std::vector<uint32_t>> readIds(index.numTargets);

    // Producers: one parser thread per input file (paired reads are two
    // independent streams), filling a bounded queue of batches. Consumers: the
    // OpenMP team. A batch marked `last` is the final one of its file; the
    // consumers stop once every file has delivered its last batch and the queue
    // is empty.
    std::mutex mu;
    std::condition_variable cvFull, cvEmpty;
    std::deque<std::unique_ptr<Batch>> ready, spare;
    size_t filesDone = 0;           // files whose last batch has been taken (under mu)
    // One batch in hand per consumer, plus QUEUE_DEPTH per reader to fill ahead.
    const size_t queueDepth = static_cast<size_t>(maxThreads) + QUEUE_DEPTH * sampleFiles.size();
    for (size_t i = 0; i < queueDepth; i++) spare.emplace_back(new Batch);

    std::vector<std::thread> readers;
    for (const std::string& path : sampleFiles) {
        readers.emplace_back([&, path] {
            gzFile fp = gzopen(path.c_str(), "r");
            if (fp) gzbuffer(fp, 1 << 20);
            kseq_t* ks = kseq_init(fp);
            bool eof = false;
            while (!eof) {
                std::unique_ptr<Batch> b;
                {
                    std::unique_lock<std::mutex> lk(mu);
                    cvEmpty.wait(lk, [&] { return !spare.empty(); });
                    b = std::move(spare.front());
                    spare.pop_front();
                }
                b->seq.clear();
                b->name.clear();
                for (size_t i = 0; i < BATCH_SIZE; ++i) {
                    if (kseq_read(ks) < 0) {
                        eof = true;
                        break;
                    }
                    b->seq.emplace_back(ks->seq.s, ks->seq.l);
                    b->name.emplace_back(ks->name.s, ks->name.l);
                }
                b->last = eof;
                {
                    std::lock_guard<std::mutex> lk(mu);
                    ready.push_back(std::move(b));
                }
                cvFull.notify_one();
            }
            kseq_destroy(ks);
            gzclose(fp);
        });
    }

    std::mutex mergeMu;
    #pragma omp parallel
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        std::vector<Hit>& mine = hits[tid];
        std::vector<uint32_t> nameId;   // per read of the current batch: its hitNames index
        for (;;) {
            std::unique_ptr<Batch> b;
            {
                std::unique_lock<std::mutex> lk(mu);
                cvFull.wait(lk, [&] { return !ready.empty() || filesDone == sampleFiles.size(); });
                if (ready.empty()) break;
                b = std::move(ready.front());
                ready.pop_front();
                if (b->last && ++filesDone == sampleFiles.size()) cvFull.notify_all();
            }
            const std::vector<std::string>& seqBatch = b->seq;
            mine.clear();

            for (size_t r = 0; r < seqBatch.size(); ++r) {
                const std::string& seq = seqBatch[r];
                if (seq.length() < cfg.kmerLen) continue;

                ntHashIterator itr(seq, cfg.nhash, cfg.kmerLen);
                while (itr != itr.end()) {
                    if (solidFilter && solidFilter->insert_make_change(*itr)) {
                        ++itr;
                        continue;
                    }
                    // One 64-bit probe decides the (overwhelmingly common) miss.
                    uint32_t n = 0;
                    if (const QueryEntry* e = index.find((*itr)[0], n)) {
                        const char* kp = seq.c_str() + itr.get_pos();
                        for (uint32_t i = 0; i < n; i++) {
                            if (kmerMatches(kp, *e[i].seq)) {
                                mine.push_back({e[i].target, e[i].slot, static_cast<unsigned>(r)});
                            }
                        }
                    }
                    ++itr;
                }
            }

            // Merge this batch: counts, and the read behind each hit, named once.
            if (!mine.empty()) {
                nameId.assign(seqBatch.size(), UINT32_MAX);
                std::lock_guard<std::mutex> lk(mergeMu);
                unsigned prevTarget = UINT32_MAX, prevRead = UINT32_MAX;
                for (const Hit& h : mine) {
                    ++tally.count[h.target];
                    if (h.slot != NO_SLOT) ++tally.slotCount[h.slot];
                    if (h.target != prevTarget || h.read != prevRead) {
                        if (nameId[h.read] == UINT32_MAX) {
                            nameId[h.read] = static_cast<uint32_t>(hitNames.size());
                            hitNames.push_back(std::move(b->name[h.read]));
                        }
                        readIds[h.target].push_back(nameId[h.read]);
                        prevTarget = h.target;
                        prevRead = h.read;
                    }
                }
            }

            {
                std::lock_guard<std::mutex> lk(mu);
                spare.push_back(std::move(b));
            }
            cvEmpty.notify_one();
        }
    }
    for (std::thread& t : readers) t.join();

    // Distinct read names per target. Two reads with the same name (mates named
    // alike, say) are one read here, as they were when names went straight into
    // the set.
    #pragma omp parallel for schedule(dynamic)
    for (size_t t = 0; t < readIds.size(); t++) {
        std::vector<uint32_t>& ids = readIds[t];
        if (ids.empty()) continue;
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        tally.reads[t].reserve(ids.size());
        for (uint32_t id : ids) tally.reads[t].insert(hitNames[id]);
        ids = std::vector<uint32_t>();
    }
}

/**
 * Flat query: every reference of one sketch against the sample.
 *
 * @param sampleFiles  All input read file paths.
 * @param sketchHash  Signature -> refId map, as loaded by loadSketch.
 * @param sketchCount  Per-reference signature hit counts (added to).
 * @param refSigCount  Per-reference signature -> hit count maps (added to).
 * @param refRead  Per-reference distinct read id sets (added to).
 */
void querySampleBatch(const QuerySketchConfig& cfg, const std::vector<std::string>& sampleFiles,
                      const SketchHash& sketchHash,
                      std::vector<unsigned>& sketchCount,
                      std::vector<SketchHash>& refSigCount,
                      std::vector<std::unordered_set<std::string>>& refRead) {
    // One target per reference; one slot per signature, pointing at its counter
    // in refSigCount so the pass writes straight into the maps the logs print.
    QueryIndex index;
    index.numTargets = static_cast<unsigned>(refSigCount.size());
    std::vector<unsigned*> slotCounter;
    for (const auto& [sig, refId] : sketchHash) {
        auto it = refSigCount[refId].find(sig);
        // A signature present in the sketch but absent from refSigCount is still
        // counted against its reference, without a per-signature counter.
        unsigned slot = NO_SLOT;
        if (it != refSigCount[refId].end()) {
            slot = static_cast<unsigned>(slotCounter.size());
            slotCounter.push_back(&it->second);
        }
        index.add(sig, refId, slot);
    }
    index.numSlots = static_cast<unsigned>(slotCounter.size());
    index.finalize();

    QueryTally tally(index);
    queryReads(cfg, sampleFiles, index, tally);

    for (unsigned ref = 0; ref < index.numTargets; ref++) {
        sketchCount[ref] += tally.count[ref];
        if (!tally.reads[ref].empty()) {
            if (refRead[ref].empty()) refRead[ref].swap(tally.reads[ref]);
            else refRead[ref].insert(tally.reads[ref].begin(), tally.reads[ref].end());
        }
    }
    for (unsigned s = 0; s < index.numSlots; s++) *slotCounter[s] += tally.slotCount[s];

    // Adjust for the first occurrence the --solid filter swallowed.
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
 * Breadth of the evidence behind a call: the number of distinct signatures hit,
 * relative to what the hit count would spread over if the whole genome were
 * present.
 *
 * Reads of a present genome land all over it, so with lambda = count / sigs_total
 * hits per signature about sigs_total * (1 - exp(-lambda)) distinct signatures
 * are hit. Reads of one locus of it, a mobile element the sample carries in some
 * other organism, pile their hits onto the few signatures in that locus. The
 * ratio observed / expected is near 1 in the first case and far below it in the
 * second, and the count alone cannot tell them apart: 717 reads of a human gut
 * metagenome with no Escherichia in it covered a 785-bp insertion of one E. coli
 * genome at 137x, hit its 80 signatures there 2 732 times, and would have been
 * reported as that strain at abundance 1. Over the 1 978 true strain rows of the
 * four validation sweeps the ratio was never below 0.73 (1 % quantile 0.88 or
 * better); that call scored 0.035, the family-level erm(B) call that began the
 * repeat investigation would have scored 0.02. A read-limited call is safe: at
 * lambda << 1 the expectation is the count itself, and distinct reads hit
 * distinct signatures, so the ratio is about 1 however few reads there are.
 *
 * @return observed / expected distinct signatures hit, or 1 when nothing was hit
 *         (the caller's count threshold decides that case).
 */
inline double breadthRatio(size_t count, size_t sigsHit, size_t sigsTotal) {
    if (count == 0 || sigsHit == 0 || sigsTotal == 0) return 1.0;
    const double lambda = static_cast<double>(count) / static_cast<double>(sigsTotal);
    const double expected = static_cast<double>(sigsTotal) * (1.0 - std::exp(-lambda));
    return expected > 0.0 ? static_cast<double>(sigsHit) / expected : 1.0;
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

    // Which leaves clear the reporting thresholds: hits, reads, and the breadth
    // of the hits over the signatures (breadthRatio). A row that clears the first
    // two and fails the third is a single locus of the genome, not the genome;
    // it is withheld and named on stderr.
    std::vector<size_t> leafSigsHit(sketchCount.size(), 0);
    std::vector<char> leafReported(sketchCount.size(), 0);
    struct Withheld { std::string name; unsigned count; size_t reads, sigsHit, sigsTotal; double ratio; };
    std::vector<Withheld> withheld;
    for (unsigned i = 0; i < sketchCount.size(); i++) {
        if (refRead[i].size() <= static_cast<size_t>(cfg.readCutoff) ||
            sketchCount[i] < static_cast<unsigned>(cfg.sketchHit)) continue;
        leafSigsHit[i] = countSigsHit(refSigCount[i]);
        const double ratio = breadthRatio(sketchCount[i], leafSigsHit[i], refSigCount[i].size());
        if (ratio < cfg.treeBreadth) {
            withheld.push_back({sketchRef[i], sketchCount[i], refRead[i].size(), leafSigsHit[i], refSigCount[i].size(), ratio});
            continue;
        }
        leafReported[i] = 1;
    }
    std::unordered_set<std::string> reportedLeaves;
    for (unsigned i = 0; i < sketchCount.size(); i++) {
        if (leafReported[i]) reportedLeaves.insert(sketchRef[i]);
    }

    // A family is only worth reporting when it is itself confident and none of the
    // leaves it covers was resolved; otherwise the leaf rows already say more.
    // The breadth floor does not apply here: a node's signatures are its
    // representative's node-unique k-mers, which sit in that genome's accessory
    // sequence, and a member carries only a part of them (36 % for one E. coli
    // node, a median 13 % across the 5 493-genome tree), so a present member
    // gives its node a low-breadth pattern that cannot be told from one locus.
    // The repeat exclusion in the group scan is what guards family rows against
    // mobile sequence.
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

    // One row per group. A node spanning several levels is scored at each of
    // them with a different signature draw and yields one call per level; keep
    // the finest. And when a group and one of its sub-groups are both called,
    // the sub-group already says more, so the enclosing group is dropped.
    {
        std::unordered_map<std::string, const FamilyCall*> finest;
        for (const FamilyCall* f : reportedFamilies) {
            auto it = finest.find(f->name);
            if (it == finest.end() || f->level > it->second->level) finest[f->name] = f;
        }
        std::vector<const FamilyCall*> kept;
        for (const auto& [name, f] : finest) {
            bool hasReportedSubgroup = false;
            for (const auto& [otherName, g] : finest) {
                if (g == f || g->coveredLeaves.size() >= f->coveredLeaves.size()) continue;
                bool subset = true;
                for (const std::string& leaf : g->coveredLeaves)
                    if (!f->coveredLeaves.count(leaf)) { subset = false; break; }
                if (subset) { hasReportedSubgroup = true; break; }
            }
            if (!hasReportedSubgroup) kept.push_back(f);
        }
        std::sort(kept.begin(), kept.end(),
                  [](const FamilyCall* a, const FamilyCall* b) { return a->count > b->count; });
        reportedFamilies.swap(kept);
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
        if (leafReported[i]) {
            rows.push_back({sketchRef[i], sketchCount[i], refRead[i].size(),
                            leafSigsHit[i], refSigCount[i].size(), false});
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
    if (!withheld.empty()) {
        std::cerr << "Tree query: " << withheld.size() << " strain call(s) withheld, hits confined to a few signatures "
                     "(breadth below " << cfg.treeBreadth << "; --tree-breadth=0 to report them):\n";
        for (const Withheld& w : withheld) {
            std::cerr << "  " << w.name << "\tcount " << w.count << "\treads " << w.reads
                      << "\tsigs_hit " << w.sigsHit << "/" << w.sigsTotal
                      << "\tbreadth " << std::setprecision(3) << w.ratio << "\n";
        }
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

    auto clearsThresholds = [&](unsigned i) {
        return refRead[i].size() > static_cast<size_t>(cfg.readCutoff) &&
               sketchCount[i] >= static_cast<unsigned>(cfg.sketchHit);
    };

    // A reference that clears the hit and read thresholds is still only a call if
    // its hits are spread over its signatures. When they sit on a handful of them
    // the evidence is one locus of the genome, not the genome: a window the
    // reference alone appears to carry because another reference is linearized at
    // a different origin, a divergent assembly window, a mobile element. Such a
    // locus can collect hundreds of thousands of hits from a sample that does not
    // hold the genome at all, and neither the count nor the read total shows it.
    // breadthRatio compares the distinct signatures hit against the number the
    // count would have spread over were the genome present, so a call carried by
    // few reads is not penalised for covering few signatures: at low counts the
    // expectation is the count itself and the ratio stays near 1 however little
    // evidence there is. Withheld references are named on stderr.
    // Off by default (--breadth=0): a floor changes which rows a flat query
    // reports, and the validation sweeps are what justify a non-zero default.
    struct Withheld { std::string name; unsigned count; size_t reads, sigsHit, sigsTotal; double ratio; };
    std::vector<Withheld> withheld;
    std::vector<char> reportedRef(sketchCount.size(), 0);
    std::vector<size_t> refSigsHit(sketchCount.size(), 0);
    for (unsigned i = 0; i < sketchCount.size(); i++) {
        if (!clearsThresholds(i)) continue;
        refSigsHit[i] = countSigsHit(refSigCount[i]);
        if (cfg.breadth > 0.0) {
            const double ratio = breadthRatio(sketchCount[i], refSigsHit[i], refSigCount[i].size());
            if (ratio < cfg.breadth) {
                withheld.push_back({sketchRef[i], sketchCount[i], refRead[i].size(),
                                    refSigsHit[i], refSigCount[i].size(), ratio});
                continue;
            }
        }
        reportedRef[i] = 1;
    }
    auto reported = [&](unsigned i) { return reportedRef[i] != 0; };

    // Abundance is each reference's hits per signature relative to all reported
    // references; see hitsPerSignature. A withheld reference is left out, as the
    // tree writer leaves out a withheld leaf, so the column still sums to 1 over
    // the rows that reach the file.
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
                            << "\t" << refSigsHit[idx]
                            << "\t" << refSigCount[idx].size()
                            << "\tstrain\n";
                }
            }
        }
    }
    outFile.close();

    if (!withheld.empty()) {
        std::cerr << withheld.size() << " call(s) withheld, hits confined to a few signatures "
                     "(breadth below " << cfg.breadth << "; --breadth=0 to report them):\n";
        for (const Withheld& w : withheld) {
            std::cerr << "  " << w.name << "\tcount " << w.count << "\treads " << w.reads
                      << "\tsigs_hit " << w.sigsHit << "/" << w.sigsTotal
                      << "\tbreadth " << std::setprecision(3) << w.ratio << "\n";
        }
    }

    writeSignatureLog(cfg, sketchRef, refSigCount, refRead);
    writeReadLog(cfg, sketchRef, refRead);
}

#endif // QUERYSKETCHUTIL_HPP_