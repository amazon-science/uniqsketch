#ifndef NODECARRIERUTIL_HPP_
#define NODECARRIERUTIL_HPP_

/**
 * Node signatures with carriers.
 *
 * A node's signatures are its representative's k-mers absent from every genome
 * outside the node. When the clustering splits a lineage across nodes (the
 * thresholds are unique-k-mer counts, and a few insertions or a deletion put a
 * close relative over one), the node's core is carried by relatives outside it and
 * nothing of the core is node-unique; what remains is the representative's private
 * sequence, which the other members do not carry. On 5 493 complete E. coli the
 * 102-leaf K-12 node (MG1655, BW25113 and their derivatives) had 1 407 such
 * candidates and 100 of its 102 members carried none of the signatures drawn from
 * them; 56 references of that index carried no signature of any kind at any level
 * and 10x of pure MG1655 yielded an empty table (docs/v2.0-validation.md, "Dense
 * lineages").
 *
 * The within-group mechanism (GroupSketchUtil.hpp) solved the same problem for
 * leaves: signatures unique within the group, each recorded with the outside
 * genomes that carry it, used only when no carrier is present. This header is
 * that idea one level up. Beside a node's candidate pool (its representative's
 * k-mers unique to the node) a second pool is written: the representative's k-mers
 * that occur once in it, are carried by at least half the node's members (and by one
 * besides the representative) and by at most `nodeCarriers` genomes outside the node.
 * The node's signatures are selected from the first pool exactly as before; only when
 * they fall short of the quota is the remainder selected from the second
 * (buildSketch), with the same margin, entropy and spread rules. The carriers of the
 * selected signatures are then recorded in the level sketch (seventh column,
 * comma-separated leaf ids or `-`), and the query counts a carried signature for its
 * node only when every carrier is absent (querysketch.cpp, before the walk: a carrier
 * is present when anything of its own was hit, excluded when its own signature loci
 * would have been hit at the coverage the node's hits imply, and otherwise may be
 * present and blocks the signature; the parent node, which contains the node and its
 * sibling carriers alike, reports the ambiguous evidence through its own signatures).
 *
 * Filling a shortfall rather than widening every pool was decided on the 120-genome
 * index: widened in proportion, with a cap of 12, every node came out 38 to 99 %
 * carried (a 275 000-candidate node selected 88 % of its signatures from the carried
 * pool), signatures that stop counting whenever one of about five outside genomes is
 * present, for nodes that members already reach. Every unreachable leaf of the
 * 5 493-genome index sits under nodes that end below quota after the margin filter
 * (the 896-leaf level-0 node selected 2 413 of 5 000; the K-12 node's pool is 1 407),
 * and on such a node filling and widening select the same signatures, since the
 * strict margin consumes nearly all of a dense lineage's carried core (on the
 * 120-genome index one node's 297 000 carried candidates yielded 3 223 signatures; on
 * the 5 493-genome one the K-12 node's 833 yielded none, every K-12 core k-mer having a
 * substitution neighbour among the 150 K-12 genomes, so that node stays
 * representative-private and its members are reached through the level-2 node's
 * 1 860 filled signatures). What filling cannot reach is a node at quota whose
 * signatures its members do not carry; neither index has one, and the reachability
 * report would show it (a leaf carrying none of a node's signatures).
 *
 * Counting the outside carriers of a k-mer needs, for each of a representative's
 * k-mers, the number of genomes of the index that contain it. An exact table over
 * every representative is too large (hundreds of nodes of 4.6 M k-mers per level), so
 * the count comes from a counting Bloom filter of per-genome carriage: one pass over
 * all genomes, each genome's distinct k-mers incrementing four 16-bit counters; the
 * estimate is the minimum of the four, an overestimate on collisions, so it can only
 * drop a candidate, never admit one with more carriers than counted, and it is a sum,
 * so it does not depend on the thread count. Sized at 16 counters per distinct k-mer
 * of the index (ntCard's estimate, already in hand), rounded up to a power of two so
 * the hash is masked rather than divided, it is 16 GB at 5 493 genomes and 2 GB at
 * 120, 32 to 64 bytes per distinct k-mer; measured on the 5 493-genome index at 2^31
 * counters (a quarter of that), 95 % of the K-12
 * representative's rare k-mers were counted exactly, none was over by more than 20,
 * and none with 50 or fewer carriers was estimated above 50. The in-node count (how
 * many members carry a representative k-mer) is exact: a hash of the representative's
 * k-mers, one pass over the members.
 */

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "UniqSketchUtil.hpp"
#include "TreeSketchUtil.hpp"
#include "FileUtil.hpp"

/**
 * Carriers recorded in a level sketch: seq -> leaf ids of the genomes outside the node
 * that contain it, for the lines whose seventh column is not `-`. Empty for a sketch
 * without the column (an index built with node carriers off, or before them).
 */
inline std::unordered_map<std::string, std::vector<unsigned>> loadLevelCarriers(const std::string& path) {
    std::unordered_map<std::string, std::vector<unsigned>> out;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        size_t p = 0;
        for (int f = 0; f < 6 && p != std::string::npos; f++) p = line.find('\t', p + (f ? 1 : 0));
        if (p == std::string::npos || p + 1 >= line.size() || line[p + 1] == '-') continue;
        std::vector<unsigned> c;
        const char* s = line.c_str() + p + 1;
        const char* end = line.c_str() + line.size();
        while (s < end) {
            char* q = nullptr;
            unsigned long v = std::strtoul(s, &q, 10);
            if (q == s) break;
            c.push_back(static_cast<unsigned>(v));
            s = (q < end && *q == ',') ? q + 1 : q;
        }
        if (!c.empty()) out.emplace(line.substr(0, line.find('\t')), std::move(c));
    }
    return out;
}

/** Default carrier cap for an index of n leaves: 50, or a tenth of the leaves if smaller. */
inline unsigned defaultNodeCarriers(size_t leaves) {
    return static_cast<unsigned>(std::min<size_t>(50, leaves / 10));
}

/**
 * Signature loci per leaf: the positions of a leaf's own signatures (leaf sketch, and
 * within-group signatures given as (refId, contig, pos)) clustered on each contig,
 * a new locus where the gap exceeds `gap` bases. The query uses this to ask how many
 * distinct reads a genome's own signatures would collect were it present at a given
 * coverage: a genome one to three substitutions from its relatives has its few
 * signatures packed into two or three sites (on the 5 493-genome E. coli index one
 * such strain's 146 signatures took 15 494 hits from 304 reads at 99 % abundance),
 * and the number of signatures says nothing about how often they are seen; the
 * number of loci does. The leaf sketch is read once more for the positions, which
 * loadSketch does not keep (hand-split, as there; 22 M lines at 5 493 genomes).
 */
inline std::vector<unsigned> countSignatureLoci(std::vector<std::vector<std::pair<uint64_t, uint32_t>>> pos,
                                                size_t nLeaves,
                                                const std::vector<std::tuple<unsigned, std::string, unsigned>>& extra,
                                                unsigned gap);
inline std::vector<unsigned> countSignatureLoci(const std::string& leafSketch, size_t nLeaves,
                                                const std::vector<std::tuple<unsigned, std::string, unsigned>>& extra,
                                                unsigned gap) {
    std::vector<std::vector<std::pair<uint64_t, uint32_t>>> pos(nLeaves);   // (contig hash, position)
    const std::hash<std::string> hs;
    std::ifstream in(leafSketch);
    std::string line;
    std::vector<std::pair<size_t, size_t>> field(6);
    while (std::getline(in, line)) {
        size_t n = 0, i = 0;
        while (n < 6) {
            while (i < line.size() && (line[i] == '\t' || line[i] == ' ')) ++i;
            if (i >= line.size()) break;
            size_t b = i;
            while (i < line.size() && line[i] != '\t' && line[i] != ' ') ++i;
            field[n++] = {b, i - b};
        }
        if (n < 6) continue;
        const unsigned long refId = std::strtoul(line.c_str() + field[1].first, nullptr, 10);
        if (refId >= nLeaves) continue;
        pos[refId].push_back({hs(line.substr(field[4].first, field[4].second)),
                              static_cast<uint32_t>(std::strtoul(line.c_str() + field[5].first, nullptr, 10))});
    }
    return countSignatureLoci(std::move(pos), nLeaves, extra, gap);
}

/**
 * countSignatureLoci on positions already in hand: per leaf the (contig hash,
 * position) of its leaf signatures, as loadSketch collects them when asked
 * (SignaturePositions), so the leaf sketch is not read a second time (1.5 GB and
 * 2.5 s of every query on the 5 493-genome index). The result is the same: the
 * positions are sorted before the loci are counted.
 */
inline std::vector<unsigned> countSignatureLoci(std::vector<std::vector<std::pair<uint64_t, uint32_t>>> pos,
                                                size_t nLeaves,
                                                const std::vector<std::tuple<unsigned, std::string, unsigned>>& extra,
                                                unsigned gap) {
    const std::hash<std::string> hs;
    pos.resize(std::max(pos.size(), nLeaves));
    for (const auto& [refId, contig, p] : extra)
        if (refId < nLeaves) pos[refId].push_back({hs(contig), p});
    std::vector<unsigned> loci(nLeaves, 0);
    #pragma omp parallel for schedule(dynamic, 64)
    for (size_t r = 0; r < nLeaves; r++) {
        auto& v = pos[r];
        std::sort(v.begin(), v.end());
        unsigned n = 0;
        for (size_t i = 0; i < v.size(); i++)
            if (i == 0 || v[i].first != v[i - 1].first || v[i].second - v[i - 1].second > gap) ++n;
        loci[r] = n;
    }
    return loci;
}

/**
 * Counting Bloom filter of per-genome k-mer carriage: how many genomes of the index
 * contain a k-mer, as the minimum of four saturating 16-bit counters.
 */
class CarriageCounter {
public:
    static constexpr unsigned HASHES = 4;
    /** @param elements  Distinct k-mers the filter will hold (16 counters each). */
    explicit CarriageCounter(size_t elements)
        : size_(roundUp(std::max<size_t>(elements, 1) * 16)), mask_(size_ - 1),
          counters_(new std::atomic<uint16_t>[size_]) {
        // Zeroed by every thread: the table is 16 counters per distinct k-mer of the
        // index (8.6 G of them on 5 493 E. coli), and one thread storing them held the
        // others for the first seconds of the carriage pass.
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < size_; i++) counters_[i].store(0, std::memory_order_relaxed);
    }
    size_t size() const { return size_; }
    size_t bytes() const { return size_ * sizeof(uint16_t); }
    /** Ask for the four counters add() or count() will touch for this k-mer (a hint, see forEachKmerAhead). */
    void prefetch(const uint64_t* hVal) const {
        for (unsigned h = 0; h < HASHES; h++) __builtin_prefetch(&counters_[hVal[h] & mask_]);
    }
    /** Count one genome for the k-mer with these HASHES hash values. */
    void add(const uint64_t* hVal) {
        for (unsigned h = 0; h < HASHES; h++) {
            std::atomic<uint16_t>& c = counters_[hVal[h] & mask_];
            uint16_t v = c.load(std::memory_order_relaxed);
            while (v != UINT16_MAX &&
                   !c.compare_exchange_weak(v, static_cast<uint16_t>(v + 1), std::memory_order_relaxed)) {}
        }
    }
    /** Genomes containing the k-mer, or more (never fewer). */
    unsigned count(const uint64_t* hVal) const {
        uint16_t m = UINT16_MAX;
        for (unsigned h = 0; h < HASHES; h++)
            m = std::min(m, counters_[hVal[h] & mask_].load(std::memory_order_relaxed));
        return m;
    }
private:
    static size_t roundUp(size_t n) {     // power of two, so the hash is masked, not divided
        size_t p = 1;
        while (p < n) p <<= 1;
        return p;
    }
    size_t size_, mask_;
    std::unique_ptr<std::atomic<uint16_t>[]> counters_;
};

/**
 * One pass over every genome into the carriage counter: each genome's distinct k-mers
 * once (a per-genome seen filter at 64 bits per base keeps a repeat within the genome
 * from counting twice). Genomes in parallel; the counter is a sum, so the result does
 * not depend on the order.
 *
 * Every k-mer is a block of the genome's seen filter and four counters spread over a
 * table of gigabytes, all of them out of cache, so the pass is bound by memory round
 * trips: the lines are asked for when the k-mer is hashed and it is counted sixteen
 * k-mers later, in the same order (forEachKmerAhead). The seen filter belongs to this
 * genome's thread alone, so the exclusive insert is the same test-and-set without
 * the stripe lock, which is shared by every filter and every thread of the process.
 */
inline void countCarriage(const UniqSketchConfig& cfg, const std::vector<std::string>& refFiles,
                          CarriageCounter& cc) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned g = 0; g < refFiles.size(); g++) {
        BloomFilter seen(std::max<size_t>(uncompressedSize(refFiles[g]), 1) * 64,
                         CarriageCounter::HASHES, cfg.kmerLen);
        scanFasta(refFiles[g], [&](const std::string&, const std::string& seq, size_t) {
            if (seq.size() < cfg.kmerLen) return;
            forEachKmerAhead(seq, CarriageCounter::HASHES, cfg.kmerLen,
                [&](const uint64_t* hVal) { seen.prefetch(hVal[0]); cc.prefetch(hVal); },
                [&](const uint64_t* hVal, size_t) { if (seen.insert_make_change_exclusive(hVal)) cc.add(hVal); });
        }, wrapPad(cfg));
    }
}

/**
 * A 64-bit hash -> uint32 table in one flat array: open addressing with linear
 * probing, the slot chosen by multiplying the key (mixed with an odd constant, see
 * slotOf) into the capacity, load kept at or below 0.6, an empty slot marked by the
 * value UINT32_MAX (the values are dense indices).
 *
 * std::unordered_map held the same mapping first, and freeing it was the largest
 * cost of the carried scan: a node is a malloc per entry, and freeing 7.5 G of them
 * (the 5 493-genome index hashes that many representative k-mers over its four
 * levels) ran on one thread at 85 ns each, 640 of the build's 2 080 extra seconds,
 * with inserts at 465 ns and lookups at 48 ns on top. One array is one free, and a
 * probe is one cache line (16 bytes per slot, 27 per k-mer at load 0.6, against 44).
 *
 * Reads are safe from any number of threads while no one writes; the member scans
 * only read. Writing is one thread per table (a representative is hashed by one
 * thread).
 */
class KmerTable {
public:
    static constexpr uint32_t EMPTY = UINT32_MAX;
    /** Room for `expected` keys without growing. */
    explicit KmerTable(size_t expected = 0) { reserve(expected); }
    void reserve(size_t expected) {
        size_t cap = std::max<size_t>(16, expected + expected * 2 / 3);   // load 0.6 when full
        if (cap <= slots_.size()) return;
        std::vector<Slot> old = std::move(slots_);
        slots_.assign(cap, Slot{0, EMPTY});
        size_ = 0;
        for (const Slot& t : old)
            if (t.val != EMPTY) *place(t.key) = t;
    }
    size_t size() const { return size_; }
    size_t capacity() const { return slots_.size(); }
    /** Ask for the slot find(key) or findOrInsert(key) will look at first (a hint, see forEachKmerAhead). */
    void prefetch(uint64_t key) const { __builtin_prefetch(&slots_[slotOf(key)]); }
    /** The value of `key`, or nullptr. */
    const uint32_t* find(uint64_t key) const {
        for (size_t i = slotOf(key);; i = (i + 1 == slots_.size()) ? 0 : i + 1) {
            const Slot& t = slots_[i];
            if (t.val == EMPTY) return nullptr;
            if (t.key == key) return &t.val;
        }
    }
    /**
     * The value of `key`, inserted as `value` if the key was absent (`inserted` says
     * which). The reference holds until the next insertion.
     */
    uint32_t& findOrInsert(uint64_t key, uint32_t value, bool& inserted) {
        if ((size_ + 1) * 5 > slots_.size() * 3) reserve(slots_.size());   // past load 0.6: double
        for (size_t i = slotOf(key);; i = (i + 1 == slots_.size()) ? 0 : i + 1) {
            Slot& t = slots_[i];
            if (t.val == EMPTY) { t.key = key; t.val = value; ++size_; inserted = true; return t.val; }
            if (t.key == key) { inserted = false; return t.val; }
        }
    }
private:
    struct Slot { uint64_t key; uint32_t val; };
    size_t slotOf(uint64_t key) const {
        // Multiply by an odd constant first: the slot comes from the product's high
        // bits, which then depend on every bit of the key. The keys are not always
        // uniform over 64 bits; the sampled k-mers of chooseNodeRepresentative are
        // the hashes whose top bits are zero, and multiplying those into the
        // capacity directly put every key in the first 1/64 of the table (a 20-member
        // node's choice took 6 s instead of 0.1).
        const uint64_t mixed = key * 0x9E3779B97F4A7C15ULL;
        return static_cast<size_t>((static_cast<unsigned __int128>(mixed) * slots_.size()) >> 64);
    }
    Slot* place(uint64_t key) {   // first empty slot for a key known to be absent
        size_t i = slotOf(key);
        while (slots_[i].val != EMPTY) i = (i + 1 == slots_.size()) ? 0 : i + 1;
        ++size_;
        return &slots_[i];
    }
    std::vector<Slot> slots_;
    size_t size_ = 0;
};

/**
 * The representative's k-mers of one node, with their multiplicity in the
 * representative and the number of members carrying each. Keyed by the canonical
 * ntHash value alone: at 4.6 M k-mers per genome a 64-bit collision is expected once in
 * about 10^6 nodes and would merge two counts, a deterministic error of one candidate.
 */
struct NodeKmers {
    KmerTable index;                                // hash -> position in the arrays below
    std::vector<uint16_t> repMult;                  // occurrences in the representative
    std::vector<std::atomic<uint16_t>> inNode;      // members (representative included) containing it
};

/**
 * A genome's sampled k-mers: the canonical hashes whose top `shift`-complement bits
 * are zero (1 in 2^(64 - shift); 1 in 64 at the default), one entry per occurrence,
 * sorted so that a k-mer repeated in the genome is a run. The sample is the same
 * k-mers in every genome (it is a property of the hash, not of the genome), so two
 * genomes' samples intersect exactly where the genomes do. The other hashes of a
 * sampled k-mer, for the filters, come back from the base hash (NTE64, as the
 * iterator derives them). 8 bytes per sampled k-mer: a 5 Mb genome is 0.6 MB, 5 493
 * of them 3.5 GB, held for the levels (sampleGenomes) and freed after them.
 */
struct GenomeSamples {
    unsigned shift = 58;                                         // 1 in 2^(64 - shift) k-mers
    std::unordered_map<std::string, std::vector<uint64_t>> byPath;
    double scale() const { return static_cast<double>(uint64_t(1) << (64 - shift)); }
};
inline bool sampledHash(uint64_t h, unsigned shift) { return shift >= 64 || (h >> shift) == 0; }
inline std::vector<uint64_t> sampleGenome(const UniqSketchConfig& cfg, const std::string& path, unsigned shift) {
    std::vector<uint64_t> out;
    scanFasta(path, [&](const std::string&, const std::string& seq, size_t) {
        if (seq.size() < cfg.kmerLen) return;
        ntHashIterator itr(seq, 1, cfg.kmerLen);
        while (itr != itr.end()) {
            if (sampledHash((*itr)[0], shift)) out.push_back((*itr)[0]);
            ++itr;
        }
    }, wrapPad(cfg));
    std::sort(out.begin(), out.end());
    return out;
}
inline GenomeSamples sampleGenomes(const UniqSketchConfig& cfg, const std::vector<std::string>& refFiles,
                                   unsigned shift = 58) {
    GenomeSamples gs;
    gs.shift = shift;
    std::vector<std::vector<uint64_t>> v(refFiles.size());
    #pragma omp parallel for schedule(dynamic)
    for (unsigned g = 0; g < refFiles.size(); g++) v[g] = sampleGenome(cfg, refFiles[g], shift);
    for (unsigned g = 0; g < refFiles.size(); g++) gs.byPath.emplace(refFiles[g], std::move(v[g]));
    return gs;
}

/** The result of chooseNodeRepresentative for one node. */
struct RepresentativeChoice {
    enum class Rule : char { None, Reached, Fill };
    unsigned member = 0;            // index into the group (0 = the genome the tree named)
    double reached = 0.0;           // members served by this member's signatures (saturated sum, see below)
    double carriage = 0.0;          // mean over the members of the share of this member's node-unique k-mers they contain
    double expectedCarried = 0.0;   // node-unique signatures of this member a random member would carry
    double poolUnique = 0.0;        // node-unique candidates, estimate
    double poolCarried = 0.0;       // carried candidates, estimate (reported, not scored)
    double fill = 0.0;              // carried candidates that fit in what the node-unique pool leaves of the quota
    bool decided = false;           // false: no member qualifies under either rule; keep the tree's genome
    Rule rule = Rule::None;         // which rule chose the member
    const char* ruleName() const { return rule == Rule::Reached ? "reached" : rule == Rule::Fill ? "fill" : "tree"; }
};

/**
 * Which member a node's signatures should come from.
 *
 * The tree names a node after its largest member, and the largest member is the one
 * with the most accessory sequence, which is what the others lack: its node-unique
 * k-mers are its private plasmids and prophages, and a member's reads hit few of
 * the signatures drawn from them (on 5 493 E. coli a member carried a median 20 % of
 * its node's signatures, 31 % of (member, node) pairs carried none). A MinHash
 * medoid was tried first and could not tell near-identical members apart: a 2 000-hash
 * sketch resolves about 1/2 000 of a genome, the size of the private pools that
 * decide this, so it degenerated into "the smallest genome" and picked a divergent
 * member for 169 of 1 143 nodes.
 *
 * So the choice is made on the quantity that matters, measured on the members'
 * sampled k-mers (the same 1 in 64 k-mers in every genome). For every member m: its
 * sampled k-mers that occur once in it and are node-unique (not in the solid filter,
 * exactly checkRef's rule) are its would-be signatures; for every other member i the
 * share of them i contains, times the signatures m would get (the pool scaled up from
 * the sample, capped at the quota), is the number of m's signatures i would carry.
 * The score is the sum over the members of min(carried_i, S) / S with S = `saturation`
 * (500): a member is served once it carries 500 of the node's signatures, and more
 * add little, so a representative that gives 100 members 300 signatures each beats
 * one that gives 20 members 1 800 (the saturation is what keeps the choice from
 * favouring a sub-clade's private core over a core the whole node shares). Ties,
 * which is every node whose members all carry 500 of any member's signatures, go to
 * the higher mean carriage (the member with the least private sequence among its
 * signatures: a family row's abundance is hits per signature, and signatures no
 * member carries dilute it), then to the smaller genome, then the earlier member. A
 * member whose node-unique pool is below `poolFloor` (200) is not scored, and a
 * member no other member carries (a score below `minReached`, 0.1: no member would
 * carry even 50 of its signatures) is not chosen through this rule.
 *
 * When no member qualifies, the members' node-unique k-mers are private to each of
 * them and the node's shared identity is in the carried pool (k-mers shared with at
 * most `cap` outside genomes), which fills what the node-unique selection leaves of
 * the quota (identifyGroupUniqKmersCarried). The representative is then the member
 * that leaves the fill the most room: the carried candidates that fit in
 * quota - min(pool, quota), the largest, ties to the smaller node-unique pool, then
 * the smaller genome; a member whose fit is below `poolFloor` is not chosen, and
 * when none qualifies the node keeps the genome the tree named. On the 5 493-genome
 * index 69 nodes fall to this rule, and in every one of them the member with the
 * smaller node-unique pool is the one the others carry more of: a 5-leaf node whose
 * tree genome had no node-unique k-mer kept its 1 233 filled signatures, carried by
 * 87 to 100 % of the members, instead of being handed to a member with 1 944
 * node-unique k-mers that no other member carries; another 5-leaf node went from
 * its tree genome (4 688 node-unique signatures and 295 filled, the members carry
 * 5 to 6 %) to the member with 998 and 2 094 (56 to 69 %).
 *
 * Under the first rule the carried pool is estimated and reported but not scored:
 * whether a carried k-mer survives the strict margin cannot be told from its hash,
 * and in the lineages where the carried pool is what the node has (a core shared
 * with close relatives outside the node) most of it does not (the K-12 node of the
 * 5 493-genome index: 833 carried candidates, 0 signatures). Scored against a
 * node-unique pool the members carry, it would hand such nodes to a member whose
 * signatures then do not exist; among members none of which is carried, the same
 * uncertainty applies to all of them and the ranking by room stands.
 */
inline RepresentativeChoice chooseNodeRepresentative(const UniqSketchConfig& cfg,
                                                     const std::vector<std::string>& group,
                                                     const GenomeSamples& samples,
                                                     const BloomFilter& sbFilter,
                                                     const CarriageCounter* cc, unsigned cap,
                                                     double saturation = 500.0, double poolFloor = 200.0,
                                                     double minReached = 0.1) {
    RepresentativeChoice best;
    const unsigned n = static_cast<unsigned>(group.size());
    if (n < 2) return best;
    std::vector<const std::vector<uint64_t>*> sv(n, nullptr);
    size_t total = 0;
    for (unsigned m = 0; m < n; m++) {
        auto it = samples.byPath.find(group[m]);
        if (it == samples.byPath.end()) return best;          // no sample: keep the tree's genome
        sv[m] = &it->second;
        total += it->second.size();
    }
    // The node's sampled k-mers, each with its key and the number of members
    // containing it; per member, the positions in its sample of the k-mers that
    // occur once in it.
    KmerTable table(total);
    std::vector<uint64_t> keys;
    std::vector<uint16_t> inNode;
    keys.reserve(total);
    inNode.reserve(total);
    std::vector<std::vector<uint32_t>> own(n);
    for (unsigned m = 0; m < n; m++) {
        const std::vector<uint64_t>& v = *sv[m];
        own[m].reserve(v.size());
        for (size_t i = 0; i < v.size();) {
            size_t j = i + 1;
            while (j < v.size() && v[j] == v[i]) ++j;        // a run: the k-mer's occurrences in m
            bool inserted = false;
            const uint32_t idx = table.findOrInsert(v[i], static_cast<uint32_t>(inNode.size()), inserted);
            if (inserted) { inNode.push_back(0); keys.push_back(v[i]); }
            if (inNode[idx] < UINT16_MAX) inNode[idx]++;
            if (j - i == 1) own[m].push_back(static_cast<uint32_t>(i));
            i = j;
        }
    }
    // Every k-mer classified once: 0 unusable, 1 node-unique, 2 carried.
    const unsigned halfMembers = std::max<unsigned>(2, (n + 1) / 2);
    const unsigned nh = std::max<unsigned>(cfg.nhash2, CarriageCounter::HASHES);
    std::vector<int8_t> cls(inNode.size(), 0);
    std::vector<uint64_t> hashes(std::max<unsigned>(nh, 1));
    for (size_t idx = 0; idx < keys.size(); idx++) {
        hashes[0] = keys[idx];
        for (unsigned i = 1; i < nh; i++) hashes[i] = NTE64(keys[idx], cfg.kmerLen, i);
        if (!sbFilter.contains(hashes.data())) cls[idx] = 1;
        else if (cc && cap > 0) {
            const unsigned in = inNode[idx];
            const unsigned all = cc->count(hashes.data());
            if (in >= halfMembers && all >= in && all - in <= cap) cls[idx] = 2;
        }
    }
    // Which members contain each node-unique k-mer: one bit row per such k-mer.
    std::vector<uint32_t> row(inNode.size(), UINT32_MAX);
    size_t nUnique = 0;
    for (size_t idx = 0; idx < cls.size(); idx++) if (cls[idx] == 1) row[idx] = static_cast<uint32_t>(nUnique++);
    const size_t words = (n + 63) / 64;
    std::vector<uint64_t> bits(nUnique * words, 0);
    for (unsigned m = 0; m < n; m++) {
        const std::vector<uint64_t>& v = *sv[m];
        for (size_t i = 0; i < v.size(); i++) {
            if (i && v[i] == v[i - 1]) continue;
            const uint32_t r = row[*table.find(v[i])];
            if (r != UINT32_MAX) bits[r * words + m / 64] |= uint64_t(1) << (m % 64);
        }
    }
    const double quota = static_cast<double>(std::max(cfg.sketchnum, 1));
    const double scale = samples.scale();
    // Per member: the two pools and, over its node-unique k-mers, the members it
    // reaches and the mean carriage.
    struct MemberStat { double poolU = 0.0, poolC = 0.0, reached = 0.0, carriage = 0.0, fill = 0.0; };
    std::vector<MemberStat> stat(n);
    std::vector<uint32_t> overlap(n);
    for (unsigned m = 0; m < n; m++) {
        const std::vector<uint64_t>& v = *sv[m];
        size_t uniq = 0, carr = 0;
        double uniqIn = 0.0;
        std::fill(overlap.begin(), overlap.end(), 0u);
        for (uint32_t pos : own[m]) {
            const uint32_t idx = *table.find(v[pos]);
            if (cls[idx] == 2) { ++carr; continue; }
            if (cls[idx] != 1) continue;
            ++uniq;
            uniqIn += inNode[idx];
            const uint64_t* r = &bits[row[idx] * words];
            for (size_t w = 0; w < words; w++)
                for (uint64_t word = r[w]; word; word &= word - 1)
                    overlap[w * 64 + __builtin_ctzll(word)]++;
        }
        MemberStat& s = stat[m];
        s.poolU = uniq * scale;
        s.poolC = carr * scale;
        const double takeU = std::min(s.poolU, quota);
        s.fill = std::min(quota - takeU, s.poolC);
        if (uniq == 0) continue;
        for (unsigned i = 0; i < n; i++) {
            if (i == m) continue;
            const double carried = static_cast<double>(overlap[i]) / static_cast<double>(uniq) * takeU;
            s.reached += std::min(carried, saturation) / saturation;
        }
        s.carriage = uniqIn / (static_cast<double>(uniq) * n);
    }
    auto take = [&](unsigned m, RepresentativeChoice::Rule rule) {
        const MemberStat& s = stat[m];
        best.decided = true;
        best.rule = rule;
        best.member = m;
        best.reached = s.reached;
        best.carriage = s.carriage;
        best.expectedCarried = std::min(s.poolU, quota) * s.carriage;
        best.poolUnique = s.poolU;
        best.poolCarried = s.poolC;
        best.fill = s.fill;
    };
    // The members reached through a member's node-unique k-mers.
    for (unsigned m = 0; m < n; m++) {
        const MemberStat& s = stat[m];
        if (s.poolU < poolFloor || s.reached < minReached) continue;
        const bool better = !best.decided ||
            s.reached > best.reached + 1e-9 ||
            (std::fabs(s.reached - best.reached) <= 1e-9 &&
             (s.carriage > best.carriage + 1e-9 ||
              (std::fabs(s.carriage - best.carriage) <= 1e-9 && sv[m]->size() < sv[best.member]->size())));
        if (better) take(m, RepresentativeChoice::Rule::Reached);
    }
    if (best.decided) return best;
    // No member reaches another through its node-unique k-mers: the node's shared
    // identity is in the carried pool, and the representative is the member that
    // leaves the fill the most room.
    for (unsigned m = 0; m < n; m++) {
        const MemberStat& s = stat[m];
        if (s.fill < poolFloor) continue;
        const bool better = !best.decided ||
            s.fill > best.fill + 1e-9 ||
            (std::fabs(s.fill - best.fill) <= 1e-9 &&
             (s.poolU < best.poolUnique - 1e-9 ||
              (std::fabs(s.poolU - best.poolUnique) <= 1e-9 && sv[m]->size() < sv[best.member]->size())));
        if (better) take(m, RepresentativeChoice::Rule::Fill);
    }
    return best;
}

/**
 * Candidate scan for one node: every k-mer not in the solid filter into the node's
 * candidate file, exactly as checkRef, and into a second file the representative's
 * k-mers that occur once in it, are carried by at least half the members (one of them
 * not the representative) and by at most `cap` genomes outside the node (none: no
 * second file). Returns (plain candidates, carried candidates).
 */
inline std::pair<size_t, size_t> checkRefCarried(const UniqSketchConfig& cfg, const std::string& repPath,
                                                 size_t members, const BloomFilter& sbFilter,
                                                 const CarriageCounter* cc, const NodeKmers& nk,
                                                 unsigned cap) {
    const std::string candPath = candidatePath(cfg.outdir, repPath);
    const std::string carrPath = carriedCandidatePath(candPath);
    // At least half the members, the representative included, and at least one
    // member besides it: for a node of two, "half" is the representative alone, and
    // the rule would admit its private sequence with a few outside carriers (on the
    // 120-genome index two plasmid k-mers of a two-leaf node, carried by three other
    // E. coli and by something in a human gut metagenome with no Escherichia, made
    // the node a family row of that sample: 125 hits on the two of them).
    const unsigned halfMembers = std::max<unsigned>(2, static_cast<unsigned>((members + 1) / 2));
    size_t plain = 0, carried = 0, all = 0;
    bool candOk = false, carrOk = false;
    {
        CandidateWriter cand(candPath, cfg.kmerLen, repPath);
        CandidateWriter carr(carrPath, cfg.kmerLen, repPath);
        scanFasta(repPath, [&](const std::string& header, const std::string& faSeq, size_t trueLen) {
            cand.beginContig(header, trueLen);
            carr.beginContig(header, trueLen);
            // The solid block, the table slot and the four carriage counters of a k-mer
            // are asked for when it is hashed, and it is judged sixteen k-mers later, in
            // position order as before (forEachKmerAhead).
            const bool useCarriage = cc && cap > 0;
            forEachKmerAhead(faSeq, std::max<unsigned>(cfg.nhash2, CarriageCounter::HASHES), cfg.kmerLen,
                [&](const uint64_t* hVal) {
                    sbFilter.prefetch(hVal[0]);
                    if (useCarriage) { nk.index.prefetch(hVal[0]); cc->prefetch(hVal); }
                },
                [&](const uint64_t* hVal, size_t pos) {
                    ++all;
                    if (!sbFilter.contains(hVal)) {
                        ++plain;
                        cand.add(static_cast<uint32_t>(pos));
                    } else {
                        const uint32_t* v = useCarriage ? nk.index.find(hVal[0]) : nullptr;
                        if (v && nk.repMult[*v] == 1) {
                            const unsigned in = nk.inNode[*v].load(std::memory_order_relaxed);
                            const unsigned total = cc->count(hVal);
                            if (in >= halfMembers && total >= in && total - in <= cap) {
                                ++carried;
                                carr.add(static_cast<uint32_t>(pos));
                            }
                        }
                    }
                });
            cand.endContig();
            carr.endContig();
        }, wrapPad(cfg));
        candOk = cand.good();
        carrOk = carr.good();
    }
    if (carried == 0 || !carrOk) std::remove(carrPath.c_str());
    // The messages checkRef gives, with the carried pool taken into account: a node
    // with nothing unique to it but a carried pool still gets signatures.
    if (!candOk) {
        std::cerr << "Error: could not write " + candPath + "\n";
    } else if (carried && !carrOk) {
        std::cerr << "Error: could not write " + carrPath + "\n";
    } else if (all == 0) {
        std::cerr << "Warning: " + repPath + ": no k-mers of length " + std::to_string(cfg.kmerLen) +
                     " (empty, not FASTA, or every sequence shorter than k); it will get no signatures\n";
    } else if (plain + carried == 0) {
        std::cerr << "Warning: " + repPath + ": no unique k-mers (identical to another reference?); "
                     "it will get no signatures\n";
    } else if (plain == 0) {
        std::cerr << "Note: " + repPath + ": no k-mers unique to its node; its signatures at this level "
                     "are all carried by genomes outside it (" + std::to_string(carried) + " candidates)\n";
    }
    return {plain, carried};
}

/**
 * identifyGroupUniqKmers with a carried pool for real nodes. The plain scan and its
 * filters are exactly those of identifyGroupUniqKmers, so the candidate files and the
 * stat file are byte-identical to its; the carried pool is a second file per node
 * (carriedCandidatePath), which buildSketch draws on only to fill a quota the plain
 * pool left short. For every group of two or more genomes the representative's k-mers
 * are hashed, the members are scanned once to count how many carry each (one flat list
 * of (group, member) jobs, so the 896-genome node does not leave the machine idle),
 * then the representative is scanned again with both rules. Groups are taken in
 * batches whose representatives' k-mers fit `budgetKmers`, since the table and counts cost about
 * 30 bytes per k-mer (KmerTable); every genome is scanned once per level whatever the batching.
 *
 * With `samples`, a node's representative is chosen first (chooseNodeRepresentative)
 * for every group whose `decide` flag is set, and that member is moved to the front of
 * the group; the function returns the first member of every group as it ended up. The
 * counter may be null (no carried pool, cap taken as 0): then the candidate files are
 * exactly identifyGroupUniqKmers's.
 */
inline std::vector<std::string> identifyGroupUniqKmersCarried(
        const UniqSketchConfig& cfg, std::vector<std::vector<std::string>> groups,
        BloomFilter& dbFilter, const CarriageCounter* cc, unsigned cap,
        const GenomeSamples* samples = nullptr, const std::vector<char>* decide = nullptr,
        std::vector<RepresentativeChoice>* choices = nullptr,
        size_t budgetKmers = size_t(150) * 1000 * 1000) {
    if (!cc) cap = 0;
    size_t members = 0;
    for (const auto& g : groups) members += g.size();
    // The members containing a k-mer are counted in 16 bits (NodeKmers::inNode, a
    // plain fetch_add in the member scan below). In a node of more members the
    // count would wrap for the k-mers 65 536 or more of them hold: those would fail
    // the half-of-the-members test and the node would lose that part of its carried
    // pool, the part it shares most widely, without a word; from 131 071 members
    // on, no count reaches half the members and that part is the whole pool. The
    // carriage counter beside it saturates at 65 535.
    if (cap > 0)
        for (const auto& g : groups)
            if (g.size() > UINT16_MAX) {
                std::cerr << "Error: a tree node of " << g.size() << " genomes; carried node signatures count a "
                             "node's members in 16 bits (65 535 at most). Use finer --tree thresholds or "
                             "--node-carriers=0\n";
                exit(EXIT_FAILURE);
            }
    std::cout << "k-mer length: " << cfg.kmerLen << "\n"
              << "Groups: " << groups.size() << " over " << members << " genomes\n";
    using Clock = std::chrono::steady_clock;
    auto secs = [](Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double>(b - a).count(); };
    auto t0 = Clock::now();
    BloomFilter sbFilter(cfg.m2(), cfg.nhash2, cfg.kmerLen);
    loadBFfaGroups(cfg, groups, dbFilter, sbFilter);
    std::vector<int> refCount(groups.size());
    // Singles: the plain rule, as before.
    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < groups.size(); ++i)
        if (groups[i].size() == 1) checkRef(cfg, groups[i][0], i, sbFilter, refCount);
    const double secFilters = secs(t0, Clock::now());
    double secReps = 0, secMembers = 0, secScan = 0;
    size_t repKmers = 0, memberGenomes = 0;
    // Nodes, in batches.
    std::vector<unsigned> nodes;
    for (unsigned i = 0; i < groups.size(); ++i) if (groups[i].size() >= 2) nodes.push_back(i);
    // The representative of each node still to be decided, chosen on the members'
    // samples and moved to the front of its group.
    double secChoice = 0;
    size_t chosen = 0, changed = 0, undecided = 0, byFill = 0;
    if (samples) {
        auto tc = Clock::now();
        std::vector<RepresentativeChoice> choice(nodes.size());
        #pragma omp parallel for schedule(dynamic)
        for (size_t j = 0; j < nodes.size(); j++) {
            const unsigned i = nodes[j];
            if (decide && !(*decide)[i]) continue;
            choice[j] = chooseNodeRepresentative(cfg, groups[i], *samples, sbFilter, cc, cap);
        }
        if (choices) choices->assign(groups.size(), RepresentativeChoice());
        for (size_t j = 0; j < nodes.size(); j++) {
            const unsigned i = nodes[j];
            if (decide && !(*decide)[i]) continue;
            ++chosen;
            if (choices) (*choices)[i] = choice[j];
            if (!choice[j].decided) { ++undecided; continue; }
            if (choice[j].rule == RepresentativeChoice::Rule::Fill) ++byFill;
            if (choice[j].member != 0) {
                std::swap(groups[i][0], groups[i][choice[j].member]);
                ++changed;
            }
        }
        secChoice = secs(tc, Clock::now());
    }
    size_t totalPlain = 0, totalCarried = 0, nodesWithCarried = 0;
    for (size_t b = 0; b < nodes.size();) {
        std::vector<unsigned> batch;
        size_t kmers = 0;
        while (b < nodes.size()) {
            const size_t sz = uncompressedSize(groups[nodes[b]][0]);
            if (!batch.empty() && kmers + sz > budgetKmers) break;
            batch.push_back(nodes[b++]);
            kmers += sz;
        }
        // Representatives' k-mers, in parallel over the batch.
        auto tb = Clock::now();
        std::vector<NodeKmers> nk(batch.size());
        #pragma omp parallel for schedule(dynamic)
        for (unsigned j = 0; j < batch.size(); j++) {
            NodeKmers& n = nk[j];
            n.index.reserve(uncompressedSize(groups[batch[j]][0]));   // distinct k-mers <= bases <= file size
            std::vector<uint16_t> mult;
            scanFasta(groups[batch[j]][0], [&](const std::string&, const std::string& seq, size_t) {
                if (seq.size() < cfg.kmerLen) return;
                // The slot is asked for when the k-mer is hashed and filled sixteen k-mers
                // later, in the same order: the order numbers the k-mers (forEachKmerAhead).
                forEachKmerAhead(seq, 1, cfg.kmerLen,
                    [&](const uint64_t* hVal) { n.index.prefetch(hVal[0]); },
                    [&](const uint64_t* hVal, size_t) {
                        bool inserted = false;
                        const uint32_t v = n.index.findOrInsert(hVal[0], static_cast<uint32_t>(mult.size()), inserted);
                        if (inserted) mult.push_back(1);
                        else if (mult[v] < UINT16_MAX) mult[v]++;
                    });
            }, wrapPad(cfg));
            n.repMult = std::move(mult);
            std::vector<std::atomic<uint16_t>> in(n.repMult.size());
            for (auto& a : in) a.store(0, std::memory_order_relaxed);
            n.inNode = std::move(in);
        }
        for (const NodeKmers& n : nk) repKmers += n.repMult.size();
        auto tm = Clock::now();
        secReps += secs(tb, tm);
        // Members, one flat job list: each member marks the representative k-mers it
        // contains (once per member, however often it has them) and adds to the count.
        std::vector<std::pair<unsigned, unsigned>> jobs;   // (batch index, member index)
        for (unsigned j = 0; j < batch.size(); j++)
            for (unsigned m = 0; m < groups[batch[j]].size(); m++) jobs.push_back({j, m});
        #pragma omp parallel for schedule(dynamic)
        for (size_t t = 0; t < jobs.size(); t++) {
            NodeKmers& n = nk[jobs[t].first];
            std::vector<bool> seen(n.repMult.size(), false);
            scanFasta(groups[batch[jobs[t].first]][jobs[t].second], [&](const std::string&, const std::string& seq, size_t) {
                if (seq.size() < cfg.kmerLen) return;
                // A member k-mer is a probe of the representative's table (a hundred
                // megabytes) and, when found there, an add to a counter the other
                // members' threads are adding to as well. Two steps of look-ahead: the
                // table slot is asked for when the k-mer is hashed and looked up sixteen
                // k-mers later (forEachKmerAhead), which asks for the counter; the count
                // is made sixteen found k-mers after that. The k-mers are still counted
                // in genome order (the count is a sum over members, each marking a k-mer
                // once, so the order would not matter either).
                uint32_t found[KMER_LOOKAHEAD];
                size_t nFound = 0;
                auto count = [&](uint32_t v) {
                    if (!seen[v]) {
                        seen[v] = true;
                        n.inNode[v].fetch_add(1, std::memory_order_relaxed);
                    }
                };
                forEachKmerAhead(seq, 1, cfg.kmerLen,
                    [&](const uint64_t* hVal) { n.index.prefetch(hVal[0]); },
                    [&](const uint64_t* hVal, size_t) {
                        const uint32_t* v = n.index.find(hVal[0]);
                        if (!v) return;
                        uint32_t& slot = found[nFound % KMER_LOOKAHEAD];
                        if (nFound >= KMER_LOOKAHEAD) count(slot);
                        slot = *v;
                        __builtin_prefetch(&n.inNode[*v]);
                        ++nFound;
                    });
                for (size_t i = nFound > KMER_LOOKAHEAD ? nFound - KMER_LOOKAHEAD : 0; i < nFound; i++)
                    count(found[i % KMER_LOOKAHEAD]);
            }, wrapPad(cfg));
        }
        memberGenomes += jobs.size();
        auto ts = Clock::now();
        secMembers += secs(tm, ts);
        // The scan of each representative. The stat file records the plain count, as
        // identifyGroupUniqKmers does; the carried pool is found by its file.
        std::vector<std::pair<size_t, size_t>> res(batch.size());
        #pragma omp parallel for schedule(dynamic)
        for (unsigned j = 0; j < batch.size(); j++) {
            const unsigned i = batch[j];
            res[j] = checkRefCarried(cfg, groups[i][0], groups[i].size(), sbFilter, cc, nk[j], cap);
            refCount[i] = static_cast<int>(res[j].first);
        }
        for (unsigned j = 0; j < batch.size(); j++) {
            totalPlain += res[j].first;
            totalCarried += res[j].second;
            if (res[j].second) ++nodesWithCarried;
        }
        secScan += secs(ts, Clock::now());
    }
    std::cerr << "Node carriers scan time (s): filters and singles " << std::setprecision(1) << std::fixed
              << secFilters << ", representatives " << secReps << " (" << repKmers / 1000000 << " M k-mers hashed), members "
              << secMembers << " (" << memberGenomes << " genomes), scan " << secScan << "\n";
    if (samples)
        std::cout << "Node representatives: " << chosen << " nodes decided on the members' samples in "
                  << std::setprecision(1) << std::fixed << secChoice << " s, " << changed
                  << " moved from the genome the tree named" << (byFill ? ", " + std::to_string(byFill) +
                     " chosen for the room they leave the carried pool (no member's node-unique k-mers are carried)" : std::string())
                  << (undecided ? ", " + std::to_string(undecided) + " kept (no member qualifies)" : std::string()) << "\n";
    std::cout << "Node candidates: " << totalPlain << " unique to their node; " << totalCarried
              << " carried by at most " << cap << " outside genomes, held for " << nodesWithCarried
              << " of " << nodes.size() << " nodes to fill a quota their unique ones leave short\n";
    std::ofstream tsvOut(cfg.refstat);
    for (size_t i = 0; i < refCount.size(); i++)
        tsvOut << candidatePath(cfg.outdir, groups[i][0]) << "\t" << refCount[i] << "\t" << i << "\n";
    std::vector<std::string> reps;
    reps.reserve(groups.size());
    for (const auto& g : groups) reps.push_back(g[0]);
    return reps;
}

#endif // NODECARRIERUTIL_HPP_
