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
        for (size_t i = 0; i < size_; i++) counters_[i].store(0, std::memory_order_relaxed);
    }
    size_t size() const { return size_; }
    size_t bytes() const { return size_ * sizeof(uint16_t); }
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
 */
inline void countCarriage(const UniqSketchConfig& cfg, const std::vector<std::string>& refFiles,
                          CarriageCounter& cc) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned g = 0; g < refFiles.size(); g++) {
        BloomFilter seen(std::max<size_t>(uncompressedSize(refFiles[g]), 1) * 64,
                         CarriageCounter::HASHES, cfg.kmerLen);
        scanFasta(refFiles[g], [&](const std::string&, const std::string& seq, size_t) {
            if (seq.size() < cfg.kmerLen) return;
            ntHashIterator itr(seq, CarriageCounter::HASHES, cfg.kmerLen);
            while (itr != itr.end()) {
                if (seen.insert_make_change(*itr)) cc.add(*itr);
                ++itr;
            }
        }, wrapPad(cfg));
    }
}

/**
 * The representative's k-mers of one node, with their multiplicity in the
 * representative and the number of members carrying each. Keyed by the canonical
 * ntHash value alone: at 4.6 M k-mers per genome a 64-bit collision is expected once in
 * about 10^6 nodes and would merge two counts, a deterministic error of one candidate.
 */
struct NodeKmers {
    std::unordered_map<uint64_t, uint32_t> index;   // hash -> position in the arrays below
    std::vector<uint16_t> repMult;                  // occurrences in the representative
    std::vector<std::atomic<uint16_t>> inNode;      // members (representative included) containing it
};

/**
 * Candidate scan for one node: every k-mer not in the solid filter into the node's
 * candidate file, exactly as checkRef, and into a second file the representative's
 * k-mers that occur once in it, are carried by at least half the members (one of them
 * not the representative) and by at most `cap` genomes outside the node (none: no
 * second file). Returns (plain candidates, carried candidates).
 */
inline std::pair<size_t, size_t> checkRefCarried(const UniqSketchConfig& cfg, const std::string& repPath,
                                                 size_t members, const BloomFilter& sbFilter,
                                                 const CarriageCounter& cc, const NodeKmers& nk,
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
            ntHashIterator itr(faSeq, std::max<unsigned>(cfg.nhash2, CarriageCounter::HASHES), cfg.kmerLen);
            while (itr != itr.end()) {
                ++all;
                if (!sbFilter.contains(*itr)) {
                    ++plain;
                    cand.add(static_cast<uint32_t>(itr.get_pos()));
                } else {
                    auto it = nk.index.find((*itr)[0]);
                    if (it != nk.index.end() && nk.repMult[it->second] == 1) {
                        const unsigned in = nk.inNode[it->second].load(std::memory_order_relaxed);
                        const unsigned total = cc.count(*itr);
                        if (in >= halfMembers && total >= in && total - in <= cap) {
                            ++carried;
                            carr.add(static_cast<uint32_t>(itr.get_pos()));
                        }
                    }
                }
                ++itr;
            }
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
 * batches whose representatives' k-mers fit `budgetKmers`, since the hash costs about
 * 40 bytes per k-mer; every genome is scanned once per level whatever the batching.
 */
inline void identifyGroupUniqKmersCarried(const UniqSketchConfig& cfg,
                                          const std::vector<std::vector<std::string>>& groups,
                                          BloomFilter& dbFilter, const CarriageCounter& cc,
                                          unsigned cap, size_t budgetKmers = size_t(150) * 1000 * 1000) {
    size_t members = 0;
    for (const auto& g : groups) members += g.size();
    std::cout << "k-mer length: " << cfg.kmerLen << "\n"
              << "Groups: " << groups.size() << " over " << members << " genomes\n";
    BloomFilter sbFilter(cfg.m2(), cfg.nhash2, cfg.kmerLen);
    loadBFfaGroups(cfg, groups, dbFilter, sbFilter);
    std::vector<int> refCount(groups.size());
    // Singles: the plain rule, as before.
    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < groups.size(); ++i)
        if (groups[i].size() == 1) checkRef(cfg, groups[i][0], i, sbFilter, refCount);
    // Nodes, in batches.
    std::vector<unsigned> nodes;
    for (unsigned i = 0; i < groups.size(); ++i) if (groups[i].size() >= 2) nodes.push_back(i);
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
        std::vector<NodeKmers> nk(batch.size());
        #pragma omp parallel for schedule(dynamic)
        for (unsigned j = 0; j < batch.size(); j++) {
            NodeKmers& n = nk[j];
            std::vector<uint16_t> mult;
            scanFasta(groups[batch[j]][0], [&](const std::string&, const std::string& seq, size_t) {
                if (seq.size() < cfg.kmerLen) return;
                ntHashIterator itr(seq, 1, cfg.kmerLen);
                while (itr != itr.end()) {
                    auto it = n.index.find((*itr)[0]);
                    if (it == n.index.end()) {
                        n.index.emplace((*itr)[0], static_cast<uint32_t>(mult.size()));
                        mult.push_back(1);
                    } else if (mult[it->second] < UINT16_MAX) {
                        mult[it->second]++;
                    }
                    ++itr;
                }
            }, wrapPad(cfg));
            n.repMult = std::move(mult);
            std::vector<std::atomic<uint16_t>> in(n.repMult.size());
            for (auto& a : in) a.store(0, std::memory_order_relaxed);
            n.inNode = std::move(in);
        }
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
                ntHashIterator itr(seq, 1, cfg.kmerLen);
                while (itr != itr.end()) {
                    auto it = n.index.find((*itr)[0]);
                    if (it != n.index.end() && !seen[it->second]) {
                        seen[it->second] = true;
                        n.inNode[it->second].fetch_add(1, std::memory_order_relaxed);
                    }
                    ++itr;
                }
            }, wrapPad(cfg));
        }
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
    }
    std::cout << "Node candidates: " << totalPlain << " unique to their node; " << totalCarried
              << " carried by at most " << cap << " outside genomes, held for " << nodesWithCarried
              << " of " << nodes.size() << " nodes to fill a quota their unique ones leave short\n";
    std::ofstream tsvOut(cfg.refstat);
    for (size_t i = 0; i < refCount.size(); i++)
        tsvOut << candidatePath(cfg.outdir, groups[i][0]) << "\t" << refCount[i] << "\t" << i << "\n";
}

#endif // NODECARRIERUTIL_HPP_
