#include <fstream>
#include <iterator>
#include <algorithm>
#include <cstdio>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <zlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "SequenceUtil.hpp"
#include "TestUtil.hpp"
#include "UniqSketchUtil.hpp"

// The candidate fixtures (expected_ref_test_*.tsv) are in the text form the
// candidate files had before they went binary: one line per unique k-mer,
// `canonical k-mer <TAB> position <TAB> contig`, in scan order. They stay the
// ground truth; the helpers below translate between them and the binary files.

static std::vector<std::string> fileLines(const std::string& path) {
    std::vector<std::string> lines;
    std::ifstream in(path);
    for (std::string l; std::getline(in, l);) lines.push_back(l);
    return lines;
}

// Spell a candidate file out as the text fixture would: the canonical k-mer at
// each recorded position of the reference, with the position and contig.
static std::vector<std::string> candidateLines(const std::string& candPath) {
    CandidateSet cand;
    bool ok = readCandidates(candPath, cand);
    assert(ok);
    std::vector<std::string> contigs;
    scanFasta(cand.refPath, [&](const std::string&, std::string& seq, size_t) { contigs.push_back(std::move(seq)); });
    assert(contigs.size() == cand.runs.size());
    std::vector<std::string> lines;
    for (size_t r = 0; r < cand.runs.size(); r++) {
        assert(contigs[r].size() == cand.runs[r].seqLen);
        for (uint32_t p : cand.runs[r].pos)
            lines.push_back(get_canonical(contigs[r].substr(p, cand.k)) + "\t" + std::to_string(p) +
                            "\t" + cand.runs[r].name);
    }
    return lines;
}

// The reverse: a binary candidate file from a text fixture, for the reference
// the fixture describes.
static void candidatesFromText(const std::string& textPath, const std::string& refPath,
                               const std::string& candPath, unsigned k) {
    std::map<std::string, std::vector<uint32_t>> posByContig;
    for (const std::string& l : fileLines(textPath)) {
        std::istringstream ls(l);
        std::string seq, contig; uint32_t pos;
        ls >> seq >> pos >> contig;
        posByContig[contig].push_back(pos);
    }
    CandidateWriter w(candPath, k, refPath);
    scanFasta(refPath, [&](const std::string& header, const std::string& seq, size_t) {
        w.beginContig(header, seq.size());
        for (uint32_t p : posByContig[header]) w.add(p);
        w.endContig();
    });
    assert(w.good());
}

void testIdentifyUniqKmers(const UniqSketchConfig& cfg,
                           const std::vector<std::string>& refFiles) {
    std::cerr << "START: Test identifyUniqKmers... \t";

    /*
    ref_test_1: {CAGGCTACAAA, GAGACTGAACT} -> 4 kmers of 10bp
    ref_test_2: {CTTAGCTGAGG, CAGGCTACAAC} -> 4 kmers of 10bp

    CAGGCTACAA is shared, so canonical unique sets are:
    ref_test_1: {AGGCTACAAA, GAGACTGAAC, AGACTGAACT}
    ref_test_2: {CTCAGCTAAG, CCTCAGCTAA, AGGCTACAAC}
    */
    BloomFilter dbFilter(cfg.m1(), cfg.nhash1, cfg.kmerLen);
    identifyUniqKmers(cfg, refFiles, dbFilter);

    assert(candidateLines("ref_test_1.cand") == fileLines("expected_ref_test_1.tsv"));
    assert(candidateLines("ref_test_2.cand") == fileLines("expected_ref_test_2.tsv"));
    // The stat file names the candidate files with their counts.
    std::vector<std::string> stat = fileLines(cfg.refstat);
    assert(stat.size() == 2);
    assert(stat[0] == "./ref_test_1.cand\t3\t0" && stat[1] == "./ref_test_2.cand\t3\t1");

    std::cerr << "PASSED: Test identifyUniqKmers\n";
}

// A k-mer that occurs twice inside the representative is not a group signature.
//
// loadBFfa marks such a k-mer solid on its second sighting, so a leaf never gets
// one. loadBFfaGroup counts each k-mer once per group, which is what makes a
// k-mer shared by every member the best kind of group signature, but left at
// that the second copy inside a single genome never reaches the solid filter
// either. Multi-copy k-mers sit in mobile elements, so they are the sequence
// most likely to be present in organisms outside the index entirely: a two-copy
// erm(B) in one representative turned ~250 reads of a human gut metagenome,
// which held no Escherichia at all, into a family call.
//
// One k-mer per contig keeps the fixture unambiguous.
void testGroupRepeatExclusion(UniqSketchConfig cfg) {
    std::cerr << "START: Test group repeat exclusion... \t";
    const std::string dup     = "ACGTACGTAC";   // twice in the representative
    const std::string shared  = "TTGGTTGGTT";   // once in the representative, once in the member
    const std::string solo    = "AAACCCGGGT";   // once, representative only
    const std::string memOnly = "GGGTTTAAAC";   // once, member only
    {
        std::ofstream rep("group_rep.fa");
        rep << ">dup_a\n" << dup << "\n>dup_b\n" << dup << "\n"
            << ">shared\n" << shared << "\n>solo\n" << solo << "\n";
        std::ofstream mem("group_mem.fa");
        mem << ">shared\n" << shared << "\n>mem_only\n" << memOnly << "\n";
    }
    // Sized well clear of a solid-filter false positive, which would show up
    // here as a missing signature rather than a surviving repeat.
    cfg.dbfSize = 64;
    cfg.sbfSize = 64;
    cfg.refstat = "db_uniq_count_group.tsv";
    std::vector<std::vector<std::string>> groups = {{"group_rep.fa", "group_mem.fa"}};
    BloomFilter dbFilter(cfg.m1(), cfg.nhash1, cfg.kmerLen);
    identifyGroupUniqKmers(cfg, groups, dbFilter);

    std::vector<std::string> seqs;
    for (const std::string& l : candidateLines("group_rep.cand"))
        seqs.push_back(l.substr(0, l.find('\t')));
    auto has = [&](const std::string& kmer) {
        return std::find(seqs.begin(), seqs.end(), get_canonical(kmer)) != seqs.end();
    };
    assert(!has(dup));      // the regression: neither copy may become a candidate
    assert(has(shared));    // a k-mer the whole group carries is still a signature
    assert(has(solo));
    assert(seqs.size() == 2);
    // Only the representative is scanned for candidates.
    std::vector<std::string> stat = fileLines(cfg.refstat);
    assert(stat.size() == 1 && stat[0] == "./group_rep.cand\t2\t0");

    std::remove("group_rep.fa");
    std::remove("group_mem.fa");
    std::remove("group_rep.cand");
    std::remove(cfg.refstat.c_str());
    std::cerr << "PASSED: Test group repeat exclusion\n";
}

// loadBFfaGroups must leave the distinct and solid filters bit for bit as the
// serial scan does (loadBFfaGroup on each group in turn), on any thread count:
// the partitioned scan replays the serial order of operations on every block, so
// the filters' false positives, and through them the candidates, cannot depend
// on scheduling. The groups exercise every path: a five-genome group with a
// contig longer than one hashing piece, the representative's internal repeat, a
// k-mer shared by two members but not the representative, k-mers shared by the
// whole group, runs of N, a contig shorter than k and a gzip member; a group of
// one (loadBFfa's rule, no group-private filters); a group of two; another
// single. A small flush size forces several flushes within and across files.
//
// The filters are thin (16 bits per element) so that false positives are
// common, and the test first checks that it can see order dependence at all:
// the serial scan with the groups in reverse order must give different filters.
// Without that, a replay that silently fell back to some other order could pass.
void testGroupScanIdentical() {
    std::cerr << "START: Test group scan identical to serial on any thread count... \t";
    uint64_t s = 0x9E3779B97F4A7C15ull;
    auto rnd = [&s]() { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; };
    auto dna = [&](size_t n) {
        std::string d(n, 'A');
        for (char& c : d) c = "ACGT"[rnd() & 3];
        return d;
    };
    auto mutate = [&](std::string d, size_t every) {
        for (size_t i = every / 2; i < d.size(); i += every) d[i] = "ACGT"[(rnd() & 3)];
        return d;
    };
    const std::string core = dna(600000);                 // > 2^18 k-mers: several hashing pieces
    const std::string repeat = dna(5000);
    const std::string pairShared = dna(2000);             // in two members, not the representative
    const std::vector<std::vector<std::string>> groups = {
        {"pgroup_rep.fa", "pgroup_m1.fa", "pgroup_m2.fa", "pgroup_m3.fa", "pgroup_m4.fa.gz"},
        {"pgroup_s1.fa"},
        {"pgroup_p1.fa", "pgroup_p2.fa"},
        {"pgroup_s2.fa"},
    };
    {
        std::ofstream rep(groups[0][0]);
        rep << ">core\n" << core << "\n"
            << ">rep_dup\n" << dna(20000) << repeat << dna(20000) << repeat << dna(3000) << "\n"
            << ">with_n\n" << dna(4000) << std::string(50, 'N') << dna(4000) << "\n"
            << ">short\nACGTACGT\n"
            << ">tail\n" << dna(30000) << "\n";
        std::ofstream m1(groups[0][1]);
        m1 << ">m1_core\n" << mutate(core, 997) << "\n>m1_own\n" << dna(100000) << pairShared << "\n";
        std::ofstream m2(groups[0][2]);
        m2 << ">m2_core\n" << mutate(core, 1009) << "\n>m2_own\n" << pairShared << dna(90000) << "\n";
        std::ofstream m3(groups[0][3]);
        m3 << ">m3_core\n" << mutate(core, 1013) << "\n>m3_own\n" << dna(110000) << "\n";
        std::ostringstream m4;
        m4 << ">m4_core\n" << mutate(core, 1019) << "\n>m4_own\n" << dna(70000) << "\n";
        gzFile g = gzopen(groups[0][4].c_str(), "wb");
        assert(g);
        const std::string m4s = m4.str();
        assert(gzwrite(g, m4s.data(), static_cast<unsigned>(m4s.size())) == static_cast<int>(m4s.size()));
        gzclose(g);
        // A single with its own internal repeat (solid for it, as loadBFfa rules),
        // a pair sharing most of their sequence, and a plain single.
        const std::string dup = dna(3000);
        std::ofstream s1(groups[1][0]);
        s1 << ">s1\n" << dna(150000) << dup << dna(50000) << dup << dna(20000) << "\n";
        const std::string pairCore = dna(200000);
        std::ofstream p1(groups[2][0]);
        p1 << ">p1\n" << pairCore << dna(30000) << "\n";
        std::ofstream p2(groups[2][1]);
        p2 << ">p2_a\n" << mutate(pairCore, 503) << "\n>p2_b\n" << dna(40000) << "\n";
        std::ofstream s2(groups[3][0]);
        s2 << ">s2\n" << dna(180000) << "\n";
    }
    UniqSketchConfig cfg;
    cfg.kmerLen = 31;
    cfg.bits = 16;                                         // thin filters: false positives are common
    cfg.dbfSize = 1800000;
    cfg.sbfSize = 900000;
    auto bitsOf = [&](const BloomFilter& db, const BloomFilter& sb) {
        std::string bits(reinterpret_cast<const char*>(const_cast<BloomFilter&>(db).getFilter()), db.sizeBytes());
        bits.append(reinterpret_cast<const char*>(const_cast<BloomFilter&>(sb).getFilter()), sb.sizeBytes());
        return bits;
    };
    auto serial = [&](const std::vector<std::vector<std::string>>& gs) {
        BloomFilter db(cfg.m1(), cfg.nhash1, cfg.kmerLen);
        BloomFilter sb(cfg.m2(), cfg.nhash2, cfg.kmerLen);
        for (const auto& g : gs) loadBFfaGroup(cfg, g, db, sb);
        assert(db.get_pop() > 0 && sb.get_pop() > 0);
        return bitsOf(db, sb);
    };
    auto scanned = [&](size_t flushBases) {
        BloomFilter db(cfg.m1(), cfg.nhash1, cfg.kmerLen);
        BloomFilter sb(cfg.m2(), cfg.nhash2, cfg.kmerLen);
        loadBFfaGroups(cfg, groups, db, sb, flushBases);
        return bitsOf(db, sb);
    };
    const std::string reference = serial(groups);
    // The test has teeth: order dependence is visible in these filters. It shows
    // in the solid filter only; the distinct filter is a union of the k-mers each
    // group admits (a group-private decision), so it is the same in any order.
    std::vector<std::vector<std::string>> reversed(groups.rbegin(), groups.rend());
    const std::string reorderedRef = serial(reversed);
    const size_t dbBytes = BloomFilter(cfg.m1(), cfg.nhash1, cfg.kmerLen).sizeBytes();
    assert(reorderedRef != reference);
    assert(reorderedRef.compare(0, dbBytes, reference, 0, dbBytes) == 0);
    auto check = [&](const std::string& got) {
        if (BloomFilter::kPartitioned) {
            assert(got == reference);
        } else {
            // Standard filter: no partitions, groups go one per thread, and the
            // solid filter depends on how they interleave. The distinct filter
            // must still come out the same.
            assert(got.compare(0, dbBytes, reference, 0, dbBytes) == 0);
        }
    };
#ifdef _OPENMP
    const int savedThreads = omp_get_max_threads();
    for (int threads : {1, 8}) {
        omp_set_num_threads(threads);
        check(scanned(100000));                           // several flushes within and across files
        check(scanned(size_t(32) << 20));                 // the default: one flush per file side
    }
    omp_set_num_threads(savedThreads);
#else
    check(scanned(100000));
    check(scanned(size_t(32) << 20));
#endif
    for (const auto& g : groups)
        for (const std::string& f : g) std::remove(f.c_str());
    std::cerr << "PASSED: Test group scan identical to serial on any thread count\n";
}

// Write `content` to `path` gzip-compressed.
static void writeGzip(const std::string& path, const std::string& content) {
    gzFile g = gzopen(path.c_str(), "wb");
    assert(g);
    int n = gzwrite(g, content.data(), static_cast<unsigned>(content.size()));
    assert(n == static_cast<int>(content.size()));
    gzclose(g);
}

static std::string slurp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

// GzLineReader has std::getline's contract on plain and gzip files alike:
// empty lines are lines, the last line needs no newline, a trailing newline
// does not make an extra empty line, and a line longer than the read buffer
// comes back whole.
void testGzLineReader() {
    std::cerr << "START: Test GzLineReader... \t";
    std::string longLine(200000, 'A');
    std::string content = "first\n\nthird\r\n" + longLine + "\nlast";
    { std::ofstream o("gzline_plain.txt", std::ios::binary); o << content; }
    writeGzip("gzline.txt.gz", content);
    std::vector<std::string> expect = {"first", "", "third\r", longLine, "last"};
    for (const char* p : {"gzline_plain.txt", "gzline.txt.gz"}) {
        GzLineReader r(p);
        assert(r);
        std::vector<std::string> got;
        for (std::string l; r.getline(l);) got.push_back(l);
        assert(got == expect);
        std::string scratch;
        assert(!r.getline(scratch));             // stays at end
    }
    { std::ofstream o("gzline_nl.txt"); o << "a\nb\n"; }
    {
        GzLineReader r("gzline_nl.txt");
        std::vector<std::string> got;
        for (std::string l; r.getline(l);) got.push_back(l);
        assert(got == std::vector<std::string>({"a", "b"}));
    }
    GzLineReader missing("gzline_missing.txt");
    assert(!missing);
    std::string l;
    assert(!missing.getline(l));

    // uncompressedSize: the plain size, or the gzip trailer's.
    assert(uncompressedSize("gzline_plain.txt") == content.size());
    assert(uncompressedSize("gzline.txt.gz") == content.size());
    assert(uncompressedSize("gzline_missing.txt") == 0);

    for (const char* f : {"gzline_plain.txt", "gzline.txt.gz", "gzline_nl.txt"}) std::remove(f);
    std::cerr << "PASSED: Test GzLineReader\n";
}

// A gzipped reference is scanned like the plain one: same candidates, same
// name, same count in the stat file. Before this was true the gzipped genome
// silently got no candidates and the others' pools were inflated by everything
// they shared with it.
void testGzipReference(UniqSketchConfig cfg) {
    std::cerr << "START: Test gzipped reference... \t";
    writeGzip("ref_test_2.fa.gz", slurp("ref_test_2.fa"));
    cfg.refstat = "db_uniq_count_gz.tsv";
    std::vector<std::string> refs = {"ref_test_1.fa", "ref_test_2.fa.gz"};
    BloomFilter dbFilter(cfg.m1(), cfg.nhash1, cfg.kmerLen);
    identifyUniqKmers(cfg, refs, dbFilter);
    assert(candidateLines("ref_test_1.cand") == fileLines("expected_ref_test_1.tsv"));
    assert(candidateLines("ref_test_2.cand") == fileLines("expected_ref_test_2.tsv"));   // named ref_test_2, not ref_test_2.fa
    std::vector<std::string> stat = fileLines(cfg.refstat);
    assert(stat.size() == 2 && stat[1] == "./ref_test_2.cand\t3\t1");
    std::remove("ref_test_2.fa.gz");
    std::remove(cfg.refstat.c_str());
    std::cerr << "PASSED: Test gzipped reference\n";
}

// Candidate files round-trip, including a contig with no candidates and one
// long enough to cross the writer's buffer boundary.
void testCandidateFile() {
    std::cerr << "START: Test candidate file round trip... \t";
    const std::string path = "cand_roundtrip.cand";
    {
        CandidateWriter w(path, 21, "some/ref.fa");
        w.beginContig("chr1_first_contig", 500000); for (uint32_t p = 0; p < 20000; p++) w.add(p * 7); w.endContig();
        w.beginContig("chr2", 12); w.endContig();
        w.beginContig("plasmid", 99999); w.add(5); w.add(6); w.endContig();
        assert(w.good());
    }
    CandidateSet c;
    bool ok = readCandidates(path, c);
    assert(ok);
    assert(c.k == 21 && c.refPath == "some/ref.fa" && c.total == 20002);
    assert(c.runs.size() == 3);
    assert(c.runs[0].name == "chr1_first_contig" && c.runs[0].seqLen == 500000 && c.runs[0].pos.size() == 20000);
    assert(c.runs[0].pos[19999] == 19999 * 7);
    assert(c.runs[1].name == "chr2" && c.runs[1].pos.empty());
    assert(c.runs[2].pos == std::vector<uint32_t>({5, 6}));
    // Not a candidate file: refused, not misread.
    { std::ofstream junk("cand_junk.cand"); junk << "AGGCTACAAA\t1\t0\n"; }
    CandidateSet j;
    bool junkOk = readCandidates("cand_junk.cand", j);
    assert(!junkOk);
    std::remove(path.c_str()); std::remove("cand_junk.cand");
    std::cerr << "PASSED: Test candidate file round trip\n";
}

void testBuildSketch(const UniqSketchConfig& cfg) {
    std::cerr << "START: Test buildSketch... \t";

    /*
    With sketchnum = 2, each ref gets 2 unique k-mers spread along its genome:
    the two 11-bp contigs make two 11-bp windows, and the first candidate of
    each is taken.
    ref_test_1: {AGGCTACAAA (contig 0), GAGACTGAAC (contig 1)}
    ref_test_2: {CTCAGCTAAG (contig 3), AGGCTACAAC (contig 4)}
    Total sketch: 4 unique k-mers. (Slots spaced over the candidate list, as
    selection used to do, took two from contig 3 and none from contig 4.)
    */
    candidatesFromText("expected_ref_test_1.tsv", "ref_test_1.fa", "expected_ref_test_1.cand", cfg.kmerLen);
    candidatesFromText("expected_ref_test_2.tsv", "ref_test_2.fa", "expected_ref_test_2.cand", cfg.kmerLen);
    BloomFilter dbFilter(cfg.m1(), cfg.nhash1, cfg.kmerLen);
    buildSketch(cfg, "expected_db_uniq_count.tsv", dbFilter);
    assert(compareFiles(cfg.outfile, "expected_out_uniqsketch.txt") == true);

    std::cerr << "PASSED: Test buildSketch\n";
}

void testLowComplexity() {
    std::cerr << "START: Test lowComplexity... \t";

    // These are 81-mers, so the entropy ceiling is that of the default k=81
    // configuration. Previously the test had to hard-code the ceiling (24.56)
    // because buildSketch computed it as a side effect on a global; now it
    // follows from the config.
    UniqSketchConfig cfg;   // kmerLen = 81
    assert(std::abs(cfg.maxEntropy() - 24.56) < 0.01);

    std::string sig1("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCT");
    std::string sig2("CGTTGCATGGCATGCGGTTCGTCTGCGCGTGCGTGCTGATGCAGTAAAAAAGCTGCATAAAATCGATCGAAAGTACATGCG");

    assert(lowComplexity(cfg, sig1) == true);
    assert(lowComplexity(cfg, sig2) == false);

    std::cerr << "PASSED: Test lowComplexity\n";
}

void testHasLongHomopolymer() {
    std::cerr << "START: Test hasLongHomopolymer... \t";

    // Filter disabled (0 = off): even a pure homopolymer must pass.
    UniqSketchConfig off;
    off.maxHomopolymer = 0;
    assert(hasLongHomopolymer(off, "AAAAAAAAAAAA") == false);

    // Two configs side by side: no shared state, no restore step needed.
    UniqSketchConfig cfg;
    cfg.maxHomopolymer = 5;

    // Boundary: a run of exactly N is allowed, N+1 is rejected.
    assert(hasLongHomopolymer(cfg, "AAAAA") == false);
    assert(hasLongHomopolymer(cfg, "AAAAAA") == true);

    // Run position must not matter (start, middle, end).
    assert(hasLongHomopolymer(cfg, "GGGGGGCTACTA") == true);
    assert(hasLongHomopolymer(cfg, "CTAGGGGGGCTA") == true);
    assert(hasLongHomopolymer(cfg, "CTACTAGGGGGG") == true);

    // Run length resets on a base change: two runs of 5 are still acceptable.
    assert(hasLongHomopolymer(cfg, "AAAAACAAAAA") == false);

    // Mixed sequence with only short runs.
    assert(hasLongHomopolymer(cfg, "AACCGGTTAACCGGTT") == false);

    // Degenerate inputs.
    assert(hasLongHomopolymer(cfg, "") == false);
    assert(hasLongHomopolymer(cfg, "A") == false);

    std::cerr << "PASSED: Test hasLongHomopolymer\n";
}

void testGetCanonical() {
    std::cerr << "START: Test get_canonical... \t";

    assert(get_canonical("GTTTATTTACC") == "GGTAAATAAAC");
    assert(get_canonical("AGGTTAGCTTTC") == "AGGTTAGCTTTC");

    std::cerr << "PASSED: Test get_canonical\n";
}

void testGetBaseId() {
    std::cerr << "START: Test getBaseId... \t";

    assert(getBaseId("/home/user/workspace/test/ref1.fasta") == "ref1");
    assert(getBaseId("ref2.fa") == "ref2");
    // Multiple dots — only strip last extension
    assert(getBaseId("/path/to/ref.test.fa") == "ref.test");
    // No directory separator
    assert(getBaseId("genome.fna") == "genome");
    // No extension
    assert(getBaseId("refonly") == "refonly");
    assert(getBaseId("/path/to/refonly") == "refonly");
    // A compression suffix is not the extension
    assert(getBaseId("ref.fa.gz") == "ref");
    assert(getBaseId("/path/to/genome.fna.gz") == "genome");
    assert(getBaseId("ref.fa.bgz") == "ref");
    assert(getBaseId("ref.gz") == "ref");

    std::cerr << "PASSED: Test getBaseId\n";
}

// NTMC64sub must reproduce NTMC64 bit for bit on every single-base substitution.
// isTwoSafe relies on that to test a candidate's 3k neighbours without
// rehashing, so an error here would change which signatures are selected.
// Exercised across k values with different remainders mod 31 and 33 (the
// rotation tables are indexed by those), at every position and for every
// base, including a no-op substitution, lowercase bases and N.
void testSubstitutionHash() {
    std::cerr << "START: Test NTMC64sub matches NTMC64... \t";
    const unsigned m = 5;
    const char* alphabet = "ACGTACGTACGTACGTacgtN";   // mostly upper case
    uint64_t x = 20260923;
    auto next = [&x]() { x = x * 6364136223846793005ULL + 1442695040888963407ULL; return x >> 33; };
    size_t checked = 0;
    for (unsigned k : {21u, 31u, 33u, 62u, 64u, 81u, 100u}) {
        for (int rep = 0; rep < 40; rep++) {
            std::string s(k, 'A');
            for (unsigned i = 0; i < k; i++) s[i] = alphabet[next() % 21];
            uint64_t fh = 0, rh = 0;
            NTC64(s.c_str(), k, fh, rh);
            uint64_t expect[m], got[m];
            for (unsigned i = 0; i < k; i++) {
                for (char b : {'A', 'C', 'G', 'T', 'N', 'a'}) {
                    std::string t = s;
                    t[i] = b;
                    NTMC64(t.c_str(), k, m, expect);
                    NTMC64sub(fh, rh, i, s[i], b, k, m, got);
                    for (unsigned j = 0; j < m; j++) assert(got[j] == expect[j]);
                    ++checked;
                }
            }
        }
    }
    std::cerr << "PASSED: Test NTMC64sub matches NTMC64 (" << checked << " substitutions)\n";
}

// forEachKmerAhead must hand use() exactly what a plain ntHashIterator loop
// yields: the same k-mers, hash values and positions, in the same order, each one
// requested before it is used and never more than KMER_LOOKAHEAD requests ahead.
// The passes rewritten with it rely on that to perform the same operations in the
// same order. Sequences shorter than k, of exactly k, around the look-ahead
// distance, with runs of N (the iterator skips them) and lower case.
void testKmerLookAhead() {
    std::cerr << "START: Test k-mer look-ahead keeps the order of a plain scan... \t";
    uint64_t x = 99;
    auto next = [&x]() { x = x * 6364136223846793005ULL + 1442695040888963407ULL; return x >> 33; };
    auto dna = [&](size_t n) { std::string d(n, 'A'); for (char& c : d) c = "ACGT"[next() % 4]; return d; };
    const unsigned k = 21;
    std::vector<std::string> seqs = {"", dna(k - 1), dna(k), dna(k + 1)};
    for (size_t extra : {KMER_LOOKAHEAD - 1, KMER_LOOKAHEAD, KMER_LOOKAHEAD + 1, 2 * KMER_LOOKAHEAD + 3, 5000u})
        seqs.push_back(dna(k + extra - 1));                // `extra` k-mers
    seqs.push_back(dna(300) + std::string(40, 'N') + dna(10) + "N" + dna(200));
    seqs.push_back(std::string(k, 'N'));
    std::string lower = dna(400);
    for (size_t i = 100; i < 200; i++) lower[i] = static_cast<char>(std::tolower(lower[i]));
    seqs.push_back(lower);
    size_t kmers = 0;
    for (const std::string& seq : seqs) {
        for (unsigned nh : {1u, 4u, 5u}) {
            std::vector<std::pair<size_t, std::vector<uint64_t>>> expect, requested, used;
            for (ntHashIterator itr(seq, nh, k); itr != itr.end(); ++itr)
                expect.push_back({itr.get_pos(), std::vector<uint64_t>(*itr, *itr + nh)});
            forEachKmerAhead(seq, nh, k,
                [&](const uint64_t* hVal) {
                    requested.push_back({0, std::vector<uint64_t>(hVal, hVal + nh)});
                    assert(requested.size() - used.size() <= KMER_LOOKAHEAD);
                },
                [&](const uint64_t* hVal, size_t pos) {
                    used.push_back({pos, std::vector<uint64_t>(hVal, hVal + nh)});
                    assert(used.size() <= requested.size());                      // requested first
                    assert(requested[used.size() - 1].second == used.back().second);
                });
            assert(used == expect);
            assert(requested.size() == expect.size());
            kmers += expect.size();
        }
    }
    assert(kmers > 15000);
    std::cerr << "PASSED: Test k-mer look-ahead keeps the order of a plain scan\n";
}

// isTwoSafe against a brute-force definition on a small universe: a candidate
// is 2-safe iff none of its Hamming-1 neighbours is in the universe. The
// filter is sized so that false positives cannot occur at this scale.
void testIsTwoSafe() {
    std::cerr << "START: Test isTwoSafe... \t";
    UniqSketchConfig cfg;
    cfg.kmerLen = 31;
    cfg.dbfSize = 100000;
    BloomFilter dbf(cfg.m1(), cfg.nhash1, cfg.kmerLen);

    uint64_t x = 777;
    auto next = [&x]() { x = x * 6364136223846793005ULL + 1442695040888963407ULL; return x >> 33; };
    auto randKmer = [&]() { std::string s(cfg.kmerLen, 'A'); for (auto& c : s) c = "ACGT"[next() % 4]; return s; };

    // Universe: 200 random k-mers plus, for 50 further candidates, exactly one
    // Hamming-1 neighbour each (at a random position, random other base).
    std::vector<std::string> safe, unsafe;
    std::vector<unsigned> plantedAt;                       // the neighbour's position, per unsafe candidate
    for (int i = 0; i < 200; i++) { std::string u = randKmer(); dbf.insert(u.c_str()); safe.push_back(u); }
    for (int i = 0; i < 50; i++) {
        std::string c = randKmer(), n = c;
        unsigned pos = next() % cfg.kmerLen;
        char b; do { b = "ACGT"[next() % 4]; } while (b == c[pos]);
        n[pos] = b;
        dbf.insert(n.c_str());
        unsafe.push_back(c);
        plantedAt.push_back(pos);
    }
    // A universe k-mer's own presence is not a neighbour; a random 31-mer has
    // no neighbour in a 250-k-mer universe, so the first set is 2-safe ...
    for (const std::string& s : safe) assert(isTwoSafe(cfg, s, dbf));
    // ... and every candidate with a planted neighbour is not.
    for (const std::string& s : unsafe) assert(!isTwoSafe(cfg, s, dbf));
    // Reverse complement of a planted neighbour is the same k-mer to ntHash.
    for (const std::string& s : unsafe) {
        std::string rc(s.rbegin(), s.rend());
        for (auto& c : rc) c = (c == 'A') ? 'T' : (c == 'T') ? 'A' : (c == 'C') ? 'G' : 'C';
        assert(!isTwoSafe(cfg, rc, dbf));
    }
    // The hint says where to probe first and never changes the answer: with every
    // offset, with none (-1) and with one past the k-mer, a candidate is judged as
    // without a hint. On a failure the hint comes back as the offset whose neighbour
    // is present (there is exactly one, the planted one), on success as -1.
    for (size_t i = 0; i < unsafe.size(); i++) {
        for (int h = -1; h <= static_cast<int>(cfg.kmerLen); h++) {
            int at = h;
            assert(!isTwoSafe(cfg, unsafe[i], dbf, &at));
            assert(at == static_cast<int>(plantedAt[i]));
        }
    }
    for (size_t i = 0; i < 20; i++) {
        for (int h = -1; h <= static_cast<int>(cfg.kmerLen); h++) {
            int at = h;
            assert(isTwoSafe(cfg, safe[i], dbf, &at));
            assert(at == -1);
        }
    }
    std::cerr << "PASSED: Test isTwoSafe\n";
}

// The order in which selection tests a candidate (margin before entropy, see
// examine in getUniqSet) must not change what it picks. With a quota of one the
// genome is a single window and the rule is short enough to state here: the first
// candidate that is usable (entropy, homopolymer) and 2-safe; failing that, in
// soft-margin mode, the first usable one; without a margin, the first usable one.
// The genomes put the cases selection could confuse ahead of the candidate it must
// take: a low-complexity stretch that is 2-safe (unusable, and not a fallback
// either), a stretch a substitution away from the other reference (usable, fails
// the margin), then ordinary sequence. Each expectation is computed from the
// candidate file with the filters called one at a time.
void testSelectionFilterOrder() {
    std::cerr << "START: Test selection picks the same candidate whatever the filter order... \t";
    uint64_t x = 424242;
    auto next = [&x]() { x = x * 6364136223846793005ULL + 1442695040888963407ULL; return x >> 33; };
    auto dna = [&](size_t n) { std::string d(n, 'A'); for (char& c : d) c = "ACGT"[next() % 4]; return d; };
    // Low entropy but not periodic, so its k-mers are distinct: A with a C, G or T
    // at irregular intervals.
    auto lowEntropy = [&](size_t n) {
        std::string d(n, 'A');
        for (size_t i = 5; i < n; i += 9 + next() % 7) d[i] = "CGT"[next() % 3];
        return d;
    };
    UniqSketchConfig cfg;
    cfg.kmerLen = 31;
    cfg.bits = 128;
    cfg.dbfSize = 20000;
    cfg.sbfSize = 2000;
    cfg.sketchnum = 1;                                      // one window: the whole genome
    cfg.outdir = ".";
    cfg.refstat = "sel_db_uniq_count.tsv";
    cfg.outfile = "sel_sketch.tsv";

    // A copy with one substitution every k bases: every k-mer of the original then
    // has exactly one base different from the copy's k-mer at the same place, so
    // none of them is 2-safe once the copy is in the universe.
    auto relative = [&](std::string d) {
        for (size_t i = cfg.kmerLen / 2; i < d.size(); i += cfg.kmerLen) d[i] = (d[i] == 'A') ? 'C' : 'A';
        return d;
    };
    // sel_A: a low-complexity stretch (2-safe), a stretch the relative makes fail
    //        the margin, then ordinary sequence.
    // sel_B: a low-complexity stretch and ordinary sequence, all of it one
    //        substitution from the relative: nothing 2-safe, and the first
    //        candidates are unusable as well as margin-failing.
    // sel_C: the relative.
    const std::string low = lowEntropy(200), near = dna(300), plain = dna(400);
    const std::string wholeB = lowEntropy(200) + dna(300);
    const std::vector<std::string> refs = {"sel_A.fa", "sel_B.fa", "sel_C.fa"};
    { std::ofstream o(refs[0]); o << ">a\n" << low << near << plain << "\n"; }
    { std::ofstream o(refs[1]); o << ">b\n" << wholeB << "\n"; }
    { std::ofstream o(refs[2]); o << ">c\n" << relative(near) << dna(60) << relative(wholeB) << "\n"; }
    cfg.numRef = refs.size();

    for (int mode = 0; mode < 3; mode++) {                  // 0: no margin, 1: soft margin 2, 2: strict margin 2
        cfg.minMargin = mode == 0 ? 1 : 2;
        cfg.strictMargin = mode == 2;
        BloomFilter dbFilter(cfg.m1(), cfg.nhash1, cfg.kmerLen);
        identifyUniqKmers(cfg, refs, dbFilter);
        buildSketch(cfg, cfg.refstat, dbFilter);
        std::map<std::string, std::string> got;             // reference -> its one signature
        for (const std::string& line : fileLines(cfg.outfile)) {
            std::istringstream ls(line);
            std::string seq, id, name;
            ls >> seq >> id >> name;
            assert(!got.count(name));
            got[name] = seq;
        }
        size_t unusableSafe = 0, usableUnsafe = 0;
        for (const std::string& ref : refs) {
            CandidateSet cand;
            assert(readCandidates(candidatePath(cfg.outdir, ref), cand));
            std::string contig;
            scanFasta(ref, [&](const std::string&, std::string& seq, size_t) { contig = seq; });
            std::string firstSafe, firstUsable;
            for (uint32_t pos : cand.runs.at(0).pos) {
                const std::string kmer = get_canonical(contig.substr(pos, cfg.kmerLen));
                const bool usable = !lowComplexity(cfg, kmer) && !hasLongHomopolymer(cfg, kmer);
                const bool safeKmer = isTwoSafe(cfg, kmer, dbFilter);
                unusableSafe += !usable && safeKmer;
                usableUnsafe += usable && !safeKmer;
                if (usable && firstUsable.empty()) firstUsable = kmer;
                if (usable && safeKmer && firstSafe.empty()) firstSafe = kmer;
            }
            const std::string expect = mode == 0 ? firstUsable
                                     : mode == 1 ? (firstSafe.empty() ? firstUsable : firstSafe)
                                     : firstSafe;
            const std::string name = getBaseId(ref);
            assert((expect.empty() && !got.count(name)) || (got.count(name) && got[name] == expect));
        }
        // The fixture holds what it is meant to: candidates each filter alone rejects.
        assert(unusableSafe > 0 && usableUnsafe > 0);
        // sel_B has no 2-safe k-mer: nothing in strict mode, and in soft mode its
        // first usable candidate, which is not its first candidate.
        if (mode == 2) assert(!got.count("sel_B") && got.count("sel_A"));
        if (mode == 1) assert(got.count("sel_B") && !lowComplexity(cfg, got["sel_B"]) &&
                              got["sel_B"] != get_canonical(wholeB.substr(0, cfg.kmerLen)));
    }
    for (const std::string& ref : refs) {
        std::remove(ref.c_str());
        std::remove(candidatePath(cfg.outdir, ref).c_str());
    }
    std::remove(cfg.refstat.c_str());
    std::remove(cfg.outfile.c_str());
    std::cerr << "PASSED: Test selection picks the same candidate whatever the filter order\n";
}

int main() {
    // Shared fixture for the two-reference, k=10 index tests. m1/m2 are now
    // derived (bits * size) rather than set by hand.
    UniqSketchConfig cfg;
    cfg.kmerLen = 10;
    cfg.bits = 128;
    cfg.dbfSize = 7;
    cfg.sbfSize = 1;
    cfg.numRef = 2;
    cfg.sketchnum = 2;
    cfg.outfile = "test_out_uniqsketch.txt";
    cfg.outdir = ".";

    std::vector<std::string> refFiles = {"ref_test_1.fa", "ref_test_2.fa"};

    testIdentifyUniqKmers(cfg, refFiles);
    testGroupRepeatExclusion(cfg);
    testGroupScanIdentical();
    testKmerLookAhead();
    testGzLineReader();
    testGzipReference(cfg);
    testCandidateFile();
    testBuildSketch(cfg);
    testLowComplexity();
    testHasLongHomopolymer();
    testGetCanonical();
    testGetBaseId();
    testSubstitutionHash();
    testIsTwoSafe();
    testSelectionFilterOrder();

    std::cerr << "UniqSketch: All tests PASSED!\n\n";
    return 0;
}