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
    scanFasta(cand.refPath, [&](const std::string&, std::string& seq) { contigs.push_back(std::move(seq)); });
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
    scanFasta(refPath, [&](const std::string& header, const std::string& seq) {
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

// loadBFfaGroupParallel must leave the distinct and solid filters bit for bit as
// loadBFfaGroup does: the partitioned scan replays the serial order of
// operations on every block, so the filters' false positives, and through them
// the candidates, cannot depend on the thread count. The group here exercises
// every path: a contig longer than one hashing piece, several flushes (a small
// flush size), the representative's internal repeat, a k-mer shared by two
// members but not the representative, k-mers shared by the whole group, runs
// of N, a contig shorter than k, and a gzip member.
void testGroupParallelScanIdentical() {
    std::cerr << "START: Test parallel group scan identical to serial... \t";
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
    const std::vector<std::string> files = {"pgroup_rep.fa", "pgroup_m1.fa", "pgroup_m2.fa",
                                            "pgroup_m3.fa", "pgroup_m4.fa.gz"};
    {
        std::ofstream rep(files[0]);
        rep << ">core\n" << core << "\n"
            << ">rep_dup\n" << dna(20000) << repeat << dna(20000) << repeat << dna(3000) << "\n"
            << ">with_n\n" << dna(4000) << std::string(50, 'N') << dna(4000) << "\n"
            << ">short\nACGTACGT\n"
            << ">tail\n" << dna(30000) << "\n";
        std::ofstream m1(files[1]);
        m1 << ">m1_core\n" << mutate(core, 997) << "\n>m1_own\n" << dna(100000) << pairShared << "\n";
        std::ofstream m2(files[2]);
        m2 << ">m2_core\n" << mutate(core, 1009) << "\n>m2_own\n" << pairShared << dna(90000) << "\n";
        std::ofstream m3(files[3]);
        m3 << ">m3_core\n" << mutate(core, 1013) << "\n>m3_own\n" << dna(110000) << "\n";
        std::ostringstream m4;
        m4 << ">m4_core\n" << mutate(core, 1019) << "\n>m4_own\n" << dna(70000) << "\n";
        gzFile g = gzopen(files[4].c_str(), "wb");
        assert(g);
        const std::string m4s = m4.str();
        assert(gzwrite(g, m4s.data(), static_cast<unsigned>(m4s.size())) == static_cast<int>(m4s.size()));
        gzclose(g);
    }
    UniqSketchConfig cfg;
    cfg.kmerLen = 31;
    cfg.bits = 32;                                         // a filled filter, so false positives occur
    cfg.dbfSize = 1200000;
    cfg.sbfSize = 700000;
    auto scan = [&](bool parallel, size_t flushBases) {
        BloomFilter db(cfg.m1(), cfg.nhash1, cfg.kmerLen);
        BloomFilter sb(cfg.m2(), cfg.nhash2, cfg.kmerLen);
        if (parallel) loadBFfaGroupParallel(cfg, files, db, sb, flushBases);
        else loadBFfaGroup(cfg, files, db, sb);
        std::string bits(reinterpret_cast<const char*>(db.getFilter()), db.sizeBytes());
        bits.append(reinterpret_cast<const char*>(sb.getFilter()), sb.sizeBytes());
        assert(db.get_pop() > 0 && sb.get_pop() > 0);
        return bits;
    };
    const std::string serial = scan(false, 0);
#ifdef _OPENMP
    const int savedThreads = omp_get_max_threads();
    omp_set_num_threads(8);
#endif
    assert(scan(true, 100000) == serial);                 // several flushes within and across files
    assert(scan(true, size_t(32) << 20) == serial);       // the default: one flush per side
#ifdef _OPENMP
    omp_set_num_threads(savedThreads);
#endif
    for (const std::string& f : files) std::remove(f.c_str());
    std::cerr << "PASSED: Test parallel group scan identical to serial\n";
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
    With sketchnum = 2, each ref gets 2 unique kmers uniformly selected:
    ref_test_1: {AGACTGAACT, AGGCTACAAA}
    ref_test_2: {AGGCTACAAC, CTCAGCTAAG}
    Total sketch: 4 unique kmers
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
    for (int i = 0; i < 200; i++) { std::string u = randKmer(); dbf.insert(u.c_str()); safe.push_back(u); }
    for (int i = 0; i < 50; i++) {
        std::string c = randKmer(), n = c;
        unsigned pos = next() % cfg.kmerLen;
        char b; do { b = "ACGT"[next() % 4]; } while (b == c[pos]);
        n[pos] = b;
        dbf.insert(n.c_str());
        unsafe.push_back(c);
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
    std::cerr << "PASSED: Test isTwoSafe\n";
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
    testGroupParallelScanIdentical();
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

    std::cerr << "UniqSketch: All tests PASSED!\n\n";
    return 0;
}