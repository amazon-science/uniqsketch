#include <getopt.h>
#include <sys/stat.h>

#include <algorithm>
#include <random>
#include <unordered_map>

#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"
#include "ntcard.hpp"
#include "UniqSketchUtil.hpp"
#include "ClusterUtil.hpp"
#include "TreeSketchUtil.hpp"
#include "TreeFileUtil.hpp"
#include "GroupSketchUtil.hpp"
#include "NodeCarrierUtil.hpp"
#include "FileUtil.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "uniqsketch"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 2.0.0\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION] @LIST_FILES (or FILES)\n"
    "Creates unique sketch from a list of fasta reference files (plain or gzip).\n"
    "A list of files containing file names in each row can be passed with @ prefix.\n"
    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N\tuse N parallel threads [1]\n"
    "  -k, --kmer=N\t\tthe length of kmer [81]\n"
    "  -b, --bit=N\t\tuse N bits per element in Bloom filter [128]\n"
    "  -d, --hash1=N\t\tdistinct Bloom filter hash number [5]\n"
    "  -s, --hash2=N\t\trepeat Bloom filter hash number [5]\n"
    "  -c, --cov=N\t\tnumber of unique k-mers to represent a reference [100]\n"
    "  -f, --outdir=STRING\tdir for universe unique k-mers [outdir_uniqsketch]\n"
    "  -o, --out=STRING\tthe output sketch file name [sketch_uniq.tsv]\n"
    "  -r, --stat=STRING\tthe output unique kmer stat file name [db_uniq_count.tsv]\n"
    "  -e, --entropy=F\tsets the aggregate entropy rate threshold [0.65]\n"
    "      --cluster=N\tcluster similar references with N unique k-mer threshold [0=off]\n"
    "      --tree=T1,T2,..\tbuild hierarchical tree index with comma-separated thresholds\n"
    "\t\t\t(k-mer counts as plain integers, e.g. 1000000,300000; the list is required)\n"
    "      --tree-file=FILE\ttake the tree from FILE, in the format of <outdir>/tree_index.tsv, instead of\n"
    "\t\t\tclustering the references (with --tree, whose thresholds then only name the levels);\n"
    "\t\t\ta node's name and path in FILE are not read, a node with one child is refused (README)\n"
    "      --exact\t\tuse exact O(n^2) clustering instead of MinHash (with --cluster)\n"
    "      --minhash-size=N\thashes in each reference's MinHash sketch, from which --tree and --cluster (without\n"
    "\t\t\t--exact) estimate distances [2000]. A hash stands for 1/N of a genome's k-mers (2 500\n"
    "\t\t\tof 5 Mb at 2000, 250 at 20000). Another N can give another tree, hence another\n"
    "\t\t\tindex (README); keep N below the distinct k-mers of the smallest reference\n"
    "      --min-margin=N\tmin Hamming distance of signatures to other refs (1=off, 2) [1; 2 with --tree]\n"
    "      --max-homopolymer=N\tdrop signatures with a homopolymer run > N bp (0=off) [0]\n"
    "      --strict-margin\tenforce --min-margin as a hard requirement (no margin-1 fallback; may yield fewer signatures) [on with --tree unless --min-margin=1]\n"
    "      --circular\treferences are circular: also take the k-1 k-mers spanning each contig's end-to-start join\n"
    "      --no-group-signatures\tdo not build within-group leaf signatures (with --tree; see README)\n"
    "      --node-carriers=N\tfill a tree node's signature quota, when its node-unique k-mers leave it short,\n"
    "\t\t\tfrom k-mers carried by up to N genomes outside the node (recorded; the query\n"
    "\t\t\tsets them aside when a carrier is present) [min(50, leaves/10); 0=off]\n"
    "      --node-representative=WHICH\tthe member a tree node's signatures come from: central (the one whose\n"
    "\t\t\tnode-unique k-mers the other members carry) or largest [central]\n"
    "      --discard-candidates\tremove each reference's candidate file (<name>.cand, .carried under -f) once its\n"
    "\t\t\tsignatures are selected; they are scratch for the build and no query reads them\n"
    "      --sensitive\tsets sensitivity parameter c to 100 (an explicit -c wins)\n"
    "      --very-sensitive\tsets sensitivity parameter c to 1000 (an explicit -c wins)\n"
    "      --help\t\tdisplay this help and exit\n"
    "      --version\t\toutput version information and exit\n"
    "\n";

static const char shortopts[] = "t:k:b:d:s:c:f:o:r:e:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_CLUSTER, OPT_MINMARGIN, OPT_EXACT, OPT_TREE,
       OPT_MAXHOMOP, OPT_STRICTMARGIN, OPT_NOGROUPSIGS , OPT_CIRCULAR, OPT_NODECARRIERS, OPT_NODEREP, OPT_DISCARDCAND,
       OPT_TREEFILE, OPT_MINHASHSIZE };

// --sensitive / --very-sensitive are getopt flag entries, which store through a
// pointer at parse time and so cannot target a struct member. They write here and
// the value is copied into the config once parsing is complete (0: neither given).
// -c is kept apart from them: when it shared this variable the last of -c and the
// flags on the command line won, and "-c500 --sensitive" built with c = 100.
static int g_sketchnum = 0;

static const struct option longopts[] = {
    {"threads",        required_argument, nullptr, 't'},
    {"kmer",           required_argument, nullptr, 'k'},
    {"bit",            required_argument, nullptr, 'b'},
    {"bits",           required_argument, nullptr, 'b'},
    {"hash1",          required_argument, nullptr, 'd'},
    {"hash2",          required_argument, nullptr, 's'},
    {"cov",            required_argument, nullptr, 'c'},
    {"outdir",         required_argument, nullptr, 'f'},
    {"out",            required_argument, nullptr, 'o'},
    {"stat",           required_argument, nullptr, 'r'},
    // --entropy was declared without an argument while -e takes one: --entropy=0.9
    // was refused, "--entropy 0.9" read 0.9 as a reference file, and a bare
    // --entropy did nothing.
    {"entropy",        required_argument, nullptr, 'e'},
    {"cluster",        required_argument, nullptr, OPT_CLUSTER},
    {"tree",           required_argument, nullptr, OPT_TREE},
    {"tree-file",      required_argument, nullptr, OPT_TREEFILE},
    {"min-margin",     required_argument, nullptr, OPT_MINMARGIN},
    {"exact",          no_argument,       nullptr, OPT_EXACT},
    {"minhash-size",   required_argument, nullptr, OPT_MINHASHSIZE},
    {"max-homopolymer",required_argument, nullptr, OPT_MAXHOMOP},
    {"strict-margin",  no_argument,       nullptr, OPT_STRICTMARGIN},
    {"circular",       no_argument,       nullptr, OPT_CIRCULAR},
    {"no-group-signatures", no_argument,  nullptr, OPT_NOGROUPSIGS},
    {"node-carriers",  required_argument, nullptr, OPT_NODECARRIERS},
    {"node-representative", required_argument, nullptr, OPT_NODEREP},
    {"discard-candidates", no_argument,   nullptr, OPT_DISCARDCAND},
    {"sensitive",      no_argument, &g_sketchnum, 100},
    {"very-sensitive", no_argument, &g_sketchnum, 1000},
    {"help",           no_argument,       nullptr, OPT_HELP},
    {"version",        no_argument,       nullptr, OPT_VERSION},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char** argv) {
#ifdef _OPENMP
    double sTime = omp_get_wtime();
#endif

    UniqSketchConfig cfg;
    unsigned clusterThreshold = 0;  // 0 = clustering disabled
    bool useExactCluster = false;   // --exact: use exact O(n^2) clustering instead of MinHash
    std::vector<unsigned> treeThresholds;  // empty = tree disabled
    std::vector<TreeNode> tree;            // the hierarchy, when --tree is used
    std::string treeFileArg;               // --tree-file: the hierarchy is read from this file, not clustered
    bool treeFileGiven = false;
    bool marginGiven = false;       // --min-margin passed explicitly (tree default otherwise)
    bool groupSignatures = true;    // within-group leaf signatures for a tree index (--no-group-signatures)
    int nodeCarriers = -1;          // --node-carriers (-1: min(50, leaves/10) once the leaf count is known)
    bool centralRepresentative = true;   // --node-representative=central|largest
    int explicitCov = 0;            // -c
    bool covGiven = false;
    bool clusterGiven = false;      // --cluster
    unsigned minhashSize = MINHASH_SKETCH_SIZE;   // --minhash-size: hashes in a reference's MinHash sketch
    bool minhashGiven = false;
    std::vector<std::string> treeOnly;   // options given that only a --tree build reads

    bool die = false;
    int longIndex = -1;     // set by getopt_long for a long option only, so reset for every option
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, &longIndex)) != -1; longIndex = -1) {
        std::istringstream arg(optarg != nullptr ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 't': arg >> cfg.threads;          break;
        case 'k': arg >> cfg.kmerLen;          break;
        case 'b': arg >> cfg.bits;             break;
        case 'd': arg >> cfg.nhash1;           break;
        case 's': arg >> cfg.nhash2;           break;
        case 'c': arg >> explicitCov; covGiven = true; break;
        case 'o': arg >> cfg.outfile;          break;
        case 'f': arg >> cfg.outdir;           break;
        case 'r': arg >> cfg.refstat;          break;
        case 'e': arg >> cfg.entropyThreshold; break;
        case OPT_CLUSTER: arg >> clusterThreshold; clusterGiven = true; break;
        case OPT_MINMARGIN: arg >> cfg.minMargin; marginGiven = true; break;
        case OPT_EXACT: useExactCluster = true; break;
        case OPT_TREE: {
            std::string treeArg;
            arg >> treeArg;
            const std::string problem = parseTreeThresholds(treeArg, treeThresholds);
            if (!problem.empty()) {
                std::cerr << PROGRAM ": --tree=" << treeArg << ": " << problem
                          << "\n  --tree takes a list of thresholds (there is no default list): k-mer counts as plain "
                             "integers, comma-separated, e.g. --tree=1000000,300000,100000,10000\n";
                exit(EXIT_FAILURE);
            }
            // Sort descending (coarsest first)
            std::sort(treeThresholds.begin(), treeThresholds.end(), std::greater<unsigned>());
            break;
        }
        case OPT_TREEFILE: arg >> treeFileArg; treeFileGiven = true; break;
        case OPT_MINHASHSIZE: {
            std::string sizeArg;
            arg >> sizeArg;
            const std::string problem = parseMinhashSize(sizeArg, minhashSize);
            if (!problem.empty()) {
                std::cerr << PROGRAM ": --minhash-size=" << sizeArg << ": " << problem
                          << "\n  --minhash-size takes the number of hashes in a reference's MinHash sketch as a plain "
                             "positive integer, e.g. --minhash-size=20000 (default " << MINHASH_SKETCH_SIZE << ")\n";
                exit(EXIT_FAILURE);
            }
            minhashGiven = true;
            break;
        }
        case OPT_CIRCULAR: cfg.circular = true; break;
        case OPT_DISCARDCAND: cfg.keepCandidates = false; break;
        case OPT_MAXHOMOP: arg >> cfg.maxHomopolymer; break;
        case OPT_STRICTMARGIN: cfg.strictMargin = true; break;
        case OPT_NOGROUPSIGS: groupSignatures = false; treeOnly.push_back("--no-group-signatures"); break;
        case OPT_NODECARRIERS: arg >> nodeCarriers; treeOnly.push_back("--node-carriers"); break;
        case OPT_NODEREP: {
            treeOnly.push_back("--node-representative");
            std::string which;
            arg >> which;
            if (which == "central") centralRepresentative = true;
            else if (which == "largest") centralRepresentative = false;
            else { std::cerr << "Error: --node-representative must be central or largest\n"; die = true; }
            break;
        }
        case OPT_HELP:
            std::cerr << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cerr << VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        }
        if (optarg != nullptr && !arg.eof()) {
            // By name: a long option without a short form has a code below 32, which
            // printed as a control character.
            if (longIndex >= 0) std::cerr << PROGRAM ": invalid value for --" << longopts[longIndex].name;
            else std::cerr << PROGRAM ": invalid value for -" << static_cast<char>(c);
            std::cerr << ": `" << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }

    // An explicit -c wins over --sensitive / --very-sensitive wherever it stands.
    cfg.sketchnum = covGiven ? explicitCov : (g_sketchnum ? g_sketchnum : 100);
    if (covGiven && g_sketchnum && g_sketchnum != explicitCov)
        std::cerr << PROGRAM ": -c " << explicitCov << " given; --sensitive / --very-sensitive ignored\n";

    // A tree index exists to separate close relatives, and that is exactly where
    // margin-1 signatures fail: a signature one substitution away from a k-mer of
    // an abundant relative is hit by every read of that relative carrying a
    // sequencing error at that base. At a few thousand-fold coverage those reads
    // outnumber the reads of a 0.01% member, so its calls rest on noise and
    // absent members get called too. Requiring Hamming distance >= 2 makes a
    // false hit need two coincident errors, which removes the noise (on the
    // 8-phage benchmark, reads landing on absent references fell from 2093 to 2
    // over 28 samples) while keeping most of the candidate pool. Strict mode is
    // needed for that guarantee; the soft fallback would refill empty slots with
    // the very signatures being excluded. An explicit --min-margin=1 wins.
    if (!treeThresholds.empty() && !marginGiven) cfg.minMargin = 2;

    if (cfg.minMargin < 1) cfg.minMargin = 1;
    if (cfg.minMargin > 2) {
        std::cerr << PROGRAM ": --min-margin > 2 is not supported; using 2\n";
        cfg.minMargin = 2;
    }

    // A margin of 2 is strict under --tree whether it is the default or spelled
    // out. Only the absence of --min-margin used to switch the strict margin on,
    // so --tree with --min-margin=2, the value the usage text gives as the tree
    // default, built a soft-margin index without a word: margin-1 signatures in
    // the slots no 2-safe one filled (228 of 3 004 on the 8-phage set), and a
    // query kept at --rcutoff=2 (236 false rows over the 280 phage mixtures
    // against 2, and 36 contaminants found at 0.01 % against 44).
    if (!treeThresholds.empty() && cfg.minMargin >= 2) {
        if (!marginGiven)
            std::cerr << PROGRAM ": --tree defaults to --min-margin=2 --strict-margin "
                         "(pass --min-margin=1 to select signatures without a margin)\n";
        else if (!cfg.strictMargin)
            std::cerr << PROGRAM ": --tree enforces --min-margin=2 strictly, as its default does "
                         "(pass --min-margin=1 to select signatures without a margin)\n";
        cfg.strictMargin = true;
    }

    // Values no build can use. They were taken as given: -k 0, -b 0 or -d 0 ran
    // to the end and wrote an empty sketch, -c 0 and below one signature per
    // reference, all with exit status 0.
    if (cfg.kmerLen < 1 || cfg.bits < 1 || cfg.nhash1 < 1 || cfg.nhash2 < 1 || cfg.sketchnum < 1) {
        std::cerr << PROGRAM ": -k, -b, -d, -s and -c must be at least 1\n";
        die = true;
    }
    // --cluster reduces the reference list to cluster representatives after the tree
    // is built, so the leaf sketch named fewer references than the tree has leaves
    // and numbered them differently: a query of such an index reported a member's
    // reads under another member's name.
    if (!treeThresholds.empty() && clusterGiven && clusterThreshold > 0) {
        std::cerr << PROGRAM ": --cluster cannot be combined with --tree (the tree keeps every reference; "
                     "its nodes are the clusters)\n";
        die = true;
    }
    // --tree-file replaces the clustering step of a tree build and nothing else.
    // The thresholds of --tree still say how many levels the index has and are
    // recorded with it, so the option needs them; a flat or --cluster build has
    // no tree to read.
    if (treeFileGiven) {
        if (treeFileArg.empty()) {
            std::cerr << PROGRAM ": --tree-file needs the name of a file\n";
            die = true;
        }
        if (treeThresholds.empty()) {
            std::cerr << PROGRAM ": --tree-file needs --tree=T1,T2,..: one threshold for each level of the tree "
                         "in the file (they name the levels and are not used to merge)\n";
            die = true;
        }
        if (treeThresholds.empty() && clusterGiven && clusterThreshold > 0) {
            std::cerr << PROGRAM ": --tree-file cannot be combined with --cluster\n";
            die = true;
        }
    }
    // Options that do nothing in this mode were ignored without a word.
    if (treeThresholds.empty())
        for (const std::string& o : treeOnly)
            std::cerr << PROGRAM ": warning: " << o << " has no effect without --tree\n";
    if (useExactCluster && clusterThreshold == 0)
        std::cerr << PROGRAM ": warning: --exact has no effect without --cluster\n";
    // --minhash-size is read where references are sketched: by the clustering of a
    // tree build, unless the tree comes from a file, and by --cluster without --exact.
    const bool sketchesReferences = treeThresholds.empty() ? (clusterThreshold > 0 && !useExactCluster) : !treeFileGiven;
    if (minhashGiven && !sketchesReferences) {
        if (treeFileGiven && !treeThresholds.empty())
            std::cerr << PROGRAM ": warning: --minhash-size has no effect with --tree-file (the tree is read from the "
                         "file; no reference is sketched)\n";
        else if (clusterThreshold > 0)
            std::cerr << PROGRAM ": warning: --minhash-size has no effect with --exact (the exact clustering takes "
                         "no sketch)\n";
        else
            std::cerr << PROGRAM ": warning: --minhash-size has no effect without --tree or --cluster\n";
    }

    if (argc - optind < 1) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }
    if (cfg.outfile.empty()) {
        std::cerr << PROGRAM ": missing outfile argument\n";
        die = true;
    }
    if (die) {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }
    // Said for another size only, so that the stderr of a default build is what it was.
    if (sketchesReferences && minhashSize != MINHASH_SKETCH_SIZE)
        std::cerr << PROGRAM ": MinHash sketch size " << minhashSize << " (--minhash-size; default "
                  << MINHASH_SKETCH_SIZE << "): the " << (treeThresholds.empty() ? "clusters" : "tree")
                  << " may differ from a default build's, and the index with " << (treeThresholds.empty() ? "them" : "it") << "\n";

#ifdef _OPENMP
    // Before the tree is built: the tree levels run ahead of the leaf sketch
    // and used to pick up -t only as a side effect of the cardinality estimate.
    omp_set_num_threads(cfg.threads);
#endif

    std::vector<std::string> refFiles;
    for (int i = optind; i < argc; i++) {
        std::string file(argv[i]);
        if (file[0] == '@') {
            std::string listPath = file.substr(1);
            if (!validateFile(listPath, "list file")) {
                exit(EXIT_FAILURE);
            }
            std::string inName;
            std::ifstream inList(listPath);
            while (std::getline(inList, inName)) {
                // A list saved with CRLF line endings: the carriage return became
                // part of the path, and every file was "not found" under a name that
                // looks right on the terminal.
                if (!inName.empty() && inName.back() == '\r') inName.pop_back();
                if (!inName.empty()) {
                    refFiles.push_back(inName);
                }
            }
        } else {
            refFiles.push_back(file);
        }
    }

    if (refFiles.empty()) {
        std::cerr << "Error: no reference files provided.\n";
        exit(EXIT_FAILURE);
    }

    if (!validateFiles(refFiles, "reference file")) {
        exit(EXIT_FAILURE);
    }
    {
        bool bad = false;
        for (const std::string& f : refFiles) {
            if (const char* format = unsupportedCompression(f)) {
                std::cerr << "Error: " << f << " is " << format << "-compressed; references must be plain or gzip-compressed\n";
                bad = true;
            }
        }
        if (bad) exit(EXIT_FAILURE);
    }

    // A reference is named by its file's base name, in the sketch and in the
    // candidate file it gets under -f. Two references with the same base name
    // would overwrite each other's candidates and be one name in every output.
    {
        std::unordered_map<std::string, std::string> byName;
        bool clash = false;
        for (const std::string& f : refFiles) {
            auto [it, fresh] = byName.emplace(getBaseId(f), f);
            if (!fresh) {
                std::cerr << "Error: reference files " << it->second << " and " << f
                          << " have the same name (" << it->first << ")\n";
                clash = true;
            }
        }
        if (clash) exit(EXIT_FAILURE);
    }
    // The stat file names one candidate file per reference and is read back split
    // at white space, like the sketch itself: a reference whose file name holds a
    // blank was cut there, reported as "not a uniqsketch candidate file" and left
    // out of the sketch, with exit status 0. (-f cannot hold one: option values
    // are read up to the first blank.)
    {
        bool bad = false;
        for (const std::string& f : refFiles) {
            if (getBaseId(f).find_first_of(" \t") != std::string::npos) {
                std::cerr << "Error: the name of a reference file cannot contain white space: " << f << "\n";
                bad = true;
            }
        }
        if (bad) exit(EXIT_FAILURE);
    }

    // The tree of a --tree-file build, read and checked against the references
    // before anything is written or any genome is read.
    if (treeFileGiven) {
        const std::string problem = loadTreeFile(treeFileArg, refFiles, treeThresholds.size(), tree);
        if (!problem.empty()) {
            std::cerr << "Error: --tree-file " << treeFileArg << ": " << problem
                      << "\n  --tree-file takes a tree in the format of the tree_index.tsv a --tree build writes: "
                         "its leaves the references of this build, numbered in the order of the list, "
                         "at level " << treeThresholds.size() << " (one level for each --tree threshold)\n";
            exit(EXIT_FAILURE);
        }
    }

    mkdir(cfg.outdir.c_str(), 0777);

    // Everything the build writes goes under -f, to -o and to -r, and none of them
    // was checked: with -o in a missing directory the build ran to the end, wrote
    // no sketch and returned 0; with -f or -r unwritable it left an empty sketch.
    {
        bool bad = false;
        if (!canWriteDir(cfg.outdir)) {
            std::cerr << "Error: cannot create or write the output directory (-f): " << cfg.outdir << "\n";
            bad = true;
        }
        if (!canWriteFile(cfg.outfile)) {
            std::cerr << "Error: cannot write the sketch file (-o): " << cfg.outfile << "\n";
            bad = true;
        }
        if (!canWriteFile(cfg.refstat)) {
            std::cerr << "Error: cannot write the stat file (-r): " << cfg.refstat << "\n";
            bad = true;
        }
        if (bad) exit(EXIT_FAILURE);
    }

    // Optional: build hierarchical tree index with per-level sketches
    if (!treeThresholds.empty()) {
        std::string treeFile = cfg.outdir + "/tree_index.tsv";
        if (treeFileGiven) {
            // The grouping is the file's. Which member stands for a node is decided
            // as the merges decide it, whatever the file names (treeFileRepresentatives
            // says why), so from here on the build is that of a tree it clustered
            // itself, and it writes the tree to the index as that build does.
            treeFileRepresentatives(tree, cfg.kmerLen);
            writeTreeIndex(tree, treeFile);
            std::cerr << "Hierarchical tree read from " << treeFileArg << ": " << refFiles.size() << " leaves, "
                      << tree.size() - refFiles.size() << " internal nodes, " << treeThresholds.size()
                      << " levels (no reference was clustered; the --tree thresholds name the levels; "
                         "node names and paths in the file are not read)\n";
        } else {
            tree = buildHierarchicalTree(refFiles, treeThresholds, cfg.kmerLen,
                                              cfg.nhash1, cfg.bits,
                                              estimateMaxGenomeSize(refFiles), treeFile, minhashSize);
        }
        // No two references merged at any threshold: the tree is its leaves and the
        // index will have no node level. That is what a list written with thousands
        // separators (1,000,000 is the thresholds 1, 0 and 0) or in another unit
        // builds, and it was said only in passing ("0 internal nodes").
        if (std::none_of(tree.begin(), tree.end(), [](const TreeNode& n) { return !n.isLeaf; })) {
            if (treeFileGiven)
                std::cerr << PROGRAM ": warning: the tree in " << treeFileArg << " has no node above its leaves; "
                             "the index has no node level, so no family rows\n";
            else
                std::cerr << PROGRAM ": warning: no two references merged at any --tree threshold (largest: "
                          << treeThresholds.front() << "); the index has no node level, so no family rows. "
                             "Thresholds are counts of distinguishing k-mers between genomes\n";
        }

        // Every level's sketch is built over the complete partition of the leaves
        // at that level (see treeLevelPartition), so its signatures are unique
        // against every genome present at the level. A level that partitions the
        // leaves exactly like the next finer level (or like the leaf level itself,
        // when nothing merged there) would yield the same sketch and is skipped;
        // the query walks whatever level files exist.
        // Each entry becomes a group: the representative first, then the node's
        // other members, so identifyGroupUniqKmers can make the representative's
        // signatures unique against every genome outside the node rather than
        // against the other representatives only.
        // The levels to build. Their groups are made from the tree when each level
        // is reached, because the level scan may change which member a node's
        // signatures come from (chooseNodeRepresentative), and a node that spans
        // several levels must keep one representative for all of them.
        auto partition = treeLevelPartition(tree);
        std::vector<unsigned> levelsToBuild;
        for (size_t li = 0; li < partition.size(); li++) {
            const auto& [level, entries] = partition[li];
            bool nodes = false;
            for (const auto& e : entries) nodes = nodes || e.isNode;
            if (!nodes) {
                std::cerr << "Tree level " << level << ": no clusters (same as leaf level), skipped\n";
                continue;
            }
            if (li + 1 < partition.size() && samePartition(entries, partition[li + 1].second)) {
                std::cerr << "Tree level " << level << ": same partition as level "
                          << partition[li + 1].first << ", skipped\n";
                continue;
            }
            levelsToBuild.push_back(level);
        }
        auto groupsAtLevel = [&](unsigned level, std::vector<std::string>& nodeNames, std::vector<char>& isNode) {
            std::vector<std::vector<std::string>> groups;
            nodeNames.clear(); isNode.clear();
            for (const auto& [lvl, entries] : treeLevelPartition(tree)) {
                if (lvl != level) continue;
                for (const auto& e : entries) {
                    std::vector<std::string> g{e.repPath};
                    for (const auto& p : e.leafPaths)
                        if (p != e.repPath) g.push_back(p);
                    groups.push_back(std::move(g));
                    nodeNames.push_back(e.name);
                    isNode.push_back(e.isNode);
                }
            }
            return groups;
        };

        // Record how signatures were selected so querysketch can pick matching
        // reporting defaults (see TreeBuildParams).
        TreeBuildParams bp;
        bp.kmerLen = cfg.kmerLen;
        bp.sketchnum = cfg.sketchnum;
        bp.minMargin = cfg.minMargin;
        bp.strictMargin = cfg.strictMargin;
        bp.maxHomopolymer = cfg.maxHomopolymer;
        bp.circular = cfg.circular;
        bp.groupSignatures = groupSignatures;
        if (nodeCarriers < 0) nodeCarriers = static_cast<int>(defaultNodeCarriers(refFiles.size()));
        bp.nodeCarriers = static_cast<unsigned>(nodeCarriers);
        bp.centralRepresentative = centralRepresentative;
        bp.bloomFilter = bloomFilterName();
        for (size_t i = 0; i < treeThresholds.size(); i++) {
            if (i) bp.thresholds += ",";
            bp.thresholds += std::to_string(treeThresholds[i]);
        }
        writeTreeBuildParams(treeBuildParamsPath(cfg.outdir), bp);
        if (treeFileGiven) recordTreeFile(treeBuildParamsPath(cfg.outdir), treeFileArg);
        // A line for the sketch size only when the tree was clustered from sketches of
        // another size than the default: the file of a default build is compared byte
        // for byte with those of earlier builds (cmp_index.sh).
        if (sketchesReferences && minhashSize != MINHASH_SKETCH_SIZE)
            recordMinhashSize(treeBuildParamsPath(cfg.outdir), minhashSize);

        // The main sketch is built from the leaf references (full-resolution level)
        std::unordered_set<std::string> uniqueRefs;
        for (const auto& node : tree) {
            if (node.isLeaf) uniqueRefs.insert(node.filePath);
        }
        refFiles.assign(uniqueRefs.begin(), uniqueRefs.end());
        // The leaf numbering, for leaves the sketches cannot name (see
        // treeLeavesPath). Written before any sketch so that it is there even if
        // the build stops later.
        writeTreeLeaves(treeLeavesPath(cfg.outdir), refFiles);
        std::cerr << "Building sketch from " << refFiles.size() << " tree leaf references\n";

        // Build a coarse-level sketch for each internal level. Each level gets its
        // own copy of the configuration, so nothing here can leak into the
        // leaf-level build below. The universe at every level is the full set of
        // leaves (grouped differently), so its k-mer cardinality is measured once.
        size_t universeKmers = 0, universeDistinct = 0, universeSolid = 0;
        if (!levelsToBuild.empty()) {
            getCardinality(universeKmers, universeDistinct, universeSolid,
                           cfg.kmerLen, cfg.threads, refFiles);
        }
        cfg.treeBuild = true;   // messages: a leaf without leaf signatures is not lost
        // Per-genome carriage of every k-mer, for the nodes' carried pools
        // (NodeCarrierUtil.hpp): one pass over all genomes, sized by the index's
        // distinct k-mers. Freed once the levels are built.
        std::unique_ptr<CarriageCounter> carriage;
        if (nodeCarriers > 0 && !levelsToBuild.empty()) {
#ifdef _OPENMP
            double tc = omp_get_wtime();
#endif
            carriage.reset(new CarriageCounter(universeDistinct));
            countCarriage(cfg, refFiles, *carriage);
            std::cerr << "Node carriers: cap " << nodeCarriers << ", carriage counted over "
                      << refFiles.size() << " genomes (" << carriage->bytes() / (1 << 20) << " MB)";
#ifdef _OPENMP
            std::cerr << " in " << std::setprecision(1) << std::fixed << omp_get_wtime() - tc << " s";
#endif
            std::cerr << "\n";
        }
        // The members' sampled k-mers, for choosing each node's representative
        // (chooseNodeRepresentative): one pass over all genomes, 1 in 64 k-mers,
        // freed with the counter.
        std::unique_ptr<GenomeSamples> samples;
        if (centralRepresentative && !levelsToBuild.empty()) {
#ifdef _OPENMP
            double ts = omp_get_wtime();
#endif
            samples.reset(new GenomeSamples(sampleGenomes(cfg, refFiles)));
            size_t sampled = 0;
            for (const auto& kv : samples->byPath) sampled += kv.second.size();
            std::cerr << "Node representatives: " << refFiles.size() << " genomes sampled, 1 in "
                      << static_cast<unsigned>(samples->scale()) << " k-mers (" << sampled * 8 / (1 << 20) << " MB)";
#ifdef _OPENMP
            std::cerr << " in " << std::setprecision(1) << std::fixed << omp_get_wtime() - ts << " s";
#endif
            std::cerr << "\n";
        }
        std::unordered_map<std::string, size_t> treeIndexByName;
        for (size_t i = 0; i < tree.size(); i++) treeIndexByName[tree[i].name] = i;
        std::unordered_set<std::string> representativeDecided;   // node names whose representative is settled
        size_t representativesMoved = 0;
        // Per node, the choice made: <outdir>/node_representatives.tsv.
        std::ofstream repTable;
        if (samples) {
            repTable.open(cfg.outdir + "/node_representatives.tsv");
            repTable << "level\tnode\tleaves\ttree_genome\trepresentative\tmoved\trule\treached\tcarriage\tpool_unique\tpool_carried\tfill\n";
        }
        for (unsigned level : levelsToBuild) {
            std::vector<std::string> nodeNames;
            std::vector<char> isNode;
            std::vector<std::vector<std::string>> groups = groupsAtLevel(level, nodeNames, isNode);
            if (groups.empty()) continue;
            {
                size_t clusters = 0;
                for (char c : isNode) clusters += c;
                std::cerr << "Tree level " << level << ": " << groups.size() << " references ("
                          << clusters << " clusters, " << groups.size() - clusters << " single)\n";
            }

            UniqSketchConfig lvl = cfg;
            lvl.levelLabel = "Tree level " + std::to_string(level);
            lvl.outdir = cfg.outdir + "/level_" + std::to_string(level);
            lvl.outfile = cfg.outdir + "/sketch_level_" + std::to_string(level) + ".tsv";
            lvl.refstat = lvl.outdir + "/db_uniq_count.tsv";
            lvl.numRef = groups.size();
            lvl.dbfSize = universeDistinct;
            lvl.sbfSize = universeSolid;
            mkdir(lvl.outdir.c_str(), 0777);

            if (resolveBloomSizes(universeKmers, lvl.dbfSize, lvl.sbfSize)) {
#ifdef _OPENMP
                double t0 = omp_get_wtime();
#endif
                BloomFilter levelDbf(lvl.m1(), lvl.nhash1, lvl.kmerLen);
                if (carriage || samples) {
                    std::vector<char> decide(groups.size(), 0);
                    if (samples)
                        for (size_t i = 0; i < groups.size(); i++)
                            decide[i] = isNode[i] && !representativeDecided.count(nodeNames[i]);
                    std::vector<RepresentativeChoice> choices;
                    std::vector<std::string> reps = identifyGroupUniqKmersCarried(
                        lvl, groups, levelDbf, carriage.get(), static_cast<unsigned>(nodeCarriers),
                        samples.get(), samples ? &decide : nullptr, samples ? &choices : nullptr);
                    // The chosen representatives into the tree: the node is named after the
                    // genome its signatures come from, at every level it spans.
                    for (size_t i = 0; i < groups.size(); i++) {
                        if (!isNode[i]) continue;
                        if (samples && decide[i]) {
                            const RepresentativeChoice& ch = choices[i];
                            repTable << level << "\t" << nodeNames[i] << "\t" << groups[i].size() << "\t"
                                     << getBaseId(groups[i][0]) << "\t" << getBaseId(reps[i]) << "\t"
                                     << (reps[i] != groups[i][0] ? 1 : 0) << "\t" << ch.ruleName() << "\t";
                            if (ch.decided) repTable << std::setprecision(2) << std::fixed << ch.reached << "\t" << ch.carriage
                                                     << "\t" << std::setprecision(0) << ch.poolUnique << "\t" << ch.poolCarried
                                                     << "\t" << ch.fill << "\n";
                            else repTable << "-\t-\t-\t-\t-\n";
                        }
                        representativeDecided.insert(nodeNames[i]);
                        if (reps[i] == groups[i][0]) continue;
                        auto it = treeIndexByName.find(nodeNames[i]);
                        if (it == treeIndexByName.end()) continue;
                        TreeNode& node = tree[it->second];
                        node.filePath = reps[i];
                        node.name = "cluster_L" + std::to_string(node.level) + "_" + getBaseId(reps[i]);
                        // `it` dies with the entry it points at: take the index first.
                        const size_t nodeIndex = it->second;
                        treeIndexByName.erase(it);
                        treeIndexByName[node.name] = nodeIndex;
                        representativeDecided.insert(node.name);
                        ++representativesMoved;
                    }
                } else {
                    identifyGroupUniqKmers(lvl, groups, levelDbf);
                }
#ifdef _OPENMP
                double t1 = omp_get_wtime();
#endif
                buildSketch(lvl, lvl.refstat, levelDbf);
                std::cerr << "Built sketch for tree level " << level
                          << ": " << lvl.outfile << "\n";
#ifdef _OPENMP
                std::cerr << "Tree level " << level << " build time (s): candidates "
                          << std::setprecision(1) << std::fixed << t1 - t0
                          << ", selection " << omp_get_wtime() - t1 << "\n";
#endif
            } else {
                std::cerr << "Warning: skipping sketch for tree level " << level
                          << " (no k-mer content)\n";
            }
        }
        carriage.reset();
        samples.reset();
        if (representativesMoved) {
            writeTreeIndex(tree, treeFile);
            std::cerr << "Node representatives: " << representativesMoved
                      << " node(s) renamed after the member their signatures come from; "
                      << treeFile << " rewritten\n";
        }
    }

    // Optional clustering: merge nearly identical references
    if (clusterThreshold > 0) {
        std::string clusterFile = cfg.outdir + "/clusters.tsv";
        refFiles = clusterReferences(refFiles, clusterThreshold, cfg.kmerLen,
                                     cfg.nhash1, cfg.bits, 5000000, clusterFile,
                                     useExactCluster, minhashSize);
    }

    cfg.numRef = refFiles.size();

    // Get cardinality of input reference universe data to set Bloom filter parameter
    size_t totalKmers = 0;
    getCardinality(totalKmers, cfg.dbfSize, cfg.sbfSize, cfg.kmerLen, cfg.threads, refFiles);

    if (!resolveBloomSizes(totalKmers, cfg.dbfSize, cfg.sbfSize)) {
        std::cerr << "Error: insufficient k-mer content for Bloom filter construction.\n"
                  << "  Total k-mers found: " << totalKmers << "\n"
                  << "  Distinct k-mers (F0): " << cfg.dbfSize << "\n"
                  << "  Solid k-mers: " << cfg.sbfSize << "\n"
                  << "  k-mer length: " << cfg.kmerLen << "\n"
                  << "  Number of input files: " << refFiles.size() << "\n"
                  << "Possible causes:\n"
                  << "  - Input sequences are shorter than k (" << cfg.kmerLen << ")\n"
                  << "  - Input files are empty or contain too few sequences\n"
                  << "  - Try a smaller k value with -k\n";
        exit(EXIT_FAILURE);
    }

    BloomFilter dbFilter(cfg.m1(), cfg.nhash1, cfg.kmerLen);
    identifyUniqKmers(cfg, refFiles, dbFilter);
    buildSketch(cfg, cfg.refstat, dbFilter);
    // A flat index without a single signature can report nothing (a k too small
    // for any k-mer to be unique, references that are all copies of one another):
    // the warnings above say why, and the build used to end with exit status 0.
    // A tree index may still report every reference through its nodes.
    if (tree.empty()) {
        struct stat st;
        if (stat(cfg.outfile.c_str(), &st) == 0 && st.st_size == 0) {
            std::cerr << "Error: no reference got a signature; the sketch " << cfg.outfile << " is empty\n";
            exit(EXIT_FAILURE);
        }
    }

    // Within-group leaf signatures (see GroupSketchUtil.hpp). Built after the
    // leaf sketch so that refIds match it and its signatures are not repeated,
    // with the universe filter for the margin check.
    if (!tree.empty()) {
        if (groupSignatures) {
            buildGroupSketch(cfg, tree, refFiles, dbFilter);
        } else {
            std::cerr << "Within-group signatures: skipped (--no-group-signatures)\n";
        }
        // What each leaf can be reported through, and which leaves nothing can
        // report (reportReachability). Only the levels whose sketch was built.
        std::vector<unsigned> builtLevels;
        for (const auto& [level, entries] : treeLevelPartition(tree)) {
            (void)entries;
            std::ifstream probe(cfg.outdir + "/sketch_level_" + std::to_string(level) + ".tsv");
            if (probe.good()) builtLevels.push_back(level);
        }
        reportReachability(cfg, tree, refFiles, builtLevels, nodeCarriers > 0);
    }

#ifdef _OPENMP
    std::cerr << "Runtime(sec): " << std::setprecision(4) << std::fixed
              << omp_get_wtime() - sTime << "\n";
#endif
    return 0;
}