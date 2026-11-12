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
    "  -e, --entropy\t\tsets the aggregate entropy rate threshold [0.65]\n"
    "      --cluster=N\tcluster similar references with N unique k-mer threshold [0=off]\n"
    "      --tree=T1,T2,..\tbuild hierarchical tree index with comma-separated thresholds\n"
    "      --exact\t\tuse exact O(n^2) clustering instead of MinHash (with --cluster)\n"
    "      --min-margin=N\tmin Hamming distance of signatures to other refs (1=off, 2) [1; 2 with --tree]\n"
    "      --max-homopolymer=N\tdrop signatures with a homopolymer run > N bp (0=off) [0]\n"
    "      --strict-margin\tenforce --min-margin as a hard requirement (no margin-1 fallback; may yield fewer signatures) [on with --tree]\n"
    "      --circular\treferences are circular: also take the k-1 k-mers spanning each contig's end-to-start join\n"
    "      --no-group-signatures\tdo not build within-group leaf signatures (with --tree; see README)\n"
    "      --node-carriers=N\tfill a tree node's signature quota, when its node-unique k-mers leave it short,\n"
    "\t\t\tfrom k-mers carried by up to N genomes outside the node (recorded; the query\n"
    "\t\t\tsets them aside when a carrier is present) [min(50, leaves/10); 0=off]\n"
    "      --sensitive\tsets sensitivity parameter c to 100\n"
    "      --very-sensitive\tsets sensitivity parameter c to 1000\n"
    "      --help\t\tdisplay this help and exit\n"
    "      --version\t\toutput version information and exit\n"
    "\n";

static const char shortopts[] = "t:k:b:d:s:c:f:o:r:e:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_CLUSTER, OPT_MINMARGIN, OPT_EXACT, OPT_TREE,
       OPT_MAXHOMOP, OPT_STRICTMARGIN, OPT_NOGROUPSIGS , OPT_CIRCULAR, OPT_NODECARRIERS };

// --sensitive / --very-sensitive are getopt flag entries, which store through a
// pointer at parse time and so cannot target a struct member. They write here and
// the value is copied into the config once parsing is complete.
static int g_sketchnum = 100;

static const struct option longopts[] = {
    {"threads",        required_argument, nullptr, 't'},
    {"kmer",           required_argument, nullptr, 'k'},
    {"bits",           required_argument, nullptr, 'b'},
    {"hash1",          required_argument, nullptr, 'd'},
    {"hash2",          required_argument, nullptr, 's'},
    {"cov",            required_argument, nullptr, 'c'},
    {"outdir",         required_argument, nullptr, 'f'},
    {"out",            required_argument, nullptr, 'o'},
    {"stat",           required_argument, nullptr, 'r'},
    {"entropy",        no_argument,       nullptr, 'e'},
    {"cluster",        required_argument, nullptr, OPT_CLUSTER},
    {"tree",           required_argument, nullptr, OPT_TREE},
    {"min-margin",     required_argument, nullptr, OPT_MINMARGIN},
    {"exact",          no_argument,       nullptr, OPT_EXACT},
    {"max-homopolymer",required_argument, nullptr, OPT_MAXHOMOP},
    {"strict-margin",  no_argument,       nullptr, OPT_STRICTMARGIN},
    {"circular",       no_argument,       nullptr, OPT_CIRCULAR},
    {"no-group-signatures", no_argument,  nullptr, OPT_NOGROUPSIGS},
    {"node-carriers",  required_argument, nullptr, OPT_NODECARRIERS},
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
    bool marginGiven = false;       // --min-margin passed explicitly (tree default otherwise)
    bool groupSignatures = true;    // within-group leaf signatures for a tree index (--no-group-signatures)
    int nodeCarriers = -1;          // --node-carriers (-1: min(50, leaves/10) once the leaf count is known)

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) {
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
        case 'c': arg >> g_sketchnum;          break;
        case 'o': arg >> cfg.outfile;          break;
        case 'f': arg >> cfg.outdir;           break;
        case 'r': arg >> cfg.refstat;          break;
        case 'e': arg >> cfg.entropyThreshold; break;
        case OPT_CLUSTER: arg >> clusterThreshold; break;
        case OPT_MINMARGIN: arg >> cfg.minMargin; marginGiven = true; break;
        case OPT_EXACT: useExactCluster = true; break;
        case OPT_TREE: {
            std::string treeArg;
            arg >> treeArg;
            std::istringstream ts(treeArg);
            std::string token;
            while (std::getline(ts, token, ',')) {
                if (!token.empty()) treeThresholds.push_back(std::stoi(token));
            }
            // Sort descending (coarsest first)
            std::sort(treeThresholds.begin(), treeThresholds.end(), std::greater<unsigned>());
            break;
        }
        case OPT_CIRCULAR: cfg.circular = true; break;
        case OPT_MAXHOMOP: arg >> cfg.maxHomopolymer; break;
        case OPT_STRICTMARGIN: cfg.strictMargin = true; break;
        case OPT_NOGROUPSIGS: groupSignatures = false; break;
        case OPT_NODECARRIERS: arg >> nodeCarriers; break;
        case OPT_HELP:
            std::cerr << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cerr << VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        }
        if (optarg != nullptr && !arg.eof()) {
            std::cerr << PROGRAM ": invalid option: `-"
                      << static_cast<char>(c) << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }

    cfg.sketchnum = g_sketchnum;

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
    // the very signatures being excluded. An explicit --min-margin wins.
    if (!treeThresholds.empty() && !marginGiven) {
        cfg.minMargin = 2;
        cfg.strictMargin = true;
        std::cerr << PROGRAM ": --tree defaults to --min-margin=2 --strict-margin "
                     "(pass --min-margin=1 to select signatures without a margin)\n";
    }

    if (cfg.minMargin < 1) cfg.minMargin = 1;
    if (cfg.minMargin > 2) {
        std::cerr << PROGRAM ": --min-margin > 2 is not supported; using 2\n";
        cfg.minMargin = 2;
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

    mkdir(cfg.outdir.c_str(), 0777);

    // Optional: build hierarchical tree index with per-level sketches
    if (!treeThresholds.empty()) {
        std::string treeFile = cfg.outdir + "/tree_index.tsv";
        tree = buildHierarchicalTree(refFiles, treeThresholds, cfg.kmerLen,
                                          cfg.nhash1, cfg.bits,
                                          estimateMaxGenomeSize(refFiles), treeFile);

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
        auto partition = treeLevelPartition(tree);
        std::vector<std::pair<unsigned, std::vector<std::vector<std::string>>>> levelGroups;
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
            std::vector<std::vector<std::string>> groups;
            size_t clusters = 0;
            for (const auto& e : entries) {
                std::vector<std::string> g{e.repPath};
                for (const auto& p : e.leafPaths)
                    if (p != e.repPath) g.push_back(p);
                groups.push_back(std::move(g));
                clusters += e.isNode;
            }
            std::cerr << "Tree level " << level << ": " << groups.size() << " references ("
                      << clusters << " clusters, " << groups.size() - clusters << " single)\n";
            levelGroups.emplace_back(level, std::move(groups));
        }

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
        bp.bloomFilter = bloomFilterName();
        for (size_t i = 0; i < treeThresholds.size(); i++) {
            if (i) bp.thresholds += ",";
            bp.thresholds += std::to_string(treeThresholds[i]);
        }
        writeTreeBuildParams(treeBuildParamsPath(cfg.outdir), bp);

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
        if (!levelGroups.empty()) {
            getCardinality(universeKmers, universeDistinct, universeSolid,
                           cfg.kmerLen, cfg.threads, refFiles);
        }
        cfg.treeBuild = true;   // messages: a leaf without leaf signatures is not lost
        // Per-genome carriage of every k-mer, for the nodes' carried pools
        // (NodeCarrierUtil.hpp): one pass over all genomes, sized by the index's
        // distinct k-mers. Freed once the levels are built.
        std::unique_ptr<CarriageCounter> carriage;
        if (nodeCarriers > 0 && !levelGroups.empty()) {
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
        for (const auto& [level, groups] : levelGroups) {
            if (groups.empty()) continue;

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
                if (carriage) identifyGroupUniqKmersCarried(lvl, groups, levelDbf, *carriage,
                                                            static_cast<unsigned>(nodeCarriers));
                else identifyGroupUniqKmers(lvl, groups, levelDbf);
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
    }

    // Optional clustering: merge nearly identical references
    if (clusterThreshold > 0) {
        std::string clusterFile = cfg.outdir + "/clusters.tsv";
        refFiles = clusterReferences(refFiles, clusterThreshold, cfg.kmerLen,
                                     cfg.nhash1, cfg.bits, 5000000, clusterFile,
                                     useExactCluster);
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