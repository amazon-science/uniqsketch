#include <getopt.h>

#include <algorithm>
#include <numeric>
#include <unordered_map>

#include "ntHashIterator.hpp"
#include "ntcard.hpp"
#include "QuerySketchUtil.hpp"
#include "TreeSketchUtil.hpp"
#include "FileUtil.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "querysketch"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.5.0\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTIONS] [ARGS]\n"
    "Identify references and their abundance in FILE(S).\n\n"
    "Acceptable file formats: fastq in compressed formats gz, bz, zip, xz.\n"
    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N\tuse N parallel threads [1]\n"
    "  -b, --bit=N\t\tuse N bits per element in Bloom filter [64]\n"
    "  -o, --out=STRING\tthe output file name\n"
    "  -l, --r1=STRING\tinput read 1\n"
    "  -r, --r2=STRING\tinput read 2\n"
    "  -g, --ref=STRING\tinput uniqsketch reference\n"
    "  -h, --hit=N\t\tnumber of uniqsketch hits to call a reference [10]\n"
    "  -a, --acutoff=N\tabundance cutoff to report [0.0]\n"
    "  -s, --rcutoff=N\tread cutoff to report [2]\n"
    "      --sensitive\tsets sensitivity parameter h=10\n"
    "      --very-sensitive\tsets sensitivity parameter h=5\n"
    "      --solid\t\tonly use solid k-mers in reads\n"
    "      --tree=DIR\thierarchical query using a --tree index directory\n"
    "      --help\t\tdisplay this help and exit\n"
    "      --version\t\tversion information and exit\n"
    "\n";

static const char shortopts[] = "t:g:l:r:o:h:a:s:b:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_TREE };

static const struct option longopts[] = {
    {"threads",        required_argument, nullptr, 't'},
    {"out",            required_argument, nullptr, 'o'},
    {"r1",             required_argument, nullptr, 'l'},
    {"r2",             required_argument, nullptr, 'r'},
    {"ref",            required_argument, nullptr, 'g'},
    {"hit",            required_argument, nullptr, 'h'},
    {"acutoff",        required_argument, nullptr, 'a'},
    {"rcutoff",        required_argument, nullptr, 's'},
    {"bit",            required_argument, nullptr, 'b'},
    {"sensitive",      no_argument, &opt::sketchHit, 10},
    {"very-sensitive", no_argument, &opt::sketchHit, 5},
    {"solid",          no_argument, &opt::solid, 1},
    {"tree",           required_argument, nullptr, OPT_TREE},
    {"help",           no_argument, nullptr, OPT_HELP},
    {"version",        no_argument, nullptr, OPT_VERSION},
    {nullptr, 0, nullptr, 0}
};

/**
 * Hierarchical (tree) query: narrow candidate references from coarse to fine.
 *
 * The coarsest level sketch is queried first; cluster representatives that clear the
 * hit/read thresholds mark their descendant leaves as candidates. Each finer level is
 * then queried with its sketch restricted to signatures of references that are still
 * candidates, so the search narrows down the tree. The final (leaf) level sketch
 * (opt::ref) is queried restricted to the surviving candidate leaves and its result is
 * written. This applies the candidate narrowing that the flat port computed but never
 * used, so the tree actually focuses the query on promising sub-trees.
 */
static void runTreeQuery(const std::string& treeDir,
                         const std::vector<std::string>& sampleFiles) {
    std::string treeFile = treeDir + "/tree_index.tsv";
    if (!validateFile(treeFile, "tree index file")) exit(EXIT_FAILURE);
    if (!validateFile(opt::ref, "leaf sketch (--ref)")) exit(EXIT_FAILURE);

    std::vector<TreeNode> tree = loadTree(treeFile);

    // Descendant leaf-name set for every node. Children always have smaller ids than
    // their parent (nodes are created bottom-up), so a single ascending pass suffices.
    std::vector<std::unordered_set<std::string>> descLeaves(tree.size());
    unsigned maxInternalLevel = 0;
    bool anyInternal = false;
    for (size_t idx = 0; idx < tree.size(); idx++) {
        const TreeNode& node = tree[idx];
        if (node.isLeaf) {
            descLeaves[idx].insert(getBaseId(node.filePath));
        } else {
            for (unsigned c : node.children)
                if (c < tree.size())
                    descLeaves[idx].insert(descLeaves[c].begin(), descLeaves[c].end());
            maxInternalLevel = std::max(maxInternalLevel, node.level);
            anyInternal = true;
        }
    }

    // level -> (representative basename -> covered leaf names). Mirrors how uniqsketch
    // emits level sketches: internal-node representatives at their level, plus any
    // never-clustered leaves placed at the deepest internal level.
    std::unordered_map<unsigned,
        std::unordered_map<std::string, std::unordered_set<std::string>>> levelMap;
    for (size_t idx = 0; idx < tree.size(); idx++) {
        const TreeNode& node = tree[idx];
        if (!node.isLeaf)
            levelMap[node.level][getBaseId(node.filePath)] = descLeaves[idx];
    }
    if (anyInternal) {
        for (size_t idx = 0; idx < tree.size(); idx++) {
            const TreeNode& node = tree[idx];
            if (node.isLeaf && node.parentId == -1)
                levelMap[maxInternalLevel][getBaseId(node.filePath)] =
                    std::unordered_set<std::string>{ getBaseId(node.filePath) };
        }
    }

    // Level sketches coarsest-first, then the leaf/full-resolution sketch (opt::ref).
    std::vector<std::pair<unsigned, std::string>> levelSketches;
    for (unsigned lvl = 0; lvl < 100; lvl++) {
        std::string p = treeDir + "/sketch_level_" + std::to_string(lvl) + ".tsv";
        if (fileExists(p)) levelSketches.emplace_back(lvl, p);
    }
    unsigned leafLevel = 0;
    for (const TreeNode& node : tree)
        if (node.isLeaf) leafLevel = std::max(leafLevel, node.level);
    levelSketches.emplace_back(leafLevel, opt::ref);
    std::sort(levelSketches.begin(), levelSketches.end());

    // Prime opt::kmerLen from a sketch, then size read Bloom filters once from the sample.
    {
        SketchHash h; std::vector<std::string> r; std::vector<SketchHash> sc;
        loadSketch(levelSketches.front().second, h, r, sc);
    }
    size_t totalKmers = 0;
    getCardinality(totalKmers, opt::dbfSize, opt::sbfSize, opt::kmerLen, opt::threads, sampleFiles);
    if (!resolveBloomSizes(totalKmers, opt::dbfSize, opt::sbfSize)) {
        std::cerr << "Error: insufficient k-mer content in the sample for query.\n";
        exit(EXIT_FAILURE);
    }

    std::unordered_set<std::string> candidates;  // surviving candidate leaf names
    bool firstLevel = true;

    for (size_t li = 0; li < levelSketches.size(); li++) {
        unsigned level = levelSketches[li].first;
        const std::string& sketchPath = levelSketches[li].second;
        bool isLeaf = (li + 1 == levelSketches.size());

        SketchHash levelHash;
        std::vector<std::string> levelRef;
        std::vector<SketchHash> levelSig;
        loadSketch(sketchPath, levelHash, levelRef, levelSig);

        // Restrict to candidate references (all are in play on the coarsest level).
        if (!firstLevel) {
            std::unordered_set<unsigned> keep;
            const auto& lm = levelMap[level];
            for (unsigned i = 0; i < levelRef.size(); i++) {
                bool inPlay = false;
                if (isLeaf) {
                    inPlay = candidates.count(levelRef[i]) > 0;
                } else {
                    auto it = lm.find(levelRef[i]);
                    if (it != lm.end())
                        for (const auto& leaf : it->second)
                            if (candidates.count(leaf)) { inPlay = true; break; }
                }
                if (inPlay) keep.insert(i);
            }
            SketchHash filtered;
            for (const auto& kv : levelHash)
                if (keep.count(kv.second)) filtered.emplace(kv.first, kv.second);
            levelHash.swap(filtered);
        }

        std::vector<unsigned> levelCount(levelRef.size(), 0);
        std::vector<std::unordered_set<std::string>> levelRead(levelRef.size());
        querySampleBatch(sampleFiles, levelHash, levelCount, levelSig, levelRead);

        if (isLeaf) {
            generateQueryResult(levelCount, levelRef, levelSig, levelRead, opt::out);
            std::cerr << "Tree query: leaf level " << level << " ("
                      << levelRef.size() << " leaves) -> " << opt::out << "\n";
            break;
        }

        // Narrow: descendants of representatives that cleared the thresholds.
        std::unordered_set<std::string> nextCandidates;
        const auto& lm = levelMap[level];
        for (unsigned i = 0; i < levelRef.size(); i++) {
            if (levelCount[i] >= static_cast<unsigned>(opt::sketchHit) &&
                levelRead[i].size() > static_cast<size_t>(opt::readCutoff)) {
                auto it = lm.find(levelRef[i]);
                if (it != lm.end())
                    nextCandidates.insert(it->second.begin(), it->second.end());
            }
        }
        candidates.swap(nextCandidates);
        firstLevel = false;
        std::cerr << "Tree query: level " << level << " -> "
                  << candidates.size() << " candidate leaves\n";
    }
}

int main(int argc, char* const argv[]) {
#ifdef _OPENMP
    double sTime = omp_get_wtime();
#endif

    std::string treeDir;  // --tree: hierarchical query index directory (empty = flat)

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) {
        std::istringstream arg(optarg != nullptr ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 't': arg >> opt::threads;         break;
        case 'o': arg >> opt::out;             break;
        case 'l': arg >> opt::r1;              break;
        case 'r': arg >> opt::r2;              break;
        case 'g': arg >> opt::ref;             break;
        case 'h': arg >> opt::sketchHit;       break;
        case 'a': arg >> opt::abundanceCutoff; break;
        case 's': arg >> opt::readCutoff;      break;
        case 'b': arg >> opt::bits;            break;
        case OPT_TREE: arg >> treeDir;         break;
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

    if (opt::r1.empty()) {
        std::cerr << "Error: no input file r1 provided.\n";
        die = true;
    }
    if (opt::ref.empty()) {
        std::cerr << "Error: no sketch file provided.\n";
        die = true;
    }
    if (die) {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

    // Validate input files exist
    if (!validateFile(opt::r1, "read file r1")) die = true;
    if (!opt::r2.empty() && !validateFile(opt::r2, "read file r2")) die = true;
    if (!validateFile(opt::ref, "sketch reference file")) die = true;
    if (die) {
        exit(EXIT_FAILURE);
    }

#ifdef _OPENMP
    omp_set_num_threads(opt::threads);
#endif

    std::vector<std::string> sampleFiles;
    sampleFiles.push_back(opt::r1);
    if (!opt::r2.empty()) {
        sampleFiles.push_back(opt::r2);
    }

    if (!treeDir.empty()) {
        // Hierarchical query: narrow candidates coarse -> fine down the tree.
        runTreeQuery(treeDir, sampleFiles);
    } else {
        // Flat query against a single sketch.
        SketchHash sketchHash;
        std::vector<std::string> sketchRef;
        std::vector<SketchHash> refSigCount;

        loadSketch(opt::ref, sketchHash, sketchRef, refSigCount);

        std::cout << "k: " << opt::kmerLen << "\n";
        std::cout << "Number of references: " << sketchRef.size() << "\n";
        std::cout << "Size of sketch: " << sketchHash.size() << "\n";

        // Get cardinality of input raw read data to set Bloom filter parameter
        size_t totalKmers = 0;
        getCardinality(totalKmers, opt::dbfSize, opt::sbfSize, opt::kmerLen, opt::threads, sampleFiles);

        if (opt::dbfSize == 0 || opt::sbfSize == 0 || opt::dbfSize > totalKmers || opt::sbfSize > totalKmers) {
            std::cerr << "Error: insufficient k-mer content for Bloom filter construction.\n"
                      << "  Total k-mers found: " << totalKmers << "\n"
                      << "  Distinct k-mers (F0): " << opt::dbfSize << "\n"
                      << "  Solid k-mers: " << opt::sbfSize << "\n"
                      << "  k-mer length: " << opt::kmerLen << "\n"
                      << "Possible causes:\n"
                      << "  - Input reads are shorter than k (" << opt::kmerLen << ")\n"
                      << "  - Input files are empty or contain too few reads\n"
                      << "  - Try a smaller k value in the sketch\n";
            exit(EXIT_FAILURE);
        }

        std::cout << "Number of weak k-mers: " << opt::dbfSize - opt::sbfSize << "\n";
        std::cout << "Number of solid k-mers: " << opt::sbfSize << "\n";

        std::vector<unsigned> sketchCount(sketchRef.size());
        std::vector<std::unordered_set<std::string>> refRead(sketchRef.size());

        querySampleBatch(sampleFiles, sketchHash, sketchCount, refSigCount, refRead);
        generateQueryResult(sketchCount, sketchRef, refSigCount, refRead, opt::out);
    }

#ifdef _OPENMP
    std::cout << "Runtime(sec): " << std::setprecision(4) << std::fixed
              << omp_get_wtime() - sTime << "\n";
#endif
    return 0;
}