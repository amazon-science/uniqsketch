#include <getopt.h>

#include <algorithm>
#include <numeric>
#include <unordered_map>

#include "ntHashIterator.hpp"
#include "ntcard.hpp"
#include "QuerySketchUtil.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "querysketch"

static const char VERSION_MESSAGE[] =
    PROGRAM "1.1.0 \n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTIONS] [ARGS]\n"
    "Identify references and their abundance in FILE(S).\n\n"
    "Acceptable file formats: fastq in compressed formats gz, bz, zip, xz.\n"

    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N	use N parallel threads [1]\n"
    "  -b, --bit=N		use N bits per element in Bloom filter [64]\n"
    "  -o, --out=STRING	the output file name\n"
    "  -l, --r1=STRING	input read 1\n"
    "  -r, --r2=STRING	input read 2\n"
    "  -g, --ref=STRING	input uniqsketch reference\n"
    "  -h, --hit=N		number of uniqsketch hits to call a reference [10]\n"
    "  -a, --acutoff=N	abundance cutoff to report [0.0]\n"
    "  -s, --rcutoff=N	read cutoff to report [2]\n"
    "      --sensitive	sets sensitivity parameter h=10\n"
    "      --very-sensitive	sets sensitivity parameter h=5\n"
    "      --solid		only use solid k-mers in reads\n"
    "      --help		display this help and exit\n"
    "      --versions	version information and exit\n"
    "\n";


static const char shortopts[] = "t:g:l:r:o:h:a:s:b:";

enum {
    OPT_HELP = 1,
    OPT_VERSION
};

static const struct option longopts[] = {
    {"threads", required_argument, NULL, 't'},
    {"out", required_argument, NULL, 'o'},
    {"r1", required_argument, NULL, 'l'},
    {"r2", required_argument, NULL, 'r'},
    {"ref", required_argument, NULL, 'g'},
    {"hit", required_argument, NULL, 'h'},
    {"acutoff", required_argument, NULL, 'a'},
    {"rcutoff", required_argument, NULL, 's'},
    {"bit", required_argument, NULL, 'b'},
    {"sensitive", no_argument, &opt::sketchHit, 10},
    {"very-sensitive", no_argument, &opt::sketchHit, 5},
    {"solid", no_argument, &opt::solid, 1},
    {"help", no_argument, NULL, OPT_HELP},
    {"version", no_argument, NULL, OPT_VERSION},
    {NULL, 0, NULL, 0}
};

int main(int argc, char *const argv[]) {
#ifdef _OPENMP
    double sTime = omp_get_wtime();
#endif

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 't':
            arg >> opt::threads;
            break;
        case 'o':
            arg >> opt::out;
            break;
        case 'l':
            arg >> opt::r1;
            break;
        case 'r':
            arg >> opt::r2;
            break;
        case 'g':
            arg >> opt::ref;
            break;
        case 'h':
            arg >> opt::sketchHit;
            break;
        case 'a':
            arg >> opt::abundanceCutoff;
            break;
        case 's':
            arg >> opt::readCutoff;
            break;
        case 'b':
            arg >> opt::bits;
            break;
        case OPT_HELP:
            std::cerr << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cerr << VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        }
        if (optarg != NULL && !arg.eof()) {
            std::cerr << PROGRAM ": invalid option: `-"
                      << (char)c << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }
    if (opt::r1.empty()) {
        std::cerr << "Error! No input file r1 provided.\n";
        die = true;
    }
    if (opt::ref.empty()) {
        std::cerr << "Error! No sketch file provided.\n";
        die = true;
    }
    if (die) {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

#ifdef _OPENMP
    omp_set_num_threads(opt::threads);
#endif

    SketchHash sketchHash;
    // vector for mapping integer reference id to real reference name
    std::vector<std::string> sketchRef;
    // vector for keeping the count of each uniq k-mer for each reference
    std::vector<SketchHash> refSigCount;

    // load uniqsketch into hash 
    loadSketch(opt::ref, sketchHash, sketchRef, refSigCount);
    
    std::cout << "k: " << opt::kmerLen << std::endl;
    std::cout << "Number of references: " << sketchRef.size() << std::endl;
    std::cout << "Size of sketch:" << sketchHash.size() << std::endl;

    std::vector<std::string> sampleFiles;
    sampleFiles.push_back(opt::r1);
    if (!opt::r2.empty()) {
        sampleFiles.push_back(opt::r2);
    } 

    // Get cardinality of input raw read data to set Bloom filter parameter
    size_t totalKmers = 0;
    getCardinality(totalKmers, opt::dbfSize, opt::sbfSize, opt::kmerLen, opt::threads, sampleFiles);

    // terminate if there is not enough distinct or solid k-mers in input data
    if (opt::dbfSize == 0 || opt::sbfSize == 0 || opt::dbfSize > totalKmers || opt::sbfSize > totalKmers) {
        std::cerr << "Error! Bloom filter issue on too few sequences in input files.\n";
        exit(EXIT_FAILURE);
    }

    std::cout << "Number of weak k-mers: " << opt::dbfSize - opt::sbfSize << std::endl;
    std::cout << "Number of solid k-mers: " << opt::sbfSize << std::endl;

    // vector for keeping the count of uniq k-mers for each reference
    std::vector<unsigned> sketchCount(sketchRef.size());

    // vector for keeping the track of distinct reads in each reference
    std::vector<std::unordered_set<std::string> > refRead(sketchRef.size());

    // query all input read files
    querySampleBatch(sampleFiles, sketchHash, sketchCount, refSigCount, refRead);

    // generate tsv output of top reference-abundance-count for input data
    generateQueryResult(sketchCount, sketchRef, refSigCount, refRead, opt::out);

#ifdef _OPENMP
    std::cout << "Runtime(sec): " << std::setprecision(4) << std::fixed << omp_get_wtime() - sTime << std::endl;
#endif
    return 0;
}
