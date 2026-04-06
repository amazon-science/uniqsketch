#include <getopt.h>

#include <algorithm>
#include <numeric>
#include <unordered_map>

#include "ntHashIterator.hpp"
#include "ntcard.hpp"
#include "QuerySketchUtil.hpp"
#include "FileUtil.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "querysketch"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.2.1\n";

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
    "      --help\t\tdisplay this help and exit\n"
    "      --version\t\tversion information and exit\n"
    "\n";

static const char shortopts[] = "t:g:l:r:o:h:a:s:b:";

enum { OPT_HELP = 1, OPT_VERSION };

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
    {"help",           no_argument, nullptr, OPT_HELP},
    {"version",        no_argument, nullptr, OPT_VERSION},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char* const argv[]) {
#ifdef _OPENMP
    double sTime = omp_get_wtime();
#endif

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

    SketchHash sketchHash;
    std::vector<std::string> sketchRef;
    std::vector<SketchHash> refSigCount;

    loadSketch(opt::ref, sketchHash, sketchRef, refSigCount);

    std::cout << "k: " << opt::kmerLen << "\n";
    std::cout << "Number of references: " << sketchRef.size() << "\n";
    std::cout << "Size of sketch: " << sketchHash.size() << "\n";

    std::vector<std::string> sampleFiles;
    sampleFiles.push_back(opt::r1);
    if (!opt::r2.empty()) {
        sampleFiles.push_back(opt::r2);
    }

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

#ifdef _OPENMP
    std::cout << "Runtime(sec): " << std::setprecision(4) << std::fixed
              << omp_get_wtime() - sTime << "\n";
#endif
    return 0;
}