#include <getopt.h>

#include <algorithm>
#include <iomanip>
#include <sstream>

#include "CompareSketchUtil.hpp"
#include "FileUtil.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "comparesketch"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.2.1\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION] LIST1 LIST2\n"
    "Compare two sets of fasta reference files |LIST1|*|LIST2|.\n"
    "Two lists of files containing file paths in each row.\n"
    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N\tuse N parallel threads [1]\n"
    "  -k, --kmer=N\t\tthe length of kmer [81]\n"
    "  -b, --bit=N\t\tuse N bits per element in Bloom filter [16]\n"
    "  -d, --hash=N\t\tBloom filter hash number [3]\n"
    "  -g, --gsize=N\t\tapproximate size for reference sequence [5000000]\n"
    "  -o, --out=STRING\tthe output similarity file name [reference_similarity.tsv]\n"
    "      --help\t\tdisplay this help and exit\n"
    "      --version\t\toutput version information and exit\n"
    "\n";

static const char shortopts[] = "t:k:d:b:g:o:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    {"threads", required_argument, nullptr, 't'},
    {"kmer",    required_argument, nullptr, 'k'},
    {"hash",    required_argument, nullptr, 'd'},
    {"bits",    required_argument, nullptr, 'b'},
    {"gsize",   required_argument, nullptr, 'g'},
    {"out",     required_argument, nullptr, 'o'},
    {"help",    no_argument,       nullptr, OPT_HELP},
    {"version", no_argument,       nullptr, OPT_VERSION},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char** argv) {
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
        case 't': arg >> opt::threads; break;
        case 'k': arg >> opt::kmerLen; break;
        case 'b': arg >> opt::bits;    break;
        case 'd': arg >> opt::nhash;   break;
        case 'g': arg >> opt::dbfSize; break;
        case 'o': arg >> opt::outfile; break;
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

    if (argc - optind < 2) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }
    if (opt::outfile.empty()) {
        std::cerr << PROGRAM ": missing outfile argument\n";
        die = true;
    }
    if (die) {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> refSet1;
    std::vector<std::string> refSet2;

    std::string inName;

    if (!validateFile(argv[optind], "list file 1")) exit(EXIT_FAILURE);
    if (!validateFile(argv[optind + 1], "list file 2")) exit(EXIT_FAILURE);

    std::ifstream inList1(argv[optind]);
    std::ifstream inList2(argv[optind + 1]);

    while (std::getline(inList1, inName)) {
        if (!inName.empty()) refSet1.push_back(inName);
    }
    while (std::getline(inList2, inName)) {
        if (!inName.empty()) refSet2.push_back(inName);
    }

    if (refSet1.empty()) {
        std::cerr << "Error: list file 1 is empty: " << argv[optind] << "\n";
        exit(EXIT_FAILURE);
    }
    if (refSet2.empty()) {
        std::cerr << "Error: list file 2 is empty: " << argv[optind + 1] << "\n";
        exit(EXIT_FAILURE);
    }

    if (!validateFiles(refSet1, "reference file (list 1)")) exit(EXIT_FAILURE);
    if (!validateFiles(refSet2, "reference file (list 2)")) exit(EXIT_FAILURE);

#ifdef _OPENMP
    omp_set_num_threads(opt::threads);
#endif

    identifyDifference(refSet1, refSet2);

#ifdef _OPENMP
    std::cout << "Runtime(sec): " << std::setprecision(4) << std::fixed
              << omp_get_wtime() - sTime << "\n";
#endif
    return 0;
}