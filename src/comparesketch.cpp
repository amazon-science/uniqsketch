#include <getopt.h>

#include <algorithm>
#include <sstream>

#include "CompareSketchUtil.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "comparesketch"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.1.0 \n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION] LIST1 LIST2\n"
    "Compare two sets of fasta reference files |LIST1|*|LIST2|.\n"
    "fasta files can be in compressed formats such as gz, bz, zip.\n"
    "Two lists of files containing file paths in each row.\n"

    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N	use N parallel threads [1]\n"
    "  -k, --kmer=N		the length of kmer [81]\n"
    "  -b, --bit=N		use N bits per element in Bloom filter [16]\n"
    "  -d, --hash=N		Bloom filter hash number [3]\n"
    "  -g, --gsize=N		approximate size for reference sequence [5000000]\n"
    "  -o, --out=STRING	the output similarity file name [reference_similarity.tsv]\n"
    "      --help		display this help and exit\n"
    "      --version		output version information and exit\n"
    "\n";

static const char shortopts[] = "t:k:d:b:g:o:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads", required_argument, NULL, 't' },
    { "kmer", required_argument, NULL, 'k' },
    { "hash", required_argument, NULL, 'd' },
    { "bits", required_argument, NULL, 'b' },
    { "gsize", required_argument, NULL, 'g' },
    { "out", required_argument, NULL, 'o'},
    { "help", no_argument, NULL, OPT_HELP },
    { "version", no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};


int main(int argc, char** argv) {
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
        case 'k':
            arg >> opt::kmerLen;
            break;
        case 'b':
            arg >> opt::bits;
            break;
        case 'd':
            arg >> opt::nhash;
            break;  
        case 'g':
            arg >> opt::dbfSize;
            break;
        case 'o':
            arg >> opt::outfile;
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
    std::ifstream inList1(argv[optind]);
    std::ifstream inList2(argv[optind+1]);

    while (getline(inList1,inName)) {
        refSet1.push_back(inName);
    }

    while (getline(inList2,inName)) {
        refSet2.push_back(inName);
    }

#ifdef _OPENMP
    omp_set_num_threads(opt::threads);
#endif

    identifyDiffernce(refSet1, refSet2);

#ifdef _OPENMP
    std::cout << "Runtime(sec): " << std::setprecision(4) << std::fixed << omp_get_wtime() - sTime << std::endl;
#endif
    return 0;
}
