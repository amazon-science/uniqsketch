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
#include "FileUtil.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "uniqsketch"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.2.1\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION] @LIST_FILES (or FILES)\n"
    "Creates unique sketch from a list of fasta reference files.\n"
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
    "      --sensitive\tsets sensitivity parameter c to 100\n"
    "      --very-sensitive\tsets sensitivity parameter c to 1000\n"
    "      --help\t\tdisplay this help and exit\n"
    "      --version\t\toutput version information and exit\n"
    "\n";

static const char shortopts[] = "t:k:b:d:s:c:f:o:r:e:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_CLUSTER };

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
    {"sensitive",      no_argument, &opt::sketchnum, 100},
    {"very-sensitive", no_argument, &opt::sketchnum, 1000},
    {"help",           no_argument,       nullptr, OPT_HELP},
    {"version",        no_argument,       nullptr, OPT_VERSION},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char** argv) {
#ifdef _OPENMP
    double sTime = omp_get_wtime();
#endif

    unsigned clusterThreshold = 0;  // 0 = clustering disabled

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) {
        std::istringstream arg(optarg != nullptr ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 't': arg >> opt::threads;          break;
        case 'k': arg >> opt::kmerLen;          break;
        case 'b': arg >> opt::bits;             break;
        case 'd': arg >> opt::nhash1;           break;
        case 's': arg >> opt::nhash2;           break;
        case 'c': arg >> opt::sketchnum;        break;
        case 'o': arg >> opt::outfile;          break;
        case 'f': arg >> opt::outdir;           break;
        case 'r': arg >> opt::refstat;          break;
        case 'e': arg >> opt::entropyThreshold; break;
        case OPT_CLUSTER: arg >> clusterThreshold; break;
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

    if (argc - optind < 1) {
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

    mkdir(opt::outdir.c_str(), 0777);

    // Optional clustering: merge nearly identical references
    if (clusterThreshold > 0) {
        std::string clusterFile = opt::outdir + "/clusters.tsv";
        refFiles = clusterReferences(refFiles, clusterThreshold, opt::kmerLen,
                                     opt::nhash1, opt::bits, 5000000, clusterFile);
    }

    opt::numRef = refFiles.size();

#ifdef _OPENMP
    omp_set_num_threads(opt::threads);
#endif

    // Get cardinality of input reference universe data to set Bloom filter parameter
    size_t totalKmers = 0;
    getCardinality(totalKmers, opt::dbfSize, opt::sbfSize, opt::kmerLen, opt::threads, refFiles);

    if (opt::dbfSize == 0 || opt::sbfSize == 0 || opt::dbfSize > totalKmers || opt::sbfSize > totalKmers) {
        std::cerr << "Error: insufficient k-mer content for Bloom filter construction.\n"
                  << "  Total k-mers found: " << totalKmers << "\n"
                  << "  Distinct k-mers (F0): " << opt::dbfSize << "\n"
                  << "  Solid k-mers: " << opt::sbfSize << "\n"
                  << "  k-mer length: " << opt::kmerLen << "\n"
                  << "  Number of input files: " << refFiles.size() << "\n"
                  << "Possible causes:\n"
                  << "  - Input sequences are shorter than k (" << opt::kmerLen << ")\n"
                  << "  - Input files are empty or contain too few sequences\n"
                  << "  - Try a smaller k value with -k\n";
        exit(EXIT_FAILURE);
    }

    opt::m1 = opt::bits * opt::dbfSize;
    opt::m2 = opt::bits * opt::sbfSize;

    identifyUniqKmers(refFiles);
    buildSketch(opt::refstat);

#ifdef _OPENMP
    std::cerr << "Runtime(sec): " << std::setprecision(4) << std::fixed
              << omp_get_wtime() - sTime << "\n";
#endif
    return 0;
}