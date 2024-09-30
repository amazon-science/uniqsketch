#include <getopt.h>
#include <sys/stat.h>

#include <algorithm>
#include <random>
#include <unordered_map>

#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"
#include "ntcard.hpp"
#include "UniqSketchUtil.hpp"

#ifdef _OPENMP
# include <omp.h>
#endif

#define PROGRAM "uniqsketch"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.1.0 \n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION] @LIST_FILES (or FILES)\n"
    "Creates unique sketch from a list of fasta reference files.\n"
    "A list of files containing file names in each row can be passed with @ prefix.\n"

    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N	use N parallel threads [1]\n"
    "  -k, --kmer=N		the length of kmer [81]\n"
    "  -b, --bit=N		use N bits per element in Bloom filter [128]\n"
    "  -d, --hash1=N		distinct Bloom filter hash number [5]\n"
    "  -s, --hash2=N		repeat Bloom filter hash number [5]\n"
    "  -c, --cov=N		number of unique k-mers to represent a reference [100]\n"
    "  -f, --outdir=STRING	dir for universe unique k-mers [outdir_uniqsketch]\n"
    "  -o, --out=STRING	the output sketch file name [sketch_uniq.tsv]\n"
    "  -r, --stat=STRING	the output unique kmer stat file name [db_uniq_count.tsv]\n"
    "  -e, --entropy		sets the aggregate entropy rate threshold [0.65]\n"
    "      --sensitive	sets sensitivity parameter c to 100\n"
    "      --very-sensitive	sets sensitivity parameter c to 1000\n"
    "      --help		display this help and exit\n"
    "      --version		output version information and exit\n"
    "\n";

static const char shortopts[] = "t:k:b:d:s:c:f:o:r:e:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads", required_argument, NULL, 't' },
    { "kmer", required_argument, NULL, 'k' },
    { "bits", required_argument, NULL, 'b' },
    { "hash1", required_argument, NULL, 'd' },
    { "hash2", required_argument, NULL, 's' },
    { "cov", required_argument, NULL, 'c' },
    { "outdir", required_argument, NULL, 'f' },
    { "out", required_argument, NULL, 'o' },
    { "stat", required_argument, NULL, 'r' },
    { "entropy", no_argument, NULL, 'e' },
    { "sensitive", no_argument, &opt::sketchnum, 100 },
    { "very-sensitive", no_argument, &opt::sketchnum, 1000 },
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
            arg >> opt::nhash1;
            break;  
        case 's':
            arg >> opt::nhash2;
            break;
        case 'c':
            arg >> opt::sketchnum;
            break;
        case 'o':
            arg >> opt::outfile;
            break;
        case 'f':
            arg >> opt::outdir;
            break;
        case 'r':
            arg >> opt::refstat;
            break;
        case 'e':
            arg >> opt::entropyThreshold;
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
            std::string inName;
            std::ifstream inList(file.substr(1,file.length()).c_str());
            while (getline(inList,inName)) {
                refFiles.push_back(inName);
            }
        } else {
            refFiles.push_back(file);
        }
    }

    mkdir(opt::outdir.c_str(), 0777);

    opt::numRef = refFiles.size();

#ifdef _OPENMP
    omp_set_num_threads(opt::threads);
#endif

    // Get cardinality of input reference universe data to set Bloom filter parameter
    size_t totalKmers = 0;
    getCardinality(totalKmers, opt::dbfSize, opt::sbfSize, opt::kmerLen, opt::threads, refFiles);

    // terminate if there is not enough distinct k-mers in universe
    if (opt::dbfSize == 0 || opt::sbfSize == 0 || opt::dbfSize > totalKmers || opt::sbfSize > totalKmers) {
        std::cerr << "Error! Bloom filter issue on too few sequences in input files.\n";
        exit(EXIT_FAILURE);
    }

    opt::m1 = opt::bits*opt::dbfSize;
    opt::m2 = opt::bits*opt::sbfSize;

    identifyUniqKmers(refFiles);

    buildSketch(opt::refstat);

#ifdef _OPENMP
    std::cerr << "Runtime(sec): " << std::setprecision(4) << std::fixed << omp_get_wtime() - sTime << std::endl;
#endif
    return 0;
}
