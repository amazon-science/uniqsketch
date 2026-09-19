#include <getopt.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "CompareSketchUtil.hpp"
#include "ntcard.hpp"
#include "FileUtil.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "comparesketch"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 2.0.0\n";

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
    "  -m, --mem=N\t\tmemory budget in GB for reference caching (0 = stream) [64]\n"
    "      --auto-gsize\tsize the Bloom filter from the largest input genome\n"
    "      --auto\t\tsize the Bloom filter from an ntCard cardinality estimate\n"
    "      --fpr=F\t\ttarget Bloom-filter false-positive rate (sets bits; e.g. 0.001)\n"
    "      --low-mem\t\talias for --mem 0 (stream: minimal memory, more I/O)\n"
    "      --help\t\tdisplay this help and exit\n"
    "      --version\t\toutput version information and exit\n"
    "\n";

static const char shortopts[] = "t:k:d:b:g:o:m:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_AUTO, OPT_AUTOGSIZE, OPT_FPR, OPT_LOWMEM };

static const struct option longopts[] = {
    {"threads", required_argument, nullptr, 't'},
    {"kmer",    required_argument, nullptr, 'k'},
    {"hash",    required_argument, nullptr, 'd'},
    {"bits",    required_argument, nullptr, 'b'},
    {"gsize",   required_argument, nullptr, 'g'},
    {"out",     required_argument, nullptr, 'o'},
    {"mem",     required_argument, nullptr, 'm'},
    {"auto-gsize", no_argument,    nullptr, OPT_AUTOGSIZE},
    {"auto",    no_argument,       nullptr, OPT_AUTO},
    {"fpr",     required_argument, nullptr, OPT_FPR},
    {"low-mem", no_argument,       nullptr, OPT_LOWMEM},
    {"help",    no_argument,       nullptr, OPT_HELP},
    {"version", no_argument,       nullptr, OPT_VERSION},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char** argv) {
#ifdef _OPENMP
    double sTime = omp_get_wtime();
#endif

    CompareSketchConfig cfg;
    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) {
        std::istringstream arg(optarg != nullptr ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 't': arg >> cfg.threads; break;
        case 'k': arg >> cfg.kmerLen; break;
        case 'b': arg >> cfg.bits;    break;
        case 'd': arg >> cfg.nhash;   break;
        case 'g': arg >> cfg.dbfSize; break;
        case 'o': arg >> cfg.outfile; break;
        case 'm': arg >> cfg.memBudgetGB; break;
        case OPT_AUTOGSIZE: cfg.autoGsize = true; break;
        case OPT_AUTO: cfg.autoSize = true; break;
        case OPT_FPR: arg >> cfg.targetFpr; break;
        case OPT_LOWMEM: cfg.memBudgetGB = 0; break;
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
    if (cfg.outfile.empty()) {
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
    omp_set_num_threads(cfg.threads);
#endif

    // Optional: derive bits-per-element from a target false-positive rate.
    // FPR p = (1 - e^(-h/b))^h  =>  b = -h / ln(1 - p^(1/h)).
    if (cfg.targetFpr > 0.0 && cfg.targetFpr < 1.0) {
        double h = static_cast<double>(cfg.nhash);
        double b = -h / std::log(1.0 - std::pow(cfg.targetFpr, 1.0 / h));
        cfg.bits = static_cast<unsigned>(std::ceil(b));
        std::cerr << PROGRAM ": target FPR " << cfg.targetFpr
                  << " -> " << cfg.bits << " bits/element\n";
    }

    // Optional Bloom-filter sizing. --auto uses an ntCard cardinality estimate of
    // the reference set; --auto-gsize uses the largest input genome size. Both fix
    // the load factor (bits per distinct k-mer) so the false-positive rate is the
    // design value rather than depending on a manually chosen --gsize.
    if (cfg.autoSize) {
        size_t kmerNum = 0, dbsize = 0, sbsize = 0;
        getCardinality(kmerNum, dbsize, sbsize, cfg.kmerLen, cfg.threads, refSet1);
        if (dbsize > 0) cfg.dbfSize = dbsize;
        std::cerr << PROGRAM ": ntCard cardinality F0=" << dbsize
                  << " -> dbfSize=" << cfg.dbfSize << "\n";
    } else if (cfg.autoGsize) {
        std::vector<std::string> allRefs(refSet1);
        allRefs.insert(allRefs.end(), refSet2.begin(), refSet2.end());
        cfg.dbfSize = estimateMaxGenomeSize(allRefs);
        std::cerr << PROGRAM ": largest input genome -> dbfSize=" << cfg.dbfSize << "\n";
    }
    {
        double load = static_cast<double>(cfg.nhash) / cfg.bits;  // h*n/m at n=dbfSize
        double fpr = std::pow(1.0 - std::exp(-load), cfg.nhash);
        std::cerr << PROGRAM ": bits=" << cfg.bits << ", dbfSize=" << cfg.dbfSize
                  << ", design FPR=" << fpr << "\n";
    }

    identifyDifference(cfg, refSet1, refSet2);

#ifdef _OPENMP
    std::cout << "Runtime(sec): " << std::setprecision(4) << std::fixed
              << omp_get_wtime() - sTime << "\n";
#endif
    return 0;
}