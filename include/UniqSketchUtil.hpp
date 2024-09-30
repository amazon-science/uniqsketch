#ifndef UNIQSKETCHUTIL_HPP_
#define UNIQSKETCHUTIL_HPP_

#include <algorithm>
#include <random>
#include <unordered_map>
#include <sstream>
#include <iomanip>

#include "SequenceUtil.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"

namespace opt {
// number of references in the universe
unsigned numRef;
// number of thread
unsigned threads(1);
// bits per element in Bloom filter
unsigned bits(128);
// number of hashes for distinct Bloom filter
unsigned nhash1(5);
// number of hashes for solid Bloom filter
unsigned nhash2(5);
// length of k-mer
unsigned kmerLen(81);
// number of unique k-mers to represent a reference
int sketchnum(100);
// bit size of distinct Bloom filter
size_t m1;
// bit size of solid Bloom filter
size_t m2;
// distinct Bloom filter size (number of elements)
size_t dbfSize(0);
// solid Bloom filter size (number of elements)
size_t sbfSize(0);
// dir for storing universe unique k-mers
std::string outdir("outdir_uniqsketch");
// universe unique k-mers stat filename
std::string refstat("db_uniq_count.tsv");
// output sketch file name
std::string outfile("sketch_uniq.tsv");
// kmer spectrum range for low-complexity filtering
unsigned kRange(5);
// maxEntropy initial value computed using shannon entropy for k={1,2,3} is 12 
double maxEntropy(12.0);
// threshold for entropy score rate
double entropyThreshold(0.65);
}

typedef std::unordered_map<std::string, unsigned> SketchHash;

/**
 * Insert k-mers fasta sequences from a reference into distinct and solid Bloom filters.
 *
 * @param fPath path to a reference file.
 * @param dbFilter distinct Bloom filter for keeping track of distinct k-mers.
 * @param sbFilter solid Bloom filter for keeping track of reapeat k-mers.
 * 
 */
void loadBFfa(const std::string &fPath, BloomFilter &dbFilter, BloomFilter &sbFilter) {
    std::ifstream refFile(fPath.c_str());
    std::string seq;    
    bool good = static_cast<bool>(getline(refFile, seq));
    while (good) {
        std::string faSeq;
        good = static_cast<bool>(getline(refFile, seq));
        while (good && seq[0] != '>') {
            faSeq += seq;
            good = static_cast<bool>(getline(refFile, seq));
        }
        ntHashIterator itr(faSeq, std::max(opt::nhash1, opt::nhash2), opt::kmerLen);
        while (itr != itr.end()) {
            if (!dbFilter.insert_make_change(*itr)) {
                sbFilter.insert(*itr);
            }
            ++itr;
        }
    }
}

/**
 * Scan through all k-mers in a reference and record all unique k-mers.
 *
 * @param fPath path to a reference file.
 * @param refId integer id assigned to a reference.
 * @param sbFilter solid Bloom filter for tracking reapeat k-mers.
 * @param refCount vector of unique k-mer counts.
 * 
 */
void checkRef(const std::string &fPath, const int refId, const BloomFilter &sbFilter, std::vector<int> &refCount) {
    std::ifstream refFile(fPath.c_str());
    std::stringstream uniqstm;
    std::string sample_id = getBaseId(fPath);
    uniqstm << opt::outdir << "/" << sample_id << ".tsv";    
    std::ofstream uniqFile(uniqstm.str().c_str());
    std::string seq;
    size_t countUnique = 0;
    bool good = static_cast<bool>(getline(refFile, seq));
    std::string headerSeq = seq.substr(1, seq.size()-1);
    std::replace(headerSeq.begin(), headerSeq.end(), ' ', '_');
    while (good) {        
        std::string faSeq;
        good = static_cast<bool>(getline(refFile, seq));
        while (good && seq[0] != '>') {
            faSeq += seq;
            good = static_cast<bool>(getline(refFile, seq));
        }        
        ntHashIterator itr(faSeq, opt::nhash2, opt::kmerLen);
        while (itr != itr.end()) {
            if (!sbFilter.contains(*itr)) {
                countUnique++;
                size_t seqLoc = itr.get_pos();
                uniqFile << get_canonical(faSeq.substr(seqLoc, opt::kmerLen)) << "\t" << seqLoc << "\t" << headerSeq << std::endl;
            }
            ++itr;
        }
        if (seq.size() > 1) {
            headerSeq = seq.substr(1, seq.size()-1);
            std::replace(headerSeq.begin(), headerSeq.end(), ' ', '_');
        }
    }
    uniqFile.close();
    refCount[refId] = countUnique;
}

/**
 * Identify all unique k-mers in the reference universe.
 *
 * @param refFile vector of all references in the universe.
 * 
 */
void identifyUniqKmers(const std::vector<std::string> &refFiles) {
    std::cout << "k-mer length: " << opt::kmerLen << std::endl;
    std::cout << "Number of distinct k-mers: " << opt::dbfSize << std::endl;
    std::cout << "Number of repeat k-mers: " << opt::sbfSize << std::endl;
    std::cout << "Number of unique k-mers: " << opt::dbfSize - opt::sbfSize << std::endl;

    // distinct and solid Bloom filters to keep track of unique/repeat k-mers
    BloomFilter dbFilter(opt::m1, opt::nhash1, opt::kmerLen);
    BloomFilter sbFilter(opt::m2, opt::nhash2, opt::kmerLen);

    // load all reference sequences into distinct and solid Bloom filters
    #pragma omp parallel for schedule(dynamic)
    for (unsigned file_i = 0; file_i < refFiles.size(); ++file_i) {
        loadBFfa(refFiles[file_i], dbFilter, sbFilter);
    }

    // check all references and output unique k-mers
    std::vector<int> refCount(opt::numRef);
    #pragma omp parallel for schedule(dynamic)
    for (unsigned file_i = 0; file_i < refFiles.size(); ++file_i) {
        checkRef(refFiles[file_i], file_i, sbFilter, refCount);
    }

    // get the unique kmer stats for all references
    std::ofstream tsvOut(opt::refstat.c_str());
    for (int i =0; i < refCount.size(); i++) {    
        tsvOut << opt::outdir << "/" << getBaseId(refFiles[i]) << ".tsv\t";    
        tsvOut << refCount[i] << "\t" << i << std::endl;
    }
    tsvOut.close();

    std::cout << "Distinct BF actual fpr: " << setprecision(4) << fixed << pow((double)dbFilter.get_pop()/opt::m1,opt::nhash1) << std::endl;
    std::cout << "Solid BF actual fpr: " << setprecision(4) << fixed << pow((double)sbFilter.get_pop()/opt::m2,opt::nhash2) << std::endl;
}

/**
 * Digitize a DNA sequence
 *
 * @param seq DNA sequence to digitize.
 * 
 * @return integer representation of the DNA sequence.
 * 
 */
unsigned digitize(const std::string &seq) {
    if (seq.length() > opt::kRange) {
        return 0;
    }
    uint64_t hVal = 0;
    for(unsigned i=0; i < seq.length(); i++) {
        hVal = ((hVal << 2) | b2f[(unsigned char)seq[i]]);
    }
    return (unsigned)hVal;
}

/**
 * Check a signature for low-complexity
 *
 * @param signature signature sequence to check
 * 
 * @return bool value 
 * 
 */
bool lowComplexity(const std::string &signature) {
    if (opt::kmerLen <= opt::kRange) {
        return true;
    }

    double aggEntropy = 0.0;
    #pragma omp parallel for
    for (unsigned k = 1; k <= opt::kRange; k++) {
        std::vector<unsigned> kmerSpec((int)(2 << 2*k-1));
        for (unsigned i = 0; i < signature.length()-k+1; i++) {
            std::string seq = signature.substr(i,k);
            unsigned dmer = digitize(seq);
            kmerSpec[dmer]++;
        }        
        double entropy = 0.0;
        for (unsigned i=0; i< kmerSpec.size(); i++) {
            if (kmerSpec[i]) {
                entropy += -(double)kmerSpec[i]/(signature.size()-k+1)*log2((double)kmerSpec[i]/(signature.size()-k+1));
            }
        }
        #pragma omp atomic
        aggEntropy += entropy;
    }
    return (aggEntropy/opt::maxEntropy <= opt::entropyThreshold);
}

/**
 * Extract the set unique k-mers for a given reference.
 *
 * @param fPath path to a reference candidate unique k-mers.
 * @param refId integer id assigned to a reference.
 * @param sketchFilter Bloom filter for unique sketch.
 * @param out output file stream for th uniqsketch.
 * 
 */
void getUniqSet(const std::string &fPath, const unsigned refId, BloomFilter &sketchFilter, std::ofstream &out) {    
    // store all unique k-mers for a reference
    std::vector<std::string> sketchCandidates;
    std::ifstream sketchFile(fPath.c_str());
    std::string line;
    while (getline(sketchFile, line)) {
        sketchCandidates.push_back(line);
    }

    // shuffle the candidate list to uniformly select a subset of them
    // TODO: explore multiple conditions (GC-content, proximity to other signatures) for selecting a uniq k-mer 
    std::default_random_engine rng(0);
    std::shuffle(std::begin(sketchCandidates), std::end(sketchCandidates), rng);

    int count = 0;
    for (unsigned i = 0; i < sketchCandidates.size() && count < opt::sketchnum; i++) { 
        istringstream seqstm(sketchCandidates[i]);
        std::string useq, contig;
        unsigned pos;
        seqstm >> useq >> pos >> contig;
        if (lowComplexity(useq)) {
            continue;   
        } 

        ntHashIterator itr(useq, opt::nhash1, opt::kmerLen);
        while (itr != itr.end()) {
            if (sketchFilter.insert_make_change(*itr)) {
                count++;
                #pragma omp critical(out)
                {
                    out << useq << "\t" << refId << "\t" << getBaseId(fPath) << "\t" << opt::numRef << "\t"  << contig << "\t" << pos << std::endl;
                }
            }
            ++itr;
        }
    }
    sketchFile.close();
}

/**
 * Get the unique k-mer count for a reference.
 *
 * @param sketchLine String of reference sketch file path and its unique k-mer count.
 *
 * @return Number of unique k-mers.
 */
inline int sketchRank(const std::string &sketchLine) {
    std::stringstream ss(sketchLine);
    std::string sketchPath;
    int sketchCount, refId;
    ss >> sketchPath >> sketchCount >> refId;
    return sketchCount;
}

/**
 * construct uniquesketch for the reference universe.
 * @param statFile path to the unique kmer stat file for all references.
 */
void buildSketch(const std::string &statFile) {
    std::ifstream tsvIn(statFile.c_str());
    std::vector<std::string> uniqStats;
    std::string line;
    while (getline(tsvIn, line)) {
        uniqStats.push_back(line);
    }

    std::sort(  uniqStats.begin(), 
                uniqStats.end(), 
                [](std::string s1, std::string s2) { return sketchRank(s1) < sketchRank(s2); });

    std::ofstream out(opt::outfile.c_str());

    BloomFilter sketchFilter(opt::bits*opt::sketchnum*uniqStats.size(), opt::nhash1, opt::kmerLen);

    for (unsigned i = 4; i <= opt::kRange; i++) {
        opt::maxEntropy += log2(opt::kmerLen - i + 1);
    }

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < uniqStats.size(); i++) {
        std::stringstream uniqstm(uniqStats[i]);
        std::string sketch_path;
        size_t sketch_count, sketch_refId;
        uniqstm >> sketch_path >> sketch_count >> sketch_refId;
        getUniqSet(sketch_path, sketch_refId, sketchFilter, out);
    }
    out.close();
}

#endif // UNIQSKETCHUTIL_HPP_
