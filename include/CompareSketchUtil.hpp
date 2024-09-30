#ifndef COMPARESKETCH_HPP_
#define COMPARESKETCH_HPP_

#include "SequenceUtil.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"
#include "ntcard.hpp"

namespace opt {
// number of thread
unsigned threads(1);
// bits per element in Bloom filter
unsigned bits(16);
// number of hashes for distinct Bloom filter
unsigned nhash(3);
// length of k-mer
unsigned kmerLen(81);
// distinct Bloom filter size (number of elements, approximately genome size)
size_t dbfSize(5000000);
// output sketch file name
std::string outfile("reference_similarity.tsv");
}


/**
 * Insert k-mers fasta sequences from a reference into a distinct Bloom filter.
 *
 * @param fPath path to a reference file.
 * @param dbFilter distinct Bloom filter for keeping track of distinct k-mers.
 * 
 */
void loadBloomFilter(const std::string &fPath, BloomFilter &dbFilter) {
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
        ntHashIterator itr(faSeq, opt::nhash, opt::kmerLen);
        while (itr != itr.end()) {
            dbFilter.insert(*itr);
            ++itr;
        }
    }
}

/**
 * Insert k-mers fasta sequences from a reference into a distinct Bloom filter.
 *
 * @param fPath path to the second reference file to be checked against first reference Bloom filter.
 * @param dbFilter distinct Bloom filter for the first reference distinct k-mers.
 * @param countAll total number of k-mers in the second reference.
 * @param countUnique number of unique k-mers for the second reference.
 * 
 * @return double representing the similarity rate between two references.
 * 
 */
double checkBloomFilter(const std::string &fPath, const BloomFilter &dbFilter, size_t &countAll, size_t &countUnique) {
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
        ntHashIterator itr(faSeq, opt::nhash, opt::kmerLen);
        while (itr != itr.end()) {
            countAll++;
            if (!dbFilter.contains(*itr)) {
                countUnique++;
            }
            ++itr;
        }
    }
    return ((double)(countAll - countUnique) / countAll);    
}

/**
 * Calculate similarity/difference between two sets of references. M*N pairwise comparison
 *
 * @param refSet1 list of all references in the first set.
 * @param refSet2 list of all references in the second set.
 * 
 */
void identifyDiffernce(const std::vector<std::string> &refSet1, const std::vector<std::string> &refSet2) {
    std::ofstream sim_out(opt::outfile.c_str());
    sim_out << "ref1\tref2\tuniq_ref1\tall_ref1\tsimilarity\n";

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < refSet1.size(); i++) {
        // distinct filter to keep track of distinct k-mers in the first reference
        BloomFilter dbFilter(opt::bits*opt::dbfSize, opt::nhash, opt::kmerLen);

        // load all reference sequences into distinct and solid Bloom filters
        loadBloomFilter(refSet1[i], dbFilter);

        // check the reference refSet1[i] against all references in refSet2
        for (unsigned f_index = 0; f_index < refSet2.size(); f_index++) {
            // check aganist all references and output similarity score
            size_t countAll = 0, countUnique = 0;
            double similarity_score = checkBloomFilter(refSet2[f_index], dbFilter, countAll, countUnique);
            #pragma omp critical(simout) 
            sim_out << getBaseId(refSet2[f_index]) << "\t" << getBaseId(refSet1[i]) << "\t" << countUnique << "\t" << countAll << "\t" << similarity_score << "\n";
        }
    }
    sim_out.close();
}

#endif // COMPARESKETCH_HPP_
