#ifndef QUERYSKETCHUTIL_HPP_
#define QUERYSKETCHUTIL_HPP_

#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>

#include "SequenceUtil.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"

#include "seqio.hpp"

using namespace klibpp;

namespace opt {
// k-mer length
unsigned kmerLen;
// number of threads
unsigned threads(1);
// output file name
std::string out("out_sketch.tsv");
// read 1 input
std::string r1;
// read 2 input
std::string r2;
// input uniqsketch reference
std::string ref;
// number of hashes for Bloom filter
unsigned nhash(3);
// distinct Bloom filter size
size_t dbfSize;
// solid Bloom filter size
size_t sbfSize;
// bits per element in Bloom filter
unsigned bits(64);
//number of unique k-mers observed to call a hit
int sketchHit(10);
// flag to only use solid k-mers from reads
int solid(0);
// abundance cutoff 
double abundanceCutoff(0.0);
// read number cutoff
int readCutoff(2);
}

typedef std::unordered_map<std::string, unsigned> SketchHash;

/**
 * Load the uniqsketch of the reference universe into a hash table.
 *
 * @param fPath path to the uniquesketch file.
 * @param sketchHash Hash table for uniquesketch.
 * @param sketchRef Vector for mapping refId--refName
 * @param refSigCount Vector for keeping the count of each uniq k-mer for each reference
 *
 */
void loadSketch(const std::string &fPath, SketchHash &sketchHash, std::vector<std::string> &sketchRef, std::vector<SketchHash> &refSigCount) {    
    std::ifstream sketchFile(fPath.c_str());
    bool init = true;
    std::string seq;
    while (getline(sketchFile, seq)) {
        std::string uniqSeq, refName, contig;
        unsigned refId, refNum, pos;
        std::istringstream seqstm(seq);
        seqstm >> uniqSeq >> refId >> refName >> refNum >> contig >> pos;
        if (seqstm.fail()) {
            std::cerr << "Error in parsing uniqsketch entry:\n" << seq << std::endl;
            exit(EXIT_FAILURE);
        }
        sketchHash.emplace(uniqSeq, refId);
        if (init) {
            opt::kmerLen = uniqSeq.length();
            sketchRef.resize(refNum);
            refSigCount.resize(refNum);
            init = false;
        }
        sketchRef[refId] = refName;
        refSigCount[refId].emplace(uniqSeq, 0);
    }
    sketchFile.close();
}

/**
 * Query an input read sequence file against uniqsketch and update the count of observed unique k-mers for references.
 *
 * @param in input file stream for a read file.
 * @param sketchHash Hash table for uniqsketch.
 * @param dbFilter Bloom filter for tracking distinct k-mers.
 * @param sbFilter Bloom filter for tracking k-mers observed at least twice.
 * @param sketchCount Vector for keeping the count of uniq k-mers for each reference
 * @param refSigCount Vector for keeping the count of each uniq k-mer for each reference
 * @param refRead Vector for keeping distinct reads for each reference
 *
 */
void querySample(const string& in, const SketchHash &sketchHash, BloomFilter &dbFilter, BloomFilter &sbFilter, std::vector<unsigned> &sketchCount, std::vector<SketchHash> &refSigCount, std::vector<std::unordered_set<std::string> > &refRead) {
    bool good = true;
    SeqStreamIn iss(in.c_str());    
    #pragma omp parallel
    for (KSeq record; good;) {
        #pragma omp critical(in)
        good = (iss >> record);
        if (good && record.seq.length() >= opt::kmerLen) {
            ntHashIterator itr(record.seq, opt::nhash, opt::kmerLen);
            while (itr != itr.end()) {
                // if solid is set, ignore weak kmers (appearing just once in data mainly due to sequencing error)
                if (opt::solid) {
                    if (dbFilter.insert_make_change(*itr)) {
                        ++itr;
                        continue;
                    }
                }

                std::string kmer = get_canonical(record.seq.substr(itr.get_pos(), opt::kmerLen));
                auto query = sketchHash.find(kmer);
                if (query != sketchHash.end()) {
                    unsigned refId = query->second;

                    #pragma omp atomic
                    ++sketchCount[refId];

                    std::string read_id = record.name;
                    #pragma omp critical(refout)
                    refRead[refId].insert(read_id);

                    #pragma omp atomic
                    refSigCount[refId][kmer]++;
                }
                ++itr;
            }
        }
    }
}

/**
 * Query batch of input read files against uniqsketch and update the count of observed unique k-mers for references.
 *
 * @param sampleFiles all input read files.
 * @param sketchHash Hash table for uniqsketch.
 * @param sketchCount Vector for keeping the count of uniq k-mers for each reference
 * @param refSigCount Vector for keeping the count of each uniq k-mer for each reference
 * @param refRead Vector for keeping distinct reads for each reference
 *
 */
void querySampleBatch(const std::vector<std::string> &sampleFiles, const SketchHash &sketchHash, std::vector<unsigned> &sketchCount, std::vector<SketchHash> &refSigCount, std::vector<std::unordered_set<std::string> > &refRead) {
    BloomFilter dbFilter(opt::dbfSize*opt::bits, opt::nhash, opt::kmerLen);
    BloomFilter sbFilter(opt::sbfSize*opt::bits, opt::nhash, opt::kmerLen);
    
    // query all input reads against uniqsketch  
    for (unsigned i = 0; i < sampleFiles.size(); i++) {
        querySample(sampleFiles[i], sketchHash, dbFilter, sbFilter, sketchCount, refSigCount, refRead);
    }

    // increase non-zero sketchCount entries by one to adjust for distinct and solid Bloom filter stages mask
    if (opt::solid) {
        for (auto &itr: sketchCount) {
            if (itr) {
                ++itr;
            }
        }
    }
}

/**
 * Generate the query result against uniqsketch and output to a tsv file.
 *
 * @param sketchCount Vector for keeping the count of uniq k-mers for each reference
 * @param sketchRef Vector for mapping refId--refName
 * @param refSigCount Vector for keeping the count of each uniq k-mer for each reference
 * @param refRead Vector for keeping distinct reads for each reference
 * @param fPath path to the output result tsv file.
 */
void generateQueryResult(std::vector<unsigned> &sketchCount, const std::vector<std::string> &sketchRef, const std::vector<SketchHash> &refSigCount, const std::vector<std::unordered_set<std::string> > &refRead, const std::string &fPath) {
    // sort the index of sketchCount vector to keep the refId-refName mapping intact
    std::vector<unsigned> sketchIndex(sketchCount.size());
    std::iota(std::begin(sketchIndex), std::end(sketchIndex), 0);
    std::sort(std::begin(sketchIndex), std::end(sketchIndex), [&](int i, int j) { return sketchCount[i] > sketchCount[j]; });

    size_t sketchSum = 0;
    for (unsigned i = 0; i < sketchCount.size(); i++) {
        if (refRead[i].size() > opt::readCutoff) {
            if (sketchCount[i] >= opt::sketchHit) {
                sketchSum += sketchCount[i];
            }
        }
    }

    // generate output
    std::ofstream outFile(opt::out.c_str());
    outFile << "ref\tabundance\tcount\n";

    if (sketchSum != 0) {
        for (unsigned i = 0; i < sketchCount.size(); i++) {
            if (refRead[sketchIndex[i]].size() > opt::readCutoff) {
                if (sketchCount[sketchIndex[i]] >= opt::sketchHit) {
                    double abundance = (double)sketchCount[sketchIndex[i]]/sketchSum;
                    if (abundance > opt::abundanceCutoff) {
                        outFile << sketchRef[sketchIndex[i]] << "\t" << abundance << "\t" << sketchCount[sketchIndex[i]] << std::endl;
                    }
                }
            }
        }
    }        
    outFile.close();

    // generate detailed log for output
    std::ofstream logRefSigCount(("log_" + opt::out).c_str());
    logRefSigCount << "ref\ttotal_reads\tsignature\tcount\n";
    for (unsigned i = 0; i < refSigCount.size(); i++) {        
        for (auto itr: refSigCount[i]) {
            if (itr.second) {
                logRefSigCount << sketchRef[i] <<  "\t" << refRead[i].size() << "\t" << itr.first << "\t" << itr.second << "\n";
            }
        }
    }
}

#endif // QUERYSKETCHUTIL_HPP_
