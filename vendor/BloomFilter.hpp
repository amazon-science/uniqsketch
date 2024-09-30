/*
  * Copyright Hamid Mohamadi.
  * SPDX-License-Identifier: MIT
  *
  * Licensed under the MIT License. See the LICENSE accompanying this file
  * for the specific language governing permissions and limitations under
  * the License.
  */

// Adapted from ntHash https://github.com/bcgsc/ntHash


#ifndef BLOOMFILTER_HPP_
#define BLOOMFILTER_HPP_

#include <stdio.h>

#include <fstream>
#include <sstream>
#include <climits>

#include "nthash.hpp"

using namespace std;

/**
 * Calculate the set bits in a byte.
 *
 * @param x input byte 
 * 
 * @return the number of bits that are set to 1.
 * 
 */
inline unsigned pop_count(uint8_t x) {
    return ((0x876543210 >>
             (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >>
            ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2))
           & 0xf;
}

class BloomFilter {
public:

    /**
     * Load Bloom filter from an existing file.
     *
     * @param filterSize byte size of the Bloom filter.
     * @param hashNum number of hashed used in Bloom filter.
     * @param kmerSize size of kmer used to build the Bloom filter.
     * @param fPath path to the file Bloom filter is stored.
     * 
     */
    BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize, const char * fPath):
        m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize) {
        m_filter = new unsigned char [(m_size + CHAR_BIT - 1) / CHAR_BIT];
        std::ifstream myFile(fPath, ios::in | ios::binary);
        myFile.seekg(0, ios::beg);
        myFile.read((char *)m_filter, (m_size + CHAR_BIT - 1) / CHAR_BIT);
        myFile.close();
    }

    /**
     * Constructor for Bloom filter 
     *
     * @param filterSize byte size of the Bloom filter.
     * @param hashNum number of hashed used in Bloom filter.
     * @param kmerSize size of kmer used to build the Bloom filter.
     * 
     */
    BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize):
        m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize) {
        m_filter = new unsigned char [(m_size + CHAR_BIT - 1) / CHAR_BIT];
        for(size_t i = 0; i < (m_size + CHAR_BIT - 1) / CHAR_BIT; i++) {
            m_filter[i]=0;
        }
    }

    /**
     * Insert a kmer into Bloom filter using precomputed kmer hash values and check if the kmer was already there.
     *
     * @param hVal the set of hash values for a kmer. 
     * 
     * @return True if the kmer is inserted for the first time, otherwise False.
     * 
     */
    bool insert_make_change(const uint64_t *hVal) {
        bool change = false;
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            unsigned char old_byte = __sync_fetch_and_or(&m_filter[hLoc / CHAR_BIT],(1 << (CHAR_BIT - 1 - hLoc % CHAR_BIT)));
            if ((old_byte & (1 << (CHAR_BIT - 1 - hLoc % CHAR_BIT))) == 0) {
                change = true;
            }
        }
        return change;
    }

    /**
     * Insert a kmer into Bloom filter using precomputed kmer hash values.
     *
     * @param hVal the set of hash values for a kmer. 
     * 
     */
    void insert(const uint64_t *hVal) {
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            __sync_or_and_fetch(&m_filter[hLoc / CHAR_BIT], (1 << (CHAR_BIT - 1 - hLoc % CHAR_BIT)));
        }
    }

    /**
     * Insert a kmer into Bloom filter.
     *
     * @param kmer the kmer sequence.  
     * 
     */
    void insert(const char* kmer) {
        uint64_t *hVal = new uint64_t[m_hashNum];
        NTMC64(kmer, m_kmerSize, m_hashNum, hVal);
        insert(hVal);
        delete [] hVal;
    }

    /**
     * Query a kmer against Bloom filter using precomputed kmer hash values.
     *
     * @param hVal the set of hash values for a kmer 
     *
     * @return True if the set of hash values are all set, otherwise False.
     * 
     */
    bool contains(const uint64_t *hVal) const {
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            if ((m_filter[hLoc / CHAR_BIT] & (1 << (CHAR_BIT - 1 - hLoc % CHAR_BIT))) == 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * Query a kmer against Bloom filter.
     *
     * @param kmer the kmer sequence
     *
     * @return True if the kmer is Bloom filter, otherwise False.
     * 
     */
    bool contains(const char* kmer) const { 
        uint64_t *hVal = new uint64_t[m_hashNum];
        NTMC64(kmer, m_kmerSize, m_hashNum, hVal);
        bool result = contains(hVal);
        delete [] hVal;
        return result;
    }

    /**
     * Store Bloom filter to file.
     *
     * @param fPath file path for Bloom filter to strore.
     * 
     */
    void storeFilter(const char * fPath) const {
        std::ofstream myFile(fPath, ios::out | ios::binary);
        myFile.write(reinterpret_cast<char*>(m_filter), (m_size + CHAR_BIT - 1) / CHAR_BIT);
        myFile.close();
    }

    /**
     * Compute the Bloom filter population count, the number of set bits.
     *
     * @return the number of set bits in Bloom filter.
     * 
     */
    size_t get_pop() const {
        size_t i, popBF = 0;
        #pragma omp parallel for reduction(+:popBF)
        for (i = 0; i < (m_size + CHAR_BIT - 1) / CHAR_BIT; i++) {
            popBF = popBF + pop_count(m_filter[i]);
        }
        return popBF;
    }

    unsigned getHashNum() const {
        return m_hashNum;
    }

    unsigned getKmerSize() const {
        return m_kmerSize;
    }

    ~BloomFilter() {
        delete[] m_filter;
    }

private:
    BloomFilter(const BloomFilter& that);
    unsigned char * m_filter;
    size_t m_size;
    unsigned m_hashNum;
    unsigned m_kmerSize;
};

#endif // BLOOMFILTER_HPP_
