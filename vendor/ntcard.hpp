/*
  * Copyright Hamid Mohamadi.
  * SPDX-License-Identifier: MIT
  *
  * Licensed under the MIT License. See the LICENSE accompanying this file
  * for the specific language governing permissions and limitations under
  * the License.
  */
 
// Adapted from ntCard https://github.com/bcgsc/ntCard


#ifndef NTCARD_HPP_
#define NTCARD_HPP_

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>

#include "ntHashIterator.hpp"
#include "seqio.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace klibpp;

namespace opt {
size_t rBuck;
unsigned rBits(27);
unsigned sBits(11);
unsigned sMask;
unsigned covMax(65536);
size_t nSamp(2);
size_t nK(1);
bool samH(true);
}

inline void ntComp(const uint64_t hVal, uint16_t *t_Counter) {
    uint64_t indBit = opt::nSamp;
    if (hVal >> (64 - opt::sBits) == 1) {
        indBit = 0;
    }
    if (hVal >> (64 - opt::sBits) == opt::sMask) {
        indBit = 1;
    }
    if (indBit < opt::nSamp) {
        size_t shVal = hVal & (opt::rBuck - 1);
        #pragma omp atomic
        ++t_Counter[indBit * opt::rBuck + shVal];
    }

}

inline void ntRead(const KSeq &record, const vector<unsigned> &kList, uint16_t *t_Counter, size_t totKmer[]) {
    for (unsigned k = 0; k < kList.size(); k++) {
        ntHashIterator itr(record.seq, 1, kList[k]);
        while (itr != itr.end()) {
            ntComp((*itr)[0], t_Counter + k * opt::nSamp * opt::rBuck);
            ++itr;
            ++totKmer[k];
        }
    }
}

void compEst(const uint16_t *t_Counter, double &F0Mean, vector<double> &fMean) {
    vector<vector<unsigned>> p(opt::nSamp, vector<unsigned>(opt::covMax, 0));
    for (size_t i = 0; i < opt::nSamp; i++) {
        for (size_t j = 0; j < opt::rBuck; j++) {
            ++p[i][t_Counter[i * opt::rBuck + j]];
        }
    }

    vector<double> pMean(opt::covMax, 0.0);
    for (size_t i=0; i < opt::covMax; i++) {
        for (size_t j=0; j < opt::nSamp; j++) {
            pMean[i]+=p[j][i];
        }
        pMean[i] /= 1.0 * opt::nSamp;
    }

    F0Mean = (ssize_t)((opt::rBits * log(2) - log(pMean[0])) * 1.0 * ((size_t)1 << (opt::sBits + opt::rBits)));

    fMean[1] = -1.0 * pMean[1] / (pMean[0] * (log(pMean[0]) - opt::rBits * log(2)));
    for (size_t i = 2; i < opt::covMax; i++) {
        double sum = 0.0;
        for (size_t j = 1; j < i; j++) {
            sum += j * pMean[i-j] * fMean[j];
        }
        fMean[i] = -1.0 * pMean[i] / (pMean[0] * (log(pMean[0]) - opt::rBits * log(2))) - sum / (i * pMean[0]);
    }
    for (size_t i = 1; i < opt::covMax; i++) {
        fMean[i] = abs((ssize_t)(fMean[i] * F0Mean));
    }
}

bool getCardinality(size_t &kmerNum, size_t &dbsize, size_t &sbsize, const unsigned kmer_len, const unsigned num_thread, const vector<string> &inFiles) {
    vector<unsigned> kList;
    kList.push_back(kmer_len);
    vector<size_t> totalKmers(kList.size(), 0);
    opt::rBuck = ((size_t)1) << opt::rBits;
    opt::sMask = (((size_t)1) << (opt::sBits - 1)) - 1;
    uint16_t *t_Counter = new uint16_t [opt::nK * opt::nSamp * opt::rBuck]();

#ifdef _OPENMP
    omp_set_num_threads(num_thread);
#endif

    #pragma omp parallel for schedule(dynamic)
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
        size_t totKmer[kList.size()];
        for (unsigned k = 0; k < kList.size(); k++) {
            totKmer[k]=0;
        }
        SeqStreamIn iss(inFiles[file_i].c_str());
        KSeq record;
        while (iss >> record) {
            ntRead(record, kList, t_Counter, totKmer);
        }
        for (unsigned k = 0; k < kList.size(); k++) {
            #pragma omp atomic
            totalKmers[k] += totKmer[k];
        }
    }
    double F0Mean = 0.0;
    vector<double> fMean(opt::covMax, 0.0);
    compEst(t_Counter, F0Mean, fMean);

    kmerNum = totalKmers[0];
    dbsize = (size_t)F0Mean;
    sbsize = dbsize - (size_t)fMean[1];

    delete [] t_Counter;
    return 0;
}

#endif // NTCARD_HPP_
