#ifndef SEQUENCEUTIL_HPP_
#define SEQUENCEUTIL_HPP_

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

/**
 * Encode nucleotide to integer.
 * A,a: 00 = 0x00
 * C,c: 01 = 0x01
 * G,g: 10 = 0x02
 * T,t: 11 = 0x03
 */
static const uint8_t b2f[256] = {
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // 0
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // 1
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // 2
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // 3
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0x00, 0xFF, 0x01, 0xFF, 0xFF, 0xFF, 0x02, // 4   'A' 'C' 'G'
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0x03, 0xFF, 0xFF, 0xFF, // 5   'T'
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0x00, 0xFF, 0x01, 0xFF, 0xFF, 0xFF, 0x02, // 6   'a' 'c' 'g'
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0x03, 0xFF, 0xFF, 0xFF, // 7   't'
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // 8
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // 9
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // 10
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // 11
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // 12
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // 13
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // 14
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // 15
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF
};

/**
 * Nucleotide to complement character table.
 * A,a → T    C,c → G    G,g → C    T,t → A
 */
static const unsigned char b2r[256] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 0
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 1
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 2
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 3
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', // 4   'A' 'C' 'G'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', // 5   'T'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', // 6   'a' 'c' 'g'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', // 7   't'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 8
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 9
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 10
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 11
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 12
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 13
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 14
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 15
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

/**
 * Return the canonical representation for a given k-mer.
 *
 * @param kmer  Input k-mer sequence.
 * @return Canonical (lexicographically smaller of forward/reverse-complement) representation.
 */
static inline std::string get_canonical(const std::string& kmer) {
    std::string canonical = kmer;
    int p = 0;
    int hLen = static_cast<int>((kmer.length() - 1) / 2);
    while (canonical[p] == b2r[static_cast<unsigned char>(canonical[kmer.length() - 1 - p])]) {
        ++p;
        if (p >= hLen) break;
    }
    if (canonical[p] > b2r[static_cast<unsigned char>(canonical[kmer.length() - 1 - p])]) {
        for (int l = p, r = static_cast<int>(kmer.length()) - 1 - p; l <= r; ++l, --r) {
            char tmp = b2r[static_cast<unsigned char>(canonical[r])];
            canonical[r] = b2r[static_cast<unsigned char>(canonical[l])];
            canonical[l] = tmp;
        }
    }
    return canonical;
}

/**
 * Return the file basename (without extension) from a full path.
 *
 * @param fPath  Full file path string.
 * @return Base name of the file without extension.
 */
static inline std::string getBaseId(const std::string& fPath) {
    size_t posPath = fPath.rfind('/');
    std::string fName = (posPath == std::string::npos)
                        ? fPath
                        : fPath.substr(posPath + 1);
    // A compression suffix is not the extension: x.fa.gz names the same
    // reference as x.fa.
    for (const char* z : {".gz", ".bgz"}) {
        std::string suffix(z);
        if (fName.size() > suffix.size() &&
            fName.compare(fName.size() - suffix.size(), suffix.size(), suffix) == 0) {
            fName.erase(fName.size() - suffix.size());
            break;
        }
    }
    size_t posFile = fName.rfind('.');
    if (posFile == std::string::npos) {
        return fName;
    }
    if (posFile == 0) {
        std::cerr << "Error in file name: " << fPath << "\n";
        exit(EXIT_FAILURE);
    }
    return fName.substr(0, posFile);
}

#endif // SEQUENCEUTIL_HPP_