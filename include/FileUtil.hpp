#ifndef FILEUTIL_HPP_
#define FILEUTIL_HPP_

#include <sys/stat.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/**
 * Check if a file exists and is readable.
 *
 * @param fPath  Path to the file.
 * @return True if the file exists and can be opened for reading.
 */
inline bool fileExists(const std::string& fPath) {
    std::ifstream f(fPath);
    return f.good();
}

/**
 * Validate that a file exists, printing an error and returning false if not.
 *
 * @param fPath  Path to the file.
 * @param context  Description of what the file is for (e.g. "reference file", "read file").
 * @return True if the file exists.
 */
inline bool validateFile(const std::string& fPath, const std::string& context) {
    if (!fileExists(fPath)) {
        std::cerr << "Error: " << context << " not found: " << fPath << "\n";
        return false;
    }
    return true;
}

/**
 * Validate that all files in a vector exist.
 *
 * @param files  Vector of file paths.
 * @param context  Description of what the files are for.
 * @return True if all files exist.
 */
inline bool validateFiles(const std::vector<std::string>& files, const std::string& context) {
    bool allGood = true;
    for (const auto& f : files) {
        if (!validateFile(f, context)) {
            allGood = false;
        }
    }
    return allGood;
}

/**
 * Prepend a prefix to the basename of a path, preserving the directory component.
 * e.g. prefixBasename("../out/x.tsv", "log_") -> "../out/log_x.tsv"
 *      prefixBasename("x.tsv", "log_")        -> "log_x.tsv"
 *
 * @param fPath  File path (may include directory components).
 * @param prefix  String to prepend to the basename.
 * @return Path with prefix applied to the basename only.
 */
inline std::string prefixBasename(const std::string& fPath, const std::string& prefix) {
    size_t slash = fPath.rfind('/');
    if (slash == std::string::npos) {
        return prefix + fPath;
    }
    return fPath.substr(0, slash + 1) + prefix + fPath.substr(slash + 1);
}

/**
 * Estimate the largest genome size (in bases) among a set of files, using file
 * size as a fast proxy. Used to size Bloom filters to the largest input genome
 * without scanning the sequences. FASTA headers/newlines overestimate by a few
 * percent, which is harmless (and conservative) for Bloom-filter sizing.
 *
 * @param files  Vector of file paths.
 * @return Largest file size in bytes, or 5000000 if none could be stat'd.
 */
inline size_t estimateMaxGenomeSize(const std::vector<std::string>& files) {
    size_t maxSize = 0;
    struct stat st;
    for (const auto& f : files) {
        if (stat(f.c_str(), &st) == 0 && static_cast<size_t>(st.st_size) > maxSize) {
            maxSize = static_cast<size_t>(st.st_size);
        }
    }
    return maxSize > 0 ? maxSize : 5000000;
}

/**
 * Resolve Bloom-filter sizing parameters, falling back to an exact bound when the
 * ntCard streaming estimate is unreliable on small inputs.
 *
 * ntCard estimates the number of distinct (dbfSize) and solid (sbfSize) k-mers by
 * sampling hashed k-mers into a fixed set of buckets. On very small inputs (for
 * example a handful of short references at a coarse tree level, or a tiny sample)
 * too few k-mers reach the sample buckets and the estimate degenerates: it can come
 * back as 0, or, through statistical noise, larger than the total number of k-mers
 * actually present. Either way the original guards rejected the input.
 *
 * The total k-mer count is known exactly (counted, not sampled, by getCardinality)
 * and is a hard upper bound on the distinct count, since distinct <= total. So when
 * an estimate is out of range we clamp it to that exact total. This never undersizes
 * the Bloom filters (an oversized filter only lowers the false-positive rate and, for
 * the small inputs that trigger this, costs trivial memory). For healthy inputs the
 * estimates already satisfy 0 < estimate <= total, so this is a no-op.
 *
 * @param totalKmers  Exact total number of k-mers processed (from getCardinality).
 * @param dbfSize  Distinct k-mer estimate to validate/clamp (in/out).
 * @param sbfSize  Solid k-mer estimate to validate/clamp (in/out).
 * @return True if there is any k-mer content to build from (totalKmers > 0).
 */
inline bool resolveBloomSizes(size_t totalKmers, size_t& dbfSize, size_t& sbfSize) {
    if (totalKmers == 0) {
        return false;  // genuinely no k-mers (e.g. all sequences shorter than k)
    }
    if (dbfSize == 0 || dbfSize > totalKmers) {
        dbfSize = totalKmers;
    }
    // Solid k-mers are a subset of the distinct k-mers, so sbfSize must satisfy
    // 0 < sbfSize <= dbfSize. ntCard's sampled estimate can violate this on small
    // inputs (returning 0, or even a value above dbfSize); clamp to dbfSize, which is
    // the correct upper bound and keeps the distinct-minus-solid diagnostic sensible.
    if (sbfSize == 0 || sbfSize > dbfSize) {
        sbfSize = dbfSize;
    }
    return true;
}

#endif // FILEUTIL_HPP_
