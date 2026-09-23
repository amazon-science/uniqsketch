#ifndef FILEUTIL_HPP_
#define FILEUTIL_HPP_

#include <sys/stat.h>
#include <zlib.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/**
 * Line reader over a plain or gzip-compressed file, with std::getline's
 * contract: a line is the bytes up to a '\n' (not included), the last line
 * needs no newline, an empty line is a line, and getline() is false only at
 * the end of the file with nothing read. zlib reads a file that is not gzip
 * as it is, so one code path serves both.
 *
 * Every reference scanner reads through this. Before, they read with
 * std::ifstream, and a gzipped FASTA given as a reference was scanned as
 * compressed bytes: the k-mer counts (ntCard, which did read gzip) sized the
 * filters correctly, the genome contributed no k-mers and got no candidates,
 * every k-mer it shared with another genome became that genome's "unique"
 * k-mer, and the build exited 0 without a word.
 */
class GzLineReader {
public:
    explicit GzLineReader(const std::string& path) : fp_(gzopen(path.c_str(), "rb")) {
        if (fp_) gzbuffer(fp_, 1 << 20);
    }
    ~GzLineReader() {
        if (fp_) gzclose(fp_);
    }
    GzLineReader(const GzLineReader&) = delete;
    GzLineReader& operator=(const GzLineReader&) = delete;

    /** False if the file could not be opened. */
    explicit operator bool() const { return fp_ != nullptr; }

    bool getline(std::string& line) {
        line.clear();
        bool any = false;
        for (;;) {
            if (pos_ == len_ && !fill()) return any;
            const char* start = buf_.data() + pos_;
            const char* nl = static_cast<const char*>(std::memchr(start, '\n', len_ - pos_));
            if (nl) {
                line.append(start, nl - start);
                pos_ = static_cast<size_t>(nl - buf_.data()) + 1;
                return true;
            }
            line.append(start, len_ - pos_);
            pos_ = len_;
            any = true;
        }
    }

private:
    bool fill() {
        if (!fp_) return false;
        int n = gzread(fp_, buf_.data(), static_cast<unsigned>(buf_.size()));
        if (n <= 0) return false;
        pos_ = 0;
        len_ = static_cast<size_t>(n);
        return true;
    }

    gzFile fp_;
    std::vector<char> buf_ = std::vector<char>(1 << 16);
    size_t pos_ = 0, len_ = 0;
};

/**
 * Bytes of sequence data in a reference file, for sizing filters: the size of a
 * plain file; for a gzip file the uncompressed size, taken from the trailer
 * (ISIZE, exact for the single-member file `gzip` writes) or, when the trailer
 * cannot be right (a multi-member bgzip file, or over 4 GB), by decompressing
 * once. Headers and newlines are counted too, which only makes a filter larger.
 * Sizing from the compressed size would give a gzipped genome a quarter of the
 * bits it needs.
 */
inline size_t uncompressedSize(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return 0;
    const size_t compressed = static_cast<size_t>(st.st_size);
    std::ifstream in(path, std::ios::binary);
    unsigned char magic[2] = {0, 0};
    in.read(reinterpret_cast<char*>(magic), 2);
    if (!in || magic[0] != 0x1f || magic[1] != 0x8b) return compressed;   // not gzip
    unsigned char t[4] = {0, 0, 0, 0};
    in.seekg(-4, std::ios::end);
    in.read(reinterpret_cast<char*>(t), 4);
    size_t isize = t[0] | (size_t(t[1]) << 8) | (size_t(t[2]) << 16) | (size_t(t[3]) << 24);
    if (in && isize >= compressed) return isize;
    gzFile fp = gzopen(path.c_str(), "rb");
    if (!fp) return compressed;
    gzbuffer(fp, 1 << 20);
    std::vector<char> buf(1 << 20);
    size_t total = 0;
    for (int n; (n = gzread(fp, buf.data(), static_cast<unsigned>(buf.size()))) > 0;) total += static_cast<size_t>(n);
    gzclose(fp);
    return total ? total : compressed;
}

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
 * Estimate the largest genome size (in bases) among a set of files, using the
 * (uncompressed) file size as a fast proxy. Used to size Bloom filters to the
 * largest input genome without scanning the sequences. FASTA headers/newlines
 * overestimate by a few percent, which is harmless (and conservative) for
 * Bloom-filter sizing.
 *
 * @param files  Vector of file paths.
 * @return Largest file size in bytes, or 5000000 if none could be stat'd.
 */
inline size_t estimateMaxGenomeSize(const std::vector<std::string>& files) {
    size_t maxSize = 0;
    for (const auto& f : files) maxSize = std::max(maxSize, uncompressedSize(f));
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
