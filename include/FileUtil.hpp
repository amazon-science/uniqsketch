#ifndef FILEUTIL_HPP_
#define FILEUTIL_HPP_

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

#endif // FILEUTIL_HPP_
