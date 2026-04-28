# Changelog

All notable changes to this project will be documented in this file.

## 1.2.2 (2026-04-28)

### Bug Fixes
* Fixed `querysketch` log files (`log_*.tsv`, `logread_*.tsv`) being silently dropped when `-o` includes a directory path. Prefixes are now applied to the basename only, so logs land next to the main output file.

## 1.2.1 (2026-04-06)
* Aligned tool `--version` output across `uniqsketch`, `querysketch`, and `comparesketch`.

## 1.2.0 (2026-03-23)

### Performance
* Replaced kseq++ (C++ wrapper) with kseq.h (Heng Li) for faster sequence I/O (~13% speedup in querysketch)
* Batch-read pattern in querySample eliminates per-read lock contention under OpenMP
* Rolling hash in lowComplexity filter replaces per-position recomputation
* Stack-allocated spectrum arrays in lowComplexity replace heap-allocated vectors
* Replaced std::endl with "\n" in hot loops to avoid unnecessary stream flushes
* BloomFilter constructor uses value-initialized allocation instead of manual zeroing
* comparesketch pre-caches all k-mer hashes to avoid redundant file I/O (M+N reads instead of M×N)
* comparesketch reuses per-thread Bloom filters with memset instead of re-allocating
* Sort lambda and range-for loops use const references to avoid copies

### Features
* Genomically spaced signature selection in buildSketch — signatures are evenly distributed across reference genome positions instead of randomly shuffled
* Added read-level log output (logread_*.tsv) in querysketch listing matched read IDs per reference
* comparesketch now supports compressed FASTA input (.fa.gz) via kseq.h/zlib
* Added CMake build system (CMakeLists.txt) alongside existing Makefile
* Added `--cluster=N` option to uniqsketch: automatically clusters nearly identical references using single-linkage union-find on pairwise unique k-mer counts, selects a representative per cluster, and outputs `clusters.tsv` mapping all members to their representatives

### Error Handling
* Added file existence validation for all input files across all three tools
* Typo'd filenames now produce clear error messages instead of core dumps
* Missing @list files are caught before processing
* Bloom filter error message now reports k-mer counts, k value, and suggests causes/fixes
* Empty lines in list files are skipped

### Code Quality
* Modern C++17 style throughout: nullptr, using aliases, static_cast, structured bindings
* Removed using namespace std from headers (BloomFilter.hpp, ntcard.hpp)
* All headers are self-contained with explicit includes
* Fixed typo: identifyDiffernce → identifyDifference
* Fixed version string bug in querysketch ("querysketch1.1.0" → "querysketch Version 1.1.0")
* Fixed signed/unsigned comparison warnings in readCutoff checks
* Consistent Doxygen-style documentation on all public functions
* Removed unused b2c lookup table from SequenceUtil.hpp
* Simplified getBaseId to handle edge cases (no slash, no extension, multiple dots)

### Tests
* Added comparesketchtest: loadBloomFilter, checkBloomFilter (identical and different references)
* Added BloomFilter tests: empty filter, get_pop, store/load roundtrip
* Extended getBaseId tests: multiple dots, no extension, no directory separator
* Total: 15 tests across 4 test suites (up from 9 across 3)

### Dependencies
* Replaced kseq++ (kseq++.hpp, seqio.hpp, config.hpp) with kseq.h from seqtk
* Updated THIRD-PARTY-LICENSES accordingly

## 1.1.0 (2024-06-01)
* First public release version 1.1.0.
