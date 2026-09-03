# Changelog

All notable changes to this project will be documented in this file.

## 1.6.1 (2026-08-24)

### Fixes
* Fixed a data race in `BloomFilter::insert_make_change` that made index construction non-deterministic and could cause a k-mer shared by two references to be used as a signature for both. Each hash bit was written atomically, but the composite "was this k-mer newly inserted?" answer was not: two threads inserting the same k-mer could each win a subset of the bits and both report success. `uniqsketch` builds its distinct and solid Bloom filters in parallel across references and relies on exactly one caller observing success, so that a k-mer occurring in more than one reference is recorded as a repeat and excluded from every reference's unique set. When the race hid that repeat, the k-mer remained a signature for both references, and reads from one inflated the reported abundance of the other.

  The impact is limited to reference sets containing near-identical genomes, but there it is significant. On an 8-genome test set with four near-identical members, repeated index builds differed from one another and a reference that was absent from the sample was reported at up to 15% abundance in two to three of every six builds. After the fix, six consecutive builds are byte-identical and no spurious call appears in any configuration. Users who previously worked around this with `--cluster` no longer need to.

  The test-then-insert sequence is now serialised on a lock striped by the leading hash value, so inserts of the same k-mer are ordered while unrelated k-mers proceed in parallel. The lock is only taken when the k-mer appears to be absent, so its cost is per distinct k-mer rather than per occurrence and index build time is unchanged. Signature selection and `querysketch --solid` are deterministic as a result. Single-threaded behaviour is unchanged, and indexes built with earlier versions remain readable — but rebuilding is recommended if your reference set contains near-identical genomes.

### Tests
* Added a concurrency regression test for `insert_make_change`: sixteen threads released from a spin barrier insert the same k-mer, and exactly one must observe success. Verified to fail against the previous implementation and pass against the fixed one. `bloomfiltertest` now links against pthreads (`-pthread` under make, `Threads::Threads` under CMake).

## 1.6.0 (2026-08-21)

### Features
* `uniqsketch`: added `--max-homopolymer=N` to drop signature candidates containing a single-base run longer than N bp. The aggregate entropy filter scores a signature over the whole k-mer, so a short local run (e.g. a poly-G tract) is diluted and can survive; this adds a direct local-run-length gate for such indel-prone tracts. Because the filter is applied while scanning each genomic selection slot, a rejected candidate is replaced by another candidate from the same slot wherever one exists, so signature counts and spacing are largely preserved. The default (`--max-homopolymer=0`) disables the filter and leaves selection unchanged.
* `uniqsketch`: added `--strict-margin` to enforce `--min-margin` as a hard requirement. By default, `--min-margin=2` falls back to a margin-1 candidate when a selection slot contains no 2-safe option; with `--strict-margin` such slots are left empty instead, so every emitted signature satisfies the requested margin at the cost of a smaller sketch. The default (off) leaves `--min-margin` behavior unchanged.

Both options are contributed features and are off by default, so 1.6.0 is backward compatible with 1.5.0 output.

## 1.5.0 (2026-07-17)

### Features
* `uniqsketch`: added `--min-margin=N` to harden signature selection against single-base sequencing errors. With `--min-margin=2`, uniqsketch prefers signatures that lie at Hamming distance ≥ 2 from every k-mer in other references, so a single substitution error in a read cannot turn a foreign k-mer into one of a reference's signatures — a source of false-positive hits. It is a tiered preference: within each genomic selection slot a 2-safe, non-low-complexity candidate is chosen when available, falling back to the previous behavior only when a slot has no 2-safe option, so per-reference signature counts and genomic spacing are preserved. A summary line reports how many selected signatures met the margin. The default (`--min-margin=1`) leaves selection unchanged. Currently supports margins 1 and 2; higher margins are planned.

## 1.4.0 (2026-07-16)

### Features
* `comparesketch`: added `--low-mem`, a streaming M×N comparison mode for large reference sets. The default mode pre-caches every reference's k-mer hashes (M+N file reads) but its memory scales with the entire reference universe — for a few thousand multi-megabase genomes this can exceed 1 TB. With `--low-mem`, only one Bloom filter and one reference's k-mers are held per thread at a time: each first-set reference is loaded into a filter and every second-set reference is re-read from disk and queried against it. Memory is O(threads × genome) instead of O(all genomes), at the cost of M×N (rather than M+N) file reads. The similarity output is identical to the default mode; only the memory/I/O profile changes.

## 1.3.0 (2026-06-11)

### Features
* `comparesketch`: added `--auto` to size the Bloom filter from an ntCard cardinality estimate of the reference set, fixing the load factor (bits per distinct k-mer) so the false-positive rate is the design value regardless of genome size.
* `comparesketch`: added `--auto-gsize` to size the Bloom filter from the largest input genome (`estimateMaxGenomeSize`), removing the need to set `--gsize` manually.
* `comparesketch`: added `--fpr=F` to set bits-per-element for a target Bloom-filter false-positive rate (e.g. `--fpr 0.001` selects 29 bits/element).
* The comparesketch false-positive rate can now be controlled precisely; the default (`--gsize`-based sizing) behavior is unchanged.

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
