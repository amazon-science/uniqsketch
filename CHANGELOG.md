# Changelog

All notable changes to this project will be documented in this file.

## 2.0.0 (unreleased)

### Features
* `uniqsketch`: added `--tree=T1,T2,...` to build a hierarchical (multi-resolution) index instead of a single flat sketch. References are clustered at each comma-separated unique-k-mer threshold — thresholds are sorted descending so the first is the coarsest level — producing a nested tree whose internal nodes are cluster representatives and whose leaves are individual references. The tree structure is written to `<outdir>/tree_index.tsv`, a coarse sketch per internal level to `<outdir>/sketch_level_<L>.tsv`, and the full-resolution leaf sketch to the usual `-o` output. Coarse levels summarize groups of near-identical references; fine levels resolve individual strains.
* `querysketch`: added `--tree=DIR` to query a hierarchical index. Pass the `uniqsketch --tree` output directory (containing `tree_index.tsv` and the `sketch_level_*.tsv` files) as `--tree` and the leaf sketch as `--ref`. The query narrows candidates coarse-to-fine: the coarsest level is queried first, and only references whose cluster representative clears the hit/read thresholds are carried to the next finer level, so the search focuses on promising sub-trees rather than scoring every reference at full resolution.
* `uniqsketch`: `--cluster` now uses a near-linear MinHash strategy (bottom-s sketch, inverted-index candidate generation, and representative-linkage union-find) in place of the previous all-pairs O(n^2 × genome) unique-k-mer computation, making clustering of large reference sets (thousands of genomes) practical. The new `--exact` flag restores the original exact O(n^2) clustering for small sets or when exact unique-k-mer counts are required.

* `querysketch --tree`: reports at the finest level it can resolve. When a coarse-level cluster representative clears the reporting thresholds but none of its member leaves can be told apart, the result now contains a `family` row for that group, named by its node in `tree_index.tsv`, instead of nothing. The coarse level is far more sensitive for near-identical strains — a representative competes only with distant genomes, so its signature pool is orders of magnitude larger than any member's (94,406 vs 151 candidates on the 8-genome set) — and previously that evidence was discarded. Coarse-level pruning is also now permissive (`--tree-prune=N`, default 1: keep any sub-tree with at least N hits and one read), because a sub-tree pruned at a coarse level could never be recovered, which made the tree strictly less sensitive than a flat index. Measured at the 0.01% detection floor over 56 mixtures: leaf-level detections 15 → 19 (from permissive pruning) plus 6 family-level calls, so samples with signal rose from 15/56 to 25/56. Output on resolvable samples is unchanged.
* `querysketch`: the primary output gains four appended columns — `reads` (distinct supporting reads), `sigs_hit` / `sigs_total` (how many of the reference's signatures were observed), and `resolution` (`strain` or `family`). The original three columns are unchanged, so existing parsers keep working. These expose the evidence behind each call, which the abundance figure alone hides: two references can share an abundance while one rests on thousands of reads and the other on three. `reads` is per-call evidence and is not comparable across references — a read overlapping several packed signatures counts once in `reads` but several times in `count`, and near-identical strains have their few unique k-mers concentrated in short regions, so at equal depth their `reads` can differ 20-fold from a distinct genome's while `count` agrees. This is exactly why `count`, not `reads`, remains the basis for abundance (a read-based estimator was measured at 100× worse error); `count / reads` is the per-read multiplicity and a useful flag for calls carried by a handful of reads on a signature cluster.

### Changes
* `comparesketch`: unified the memory controls into a single `-m, --mem=N` budget in GB for reference caching (default 64). The tool caches as many second-list references as fit within the budget and streams the remainder, so peak memory stays bounded regardless of reference-set size. `--low-mem` is now an alias for `--mem 0` (fully streaming). Similarity output is unchanged.
* `uniqsketch` / `querysketch`: harmonized Bloom-filter sizing between the flat and tree code paths. The flat paths now use the same `resolveBloomSizes` fallback as the tree paths, so low-coverage inputs that previously aborted on the flat path (but succeeded on the tree path) are handled consistently.

### Fixes
* Fixed a data race in `BloomFilter::insert_make_change` that made index construction non-deterministic and could turn a shared k-mer into a signature for two references at once. The individual bit writes were atomic, but the composite "was this k-mer newly inserted?" answer was not: two threads inserting the same k-mer could each win a subset of the hash bits and both report success. `uniqsketch` builds its distinct/solid filter pair in parallel and relies on exactly one caller observing success, so that a k-mer occurring in two references is recorded as a repeat and excluded from both. When the race hid the repeat, the k-mer stayed a signature for both references and reads from one inflated the abundance of the other. The test-then-insert sequence is now serialised on a lock striped by the leading hash, so inserts of the same k-mer are ordered while unrelated k-mers stay parallel; the lock is only taken when the k-mer looks absent, so cost is per distinct k-mer and index build time is unchanged. Index builds are now byte-identical across runs, and `querysketch --solid` and signature selection are likewise deterministic. Measured on an 8-genome set with near-identical members: before the fix, repeated builds differed and a spurious reference was reported at up to 15% abundance in 2-3 of every 6 builds; after it, six consecutive builds are identical and the spurious call is absent in all configurations.

### Internal
* Replaced the mutable `namespace opt` globals (42 variables defined at file scope across the three `*Util` headers, 267 usages) with one config struct per tool — `UniqSketchConfig`, `QuerySketchConfig`, `CompareSketchConfig` — passed explicitly to every function that needs settings. The globals meant no two `*Util` headers could be included in one translation unit (each defined the same names), no function could be tested under two settings in one process, and `buildSketch` mutated the entropy ceiling as a side effect, which forced the tree path to snapshot and restore nine globals around every level. All three problems are gone: the headers compose, tests build their own configs side by side, the ceiling is a const method, and the tree path copies the config per level. Behaviour is unchanged — indexes, sketches and query outputs rebuilt with the refactored binaries are byte-identical to the previous ones. Added `headerstest`, which compiles every `*Util` header into one translation unit so the fix cannot silently regress (it produces 18 errors against the previous headers).
* `querysketch`: the query hot loop no longer rebuilds an 81-character string at every k-mer position. The sketch is viewed through a map keyed by the signature's canonical ntHash — the value the read iterator already computes — so a position is tested with one 64-bit probe and no allocation, with the stored sequence compared on a probe hit so a hash collision can never produce a false match. Hit tallies are accumulated per thread and merged per batch, replacing a global critical section taken on every hit and an `omp atomic` applied to a map subscript (which only covered the final increment and left the lookup unsynchronised). Results are byte-identical; the per-position cost fell 2.8× (32.3 s → 11.4 s single-threaded on a 3M-pair sample). At high thread counts the gain is smaller because the single-threaded reader is now the bottleneck.
* Added a CI workflow that builds all three tools and runs the test suite on every push and pull request, across Make and CMake with GCC and Clang. Previously the only automated check was CodeQL, so the tests never ran on a PR; the Clang legs are the first time that compiler has built the codebase.

### Tests
* Added an end-to-end test (`e2etest`) that builds an index from synthetic genomes and queries synthetic reads against it in-process, then checks the result against known ground truth — the first automated protection for the abundance formula the paper describes. Fixtures (`test/e2e/`, generated by `make_fixtures.py` from a fixed seed) are two unrelated 20 kb genomes, a near-identical variant of one (12 substitutions, so its signatures are packed into short regions exactly as real near-identical strains' are), and a relative of the other that is sequenced but deliberately not indexed. The test asserts abundances against the mixture, that the unindexed organism produces no row, that its reads are absorbed by its indexed relative (signatures are unique against the index, not the world — a documented contract now pinned by a test), that the evidence columns are internally consistent, and that two index builds are byte-identical. Mutation-checked: it fails against a broken normalisation and against a build that stops recording shared k-mers (which drives a 0.67 share to 0.42); it does not, and is documented not to, catch data races, which need the dedicated concurrency test.
* `uniqsketch`: the index file is now byte-reproducible across builds. `buildSketch` wrote each reference's signatures from inside a parallel loop, so the block order followed thread completion; two builds of the same index produced the same signatures in a different line order. Each reference now writes to its own buffer and the buffers are emitted in reference order. Found by the end-to-end test's determinism check; every earlier determinism check had sorted the file first and so could not see it.
* Added a concurrency regression test for `insert_make_change`: sixteen threads released from a spin barrier insert the same k-mer and exactly one must observe success. Verified to fail against the previous implementation and pass against the fixed one.
* Added `treesketchtest`: tree TSV round-trip (`loadTree`, including leaf rows with an empty children column), per-level representative extraction (`getTreeLevelRefs`), and `buildHierarchicalTree` clustering (near-identical references grouped, distinct references kept apart). The suite now has 5 test executables.

### Inherited from 1.6.0
* The `uniqsketch --max-homopolymer=N` and `--strict-margin` options released in 1.6.0 apply to the tree index paths as well: both act inside the shared signature-selection stage, so every tree level is filtered on the same terms as a flat sketch.

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
