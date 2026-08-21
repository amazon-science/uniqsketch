# UniqSketch: Sensitive and resource-efficient strain/abundance identification in metagenomics data

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/uniqsketch/README.html)

UniqSketch is a sensitive tool for identifying strains and their relative abundances in metagenomics samples. The algorithm consists of two stages: Indexing and Querying.

* **Indexing**: takes a list of target sequences (usually reference genomes) and identifies all k-mers within each target that uniquely belong to that specific sequence. It then checks a set of criteria such as low-complexity on the k-mer candidates to refine the list and provide the final set of unique k-mers for each target. To efficiently keep track of k-mers and their frequencies, UniqSketch utilizes the Bloom filter data structure as an alternative to hash tables.

* **Querying**: raw read sequencing data are queried against the sketch index. To remove k-mers caused by sequencing errors (k-mers with very low depth), UniqSketch utilizes a cascading Bloom filter to discard weak k-mers. The solid k-mers are then queried against the index and the count for the related reference target is updated accordingly. After all solid k-mers in the raw read input data are queried, the reference count table is computed to generate the final output of most abundant references.

## Installation

### conda installation (recommended)
```bash
conda install -c bioconda -c conda-forge uniqsketch
```

### Using CMake
```bash
git clone https://github.com/amazon-science/uniqsketch
cd uniqsketch
cmake -S . -B build
cmake --build build
```
Binaries are placed in `build/bin/`.

### Using Make
```bash
git clone https://github.com/amazon-science/uniqsketch
cd uniqsketch
make
```
Binaries are placed in `bin/`.

### Verifying installation
```bash
uniqsketch --help
```

### Running Tests
```bash
# CMake
cmake --build build --target test

# Make
make test
```

## Quick Start

UniqSketch consists of three utilities:
* `uniqsketch`: build a signature index from a reference/target database
* `querysketch`: query input read files against a constructed index
* `comparesketch`: compute pairwise similarity between two sets of reference genomes

1. Construct a sketch index of size `81bp` from three references and store it in `sketch_index.tsv`:
```bash
uniqsketch --sensitive -k81 -o sketch_index.tsv ref_1.fa ref_2.fa ref_3.fa
```

2. Cluster nearly identical references and build a sketch from representatives only:
```bash
uniqsketch --cluster=2000 --sensitive -k81 -o sketch_index.tsv @refs.txt
```
References with fewer than 2000 unique k-mers between them are merged into clusters. A `clusters.tsv` file is written to the output directory mapping each reference to its cluster representative.

3. Query a metagenomics sample with read files `r1.fq.gz` and `r2.fq.gz`:
```bash
querysketch --r1 r1.fq.gz --r2 r2.fq.gz --ref sketch_index.tsv --out out_sample.tsv
```

## Output

### uniqsketch
* `sketch_uniq.tsv`: primary output — tab-separated file storing the final set of signatures for each reference.
* `outdir_uniqsketch/`: folder containing a file per reference with all signature candidates.
* `db_uniq_count.tsv`: tab-separated summary of all references and total number of candidate signatures.
* `clusters.tsv` (when `--cluster` is used): tab-separated file mapping each reference to its cluster representative, with columns: `cluster_id`, `representative`, `member`, `num_kmers`.

### querysketch
* `out_sample.tsv`: primary output reporting references and their abundance within the sample:

```
ref         abundance   count
ref_1       0.7862      3033
ref_2       0.1993      769
ref_3       0.01426     55
```

* `log_out_sample.tsv`: detailed log showing for each reference the total number of assigned reads and all `signature:count` pairs.
* `logread_out_sample.tsv`: read-level log listing all matched read IDs per reference.

## Usage

### uniqsketch
```
Usage: uniqsketch [OPTION] @LIST_FILES (or FILES)
Creates unique sketch from a list of fasta reference files.
A list of files containing file names in each row can be passed with @ prefix.

 Options:

  -t, --threads=N       use N parallel threads [1]
  -k, --kmer=N          the length of kmer [81]
  -b, --bit=N           use N bits per element in Bloom filter [128]
  -d, --hash1=N         distinct Bloom filter hash number [5]
  -s, --hash2=N         repeat Bloom filter hash number [5]
  -c, --cov=N           number of unique k-mers to represent a reference [100]
  -f, --outdir=STRING   dir for universe unique k-mers [outdir_uniqsketch]
  -o, --out=STRING      the output sketch file name [sketch_uniq.tsv]
  -r, --stat=STRING     the output unique kmer stat file name [db_uniq_count.tsv]
  -e, --entropy         sets the aggregate entropy rate threshold [0.65]
      --cluster=N       cluster similar references with N unique k-mer threshold [0=off]
      --min-margin=N    min Hamming distance of signatures to other references (1=off, 2) [1]
      --max-homopolymer=N  drop signatures with a homopolymer run > N bp (0=off) [0]
      --strict-margin   enforce --min-margin as a hard requirement (no margin-1 fallback)
      --sensitive       sets sensitivity parameter c to 100
      --very-sensitive  sets sensitivity parameter c to 1000
      --help            display this help and exit
      --version         output version information and exit
```

With `--min-margin=2`, `uniqsketch` prefers signatures that lie at least 2 bases
(Hamming distance) away from every k-mer in the other references, so a single
sequencing substitution cannot turn a foreign k-mer into one of a reference's
signatures (a source of false-positive hits). It is a soft preference with
fallback — when a genomic region has no such signature, a normal one is used —
so per-reference signature counts are preserved. The default `--min-margin=1`
leaves selection unchanged. Pass `--strict-margin` to remove that fallback: slots
with no 2-safe candidate are left empty, so every signature meets the requested
margin, at the cost of a smaller sketch.

`--max-homopolymer=N` drops signature candidates containing a single-base run
longer than N bp. The entropy filter scores a signature across the whole k-mer, so
a short local run (a poly-G tract, say) gets diluted and can slip through;
this is a direct gate on such indel-prone tracts. Rejected candidates are replaced
from the same genomic slot where possible, so signature counts and spacing are
largely preserved. The default `--max-homopolymer=0` disables the filter.

### querysketch
```
Usage: querysketch [OPTIONS] [ARGS]
Identify references and their abundance in FILE(S).

Acceptable file formats: fastq in compressed formats gz, bz, zip, xz.

 Options:

  -t, --threads=N       use N parallel threads [1]
  -b, --bit=N           use N bits per element in Bloom filter [64]
  -o, --out=STRING      the output file name
  -l, --r1=STRING       input read 1
  -r, --r2=STRING       input read 2
  -g, --ref=STRING      input uniqsketch reference
  -h, --hit=N           number of uniqsketch hits to call a reference [10]
  -a, --acutoff=N       abundance cutoff to report [0.0]
  -s, --rcutoff=N       read cutoff to report [2]
      --sensitive       sets sensitivity parameter h=10
      --very-sensitive  sets sensitivity parameter h=5
      --solid           only use solid k-mers in reads
      --help            display this help and exit
      --version         version information and exit
```

### comparesketch
```
Usage: comparesketch [OPTION] LIST1 LIST2
Compare two sets of fasta reference files |LIST1|*|LIST2|.
Two lists of files containing file paths in each row.

 Options:

  -t, --threads=N       use N parallel threads [1]
  -k, --kmer=N          the length of kmer [81]
  -b, --bit=N           use N bits per element in Bloom filter [16]
  -d, --hash=N          Bloom filter hash number [3]
  -g, --gsize=N         approximate size for reference sequence [5000000]
  -o, --out=STRING      the output similarity file name [reference_similarity.tsv]
      --auto-gsize      size the Bloom filter from the largest input genome
      --auto            size the Bloom filter from an ntCard cardinality estimate
      --fpr=F           target Bloom-filter false-positive rate (sets bits; e.g. 0.001)
      --low-mem         low-memory streaming mode for large reference sets
      --help            display this help and exit
      --version         output version information and exit
```

For large reference sets (e.g. thousands of multi-megabase genomes), the default
`comparesketch` caches every reference's k-mer hashes in memory, which can require
more RAM than is available. Pass `--low-mem` to stream the comparison instead: it
holds only one Bloom filter and one reference at a time per thread (memory scales
with the thread count and genome size rather than the whole set), at the cost of
re-reading the second list once per first-list reference. Results are identical.

## Dependencies

UniqSketch vendors the following libraries:
* [ntHash](https://github.com/bcgsc/ntHash) — recursive nucleotide hashing
* [ntCard](https://github.com/bcgsc/ntCard) — streaming k-mer cardinality estimation
* [kseq.h](https://github.com/lh3/seqtk) — fast FASTA/FASTQ parsing (Heng Li)

System requirements:
* C++17 compiler (GCC 7+, Clang 5+)
* zlib
* OpenMP (optional, for multi-threading)
