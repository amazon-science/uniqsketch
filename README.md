# UniqSketch

**Sensitive, resource-efficient strain-level detection and abundance estimation in metagenomes.**

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/uniqsketch/README.html)
[![conda version](https://img.shields.io/conda/vn/bioconda/uniqsketch.svg)](https://anaconda.org/bioconda/uniqsketch)
[![conda downloads](https://img.shields.io/conda/dn/bioconda/uniqsketch.svg)](https://anaconda.org/bioconda/uniqsketch)
[![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Give UniqSketch a set of reference genomes and a sequencing sample; it reports which
references are present and at what relative abundance. It is built for the case where
the candidates are closely related — strains of the same species that differ across a
small fraction of their sequence.

## How it works

Strain-level assignment is hard because near-identical genomes share almost all of
their sequence. Most reads are therefore consistent with several references at once,
and deciding which one they came from is guesswork.

UniqSketch avoids that ambiguity rather than trying to resolve it. It considers only
**signatures**: k-mers that occur in exactly one reference across the entire database.
A read carrying a signature is unambiguous evidence for that single reference, so
counting signature hits gives a direct read-out of which references are present. Two
choices make this work in practice — a long k (81 bp by default), which makes a
coincidental match between unrelated genomes vanishingly unlikely, and Bloom filters
in place of hash tables, which track k-mer frequencies in a fraction of the space
storing every k-mer would require.

The workflow has two stages:

* **Indexing** (`uniqsketch`) scans the reference set, finds the k-mers unique to each
  reference, and filters the candidates for qualities that make a signature
  trustworthy — discarding low-complexity sequence, and optionally long homopolymer
  tracts (`--max-homopolymer`) or k-mers lying only one substitution away from another
  reference (`--min-margin`). Surviving candidates are sampled across the genome so
  coverage is even, producing the sketch index.

* **Querying** (`querysketch`) streams the reads and matches their k-mers against the
  index. A cascading Bloom filter can first drop k-mers seen only once (`--solid`),
  which removes most sequencing-error artefacts before they reach the index. Hits are
  tallied per reference, thresholded to suppress spurious calls, and reported as a
  ranked table of references with relative abundances.

A third utility, `comparesketch`, measures pairwise similarity within a reference set
— useful for deciding which references are too close to separate and should be
clustered before indexing.

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

Reading this table:

* `count` is the number of signature k-mer hits observed for that reference.
* `abundance` is `count` divided by the summed `count` of all **reported** references.
  Rows are sorted by descending `count`.

Two consequences worth knowing. First, `abundance` is a *relative* figure: it sums
to 1.0 across the rows in the file, so it answers "of the material I could identify,
what fraction is this reference?" — not "what fraction of the sample is this
reference?" Reads from organisms absent from your index, and reads filtered out by
the thresholds below, are not in the denominator. Second, references failing any
reporting threshold are omitted entirely rather than listed with a zero, so an
absent row means "not called", which is not the same as "not present".

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
longer than N bp — a run of exactly N is kept, N+1 is rejected.

Why this is worth doing: long homopolymer tracts are the least reliable part of a
read. Sequencers lose track of run length in them, so the same genomic tract is
reported as different lengths in different reads. An insertion or deletion inside a
signature shifts every base after it, and the k-mer no longer matches, which costs
you real hits on a reference that is genuinely present. The existing entropy filter
does not catch this, because it scores complexity across the whole k-mer: at k=81,
a 7 bp poly-G tract is diluted by 74 well-behaved bases and comfortably passes.
`--max-homopolymer` is a direct gate on the local run length instead.

**A good starting point is `--max-homopolymer=6`:**

```bash
uniqsketch --sensitive -k81 --max-homopolymer=6 -o sketch_index.tsv @refs.txt
```

The reasoning behind 6 is that indel error rates stay low through short runs and
climb once a tract reaches roughly 7–8 bp, while signatures containing runs that
long are rare — so you remove the unreliable tail cheaply. In an 8-genome test set,
about 92% of selected signatures had a longest run of 5 bp or less and only ~4%
reached 7 bp or more; filtering at 6 cost a few percent of signatures. Your own
distribution will differ with genome composition and `-k`, so treat 6 as a default
to adjust rather than a universal constant.

Choosing a value:

- **`0` (default)** — filter off, behaviour identical to earlier versions.
- **`6`** — recommended; removes the indel-prone tail at little cost.
- **`4`–`5`** — stricter. Viable if your reads are homopolymer-error-prone, but it
  starts cutting into signature yield noticeably.
- **below `4`** — not advised. Runs of 3–4 bp are extremely common (they were ~72%
  of signatures in the test set above), so you would discard most candidates and
  weaken detection sensitivity.

Rejected candidates are replaced from the same genomic slot wherever another
candidate exists, so the net loss in signature count is much smaller than the
number rejected, and genomic spacing is preserved. If a reference is left with too
few signatures, raise N or turn the filter off for that run. The filter applies
wherever signatures are selected, so it covers clustered (`--cluster`) builds on the
same terms.

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

Reads may be paired or single-end. Pass both `--r1` and `--r2` for paired data; pass
`--r1` alone for single-end. `--ref` takes the sketch file produced by `uniqsketch`
(the `-o` output), and `k` is read from that file rather than supplied again.

**How a reference gets called.** A reference is reported only if it clears all three
thresholds:

| flag | default | meaning |
|---|---|---|
| `-h, --hit=N` | 10 | needs **at least** N signature k-mer hits |
| `-s, --rcutoff=N` | 2 | needs **more than** N distinct matching reads |
| `-a, --acutoff=F` | 0.0 | needs relative abundance **greater than** F |

The first two are the ones that matter in practice. `--hit` guards against a
reference being called on a handful of stray k-mer matches, and `--rcutoff` guards
against many hits that all came from a single read — which is what a chimeric or
low-quality read looks like. Lower values make detection more sensitive at low
abundance and raise the false-positive rate; higher values do the reverse.
`--sensitive` sets `h=10` (the same as the default) and `--very-sensitive` sets
`h=5`, so reach for `--very-sensitive` when hunting trace-level organisms.

**`--solid`** discards the first occurrence of every k-mer in the sample and queries
only k-mers seen at least twice. A sequencing error creates a k-mer that almost
never recurs, so this removes most error-induced k-mers before they can be matched
against the index. It is worth enabling on real data of reasonable depth. The
trade-off is coverage-dependent: at very low coverage a genuine signature may
legitimately appear only once, so `--solid` can suppress a true low-abundance call.
Prefer it for typical-depth samples, and leave it off when chasing organisms near
the detection floor.

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

## Choosing parameters

The defaults are sensible for bacterial-scale references and typical short-read
samples. The knobs below are the ones worth reaching for, roughly in order of how
often they matter.

**`-k, --kmer` (default 81).** A long k is what makes signatures specific: the
chance of an 81-mer occurring in an unrelated genome by coincidence is negligible,
so a match is strong evidence. The cost is that k must fit inside your reads, and
every sequencing error invalidates the k-mers overlapping it — with k=81, one error
knocks out up to 81 k-mer positions. Keep 81 for 100 bp+ reads. Lower it for shorter
reads, accepting that specificity drops with it. `k` is baked into the index, so the
same value is used automatically at query time.

**Number of signatures per reference (`-c`, `--sensitive`, `--very-sensitive`,
default 100).** This sets how many signatures represent each reference. More
signatures means more chances to hit a low-abundance organism, at the cost of a
larger index and slower queries. `--sensitive` selects 100 — identical to the
default — and `--very-sensitive` selects 1000, which is the setting to use when you
care about trace-level detection. Signatures are spread across the genome rather
than clustered, so a higher count also buys robustness against uneven coverage.

**`--cluster=N`.** When your reference set contains near-identical genomes, their
signatures compete and abundance gets split arbitrarily between them. Clustering
collapses references within N distinguishing k-mers of each other and indexes one
representative, which produces cleaner quantification at the cost of strain-level
resolution within a cluster. Consult `clusters.tsv` to see what was merged.

**`--max-homopolymer` and `--min-margin`.** Both harden signatures against read
errors and are described in detail above. `--max-homopolymer=6` is a reasonable
default to adopt; `--min-margin=2` is worth it when your references include close
relatives and you are seeing cross-assignment.

**`-t, --threads`.** Set this to your core count for both indexing and querying;
both stages parallelize well.

## Troubleshooting

**`Error: insufficient k-mer content for Bloom filter construction`.** The sample or
reference set does not contain enough distinct k-mers to size the Bloom filters. The
usual cause is reads shorter than `k`: with the default k=81, 75 bp reads yield no
k-mers at all. Check your read length and, if it is short, rebuild the index with a
smaller `-k`. Near-empty or truncated input files produce the same error.

**No references reported (output has only the header).** Nothing cleared the
reporting thresholds. Confirm the sample really should contain something in your
index, then retry with `--very-sensitive` (`h=5`) and without `--solid`, which is the
most permissive combination. Inspect `log_out_sample.tsv` — it lists per-signature
counts regardless of whether the reference was reported, so it shows whether there
were hits that merely fell short of the thresholds.

**Abundances do not match expectations across similar strains.** Near-identical
references split each other's signatures. Run `comparesketch` on the reference set to
quantify how similar they are, then use `--cluster=N` to collapse the ones that are
too close to separate reliably.

**A reference gets called that should not be there.** Check `log_out_sample.tsv` for
how many distinct signatures carried the hits. Many hits concentrated on one or two
signatures suggests cross-mapping or a repeat, rather than genuine presence. Raising
`-s, --rcutoff`, enabling `--solid`, or rebuilding with `--min-margin=2` and
`--max-homopolymer=6` all reduce this class of false positive.

## Dependencies

UniqSketch vendors the following libraries:
* [ntHash](https://github.com/bcgsc/ntHash) — recursive nucleotide hashing
* [ntCard](https://github.com/bcgsc/ntCard) — streaming k-mer cardinality estimation
* [kseq.h](https://github.com/lh3/seqtk) — fast FASTA/FASTQ parsing (Heng Li)

System requirements:
* C++17 compiler (GCC 7+, Clang 5+)
* zlib
* OpenMP (optional, for multi-threading)
