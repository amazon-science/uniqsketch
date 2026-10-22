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
References with fewer than 2000 unique k-mers between them are merged into clusters, and only one representative per cluster is indexed. Note that the other members of a cluster then disappear from the output: a sample containing one of them is reported under the representative's name. Use `--tree` (example 3) if you need results at individual-reference resolution. A `clusters.tsv` file is written to the output directory mapping each reference to its cluster representative. Clustering uses a near-linear MinHash strategy by default; add `--exact` for the exact O(n^2) method on small reference sets.

3. Build a hierarchical (multi-resolution) tree index and query it. The full-resolution leaf sketch goes to `-o`; the tree structure and per-level sketches go to the output directory (`-f`):
```bash
uniqsketch --tree=5000,1000 --sensitive -k81 -o sketch_index.tsv -f tree_out @refs.txt
querysketch --r1 r1.fq.gz --r2 r2.fq.gz --ref sketch_index.tsv --tree tree_out --out out_sample.tsv
```
The query narrows candidates from the coarsest level down to individual references. Tree builds select signatures with `--min-margin=2 --strict-margin` by default, and the query then reports single-read calls (`--rcutoff=0`); see [querysketch](#querysketch-1) for why. For trace-level (0.01%) detection use `-c500` or `--very-sensitive` instead of `--sensitive`.

4. Query a metagenomics sample against a flat index with read files `r1.fq.gz` and `r2.fq.gz`:
```bash
querysketch --r1 r1.fq.gz --r2 r2.fq.gz --ref sketch_index.tsv --out out_sample.tsv
```

## Output

### uniqsketch
* `sketch_uniq.tsv`: primary output — tab-separated file storing the final set of signatures for each reference.
* `outdir_uniqsketch/`: folder containing a file per reference (`<name>.cand`) with all signature candidates: the positions of the reference's unique k-mers, per contig, in a compact binary form (`UQSKCAND` header, then per contig its name, length and a `uint32` position list). Scratch for the build, about 4 bytes per candidate; a 120-genome tree index writes under 1 GB of it.
* `db_uniq_count.tsv`: tab-separated summary of all references and total number of candidate signatures.
* `clusters.tsv` (when `--cluster` is used): tab-separated file mapping each reference to its cluster representative, with columns: `cluster_id`, `representative`, `member`, `num_kmers`.
* With `--tree`, inside the `-f` output directory: `tree_index.tsv` (the hierarchy), `sketch_level_<L>.tsv` (one coarse sketch per level), `build_params.tsv` (how signatures were selected; read by `querysketch --tree`), `sketch_group.tsv` (within-group leaf signatures: `seq`, `refId`, `refName`, `numRef`, `contig`, `pos`, `carriers` as comma-separated refIds or `-`, `node`; absent with `--no-group-signatures`), and `groups/g<i>/` with each group's candidate files.

### querysketch
* `out_sample.tsv`: primary output reporting references and their abundance within the sample:

```
ref         abundance   count   reads   sigs_hit   sigs_total   resolution
ref_1       0.7862      3033    2871    100        100          strain
ref_2       0.1993      769     751     98         100          strain
ref_3       0.01426     55      7       27         100          strain
```

Reading this table:

* `count` is the number of signature k-mer hits observed for that reference.
* `abundance` is the reference's hits per signature (`count / sigs_total`) divided by
  the sum of that quantity over all **reported** references. Hits per signature is
  proportional to coverage; the raw `count` is not once references carry different
  numbers of signatures, which happens whenever a reference's whole candidate pool
  is smaller than `-c` (near-identical strains routinely have only a few hundred
  unique k-mers). When every reference has the same number of signatures the two
  definitions agree exactly, so results on such indexes are unchanged. Rows are
  sorted by descending `count`.
* `reads` is the number of distinct reads that contributed at least one hit.
* `sigs_hit` / `sigs_total` is how many of the reference's signatures were observed
  at least once, out of how many it has.
* `resolution` is `strain` for an individually resolved reference. With a `--tree`
  index it may also be `family` — see [querysketch](#querysketch-1) below.

The first three columns are unchanged from earlier versions; the rest are appended,
so existing parsers keep working.

Two consequences worth knowing. First, `abundance` is a *relative* figure: it sums
to 1.0 across the rows in the file, so it answers "of the material I could identify,
what fraction is this reference?" — not "what fraction of the sample is this
reference?" Reads from organisms absent from your index, and reads filtered out by
the thresholds below, are not in the denominator. Second, references failing any
reporting threshold are omitted entirely rather than listed with a zero, so an
absent row means "not called", which is not the same as "not present".

**Judging a call.** `abundance` says how much; the remaining columns say how much to
trust it, which the abundance figure alone hides — two references can share an
abundance while one rests on thousands of reads and the other on three. Read them
together:

* **`reads` is per-call evidence. Do not compare it across references.** A read
  overlapping several signatures at once counts once in `reads` but several times in
  `count`, and how often that happens depends on how tightly a reference's signatures
  are packed along its genome. Near-identical strains have their few unique k-mers
  concentrated in the short regions where they differ, so a single read can carry 20
  of their signatures, while a distinct genome's signatures are spread out and a read
  carries about one. At the same depth these two references show similar `count` but
  `reads` figures 20-fold apart. `count` is the right basis for abundance precisely
  because this packing effect cancels in expectation; `reads` does not, which is why
  it is reported as evidence rather than used for abundance.
* **`count / reads` is the per-read multiplicity.** A high ratio means a clustered
  signature region, so one read there contributes many hits. That is what a genuine
  read from a near-identical strain looks like, and also what an error-carrying read
  from an abundant relative looks like; which one you are seeing depends on the
  index. On a strict-margin index (the `--tree` default) a foreign read needs two
  coincident errors to hit a signature, so `count=60, reads=1` is ordinarily a real
  read on a divergent locus. On a margin-1 index a single error suffices, and at
  high coverage of a relative such reads are common, so the same row deserves
  suspicion. `count=60, reads=55` is independent evidence either way.
* **`sigs_hit` low relative to `sigs_total`, with a high `count`,** means many hits on
  few signatures — the pattern of a repeat or a cross-mapping rather than genuine
  presence. Support spread across most of a reference's signatures is what a real
  detection looks like.
* Every call has already cleared `--hit` (on `count`) and `--rcutoff` (on `reads`).
  These columns let you see how far above those floors a call sits: `reads=3` is one
  read above the flat-index cutoff; `reads=1` is the floor on a strict-margin tree
  index; `reads=30000` is unimpeachable.

* `log_out_sample.tsv`: detailed log showing for each reference the total number of assigned reads and all `signature:count` pairs.
* `logread_out_sample.tsv`: read-level log listing all matched read IDs per reference.

## Usage

### uniqsketch
```
Usage: uniqsketch [OPTION] @LIST_FILES (or FILES)
Creates unique sketch from a list of fasta reference files (plain or gzip).
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
      --tree=T1,T2,..   build hierarchical tree index with comma-separated thresholds
      --exact           use exact O(n^2) clustering instead of MinHash (with --cluster)
      --min-margin=N    min Hamming distance of signatures to other references (1=off, 2) [1; 2 with --tree]
      --max-homopolymer=N  drop signatures with a homopolymer run > N bp (0=off) [0]
      --strict-margin   enforce --min-margin as a hard requirement (no margin-1 fallback) [on with --tree]
      --no-group-signatures  do not build within-group leaf signatures (with --tree)
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

`--tree` builds default to `--min-margin=2 --strict-margin` (pass `--min-margin=1`
to opt out); flat builds keep `--min-margin=1`. The reason is measured rather than
theoretical. A tree index exists to separate close relatives, and a margin-1
signature of one relative is hit by every read of an abundant sibling that carries a
sequencing error at the distinguishing base. With the sibling at ~2700× coverage and
a 0.3% per-base error rate, each such site yields a few error reads, which is more
than the entire contribution of a member present at 0.01%: its ~300 reads leave at
most a handful on its distinguishing loci. Across 28 such mixtures on the 8-genome
test set, margin-1 signatures collected 60 error reads on the true low-abundance
member against 31 genuine ones, and 2093 reads on members that were not in the
sample at all, producing 42 false strain calls. Requiring two substitutions cut the
reads on absent members from 2093 to 2 and kept 25 of the 31 genuine reads, while
the close members retained 62–99% of their candidate pools; the tree's coarse
levels, whose representatives compete only with distant genomes, lose almost
nothing. The selection parameters are recorded in `<outdir>/build_params.tsv`, which
`querysketch --tree` reads to choose its read cutoff (see
[querysketch](#querysketch-1)).

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

With `--tree=T1,T2,...`, `uniqsketch` builds a hierarchical index instead of a flat
sketch: references are clustered at each threshold (sorted descending, so the first
is the coarsest level) into a nested tree. The tree is written to
`<outdir>/tree_index.tsv`, a coarse sketch per level to
`<outdir>/sketch_level_<L>.tsv`, the selection parameters (`k`, `-c`, margin,
homopolymer filter, thresholds) to `<outdir>/build_params.tsv`, and the
full-resolution leaf sketch to `-o`. Query it with `querysketch --tree=DIR`
(pointing at that output directory) together with `--ref` set to the leaf sketch;
the query narrows candidates from the coarsest level down to individual references.

At every level, each reference is represented by the cluster that spans that level
on its path to the root, or by itself if it never clustered, so a level's sketch
covers every reference in the index. A cluster's signatures are k-mers that occur
once in its representative and in no genome outside the cluster (members of other
clusters included), which keeps a genome hidden inside one cluster from feeding
another cluster's family call. The single-copy requirement is the one the leaf
level applies too: a k-mer repeated inside the representative sits in a mobile
element, so it is liable to turn up in an organism missing from the index and
collect that organism's reads for the cluster. Levels whose clustering is
identical to the next finer level are not written; the query walks whatever level
files exist.

**Within-group signatures.** A tree build also writes `<outdir>/sketch_group.tsv`,
a second signature set for the leaves (turn it off with `--no-group-signatures`).
A leaf's ordinary signatures are its k-mers that occur in no other genome of the
index. For a strain with close relatives that pool is small, and mostly for a
reason that has nothing to do with the relatives: most of what the strain does
not share with its siblings it does share with some unrelated genome elsewhere in
the index, so global uniqueness discards it. On 120 *E. coli* genomes the k-mers
of a strain absent from its siblings but present somewhere else outnumbered its
globally unique k-mers by a median of 3.2×, and by 26× for the strain with the
smallest global pool. Those k-mers still separate the strain from its siblings
whenever the genomes that carry them are not in the sample. So for every leaf,
`uniqsketch` also selects k-mers unique *within its group* (the leaves under its
parent node), with the same margin and filters as ordinary signatures, and
records for each one its *carriers*: the other genomes of the index that contain
it. `querysketch --tree` uses a within-group signature only when none of its
carriers is a candidate for the sample (see [querysketch](#querysketch-1)); the
globally unique signatures stay in use throughout, so the extra set can add
evidence but never remove any. Building it costs about 20 % more time on the
120-genome *E. coli* index (153 s against about 125 s at 32 threads; the group
builds run as one parallel job list, so many small groups keep the machine busy)
and about 150 MB more candidate scratch under `<outdir>/groups/`. Query time is
unchanged.

The distinction from `--cluster` matters: clustering is used here only to *guide*
the search, and the final scoring runs against the leaf sketch, so **the output
names individual references rather than cluster representatives.** That is the
reason to prefer `--tree` over `--cluster` whenever strain identity matters — you
get the narrowing without giving up resolution. See
[Choosing parameters](#choosing-parameters) for a direct comparison of the three
index layouts.

By default `--cluster` uses a near-linear MinHash method; `--exact` selects the
exact O(n^2) method for small reference sets.

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
  -s, --rcutoff=N       read cutoff to report [2; 0 with --tree on a strict-margin index]
      --sensitive       sets sensitivity parameter h=10
      --very-sensitive  sets sensitivity parameter h=5
      --solid           only use solid k-mers in reads
      --tree=DIR        hierarchical query using a --tree index directory
      --tree-prune=N    min coarse-level hits to keep a sub-tree in play (with --tree) [1]
      --no-group-signatures  leaf level uses globally unique signatures only (with --tree)
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
| `-s, --rcutoff=N` | 2 (0 with `--tree` on a strict-margin index) | needs **more than** N distinct matching reads |
| `-a, --acutoff=F` | 0.0 | needs relative abundance **greater than** F |

The first two are the ones that matter in practice. `--hit` guards against a
reference being called on a handful of stray k-mer matches, and `--rcutoff` guards
against many hits that all came from a single read — which is what a chimeric or
low-quality read looks like. Lower values make detection more sensitive at low
abundance and raise the false-positive rate; higher values do the reverse.
`--sensitive` sets `h=10` (the same as the default) and `--very-sensitive` sets
`h=5`, so reach for `--very-sensitive` when hunting trace-level organisms.

**Why the read cutoff depends on the index.** Whether one read is trustworthy
depends on how the signatures were selected. On a strict-margin index every
signature is two or more substitutions from every k-mer of every other reference,
so a read from another organism needs two coincident sequencing errors to hit one.
A single read that does hit is then ordinarily a real read, and at the detection
floor it is often all the evidence a sample holds: a member at 0.01% in 3 million
read pairs contributes ~300 reads, of which 0–4 land on its distinguishing loci.
`querysketch --tree` therefore reads the index's `build_params.tsv` and, when the
index was built with `--min-margin=2 --strict-margin` (the `--tree` default),
lowers the read cutoff to 0 so single-read calls are reported. An explicit
`--rcutoff` always wins, and indexes built without a strict margin or by an earlier
version keep the cutoff at 2, with a note on stderr saying so. On the 8-genome test
set this took strain-level detection of close members at 0.01% from 3 of 28
mixtures to 16 of 28, which is every mixture whose reads carried any evidence at
all, at the price of 2 single-read false calls in 280 mixtures (both below 0.1%
abundance, both with `count` ≤ 20 where genuine single-read calls carried 24–71).
Flat indexes keep `--rcutoff=2`: they are usually built without a margin, and at a
sibling's coverage of thousands, single-error reads on margin-1 signatures are
routine.

**`--solid`** discards the first occurrence of every k-mer in the sample and queries
only k-mers seen at least twice. A sequencing error creates a k-mer that almost
never recurs, so this removes most error-induced k-mers before they can be matched
against the index. It is worth enabling on real data of reasonable depth. The
trade-off is coverage-dependent: at very low coverage a genuine signature may
legitimately appear only once, so `--solid` can suppress a true low-abundance call.
Prefer it for typical-depth samples, and leave it off when chasing organisms near
the detection floor.

**Tree queries and family-level calls.** With `--tree=DIR`, the query walks the
hierarchy from the coarsest level down and reports at the finest level it can
resolve. Each row's `resolution` column says which:

* `strain` — an individual reference, resolved at the leaf level. Identical in
  meaning to a flat-index call.
* `family` — a group of near-identical references that the coarse level detected
  confidently but none of whose members could be told apart at the leaf level. The
  `ref` column carries the group's node name from the index's `tree_index.tsv`
  (e.g. `cluster_L1_AWS_Ec377`), where that file lists the members it stands for.

Family rows exist because the coarse level is far more sensitive than the leaf
level for near-identical strains. A cluster representative competes only with
distant genomes, so its signatures are drawn from a pool orders of magnitude larger
than any member's (94,406 vs 151 candidates in the 8-genome test set), and it stays
detectable at abundances where no individual member has enough unique sequence to
be separated. Rather than report nothing, the query reports the family. A family
row is a reliable presence flag, and its abundance is that of the group as a whole,
since the representative's signatures are shared by every member; read it as "one
or more of these is here, at about this level". On the 8-genome test set a close
member at 0.01% that could not be resolved was reported as its family at
0.009–0.012%; on 120 *E. coli* genomes, a member of a 4-strain family at 0.01% was
reported as that family in 13 of the 28 mixtures where it could not be resolved,
at 0.006–0.013%. A family row is never emitted when any of its members is reported
as a strain, a group is reported once (at the finest level that detected it), and
a group is not reported when one of its sub-groups is.

`--tree-prune=N` controls how eagerly sub-trees are dropped during the walk. A
sub-tree stays in play if its representative gathers at least N hits and one read
at the coarse level; the default of 1 is deliberately permissive, because a sub-tree
pruned here can never be recovered at a finer level. Raise it to speed up queries
against very large reference sets at some cost in sensitivity at the detection floor.

**Within-group signatures at the leaf level.** If the index has a
`sketch_group.tsv` (every `uniqsketch --tree` build writes one unless told not to;
see [uniqsketch](#uniqsketch-1)), the leaf level also scores each candidate leaf's
within-group signatures: k-mers unique among the leaves of its group but present in
some other genome of the index, each recorded with those *carriers*. A signature is
used only if none of its carriers is itself a candidate after the coarse levels, so
a read from a carrier that is present can never be credited to the leaf; the cost of
a carrier being present is a smaller supplementary set, never a wrong hit. The
candidate set at that point is the permissive one the coarse levels leave (any
sub-tree with one read), so a carrier with any trace of evidence blocks its
signatures. What remains is a genome present below that floor, whose few reads at
the shared loci then count for the leaf; on a strict-margin index that needs a
genome too rare for a single read to reach one of its group's several thousand
coarse signatures. The globally unique signatures are always scored, so this adds
evidence and cannot take any away. On the 120-genome *E. coli* set it raised
strain-level resolution of close members at 0.01 % from 8/28 to 11/28 with no new
false calls; stderr reports how many within-group signatures were used and how
many were withheld because a carrier was a candidate. `--no-group-signatures`
ignores the file, which is useful for comparing the two on a sample.

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
  -m, --mem=N           memory budget in GB for reference caching (0 = stream) [64]
      --auto-gsize      size the Bloom filter from the largest input genome
      --auto            size the Bloom filter from an ntCard cardinality estimate
      --fpr=F           target Bloom-filter false-positive rate (sets bits; e.g. 0.001)
      --low-mem         alias for --mem 0 (stream: minimal memory, more I/O)
      --help            display this help and exit
      --version         output version information and exit
```

For large reference sets (e.g. thousands of multi-megabase genomes), caching every
reference's k-mer hashes in memory can require more RAM than is available. The
`-m, --mem=N` budget (GB, default 64) bounds this: `comparesketch` caches as many
second-list references as fit within the budget and streams the rest, re-reading
them per first-list reference. `--mem 0` (equivalently `--low-mem`) streams fully,
holding only one Bloom filter and one reference per thread at a time — minimal
memory, more I/O. The similarity output is identical regardless of the budget.

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

The detection floor this sets is easy to estimate. A reference at abundance *a* in
*N* reads of length *L* has coverage *a·N·L/G* for genome size *G*, and each
signature is hit by roughly coverage × (L−k+1)/L × P(error-free k-mer) reads: about
0.1 hits per signature for a 170 kb genome at 0.01% in 3 million read pairs of
151 bp with k=81. With `-h 10` that needs ~100 signatures to reach the threshold on
average and ~500 to reach it reliably. On the 8-genome test set at 0.01%, `-c50`
detected 7 of 28 distant contaminants, the default 100 detected 20 (measured on a
flat index; for distant references the layout makes no difference), and `-c500` all
28. Raising `-c` past a near-identical strain's candidate pool (typically a few
hundred k-mers) cannot add signatures, but it does take the whole pool, and the
strain's detection is then bounded by how many reads land on those loci rather than
by the index. Index size and query memory grow linearly with `-c`.

**Index layout: flat, `--cluster=N`, or `--tree`.** This is the choice that most
affects what your results mean, so it is worth understanding before picking one.

A flat index (the default) holds every reference at full resolution and reports each
one by name. Query cost grows with the size of the whole reference set, since every
signature is in play.

`--cluster=N` merges references that lie within N distinguishing k-mers of each
other and indexes a single representative per cluster. This makes the index smaller
and the query cheaper, and it stops near-identical references from splitting each
other's abundance. The cost is easy to understate: **the members of a cluster stop
existing as far as the output is concerned.** If the organism in your sample is a
cluster member that is not the representative, its abundance is reported correctly
but attributed to the representative's name. In an 8-genome test set where four
near-identical references clustered under one representative, samples containing
each of the three non-representative members were all reported as the
representative, at the correct ~90% abundance. That is a confident wrong
identification rather than a missing one, which is worse if strain identity is what
you care about. Check `clusters.tsv` to see exactly which references were folded
together, and only use `--cluster` when resolution within a cluster genuinely does
not matter to you.

`--tree=T1,T2,...` is the option to reach for when you want both. It clusters
internally to build a hierarchy and uses it to narrow the search coarse-to-fine, but
the final query runs against the full-resolution leaf sketch, so **results come back
as individual references, not representatives.** On the same test set, tree mode
resolved all three non-representative members to their own names at the correct
abundance. At matched signature count its accuracy and low-abundance sensitivity
track a flat index closely, so the hierarchy is about keeping queries tractable as
the reference set grows rather than about accuracy.

In short: flat for full resolution on a reference set small enough to query
directly, `--tree` for full resolution with the search narrowed, and `--cluster`
only when you accept representative-level answers.

**`--max-homopolymer` and `--min-margin`.** Both harden signatures against read
errors and are described in detail above. `--max-homopolymer=6` is a reasonable
default to adopt. `--min-margin=2 --strict-margin` is the default for `--tree`
builds; add it to flat builds too when your references include close relatives,
since it removes the error-driven cross-hits that otherwise appear at the 0.01–0.1%
level whenever one relative is abundant.

**`-t, --threads`.** Set this to your core count for both indexing and querying;
both stages parallelize well.

## Troubleshooting

**`Error: insufficient k-mer content for Bloom filter construction`.** The sample or
reference set does not contain enough distinct k-mers to size the Bloom filters. The
usual cause is reads shorter than `k`: with the default k=81, 75 bp reads yield no
k-mers at all. Check your read length and, if it is short, rebuild the index with a
smaller `-k`. Near-empty or truncated input files produce the same error.

**`Warning: <file>: no k-mers of length 81 ...` or `... no unique k-mers ...`
during a build.** The first means the file yielded nothing to index: it is empty,
not FASTA, or every sequence in it is shorter than `k`. The second means every
k-mer of that reference also occurs in another reference, which is what a
duplicate genome under two names looks like. In both cases the reference gets no
signatures and can never be reported; fix the input or drop the file. Two
references whose file names differ only in directory or extension (`a/x.fa` and
`b/x.fa`) are refused outright, because the reference is named by its base name.

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

## License

UniqSketch is released under the [MIT License](LICENSE), an OSI-approved open-source
license. The same license applies to the archived releases on Zenodo.

## Citation

If you use UniqSketch, please cite:

> Mohamadi H, Turner KH, Schmieder R, Hurwitz BL. UniqSketch: sensitive and
> resource-efficient strain-level detection in metagenomes. *Genome Biology* (2026,
> in press).

The source code is archived on Zenodo at
[10.5281/zenodo.23046441](https://doi.org/10.5281/zenodo.23046441). That is the concept
DOI, so it always resolves to the latest archived release; each release also has its own
version DOI. A machine-readable citation is in [`CITATION.cff`](CITATION.cff).
