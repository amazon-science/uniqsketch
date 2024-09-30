# UniqSketch: Sensitive and resource-efficient strain/abundance identification in metagenomics data   

UniqSketch is a sensitive tool for identifying strains and their relative abundances in metagenomics samples. UniqSketch algorithm consists of two stages: Indexing and Querying. 
* **Indexing**: the algorithm takes a list of target sequences (usually reference genomes) and then identifies all k-mers within each target in the list that uniquely belong to that specific sequence. 
It then checks a set of criteria such as low-complexity on the k-mer candidates to refine the list and provide the final set of unique k-mers for each target.  
To efficiently (memory and runtime) keep track of k-mers and their frequencies at indexing and querying stages, UniqSketch utilizes the Bloom filter data structure as an alternative for hash tables due to their costly computational resource requirements. 
* **Query**: raw read sequencing data are queried against the sketch index. 
 To remove the k-mers caused by sequencing errors (k-mers with very low depth), 
 UniqSketch utilizes a cascading Bloom filter to discard weak k-mers and. 
 The solid k-mers then queried against the index and the count for related reference target will be updated accordingly. 
 After all solid k-mers in the raw read input data are queried, 
 the reference count table is calculated computed to generate the final output of most abundant references. 
# Installation
UniqSketch builds into a stand-alone executable:
```
git clone https://github.com/amazon-science/uniqsketch
cd uniqsketch
make
```
# Quick Start
UniqSketch consists of two modules: 
* `uniqsketch`: for building a signature index from reference/target database
* `querysketch`: for querying the input read files against constructed index 

1. Construct sketch index of size `81bp` from three references, `ref_1.fa`, `ref_2.fa`, and `ref_3.fa` and store it in `sketch_index`:
```
uniqsketch --sensitive -k81 -o sketch_index.tsv ref_1.fa ref_2.fa ref_3.fa
```
2. Query a metagenomics sample with read files `r1.fq.gz` and `r2.fq.gz`:  
```
querysketch --r1 r1.fq.gz --r2 r2.fq.gz --ref sketch_index.tsv --out out_sample.tsv 
```

# Output
### uniqsketch
* `sketch_uniq.tsv`: primary output file  is a tab-separated file storing the final set of signatures for each reference in the universe.
* `outdir_uniqsketch`: folder contains a file per reference storing all signature candidates.
* `db_uniq_count.tsv`: a tab-separated file summarizing all references and total number of candidate signatures.

### querysketch
* `out_sample.tsv`: primary output file reporting references and their abundance within the sample. The `count` column is the total number of uniqsketch signatures called by querysketch:

```bash
ref	abundance	count
-----------------------------
ref_1	0.7862		3033
ref_2	0.1993		769
ref_3	0.01426		55
```

* `log_out_sample.tsv`: detailed log showing for each reference the total number of assigned reads, and all `signature`:`count`.   

# Usage
### uniqsketch
```bash
Usage: uniqsketch [OPTION] @LIST_FILES (or FILES)
Creates unique sketch from a list of fasta reference files.
A list of files containing file names in each row can be passed with @ prefix.

 Options:

  -t, --threads=N	use N parallel threads [1]
  -k, --kmer=N		the length of kmer [81]
  -b, --bit=N		use N bits per element in Bloom filter [128]
  -d, --hash1=N		distinct Bloom filter hash number [5]
  -s, --hash2=N		repeat BF hash number [5]
  -c, --cov=N		number of unique k-mers to represent a reference [100]
  -f, --outdir=STRING	dir for universe unique k-mers [outdir_uniqsketch]
  -o, --out=STRING	the output sketch file name [sketch_uniq.tsv]
  -r, --stat=STRING	the output unique kmer stat file name [db_uniq_count.tsv]
  -e, --entropy		sets the aggregate entropy rate threshold [0.65]
      --sensitive	sets sensitivity parameter c to 100
      --very-sensitive	sets sensitivity parameter c to 1000
      --help		display this help and exit
      --version		output version information and exit
```
### querysketch
```bash
Usage: querysketch [OPTIONS] [ARGS]
Identify references and their abundance in FILE(S).

Acceptable file formats: fastq in compressed formats gz, bz, zip, xz.

 Options:

  -t, --threads=N	use N parallel threads [1]
  -b, --bit=N		use N bits per element in Bloom filter [64]
  -o, --out=STRING	the output file name
  -l, --r1=STRING	input read 1
  -r, --r2=STRING	input read 2
  -g, --ref=STRING	input uniqsketch reference
  -h, --hit=N		number of uniqsketch hits to call a refernce [10]
  -a, --acutoff=N	abundance cutoff to report [0.0]
  -s, --rcutoff=N	read cutoff to report [2]
      --sensitive	sets sensitivity parameter h=10
      --very-sensitive	sets sensitivity parameter h=5
      --solid		only use solid k-mers in reads
      --help		display this help and exit
      --versions	version information and exit
  ```

