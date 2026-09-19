#!/usr/bin/env python3
"""Generate the deterministic fixtures for test/e2etest.cpp.

Everything comes from a fixed-seed RNG, so re-running this script reproduces the
committed files byte for byte. Re-run only if the fixture design changes, then
commit the regenerated files.

Design
------
Four synthetic genomes, 20 kb each:

  refA        random sequence
  refB        random sequence, unrelated to refA
  refA_var    refA with a handful of point substitutions -- a near-identical
              variant. Its only unique k-mers are the ones overlapping those
              edits, which is the close-strain regime that has exposed real bugs
              (the Bloom-filter race turned shared k-mers into signatures for
              both siblings).
  refC        refB with a handful of point substitutions; sequenced into the
              sample but deliberately NOT indexed. It must never produce a row.
              Making it a relative of an indexed genome (rather than an
              unrelated random one) is what gives that check teeth: an
              unrelated genome shares no 31-mers with anything and cannot leak
              whatever the code does, whereas refC's reads carry thousands of
              k-mers refB also has. The only reason refB is not inflated by
              refC's reads is that those shared k-mers are correctly excluded
              from refB's signatures -- so this also checks that refB's
              abundance is not pushed up by an unindexed relative.

Reads: 100 bp single-end, error-free, sampled uniformly from each genome at a
known mixture. Error-free is a deliberate choice: the test is checking the
counting and abundance arithmetic, not error tolerance, and noise would only
widen the tolerances needed.

Mixture (reads per genome):
  refA      7000
  refB      2500
  refA_var   500
  refC       500   -> present in reads, absent from index

Expected abundance. Signatures are unique against the INDEX, not against the
world. refC is not indexed, so the k-mers it shares with refB (99.94% of them)
were never seen as repeats and remain legitimate refB signatures; refC's reads
therefore land on refB and are counted as refB. This is documented, intended
behaviour -- an unindexed organism is invisible as its own row, and its reads are
absorbed by its nearest indexed relative -- and it is why the paper stresses
building the database from the strains you actually expect to see. The truth
table encodes that: refB's expected share is (refB + refC) reads over all
indexed-or-absorbed reads. The test asserts it, so the contract is pinned rather
than treated as noise.
"""
import os
import random

SEED = 20260919
GENOME_LEN = 20_000
READ_LEN = 100
K = 31                       # short k keeps fixture reads small while k << read length
N_VARIANT_EDITS = 12

MIX = [("refA", 7000), ("refB", 2500), ("refA_var", 500), ("refC", 500)]

HERE = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(SEED)


def random_genome(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, n_edits):
    s = list(seq)
    # spread edits out so their k-mer footprints do not overlap
    positions = sorted(rng.sample(range(K, len(s) - K), n_edits))
    for p in positions:
        s[p] = rng.choice([b for b in "ACGT" if b != s[p]])
    return "".join(s), positions


def write_fasta(name, seq):
    with open(os.path.join(HERE, name + ".fa"), "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")


genomes = {}
genomes["refA"] = random_genome(GENOME_LEN)
genomes["refB"] = random_genome(GENOME_LEN)
genomes["refA_var"], edit_pos = mutate(genomes["refA"], N_VARIANT_EDITS)
genomes["refC"], _ = mutate(genomes["refB"], N_VARIANT_EDITS)

for name, seq in genomes.items():
    write_fasta(name, seq)

# Reads: uniform start positions, forward strand, error-free. Written gzipped
# (mtime fixed so the archive is reproducible): it is an eighth of the size to
# commit, and it exercises the compressed-input path that no other test covers.
import gzip
with gzip.GzipFile(os.path.join(HERE, "reads.fq.gz"), "wb", mtime=0) as fq:
    rid = 0
    for name, n in MIX:
        g = genomes[name]
        for _ in range(n):
            start = rng.randrange(0, len(g) - READ_LEN + 1)
            rid += 1
            fq.write(f"@r{rid}_{name}\n{g[start:start + READ_LEN]}\n+\n{'I' * READ_LEN}\n".encode())

# Ground truth for the test to read. refC's reads are attributed to refB, the
# indexed relative that absorbs them (see docstring).
counts = dict(MIX)
absorbed = {"refA": counts["refA"], "refB": counts["refB"] + counts["refC"],
            "refA_var": counts["refA_var"]}
total = sum(absorbed.values())
with open(os.path.join(HERE, "truth.tsv"), "w") as t:
    t.write("ref\treads\texpected_abundance\n")
    for n, c in absorbed.items():
        t.write(f"{n}\t{c}\t{c / total:.6f}\n")

print(f"seed={SEED} k={K} genome={GENOME_LEN}bp read={READ_LEN}bp")
print(f"refA_var edits at: {edit_pos}")
for n, c in MIX:
    print(f"  {n:<9} {c:>5} reads")
