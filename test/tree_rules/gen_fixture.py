#!/usr/bin/env python3
"""Synthetic genomes for the tree query's rule-stack test (run_tree_rules.py).

usage: gen_fixture.py WORKDIR

Writes WORKDIR/genomes/<name>.fa (16 genomes, 13 to 40 kb, 444 kb in all) and
WORKDIR/list.txt. The list holds paths relative to WORKDIR, and the index is
built from inside WORKDIR, so the path strings the build sees are the same
wherever the repository is checked out (the leaf numbering follows them). No
randomness beyond a fixed-seed xorshift64*, so the files are byte-identical on
every run; nothing generated is committed.

Blocks are independent random sequences; a genome is a concatenation of blocks
or slices of them. Each structure is there for one rule of the query (see
expected.tsv for the samples that use it). Sizes in bases.

 lineage PP  (level 0)          X3  = CXa CXb CXc                 core only: no signature of any kind
   node X    (level 1)          X2  = CXa CXb CXc P2              resolvable member (P2 private)
     node N2 (level 2)          X1  = CXa CXb CXc P1              one substitution apart: neither has
                                X1s = X1 with one substitution    a signature (strict margin)
   sibling leaf                 SB  = SBp1 CXa SBp2               carries CXa; own signatures at 2 loci
   sibling leaf, no signatures  SL  = tail(CXa) CXb head(CXc)     carries CXb; a substring of X3
 lineage Y   (level 1, orphan)  Y1  = YC CXc Y1p                  carries CXc; own signatures at 1 locus
                                Y2  = YC CXc                      carries CXc; no signatures; has a lineage
 lineage R   (level 0)
   node QA   (level 1)          A1  = RC QA a1p ; A2 = RC QA      A2 without signatures
   node QB   (level 1)          B1  = RC QB b1p ; B2 = RC QB      B2 without signatures
 decoy                          D   = Dp RC'                      RC' = RC with a substitution every 81 bases:
                                                                  every k-mer of RC then has a Hamming-1
                                                                  neighbour in the index and fails the margin,
                                                                  so node R's signatures are QB's and the A
                                                                  lineage carries none of them (rejoin)
 node WV     (level 1)          W   = WC Z                        Z and its junction also in U: no global
                                                                  signatures, within-group ones carried by U
                                V   = WC Vp Z2                    Vp private, Z2 also in U: its within-group
                                                                  selection draws Vp k-mers the global one did
                                                                  not (carrier-free) and Z2 k-mers (carrier U)
 unrelated                      U   = Up tail(WC) Z tail(Vp) Z2
 not indexed                    NEG = random                      negative control (reads only)

Node X is the settlement case: its representative X1 has 517 node-unique
signatures (P1) and the quota is filled with 483 carried ones, whose carriers
are a sibling with signatures (SB), a sibling without (SL) and a pair in
another lineage (Y1, Y2).
"""
import os
import sys

MASK = (1 << 64) - 1
NAMES = ["X1", "X1s", "X2", "X3", "SB", "SL", "Y1", "Y2", "A1", "A2", "B1", "B2", "D", "W", "V", "U"]


class Rng:
    def __init__(self, seed):
        self.x = (seed * 0x9E3779B97F4A7C15 + 0x1234567) & MASK or 1

    def next(self):
        x = self.x
        x ^= (x >> 12)
        x ^= (x << 25) & MASK
        x ^= (x >> 27)
        self.x = x
        return (x * 0x2545F4914F6CDD1D) & MASK


def rand_dna(n, seed):
    r = Rng(seed)
    return "".join("ACGT"[(r.next() >> 33) & 3] for _ in range(n))


def substitute(seq, pos):
    b = seq[pos]
    return seq[:pos] + "CGTA"["ACGT".index(b)] + seq[pos + 1:]


def build():
    """The blocks and the genomes, as {name: sequence}."""
    B = {}
    for name, n in [("CXa", 15000), ("CXb", 13000), ("CXc", 6000), ("P1", 600), ("P2", 600),
                    ("SBp1", 500), ("SBp2", 500), ("YC", 30000), ("Y1p", 1000),
                    ("RC", 20000), ("QA", 5000), ("QB", 2000), ("a1p", 1000), ("b1p", 1000),
                    ("Dp", 20000), ("WC", 16000), ("Z", 600), ("Vp", 2200), ("Up", 30000), ("NEG", 20000), ("Z2", 500)]:
        B[name] = rand_dna(n, 1000 + len(B))
    rcp = B["RC"]
    for p in range(40, len(rcp), 81):
        rcp = substitute(rcp, p)
    cx = B["CXa"] + B["CXb"] + B["CXc"]
    x1 = cx + B["P1"]
    G = {
        "X1": x1,
        "X1s": substitute(x1, len(cx) + 300),
        "X2": cx + B["P2"],
        "X3": cx,
        "SB": B["SBp1"] + B["CXa"] + B["SBp2"],
        "SL": B["CXa"][-100:] + B["CXb"] + B["CXc"][:100],
        "Y1": B["YC"] + B["CXc"] + B["Y1p"],
        "Y2": B["YC"] + B["CXc"],
        "A1": B["RC"] + B["QA"] + B["a1p"],
        "A2": B["RC"] + B["QA"],
        "B1": B["RC"] + B["QB"] + B["b1p"],
        "B2": B["RC"] + B["QB"],
        "D": B["Dp"] + rcp,
        "W": B["WC"] + B["Z"],
        "V": B["WC"] + B["Vp"] + B["Z2"],
        "U": B["Up"] + B["WC"][-100:] + B["Z"] + B["Vp"][-100:] + B["Z2"],
    }
    assert list(G) == NAMES
    return B, G


def write(work):
    """Write the genomes and the list; return {name: sequence} including NEG."""
    os.makedirs(os.path.join(work, "genomes"), exist_ok=True)
    B, G = build()
    for name, seq in G.items():
        with open(os.path.join(work, "genomes", name + ".fa"), "w") as fh:
            fh.write(">%s synthetic\n" % name)
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")
    with open(os.path.join(work, "list.txt"), "w") as fh:
        fh.write("".join("genomes/%s.fa\n" % name for name in G))
    G = dict(G)
    G["NEG"] = B["NEG"]
    return G


if __name__ == "__main__":
    genomes = write(sys.argv[1])
    print("genomes: %d, bases: %d" % (len(genomes) - 1, sum(len(s) for n, s in genomes.items() if n != "NEG")))
