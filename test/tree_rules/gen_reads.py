#!/usr/bin/env python3
"""Error-free reads for the tree query's rule-stack test (run_tree_rules.py).

Called after the index is built: a few samples are made of reads that hold
exactly one chosen signature, read from the index's own sketch files, so that a
threshold can be approached from either side without depending on which k-mers
the build happened to select.

Reads are 150 bases tiled at a fixed stride on alternating strands. The only
random numbers are the flanks of the one-signature reads (fixed-seed
xorshift64*, as in gen_fixture.py), so the files are the same on every run for
a given index.
"""
import gzip
import os

from gen_fixture import rand_dna

READ = 150
K = 81
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def tile(name, seq, cov, begin=0, end=None):
    """Reads covering seq[begin:end] at depth cov, evenly spaced, alternating strands."""
    end = len(seq) if end is None else end
    n = max(2, int(round(cov * (end - begin) / READ)))
    out = []
    for i in range(n):
        start = begin + i * (end - begin - READ) // (n - 1)
        r = seq[start:start + READ]
        if i % 2:
            r = revcomp(r)
        out.append(("%s_%d" % (name, i), r))
    return out


def window(name, seq, start, copies, length=READ):
    return [("%s_w%d_%d" % (name, start, i), seq[start:start + length]) for i in range(copies)]


def sketch(path):
    """A sketch file as (sequence, name, carriers) rows; carriers is '-' when the column is absent."""
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append((f[0], f[2], f[6] if len(f) > 6 else "-"))
    return rows


def locate(genome, sig):
    """Position of a signature in a genome, whichever strand the sketch stores."""
    p = genome.find(sig)
    return p if p >= 0 else genome.find(revcomp(sig))


def one_signature_read(name, genome, sig, seed):
    """A read holding one k-mer of the genome and nothing else of it.

    The k-mer sits between random flanks whose bases next to it differ from the
    genome's, so the neighbouring k-mers, which in these dense sketches are
    usually signatures too, are not reproduced: the read is one hit, exactly.
    """
    p = locate(genome, sig)
    assert p >= 0, "signature of %s not found in its genome" % name
    left, right = rand_dna(30, seed), rand_dna(READ - K - 30, seed + 500)
    if p > 0 and left[-1] == genome[p - 1]:
        left = left[:-1] + "CGTA"["ACGT".index(left[-1])]
    if p + K < len(genome) and right[0] == genome[p + K]:
        right = "CGTA"["ACGT".index(right[0])] + right[1:]
    return ("%s_sig%d" % (name, seed), left + genome[p:p + K] + right)


def spread(items, n):
    """n items taken evenly across a list (fewer if the list is shorter)."""
    step = max(1, len(items) // n)
    return items[::step][:n]


def write(path, reads, gz=False):
    op = gzip.open if gz else open
    with op(path, "wt") as fh:
        for name, seq in reads:
            fh.write("@%s\n%s\n+\n%s\n" % (name, seq, "I" * len(seq)))


def samples(G, idx):
    """({sample: reads}, facts). G: genomes by name; idx: the index directory (sketch_uniq.tsv and outdir/).

    facts are numbers about the index-dependent reads that the test asserts on, so
    that a sample which no longer approaches its threshold as designed is reported
    as that and not as a rule failure.
    """
    def cov(g, c):
        return tile(g, G[g], c)

    leaf = sketch(os.path.join(idx, "sketch_uniq.tsv"))
    lvl0 = sketch(os.path.join(idx, "outdir", "sketch_level_0.tsv"))
    lvl1 = sketch(os.path.join(idx, "outdir", "sketch_level_1.tsv"))
    grp = sketch(os.path.join(idx, "outdir", "sketch_group.tsv"))

    # V's within-group signatures that no other genome carries, for the within-group rescue:
    # twelve reads, each holding one of them between random flanks, so that nothing else of V is hit.
    free = [s for s, name, car in grp if name == "V" and car == "-"]
    v_rescue = [("Vgrp_%d" % i, rand_dna(30, 9000 + i) + sig + rand_dna(39, 9100 + i))
                for i, sig in enumerate(spread(free, 12))]

    # Node WV is scored at levels 0 and 1. Three hits at each (the same three
    # signatures when both levels drew them, as they do today): enough to survive
    # the permissive pruning (one hit), short of --tree-prune=5 and far from a call.
    w0 = set(s for s, name, car in lvl0 if name == "W")
    w1 = set(s for s, name, car in lvl1 if name == "W")
    both = sorted(w0 & w1)
    picks = spread(both, 3) if len(both) >= 3 else spread(sorted(w0 - w1), 3) + spread(sorted(w1 - w0), 3)
    w_trace = [one_signature_read("W", G["W"], s, 9200 + i) for i, s in enumerate(picks)]
    w_hits = (sum(1 for s in picks if s in w0), sum(1 for s in picks if s in w1))

    # The 150-base window of A2 holding the most level-1 signatures of node QA: one read
    # that clears the hit threshold for the node by itself.
    a2 = G["A2"]
    at = sorted(locate(a2, s) for s, name, car in lvl1 if name == "A2")
    best = max(at, key=lambda p: (sum(1 for q in at if p <= q <= p + READ - K), -p))
    a2_read = window("A2", a2, best, 1)
    a2_hits = sum(1 for q in at if best <= q <= best + READ - K)

    # 85 bases of X2's private block: five k-mers, every one a leaf signature of X2.
    x2 = G["X2"]
    x2sig = set(locate(x2, s) for s, name, car in leaf if name == "X2")
    start = next(p for p in range(len(x2) - 350, len(x2) - K - 4) if all(q in x2sig for q in range(p, p + 5)))
    x2_short = window("X2", x2, start, 1, length=K + 4)

    # SL at a quarter of 1x: its reads reach a few dozen of lineage PP's level-0
    # signatures, the ones SB does not carry.
    sl_low = cov("SL", 0.25)
    pp_in_sl = [p for p in (locate(G["SL"], s) for s, name, car in lvl0 if name == "X3") if p >= 0]
    sl_hits = sum(1 for p in pp_in_sl if any(G["SL"].find(r if i % 2 == 0 else revcomp(r)) <= p <=
                                             G["SL"].find(r if i % 2 == 0 else revcomp(r)) + READ - K
                                             for i, (_, r) in enumerate(sl_low)))

    b1 = G["B1"]
    S = {
        "S01_X2": cov("X2", 20),
        "S02_X3": cov("X3", 20),
        "S03_X3_low": cov("X3", 3),
        "S04_SL": cov("SL", 20),
        "S05_SB": cov("SB", 20),
        "S06_SB_X3": cov("SB", 20) + cov("X3", 20),
        "S07_Y2": cov("Y2", 20),
        "S08_A2": cov("A2", 20),
        "S09_B2": cov("B2", 20),
        "S10_A1": cov("A1", 20),
        "S11_X1": cov("X1", 20),
        "S12_W": cov("W", 20),
        "S13_W_U": cov("W", 20) + cov("U", 20),
        "S14_U": cov("U", 20),
        "S15_V_amplicon": window("V", G["V"], 16400, 60),
        "S16_PP_locus": window("X3", G["X3"], 2000, 100),
        "S17_B2_b1read": cov("B2", 20) + window("B1", b1, 22400, 1),
        "S18_D_b1read": cov("D", 20) + window("B1", b1, 22400, 1),
        "S19_D_b1read2": cov("D", 20) + window("B1", b1, 22400, 1) + window("B1", b1, 22600, 1),
        "S20_NEG": cov("NEG", 20),
        "S22_mix": cov("X2", 20) + cov("A1", 10) + cov("U", 5) + cov("D", 5),
        "S23_V_grouprescue": v_rescue,
        "S24_X2_X3": cov("X2", 20) + cov("X3", 20),
        "S26_A2_a1read": cov("A2", 20) + window("A1", G["A1"], 25400, 1),
        "S27_D_vread": cov("D", 20) + window("V", G["V"], 16400, 1),
        "S28_X3_A2": cov("X3", 20) + cov("A2", 20),
        "S29_X2_12x": cov("X2", 12),
        "S30_W_trace": tile("Wz", G["W"], 20, begin=len(G["W"]) - 700) + w_trace,
        "S31_X3_x2short": cov("X3", 20) + x2_short,
        "S32_D_a2read": cov("D", 20) + a2_read,
        "S33_SB_sl_low": cov("SB", 20) + sl_low,
        "S34_X3_sbread": cov("X3", 20) + window("SB", G["SB"], 100, 1),
        "S36_D_b2read": cov("D", 20) + window("B2", G["B2"], 20500, 1),
        "S38_X2_X3x3": cov("X2", 20) + cov("X3", 60),
        "S39_B1_a1read": cov("B1", 20) + window("A1", G["A1"], 25400, 1),
        "S40_D_b1read_qbread": cov("D", 20) + window("B1", b1, 22400, 1) + window("B2", G["B2"], 20500, 1),
        "EMPTY": [],
    }
    facts = dict(w_trace_hits=w_hits, a2_read_signatures=a2_hits, sl_low_signatures=sl_hits)
    return S, facts


def write_all(work, G, idx, files=True):
    """Write WORKDIR/reads/*.fq (unless files is False) and return the facts of samples()."""
    S, facts = samples(G, idx)
    if not files:
        return facts
    rdir = os.path.join(work, "reads")
    os.makedirs(rdir, exist_ok=True)
    for name, reads in S.items():
        write(os.path.join(rdir, name + ".fq"), reads)
    # the same mixture as two files (two reader threads), and one sample gzip-compressed
    mix = S["S22_mix"]
    write(os.path.join(rdir, "S21_mix_R1.fq"), mix[0::2])
    write(os.path.join(rdir, "S21_mix_R2.fq"), mix[1::2])
    write(os.path.join(rdir, "S25_X3.fq.gz"), S["S02_X3"], gz=True)
    # Two reads of B1 named as the mates of one fragment, beside D: R1 holds D and
    # the first mate, R2 the second.
    b1 = G["B1"]
    write(os.path.join(rdir, "S37_D_b1mates_R1.fq"), S["S18_D_b1read"][:-1] + [("frag_B1/1", b1[22400:22400 + READ])])
    write(os.path.join(rdir, "S37_D_b1mates_R2.fq"), [("frag_B1/2", revcomp(b1[22600:22600 + READ]))])
    return facts
