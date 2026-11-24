#!/usr/bin/env python3
"""Rule-stack test of the tree query, through the real binaries.

usage: run_tree_rules.py BINDIR [WORKDIR] [--jobs N] [--query-threads N] [--build-threads N]
                         [--reuse] [--out NAME] [--index-digest FILE] [-v]

runTreeQuery lives in src/querysketch.cpp and no unit test can call it, so the
rules it applies after the read pass (carried-signature settlement, pruning and
rejoin, leaf rescue and its single-read rule, the within-group mask, family rows
and their suppression, the breadth floors) were exercised only by the validation
sweeps, which take hours and live outside the repository. This test builds a
16-genome synthetic index in which each of those rules has a structure of its own
(gen_fixture.py), queries it with samples that each make one rule decide the
result (gen_reads.py), and compares the rows with expected.tsv. It also runs the
index variants an older uniqsketch would have written, so the fallbacks stay
exercised. Nothing is committed but the generators and the expectations; the
whole run takes about ten seconds.

Three stages, so that a failure says what moved:
  INDEX  the index holds the structures the samples rely on (read from the index
         files, by name; nothing here depends on the leaf numbering or on which
         k-mers were selected);
  ROWS   each query's exit status and result rows (name, resolution, members)
         are those of expected.tsv;
  RULE   counters on stderr and numbers in the table show that the rule the case
         is named for is what produced the rows.
A failed INDEX check means the build changed and the samples no longer test what
they claim; the ROWS and RULE failures that follow it are consequences.

What the cases do not pin (each survived a one-edit mutant of querysketch):
  - the family hit threshold in runTreeQuery: generateTreeQueryResult applies it
    again, so removing either copy changes nothing;
  - which of a node's calls is kept when it is scored at two levels: both levels
    hold the same signatures here, so the two calls are the same numbers;
  - the cap of a member's carried count at the signatures in use (family breadth):
    node X always has more node-unique signatures in use (517) than its
    least-carrying member carries (483);
  - "explained" as the sum over resolved members: no node here has two resolvable
    members carrying different parts of its signatures;
  - a carrier seen through its within-group signatures, a carrier's lineage
    attested on two reads rather than one, and a signature-less leaf blocking as a
    present carrier on coarse-level evidence: the genomes have no carrier of that
    kind (a carrier with within-group signatures; a signature-less carrier under a
    node with ten signatures in one read; a signature-less carrier of a
    within-group signature that the narrowing dropped).

Exit status: 0 when every check passes, 1 otherwise.
"""
import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.dont_write_bytecode = True      # nothing is written into the source tree
import gen_fixture  # noqa: E402
import gen_reads    # noqa: E402

# Thresholds sized to the fixture: clusters at 12 000, 4 000 and 300 distinguishing
# k-mers give the three internal levels; two carriers per carried node signature.
BUILD = ["--tree=12000,4000,300", "--node-carriers=2", "-c1000", "-k81"]
QUERY_FILES = ["tree_index.tsv", "leaves.tsv", "build_params.tsv", "sketch_level_0.tsv", "sketch_level_1.tsv",
               "sketch_level_2.tsv", "sketch_group.tsv", "ancestor_signatures.tsv"]
LOCUS_GAP = 500            # QuerySketchConfig::locusGap

checks = []                # (stage, case, description, ok, detail)


def check(stage, case, what, ok, detail=""):
    checks.append((stage, case, what, bool(ok), detail))


def tsv(path, skip=0):
    with open(path) as fh:
        return [line.rstrip("\n").split("\t") for line in fh][skip:]


def run(cmd, cwd, log):
    with open(log, "w") as fh:
        try:
            return subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, timeout=300).returncode
        except subprocess.TimeoutExpired:
            fh.write("\nrun_tree_rules.py: timed out after 300 s\n")
            return -1


# ---------------------------------------------------------------- the index and its variants
def build_index(bindir, work, threads):
    G = gen_fixture.write(work)
    os.makedirs(os.path.join(work, "idx"), exist_ok=True)
    rc = run([os.path.join(bindir, "uniqsketch")] + BUILD + ["-t", str(threads), "-o", "idx/sketch_uniq.tsv",
              "-f", "idx/outdir", "-r", "idx/db.tsv", "@list.txt"], work, os.path.join(work, "idx", "build.log"))
    if rc != 0:
        sys.exit("uniqsketch failed (exit %d); see %s" % (rc, os.path.join(work, "idx", "build.log")))
    return G


def make_variants(work):
    """Copies of the index with one piece removed or altered, as an index from an older uniqsketch would be."""
    idx = os.path.join(work, "idx")
    vdir = os.path.join(work, "variants")
    shutil.rmtree(vdir, ignore_errors=True)

    def mk(name):
        out = os.path.join(vdir, name, "outdir")
        os.makedirs(out)
        shutil.copy(os.path.join(idx, "sketch_uniq.tsv"), os.path.join(vdir, name))
        for f in QUERY_FILES:
            shutil.copy(os.path.join(idx, "outdir", f), out)
        return out

    def edit(path, old, new):
        with open(path) as fh:
            text = fh.read()
        assert old in text, (path, old)
        with open(path, "w") as fh:
            fh.write(text.replace(old, new))

    os.remove(os.path.join(mk("no_ancestor"), "ancestor_signatures.tsv"))
    os.remove(os.path.join(mk("no_params"), "build_params.tsv"))
    os.remove(os.path.join(mk("no_group"), "sketch_group.tsv"))
    os.remove(os.path.join(mk("no_leaves"), "leaves.tsv"))
    os.remove(os.path.join(mk("no_level0"), "sketch_level_0.tsv"))
    os.remove(os.path.join(mk("no_level1"), "sketch_level_1.tsv"))
    os.remove(os.path.join(mk("no_tree"), "tree_index.tsv"))
    out = mk("no_carriers")
    for lvl in "012":
        p = os.path.join(out, "sketch_level_%s.tsv" % lvl)
        rows = tsv(p)
        with open(p, "w") as fh:
            fh.write("".join("\t".join(r[:6]) + "\n" for r in rows))
    edit(os.path.join(mk("largest_rep"), "build_params.tsv"), "node_representative\tcentral", "node_representative\tlargest")
    edit(os.path.join(mk("margin1"), "build_params.tsv"), "strict_margin\t1", "strict_margin\t0")
    edit(os.path.join(mk("other_bf"), "build_params.tsv"), "bloom_filter\tblocked", "bloom_filter\tstandard")
    # A k-mer with two lines in a level sketch, as a carried k-mer that filled two sibling
    # nodes' quotas has: node QA's first signature again, on a line of node QB's.
    p = os.path.join(mk("dup_line"), "sketch_level_1.tsv")
    rows = tsv(p)
    a = next(r for r in rows if r[2] == "A2")
    b = next(r for r in rows if r[2] == "B2")
    with open(p, "a") as fh:
        fh.write("\t".join([a[0], b[1], "B2"] + a[3:]) + "\n")
    # leaves.tsv one leaf short of the leaf sketch
    p = os.path.join(mk("short_leaves"), "leaves.tsv")
    with open(p) as fh:
        lines = fh.readlines()
    with open(p, "w") as fh:
        fh.writelines(lines[:-1])


def index_digest(work, path):
    """One line per index file the query reads (and the build's two reports): sha256 and name.

    The index is built from relative paths, so two builds of one source are expected
    to give the same lines whatever the machine, compiler or standard library; `diff`
    of two such files names what differs.
    """
    idx = os.path.join(work, "idx")
    names = ["sketch_uniq.tsv"] + sorted(os.path.join("outdir", f) for f in os.listdir(os.path.join(idx, "outdir"))
                                         if f.endswith(".tsv"))
    with open(path, "w") as out:
        for name in names:
            with open(os.path.join(idx, name), "rb") as fh:
                out.write("%s  %s\n" % (hashlib.sha256(fh.read()).hexdigest(), name))


# ---------------------------------------------------------------- INDEX: what the samples rely on
def loci(positions):
    """Number of loci: runs of signature positions less than LOCUS_GAP apart."""
    positions = sorted(positions)
    return sum(1 for i, p in enumerate(positions) if i == 0 or p - positions[i - 1] > LOCUS_GAP)


def index_checks(work):
    """Check the planted structures; return what the RULE stage needs from the index."""
    out = os.path.join(work, "idx", "outdir")
    reach = {r[1]: r for r in tsv(os.path.join(out, "reachability.tsv"), 1)}
    status = {n: r[5] for n, r in reach.items()}
    check("INDEX", "-", "which references have signatures of their own, within-group ones only, or none",
          status == dict(W="group", X1="family", X1s="family", X3="family", SL="family", Y2="family", A2="family",
                         B2="family", X2="own", SB="own", Y1="own", A1="own", B1="own", D="own", V="own", U="own"), str(status))
    leaf_id = {r[1]: r[0] for r in tsv(os.path.join(out, "leaves.tsv"), 1)}
    nodes = {r[1]: r for r in tsv(os.path.join(out, "tree_index.tsv"), 1)}
    want = ["cluster_L0_B2", "cluster_L0_X3", "cluster_L1_A2", "cluster_L1_B2", "cluster_L1_W", "cluster_L1_X1",
            "cluster_L1_Y2", "cluster_L2_X1"]
    got = sorted(n for n in nodes if n.startswith("cluster"))
    check("INDEX", "-", "the eight nodes and their names", got == want, str(got))
    if got != want:
        return None
    check("INDEX", "-", "three nested levels over X1 (level 2 in level 1 in level 0)",
          nodes["cluster_L2_X1"][3] == nodes["cluster_L1_X1"][0] and nodes["cluster_L1_X1"][3] == nodes["cluster_L0_X3"][0])
    lvl0 = tsv(os.path.join(out, "sketch_level_0.tsv"))
    lvl1 = tsv(os.path.join(out, "sketch_level_1.tsv"))
    # Node X's level-1 signatures by carrier set, as leaf names.
    name_of = {v: k for k, v in leaf_id.items()}
    xsig = {}
    for r in lvl1:
        if r[2] == "X1":
            key = frozenset(name_of[c] for c in r[6].split(",")) if r[6] != "-" else frozenset()
            xsig[key] = xsig.get(key, 0) + 1
    shown = {",".join(sorted(k)) or "-": v for k, v in xsig.items()}
    check("INDEX", "-", "node X: node-unique signatures, and carried ones with a sibling carrier (SB), "
          "an untestable one (SL) and a pair in another lineage (Y1, Y2)",
          xsig.get(frozenset(), 0) >= 200 and xsig.get(frozenset(["SB"]), 0) >= 50 and
          xsig.get(frozenset(["SL"]), 0) >= 50 and xsig.get(frozenset(["Y1", "Y2"]), 0) >= 50, str(shown))
    anc = {(r[0], r[2]): (int(r[3]), int(r[4])) for r in tsv(os.path.join(out, "ancestor_signatures.tsv"), 1)}
    x3 = anc.get(("X3", "cluster_L1_X1"), (1, 1))
    check("INDEX", "-", "node X: X3 carries under 60 % of its signatures (the rest are private to the representative)",
          x3[1] < 0.6 * x3[0], str(x3))
    sb = anc.get(("SB", "cluster_L0_X3"), (0, 0))
    check("INDEX", "-", "lineage PP (level 0): SB carries part of its signatures, with a margin of 100 or more",
          sb[0] >= 500 and 100 <= sb[1] <= sb[0] - 250, str(sb))
    r0 = anc.get(("B2", "cluster_L0_B2"), (0, 0))
    check("INDEX", "-", "node R (level 0): A1 and A2 carry none of its signatures, B1 and B2 all",
          anc.get(("A1", "cluster_L0_B2"), (0, 1))[1] == 0 and anc.get(("A2", "cluster_L0_B2"), (0, 1))[1] == 0 and
          r0[1] == r0[0] > 50 and anc.get(("B1", "cluster_L0_B2")) == r0, str(r0))
    grp = tsv(os.path.join(out, "sketch_group.tsv"))
    wg = [r for r in grp if r[2] == "W"]
    vg = [r for r in grp if r[2] == "V"]
    check("INDEX", "-", "W: within-group signatures only, every one carried by U",
          len(wg) >= 100 and all(r[6] == leaf_id["U"] for r in wg) and reach["W"][2] == "0", "%d" % len(wg))
    check("INDEX", "-", "V: within-group signatures nobody else carries, and some U carries",
          sum(r[6] == "-" for r in vg) >= 12 and sum(r[6] == leaf_id["U"] for r in vg) >= 12, "%d" % len(vg))
    for rep, node in (("W", "WV"), ("Y2", "Y")):
        a = sum(1 for r in lvl0 if r[2] == rep)
        b = sum(1 for r in lvl1 if r[2] == rep)
        check("INDEX", "-", "node %s is scored at levels 0 and 1" % node, a >= 100 and b >= 100, "%d and %d signatures" % (a, b))
    pos = {}
    for r in tsv(os.path.join(work, "idx", "sketch_uniq.tsv")):
        pos.setdefault(r[2], []).append(int(r[5]))
    nloci = {n: loci(pos.get(n, [])) for n in ("SB", "Y1")}
    check("INDEX", "-", "carriers of node X: SB has its signatures at two loci, Y1 at one", nloci == dict(SB=2, Y1=1), str(nloci))
    return dict(xsig=xsig)


# ---------------------------------------------------------------- the queries
def load_cases():
    cases = []
    for r in tsv(os.path.join(HERE, "expected.tsv")):
        if not r or not r[0] or r[0].startswith("#"):
            continue
        assert len(r) == 7, r
        cases.append(dict(id=r[0], index=r[1], reads=r[2].split("+"), options=[] if r[3] == "-" else r[3].split(),
                          exit=int(r[4]), rows=r[5], pins=r[6]))
    return cases


def query(bindir, work, outdir, threads, case):
    idx = "idx" if case["index"] == "idx" else os.path.join("variants", case["index"])
    cmd = [os.path.join(bindir, "querysketch"), "-t", str(threads), "--tree=" + os.path.join(idx, "outdir"),
           "--r1", os.path.join("reads", case["reads"][0])]
    if len(case["reads"]) > 1:
        cmd += ["--r2", os.path.join("reads", case["reads"][1])]
    cmd += ["--ref", os.path.join(idx, "sketch_uniq.tsv"), "--out", os.path.join(outdir, case["id"] + ".tsv")] + case["options"]
    case["rc"] = run(cmd, work, os.path.join(work, outdir, case["id"] + ".stderr"))


class Result:
    """One query's table and stderr."""
    def __init__(self, work, outdir, case):
        self.table = os.path.join(work, outdir, case["id"] + ".tsv")
        self.log = os.path.join(work, outdir, "log_" + case["id"] + ".tsv")
        with open(os.path.join(work, outdir, case["id"] + ".stderr")) as fh:
            self.err = fh.read()
        self.rows = {}
        if os.path.exists(self.table):
            for r in tsv(self.table, 1):
                self.rows[r[0]] = dict(abundance=float(r[1]), count=int(r[2]), reads=int(r[3]), sigs_hit=int(r[4]),
                                       sigs_total=int(r[5]), resolution=r[6], members=r[7])

    def kinds(self):
        return set("%s:%s%s" % (n, v["resolution"], "" if v["members"] == "-" else ":" + v["members"])
                   for n, v in self.rows.items())

    def grab(self, pattern, keys):
        m = re.search(pattern, self.err)
        return dict(zip(keys, map(int, m.groups()))) if m else {}

    def settle(self):
        return self.grab(r"(\d+) in use \((\d+) with hits, for (\d+) node\(s\)\), (\d+) set aside with their hits "
                         r"\((\d+) had hits\): (\d+) for a carrier seen or attested by its lineage, (\d+) for a carrier "
                         r"that may be present; carriers of nodes with hits: (\d+) seen or attested, (\d+) excluded on "
                         r"their own loci, (\d+) possible, (\d+) without signatures \((\d+) of them possible",
                         ["used", "used_hit", "nodes", "aside", "aside_hit", "aside_present", "aside_possible", "seen",
                          "excluded", "possible", "untestable", "untestable_possible"])

    def leaf(self):
        return self.grab(r"leaf level: (\d+) leaves in play after narrowing \((\d+) with no global signatures\), (\d+) "
                         r"more by hits on their own global signatures, (\d+) by hits on their within-group signatures "
                         r"\((\d+) of these on one read", ["narrowed", "sigless", "rescued_global", "rescued_group", "one_read"])

    def group(self):
        return self.grab(r"within-group signatures: (\d+) in use for (\d+) leaves, (\d+) withheld", ["used", "leaves", "withheld"])

    def rejoined(self):
        m = re.search(r"(\d+) node\(s\) pruned at a coarser level scored on their own evidence", self.err)
        return int(m.group(1)) if m else 0

    def leaf_hits(self, name):
        """Hits on a leaf's signatures, reported or not (the per-signature log)."""
        if not os.path.exists(self.log):
            return None
        return sum(int(r[3]) for r in tsv(self.log, 1) if r[0] == name)


def in_use(xsig, absent):
    """Node X's level-1 signatures in use when exactly the carriers in `absent` are absent."""
    return sum(n for carriers, n in xsig.items() if carriers <= frozenset(absent))


def rule_checks(R, facts, ix):
    """Evidence that the rule each case is named for produced its rows."""
    def rule(case, what, ok, detail=""):
        check("RULE", case, what, ok, detail)

    xsig = ix["xsig"]
    s = R["S01_X2"].settle()
    rule("S01_X2", "carriers: 2 excluded on their own loci, 2 without signatures of which 1 possible; carried signatures in use were hit",
         s.get("excluded") == 2 and s.get("untestable") == 2 and s.get("untestable_possible") == 1 and s.get("used_hit", 0) > 0, str(s))
    r = R["S02_X3"]
    s = r.settle()
    rule("S02_X3", "carried signatures in use were hit; the level-0 call was dropped for the level-1 one",
         s.get("used_hit", 0) > 0 and "cluster_L0_X3\tfor cluster_L1_X1" in r.err, str(s))
    want = in_use(xsig, ["SB", "Y1", "Y2"])
    got = r.rows.get("cluster_L1_X1", {}).get("sigs_total")
    rule("S02_X3", "the row's sigs_total is the signatures in use: node-unique, plus carried with every carrier absent "
         "(SB and Y1 excluded, Y2 silent), minus those set aside for SL", got == want, "%s, expected %d" % (got, want))
    s = R["S03_X3_low"].settle()
    rule("S03_X3_low", "both carriers with signatures are possible; no carried signature is in use",
         s.get("possible") == 2 and s.get("used") == 0, str(s))
    s = R["S04_SL"].settle()
    rule("S04_SL", "carried signatures set aside for a carrier that may be present; none in use was hit",
         s.get("aside_possible", 0) > 0 and s.get("used_hit") == 0, str(s))
    s = R["S05_SB"].settle()
    rule("S05_SB", "one carrier seen; its carried signatures set aside", s.get("seen") == 1 and s.get("aside_present", 0) > 0, str(s))
    r = R["S06_SB_X3"]
    rule("S06_SB_X3", "the family row is reported beside a resolved member; the level-1 call is withheld by the breadth floor",
         "reported beside a resolved member" in r.err and re.search(r"family-level call\(s\) withheld.*\n  cluster_L1_X1\t", r.err) is not None)
    s = R["S07_Y2"].settle()
    rule("S07_Y2", "both carriers in lineage Y are attested by it, the one with signatures included (it is not put to the "
         "exclusion test); no carried signature in use was hit",
         s.get("seen") == 2 and s.get("excluded") == 0 and s.get("used_hit") == 0, str(s))
    rule("S08_A2", "one node rejoined", R["S08_A2"].rejoined() == 1)
    rule("S09_B2", "the level-0 call was dropped for the level-1 one; nothing rejoined",
         "cluster_L0_B2\tfor cluster_L1_B2" in R["S09_B2"].err and R["S09_B2"].rejoined() == 0)
    rule("S10_A1", "one leaf rescued on its global signatures", R["S10_A1"].leaf().get("rescued_global") == 1, str(R["S10_A1"].leaf()))
    rule("S11_X1", "two enclosing calls dropped", R["S11_X1"].err.count("\tfor cluster_") == 2)
    g = R["S12_W"].group()
    rule("S12_W", "within-group signatures in use", g.get("used", 0) >= 100, str(g))
    g = R["S13_W_U"].group()
    rule("S13_W_U", "within-group signatures withheld, none in use", g.get("withheld", 0) >= 100 and g.get("used") == 0, str(g))
    rule("S15_V_amplicon", "a strain call was withheld by the breadth floor", "strain call(s) withheld" in R["S15_V_amplicon"].err)
    rule("S16_PP_locus", "a family call was withheld by the breadth floor", "family-level call(s) withheld" in R["S16_PP_locus"].err)
    r = R["S17_B2_b1read"]
    rule("S17_B2_b1read", "the strain row stands on one read", r.rows.get("B1", {}).get("reads") == 1, str(r.rows.get("B1")))
    rule("S19_D_b1read2", "one leaf rescued on its global signatures", R["S19_D_b1read2"].leaf().get("rescued_global") == 1)
    ab = {n: v["abundance"] for n, v in R["S22_mix"].rows.items()}
    rule("S22_mix", "abundances within 0.03 of 0.5, 0.25, 0.125, 0.125",
         all(abs(ab.get(n, 0) - t) < 0.03 for n, t in (("X2", 0.5), ("A1", 0.25), ("U", 0.125), ("D", 0.125))), str(ab))
    rule("S23_V_grouprescue", "one leaf rescued on its within-group signatures", R["S23_V_grouprescue"].leaf().get("rescued_group") == 1)
    r = R["S26_A2_a1read"]
    rule("S26_A2_a1read", "one node rejoined; one rescue on a single read", r.leaf().get("one_read") == 1 and r.rejoined() == 1, str(r.leaf()))
    s = R["S29_X2_12x"].settle()
    rule("S29_X2_12x", "of the two carriers with signatures, one is excluded (two loci) and one is possible (one locus)",
         s.get("excluded") == 1 and s.get("possible") == 1, str(s))
    r = R["S30_W_trace"]
    g = r.group()
    rule("S30_W_trace", "within-group signatures in use; nothing rejoined; node WV has three hits at each level",
         g.get("used", 0) >= 100 and r.rejoined() == 0 and facts["w_trace_hits"] == (3, 3), "%s %s" % (g, facts["w_trace_hits"]))
    r = R["S31_X3_x2short"]
    rule("S31_X3_x2short", "X2 has exactly five hits, on one read", r.leaf_hits("X2") == 5, str(r.leaf_hits("X2")))
    r = R["S32_D_a2read"]
    rule("S32_D_a2read", "the read holds ten or more signatures of node QA; nothing rejoined",
         facts["a2_read_signatures"] >= 10 and r.rejoined() == 0, "%d signatures" % facts["a2_read_signatures"])
    rule("S33_SB_sl_low", "SL's reads hit 10 to 90 of the node's signatures beyond SB's; no family row beside SB",
         10 <= facts["sl_low_signatures"] <= 90 and "reported beside a resolved member" not in R["S33_SB_sl_low"].err,
         "%d signatures" % facts["sl_low_signatures"])
    r = R["S34_X3_sbread"]
    s = r.settle()
    rule("S34_X3_sbread", "SB is seen on its one read and what it carries is set aside",
         s.get("seen", 0) >= 1 and s.get("aside_present", 0) >= 50 and r.rows.get("SB", {}).get("reads") == 1, str(s))
    r = R["S36_D_b2read"]
    rule("S36_D_b2read", "the family row stands on one read; nothing rejoined",
         r.rows.get("cluster_L1_B2", {}).get("reads") == 1 and r.rejoined() == 0, str(r.rows.get("cluster_L1_B2")))
    r = R["S37_D_b1mates"]
    rule("S37_D_b1mates", "B1 has two reads, one per mate", r.rows.get("B1", {}).get("reads") == 2, str(r.rows.get("B1")))
    a, b = R["S38_X2_X3x3"].rows.get("X2", {}), R["S24_X2_X3"].rows.get("X2", {})
    rule("S38_X2_X3x3", "X2's own evidence is that of S24 (the extra depth is all the sibling's)",
         a.get("count") == b.get("count") and a.get("count", 0) > 0, "%s vs %s" % (a.get("count"), b.get("count")))
    r = R["S39_B1_a1read"]
    rule("S39_B1_a1read", "A1 is rescued on its global signatures, on one read",
         r.leaf().get("rescued_global") == 1 and r.leaf().get("one_read") == 1 and r.rows.get("A1", {}).get("reads") == 1, str(r.leaf()))
    r = R["S40_D_b1read_qbread"]
    rule("S40_D_b1read_qbread", "B1 stands on one read, was kept by the narrowing (no rescue), and nothing rejoined",
         r.rows.get("B1", {}).get("reads") == 1 and r.leaf().get("rescued_global") == 0 and r.rejoined() == 0, str(r.leaf()))
    a, b = R["V_dup_line_S09"].rows.get("cluster_L1_B2", {}), R["S09_B2"].rows.get("cluster_L1_B2", {})
    rule("V_dup_line_S09", "sigs_total is one more than on the index as built; hits and reads are the same",
         a.get("sigs_total") == b.get("sigs_total", 0) + 1 and a.get("count") == b.get("count") and a.get("reads") == b.get("reads"),
         "%s vs %s" % (a, b))
    for case, text in (("V_no_ancestor_S16", "has no ancestor_signatures.tsv"), ("V_no_params_S19", "predates build_params.tsv"),
                       ("V_margin1_S19", "built without a strict margin"), ("V_no_group_S12", "has no within-group signatures"),
                       ("V_no_leaves_S02", "unnamed and kept in play"), ("V_short_leaves_S12", "ignoring it"),
                       ("V_largest_rep_S16", "largest member as node representative"), ("V_other_bf_S12", "Bloom filter"),
                       ("V_no_tree_S12", "tree index file not found"), ("O_solid_empty", "insufficient k-mer content")):
        rule(case, "stderr says: %s" % text, text in R[case].err)
    rule("V_no_carriers_S02", "no settlement without the carriers column", "carried node signatures" not in R["V_no_carriers_S02"].err)
    rule("V_no_level0_S08", "nothing rejoined", R["V_no_level0_S08"].rejoined() == 0)


def main():
    ap = argparse.ArgumentParser(description="Rule-stack test of the tree query (see the header of this file).")
    ap.add_argument("bindir", help="directory holding uniqsketch and querysketch")
    ap.add_argument("workdir", nargs="?", default=os.path.join(HERE, "..", "tree_rules_work"),
                    help="where the genomes, index, reads and outputs go; emptied first [test/tree_rules_work]")
    ap.add_argument("--jobs", type=int, default=2, help="queries run side by side [2]")
    ap.add_argument("--query-threads", type=int, default=2, help="threads of each query [2]")
    ap.add_argument("--build-threads", type=int, default=4, help="threads of the index build [4]")
    ap.add_argument("--reuse", action="store_true", help="keep the index and reads of an earlier run in WORKDIR")
    ap.add_argument("--out", default="out", help="sub-directory of WORKDIR for the query outputs [out]")
    ap.add_argument("--index-digest", metavar="FILE", help="also write a sha256 line per index file, for comparing builds")
    ap.add_argument("-v", "--verbose", action="store_true", help="list every check, not only the failures")
    args = ap.parse_args()
    bindir = os.path.abspath(args.bindir)
    work = os.path.abspath(args.workdir)
    for tool in ("uniqsketch", "querysketch"):
        if not os.path.exists(os.path.join(bindir, tool)):
            sys.exit("%s not found in %s (build the tools first)" % (tool, bindir))

    reuse = args.reuse and os.path.exists(os.path.join(work, "reads", "EMPTY.fq"))
    if reuse:
        B, G = gen_fixture.build()
        G = dict(G, NEG=B["NEG"])
    else:
        shutil.rmtree(work, ignore_errors=True)
        os.makedirs(work)
        G = build_index(bindir, work, args.build_threads)
        make_variants(work)
    if args.index_digest:
        index_digest(work, os.path.abspath(args.index_digest))
    ix = index_checks(work)
    facts = gen_reads.write_all(work, G, os.path.join(work, "idx"), files=not reuse) if ix is not None else {}

    cases = load_cases()
    if ix is not None:
        shutil.rmtree(os.path.join(work, args.out), ignore_errors=True)
        os.makedirs(os.path.join(work, args.out))
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            list(pool.map(lambda c: query(bindir, work, args.out, args.query_threads, c), cases))
        R = {c["id"]: Result(work, args.out, c) for c in cases}
        for c in cases:
            r = R[c["id"]]
            ok = c["rc"] == c["exit"]
            detail = "exit %d, expected %d" % (c["rc"], c["exit"])
            if ok and c["rows"].startswith("="):
                other = R[c["rows"][1:]]
                ok = os.path.exists(r.table) and os.path.exists(other.table) and open(r.table).read() == open(other.table).read()
                detail = "table differs from that of %s" % c["rows"][1:]
            elif ok and c["rows"] != "*":
                want = set() if c["rows"] == "-" else set(c["rows"].split(","))
                ok = os.path.exists(r.table) and r.kinds() == want
                detail = "got %s; expected %s" % (", ".join(sorted(r.kinds())) or "no rows", ", ".join(sorted(want)) or "no rows")
            check("ROWS", c["id"], c["pins"], ok, detail)
        rule_checks(R, facts, ix)

    failed = [c for c in checks if not c[3]]
    for stage, case, what, ok, detail in checks:
        if args.verbose or not ok:
            print("%-4s %-5s %-18s %s%s" % ("ok" if ok else "FAIL", stage, case, what, "" if ok or not detail else "\n       [" + detail + "]"))
    if ix is None:
        print("the index does not hold the fixture's structures; no query was run")
    print("tree_rules: %d checks passed, %d failed (%d queries)" % (len(checks) - len(failed), len(failed), len(cases) if ix is not None else 0))
    sys.exit(1 if failed or ix is None else 0)


if __name__ == "__main__":
    main()
