#!/usr/bin/env python3
"""
us_workflow.png: the UniqSketch workflow, three stacked panels (optional
pre-clustering, indexing, querying). Regenerate with
    python3 docs/img/us_workflow.py
(matplotlib and numpy only; writes us_workflow.png beside this script).
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, PathPatch, Circle
from matplotlib.path import Path
import numpy as np

# ── Palette ──────────────────────────────────────────────────────────────────
C_BG_CMP   = "#EEF4E8"
C_BG_IDX   = "#F3F7FB"
C_BG_QRY   = "#FBF5EF"

C_PRIMARY  = "#124678"
C_ACCENT   = "#DE8028"
C_CLUSTER  = "#4A8A5E"
C_BLOOM_D  = "#4A90C2"
C_BLOOM_S  = "#23507F"
C_REFINE   = "#C7DFEE"
C_TEXT     = "#1C1C1C"
C_ARROW    = "#4A4A4A"

C_A = "#3CB371"
C_C = "#4A90C2"
C_G = "#F2A93B"
C_T = "#D84F4F"
BASE_COLOR = {"A": C_A, "C": C_C, "G": C_G, "T": C_T}

TITLE_KW = dict(fontsize=26, fontweight="bold", color=C_PRIMARY)

# ── Helpers ─────────────────────────────────────────────────────────────────
def rounded(ax, x, y, w, h, facecolor, edgecolor=None, lw=1.3,
            rounding=0.22):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle=f"round,pad=0,rounding_size={rounding}",
                                facecolor=facecolor,
                                edgecolor=edgecolor or facecolor,
                                linewidth=lw))

def arrow(ax, x1, y1, x2, y2, color=C_ARROW, lw=2.2, mutation_scale=22):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                                 arrowstyle="-|>", mutation_scale=mutation_scale,
                                 color=color, lw=lw, shrinkA=0, shrinkB=0))

def bloom_filter(ax, x, y, w, h, title, subtitle, color):
    rounded(ax, x, y, w, h, facecolor=color, edgecolor=C_PRIMARY, lw=1.5)
    ax.text(x + w/2, y + h - 0.42, title,
            ha="center", va="center", fontsize=16, fontweight="bold",
            color="white")
    ax.text(x + w/2, y + h - 0.82, subtitle,
            ha="center", va="center", fontsize=12, color="white",
            alpha=0.92, style="italic")
    pad = 0.22
    n = 14
    cell_w = (w - 2*pad) / n * 0.80
    gap = ((w - 2*pad) - n*cell_w) / (n-1)
    cy = y + 0.24
    ch = 0.34
    for i in range(n):
        cx = x + pad + i*(cell_w + gap)
        filled = (i * 7) % 5 < 3
        ax.add_patch(Rectangle((cx, cy), cell_w, ch,
                               facecolor="white" if filled else (1,1,1,0.25),
                               edgecolor="white", lw=0.4))

def dna_block(ax, x, y, w, h, title, title_color=C_PRIMARY,
              sequences=None, stroke=None):
    stroke = stroke or title_color
    rounded(ax, x, y, w, h, facecolor="white",
            edgecolor=stroke, lw=1.6, rounding=0.25)
    ax.text(x + w/2, y + h - 0.40, title,
            ha="center", va="center",
            fontsize=15, fontweight="bold", color=title_color)
    if sequences is None:
        sequences = ["GATTACAGCTA", "TAACGTTGCAT", "ACCGTATGCAC"]
    strand_area_top = y + h - 0.75
    n_rows = len(sequences)
    row_h = (strand_area_top - (y + 0.45)) / n_rows
    for row, seq in enumerate(sequences):
        row_y = strand_area_top - (row + 1) * row_h + row_h * 0.15
        xs = np.linspace(x + 0.25, x + w - 0.25, 80)
        ys = row_y + row_h * 0.25 + 0.06*np.sin(np.linspace(0, 3.5*np.pi, 80))
        ax.plot(xs, ys, color="#CBD5E0", lw=2.2, zorder=1, solid_capstyle="round")
        usable_w = w - 0.50
        step = usable_w / len(seq)
        for i, b in enumerate(seq):
            lx = x + 0.25 + step*(i + 0.5)
            ly = row_y + row_h * 0.55
            ax.text(lx, ly, b, ha="center", va="center",
                    fontsize=11, fontweight="bold",
                    color=BASE_COLOR.get(b, C_TEXT), zorder=3)

def curly_bracket(ax, x, y_top, y_bot, text, color=C_PRIMARY,
                  text_offset=0.55):
    span = y_top - y_bot
    mid  = (y_top + y_bot) / 2
    w    = 0.18
    verts = [
        (x, y_top),
        (x - w, y_top), (x - w, mid + span*0.08), (x - w*1.6, mid),
        (x - w, mid - span*0.08), (x - w, y_bot), (x, y_bot),
    ]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CURVE4, Path.CURVE4, Path.CURVE4]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor="none",
                           edgecolor=color, lw=1.8))
    ax.text(x - w*1.6 - text_offset, mid, text,
            rotation=90, ha="center", va="center",
            fontsize=13, fontweight="bold", color=color)

def cluster_graph(ax, cx, cy):
    coords = [
        (cx - 0.55, cy + 0.28),
        (cx - 0.22, cy - 0.35),
        (cx + 0.55, cy + 0.14),
    ]
    ax.add_patch(FancyBboxPatch(
        (cx - 0.80, cy - 0.58), 0.90, 0.97,
        boxstyle="round,pad=0.06,rounding_size=0.16",
        facecolor="none", edgecolor=C_CLUSTER, lw=1.8, linestyle="--"))
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
            ax.plot([coords[i][0], coords[j][0]],
                    [coords[i][1], coords[j][1]],
                    color="#A0A8B0", lw=1.0, zorder=1)
    colors = ["#E57272", "#E57272", "#6FA4D6"]
    for (x, y), col in zip(coords, colors):
        ax.add_patch(Circle((x, y), 0.13, facecolor=col,
                            edgecolor=C_PRIMARY, lw=1.3, zorder=2))

def panel_background(ax, bg, W=22, H=7.2):
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), W, H, facecolor=bg,
                           edgecolor=C_PRIMARY, lw=1.5, zorder=-1))

# ── Figure ──────────────────────────────────────────────────────────────────
W = 22
H = 7.2
fig, (ax_top, ax_mid, ax_bot) = plt.subplots(3, 1, figsize=(22, 22))
panel_background(ax_top, C_BG_CMP, W, H)
panel_background(ax_mid, C_BG_IDX, W, H)
panel_background(ax_bot, C_BG_QRY, W, H)

DNA_X, DNA_Y, DNA_W, DNA_H = 0.7, 2.1, 2.4, 2.5
BF_X,  BF_W              = 6.0, 3.0
BF_TOP_Y, BF_BOT_Y, BF_H = 3.30, 1.10, 1.65

# ═════════════════════════════════════════════════════════════════════════════
# PANEL 1 — COMPARESKETCH (unchanged — was good)
# ═════════════════════════════════════════════════════════════════════════════
ax_top.text(0.4, 6.70, "Pre-clustering (optional)",
            **TITLE_KW, ha="left", va="center")
ax_top.text(0.4, 6.10,
            "Collapse near-identical references before indexing: comparesketch measures, uniqsketch --cluster=N merges",
            fontsize=15, color="#555", style="italic", ha="left", va="center")

dna_block(ax_top, DNA_X, DNA_Y, DNA_W, DNA_H, "References",
          sequences=["GATTACAGCTA", "TAACGTTGCAT",
                     "ACCGTATGCAC", "CGATACCGGCT"])
arrow(ax_top, DNA_X + DNA_W + 0.15, DNA_Y + DNA_H/2,
      5.55, DNA_Y + DNA_H/2, lw=2.4)
ax_top.text((DNA_X + DNA_W + 5.55)/2 - 0.3, DNA_Y + DNA_H/2,
            "k-mers", fontsize=12, color="#555", style="italic",
            ha="center", va="center", rotation=90,
            bbox=dict(facecolor=C_BG_CMP, edgecolor="none", pad=1.5))

ps_x, ps_y, ps_w, ps_h = 5.60, 1.80, 4.1, 3.40
rounded(ax_top, ps_x, ps_y, ps_w, ps_h,
        facecolor="white", edgecolor=C_CLUSTER, lw=1.8)
ax_top.text(ps_x + ps_w/2, ps_y + ps_h - 0.55,
            "Pairwise unique\nk-mer distance",
            ha="center", va="center", fontsize=14, fontweight="bold",
            color=C_CLUSTER)
gs = 0.34
grid_x = ps_x + ps_w/2 - 2*gs
grid_y = ps_y + 0.55
for r in range(4):
    for c in range(4):
        val = abs(r - c)
        shade = 1.0 - min(val/3.0, 0.8)
        ax_top.add_patch(Rectangle(
            (grid_x + c*gs, grid_y + r*gs), gs, gs,
            facecolor=(shade*0.55 + 0.25, shade*0.75 + 0.15, shade*0.55 + 0.30),
            edgecolor="white", lw=0.8))
ax_top.text(ps_x + ps_w/2, ps_y + 0.30,
            "comparesketch, via Bloom filters",
            ha="center", va="center", fontsize=11, style="italic",
            color=C_CLUSTER)

arrow(ax_top, ps_x + ps_w + 0.15, ps_y + ps_h/2,
      11.20, ps_y + ps_h/2,
      color=C_CLUSTER, lw=2.2)
ax_top.text((ps_x + ps_w + 11.20)/2, ps_y + ps_h/2 + 0.50,
            "distance matrix",
            fontsize=12, color=C_CLUSTER, ha="center", va="bottom",
            style="italic")

cl_x, cl_y, cl_w, cl_h = 11.25, 1.80, 4.0, 3.40
rounded(ax_top, cl_x, cl_y, cl_w, cl_h,
        facecolor="white", edgecolor=C_CLUSTER, lw=1.8)
ax_top.text(cl_x + cl_w/2, cl_y + cl_h - 0.55,
            "Single-linkage clustering\n+ representative selection",
            ha="center", va="center", fontsize=14, fontweight="bold",
            color=C_CLUSTER)
cluster_graph(ax_top, cl_x + cl_w/2, cl_y + cl_h/2 - 0.15)
ax_top.text(cl_x + cl_w/2, cl_y + 0.35,
            "--cluster=N: pairs below N → one group",
            ha="center", va="center", fontsize=11, style="italic", color="#555")

arrow(ax_top, cl_x + cl_w + 0.15, cl_y + cl_h/2,
      17.25, cl_y + cl_h/2,
      color=C_CLUSTER, lw=2.2)

dna_block(ax_top, 17.30, DNA_Y, DNA_W, DNA_H, "Representatives",
          title_color=C_CLUSTER, stroke=C_CLUSTER,
          sequences=["GATTACAGCTA", "ACCGTATGCAC", "CGATACCGGCT"])
arrow(ax_top, 18.5, 2.05, 18.5, 0.35,
      color=C_ACCENT, lw=2.8, mutation_scale=24)
ax_top.text(19.05, 1.20, "feeds into\nindexing",
            fontsize=11, style="italic", color=C_ACCENT,
            ha="left", va="center")

# ═════════════════════════════════════════════════════════════════════════════
# PANEL 2 — UNIQSKETCH
# Layout:
#   col 1: Clustered refs (DNA)
#   col 2: Cascading BF (Distinct stacked on Solid/Repeat)
#   col 3: 3 refinement boxes stacked vertically in a column
#   col 4: UniqSketch Index
# Arrows flow left -> right; refinement column connects with a vertical
# arrow chain, and a final arrow goes across to the Index.
# ═════════════════════════════════════════════════════════════════════════════
ax_mid.text(0.4, 6.70, "UniqSketch — indexing",
            **TITLE_KW, ha="left", va="center")
ax_mid.text(0.4, 6.10,
            "Build a signature index: k-mers unique to one reference, filtered and spread along the genome",
            fontsize=15, color="#555", style="italic", ha="left", va="center")

# ── Column 1: DNA block ───────────────────────────────────────────────────
DNA2_X, DNA2_Y, DNA2_W, DNA2_H = 0.5, 2.10, 2.3, 2.50
dna_block(ax_mid, DNA2_X, DNA2_Y, DNA2_W, DNA2_H, "References",
          sequences=["GATTACAGCTA", "ACCGTATGCAC", "CGATACCGGCT"])

# ── Column 2: Cascading Bloom filter (two BFs stacked) ────────────────────
BF2_X, BF2_W = 4.70, 3.00
BF2_TOP_Y, BF2_H = 3.30, 1.55
BF2_BOT_Y = 1.10

# k-mers arrow into the top BF (clear of the curly-brace text on BF left side)
arrow(ax_mid, DNA2_X + DNA2_W + 0.10, DNA2_Y + DNA2_H/2,
      BF2_X - 0.15, DNA2_Y + DNA2_H/2, lw=2.4)
ax_mid.text((DNA2_X + DNA2_W + BF2_X)/2 - 0.3, DNA2_Y + DNA2_H/2,
            "k-mers", fontsize=12, color="#555", style="italic",
            ha="center", va="center", rotation=90,
            bbox=dict(facecolor=C_BG_IDX, edgecolor="none", pad=1.5))

bloom_filter(ax_mid, BF2_X, BF2_TOP_Y, BF2_W, BF2_H,
             "Distinct BF", "seen ≥ 1", color=C_BLOOM_D)
bloom_filter(ax_mid, BF2_X, BF2_BOT_Y, BF2_W, BF2_H,
             "Repeat BF", "seen ≥ 2", color=C_BLOOM_S)
arrow(ax_mid, BF2_X + BF2_W/2, BF2_TOP_Y + 0.05,
      BF2_X + BF2_W/2, BF2_BOT_Y + BF2_H - 0.05,
      color=C_ACCENT, lw=2.6, mutation_scale=24)
ax_mid.text(BF2_X + BF2_W + 0.25 - 1.8 + 0.3,
            (BF2_TOP_Y + BF2_BOT_Y + BF2_H) / 2,
            "if already in\ndistinct BF",
            fontsize=12, style="italic", color=C_ACCENT,
            ha="left", va="center")
curly_bracket(ax_mid, x=BF2_X - 0.10,
              y_top=BF2_TOP_Y + BF2_H, y_bot=BF2_BOT_Y,
              text="Cascading Bloom filter")

# ── Column 3: refinement boxes stacked vertically ─────────────────────────
REF_X, REF_W = 11.30, 4.10
REF_BOX_H    = 1.45
REF_GAP      = 0.35
# Top box aligned with top of Distinct BF, then shifted up one box height,
# then nudged down by one-third of that shift.
REF_TOP_Y = BF2_TOP_Y + BF2_H - REF_BOX_H + REF_BOX_H - REF_BOX_H / 3
refine = [
    ("Low-complexity filter\n(Shannon entropy)",   REF_TOP_Y),
    ("Cross-target inspection\n(--min-margin: no other reference\nwithin one substitution)", REF_TOP_Y - REF_BOX_H - REF_GAP),
    ("Select c signatures per reference\n(spread along the genome)", REF_TOP_Y - 2*REF_BOX_H - 2*REF_GAP),
]
for text, by in refine:
    rounded(ax_mid, REF_X, by, REF_W, REF_BOX_H,
            facecolor=C_REFINE, edgecolor=C_PRIMARY, lw=1.3)
    ax_mid.text(REF_X + REF_W/2, by + REF_BOX_H/2, text,
                ha="center", va="center", fontsize=13, color=C_PRIMARY,
                linespacing=1.15)

# Visible down-arrows between the refinement boxes
for (t1, y1), (t2, y2) in zip(refine[:-1], refine[1:]):
    arrow(ax_mid,
          REF_X + REF_W/2, y1,                  # bottom of upper box
          REF_X + REF_W/2, y2 + REF_BOX_H,      # top of lower box
          color=C_PRIMARY, lw=2.2, mutation_scale=20)

# "candidate unique k-mers" arrow: from right side of Solid BF straight across
# into the *middle* box of the refinement column on its left edge
mid_box_y = refine[1][1]
mid_entry_y = mid_box_y + REF_BOX_H/2
CAND_START_X = BF2_X + BF2_W + 1.85 - 1.5    # start 1.5 units further left (10 chars)
arrow(ax_mid,
      CAND_START_X, mid_entry_y,
      REF_X - 0.05, mid_entry_y,
      color=C_PRIMARY, lw=2.4, mutation_scale=22)
ax_mid.text((CAND_START_X + REF_X)/2 - 0.6, mid_entry_y + 0.45,
            "candidate unique k-mers",
            fontsize=12, color=C_PRIMARY, ha="center", va="bottom",
            style="italic")

# ── Column 4: UniqSketch Index ────────────────────────────────────────────
IDX_W, IDX_H = 3.2, 1.70
IDX_X = W - 1.3 - IDX_W
# Vertically centered against the middle refinement box so the exit arrow is short
IDX_Y = mid_box_y + REF_BOX_H/2 - IDX_H/2
rounded(ax_mid, IDX_X, IDX_Y, IDX_W, IDX_H,
        facecolor=C_PRIMARY, edgecolor=C_PRIMARY, lw=1.4)
ax_mid.text(IDX_X + IDX_W/2, IDX_Y + IDX_H/2 + 0.28,
            "UniqSketch Index", ha="center", va="center",
            fontsize=17, fontweight="bold", color="white")
ax_mid.text(IDX_X + IDX_W/2, IDX_Y + 0.40,
            "(sketch_uniq.tsv)",
            ha="center", va="center", fontsize=11, color="white",
            alpha=0.85, style="italic")

# Final signatures arrow: from right side of the MIDDLE "Cross-target inspection"
# box, elbow up (if needed) and across to the Index.
exit_sx = REF_X + REF_W
exit_sy = refine[1][1] + REF_BOX_H/2          # middle box center
elbow_x = exit_sx + 0.60
elbow_y = IDX_Y + IDX_H/2
ax_mid.add_patch(FancyArrowPatch(
    (exit_sx, exit_sy), (elbow_x, exit_sy),
    arrowstyle="-", color=C_ACCENT, lw=2.4, shrinkA=0, shrinkB=0))
ax_mid.add_patch(FancyArrowPatch(
    (elbow_x, exit_sy), (elbow_x, elbow_y),
    arrowstyle="-", color=C_ACCENT, lw=2.4, shrinkA=0, shrinkB=0))
arrow(ax_mid, elbow_x, elbow_y, IDX_X - 0.02, elbow_y,
      color=C_ACCENT, lw=2.4, mutation_scale=22)
ax_mid.text((exit_sx + elbow_x)/2 + 0.15 + 0.45, exit_sy + 0.70,
            "final signatures",
            fontsize=12, color=C_ACCENT, ha="center", va="bottom",
            style="italic")

# Flow to Query panel: from the bottom of the Index straight down off-panel
arrow(ax_mid, IDX_X + IDX_W/2, IDX_Y - 0.02,
      IDX_X + IDX_W/2, 0.10,
      color=C_ACCENT, lw=2.8, mutation_scale=24)
ax_mid.text(IDX_X + IDX_W/2 + 0.25, 0.60,
            "to query", fontsize=11, style="italic",
            color=C_ACCENT, ha="left", va="center")

# ═════════════════════════════════════════════════════════════════════════════
# PANEL 3 — QUERYSKETCH (give Output more room, separate arrow & its label)
# ═════════════════════════════════════════════════════════════════════════════
ax_bot.text(0.4, 6.70, "QuerySketch — matching & abundance",
            **TITLE_KW, ha="left", va="center")
ax_bot.text(0.4, 6.10,
            "Match metagenomic reads against the signature index",
            fontsize=15, color="#555", style="italic", ha="left", va="center")

dna_block(ax_bot, DNA_X, DNA_Y, DNA_W, DNA_H, "Reads",
          sequences=["GATTACAGCTA", "ACCGTNTGCAC", "CGATNNCGGCT"])

arrow(ax_bot, DNA_X + DNA_W + 0.15, DNA_Y + DNA_H/2,
      BF_X - 0.15, DNA_Y + DNA_H/2, lw=2.4)
ax_bot.text((DNA_X + DNA_W + BF_X)/2 - 0.3, DNA_Y + DNA_H/2,
            "k-mers", fontsize=12, color="#555", style="italic",
            ha="center", va="center", rotation=90,
            bbox=dict(facecolor=C_BG_QRY, edgecolor="none", pad=1.5))

bloom_filter(ax_bot, BF_X, BF_TOP_Y, BF_W, BF_H,
             "Distinct BF", "sample k-mers seen ≥ 1", color=C_BLOOM_D)
bloom_filter(ax_bot, BF_X, BF_BOT_Y, BF_W, BF_H,
             "Solid BF", "sample k-mers seen ≥ 2", color=C_BLOOM_S)
arrow(ax_bot, BF_X + BF_W/2, BF_TOP_Y + 0.05,
      BF_X + BF_W/2, BF_BOT_Y + BF_H - 0.05,
      color=C_ACCENT, lw=2.6, mutation_scale=24)
ax_bot.text(BF_X + BF_W + 0.25 - 2.25 + 0.6,
            (BF_TOP_Y + BF_BOT_Y + BF_H) / 2,
            "--solid: drop k-mers\nseen once (errors)",
            fontsize=12, style="italic", color=C_ACCENT,
            ha="left", va="center")
curly_bracket(ax_bot, x=BF_X - 0.10,
              y_top=BF_TOP_Y + BF_H, y_bot=BF_BOT_Y,
              text="Cascading Bloom filter (--solid)")

# Right-side chain (Index, Output, connecting arrows) all sits on the
# mid-row between the two BFs so the arrows are horizontal and line up.
row_y = (BF_TOP_Y + BF_BOT_Y + BF_H) / 2

# Index block
idx_w, idx_h = 2.8, 2.0
idx_x = 12.55
idx_y = row_y - idx_h/2
rounded(ax_bot, idx_x, idx_y, idx_w, idx_h,
        facecolor=C_PRIMARY, edgecolor=C_PRIMARY, lw=1.4)
ax_bot.text(idx_x + idx_w/2, idx_y + idx_h - 0.55,
            "UniqSketch\nIndex",
            ha="center", va="center", fontsize=16,
            fontweight="bold", color="white")
ax_bot.text(idx_x + idx_w/2, idx_y + 0.45,
            "lookup",
            ha="center", va="center", fontsize=13, color="white",
            alpha=0.9, style="italic")

# Solid k-mers -> Index: horizontal arrow at row_y
solid_sx = BF_X + BF_W
idx_enter_x = idx_x - 0.05
arrow(ax_bot, solid_sx + 0.65, row_y, idx_enter_x, row_y,
      color=C_PRIMARY, lw=2.2, mutation_scale=22)
ax_bot.text((solid_sx + 0.65 + idx_enter_x)/2, row_y + 0.45,
            "k-mers (solid ones with --solid)",
            fontsize=12, color=C_PRIMARY, ha="center",
            va="bottom", style="italic")

# Index -> Output: horizontal arrow at same row_y
out_w, out_h = 4.3, 2.2
out_x = 17.30
out_y = row_y - out_h/2
arrow(ax_bot, idx_x + idx_w + 0.10, row_y,
      out_x - 0.05, row_y,
      color=C_ACCENT, lw=2.2, mutation_scale=20)
ax_bot.text((idx_x + idx_w + out_x)/2, row_y + 0.55,
            "per-reference counts",
            fontsize=12, color=C_ACCENT, ha="center",
            va="bottom", style="italic")

# Output card
rounded(ax_bot, out_x, out_y, out_w, out_h,
        facecolor="white", edgecolor=C_PRIMARY, lw=1.4)
ax_bot.text(out_x + out_w/2, out_y + out_h - 0.35,
            "Output", ha="center", va="center",
            fontsize=16, fontweight="bold", color=C_PRIMARY)
rows = [
    ("ref",   "abundance", "count",  True),
    ("ref_A", "0.79",  "3033",   False),
    ("ref_B", "0.20",  " 769",   False),
    ("ref_C", "0.01",  "  55",   False),
]
for i, (a, b, c, hdr) in enumerate(rows):
    yy = out_y + out_h - 0.80 - i*0.30
    fw = "bold" if hdr else "normal"
    ax_bot.text(out_x + 0.40, yy, a, ha="left",  va="center",
                fontsize=12, family="monospace", fontweight=fw, color=C_TEXT)
    ax_bot.text(out_x + 2.05, yy, b, ha="right", va="center",
                fontsize=12, family="monospace", fontweight=fw, color=C_TEXT)
    ax_bot.text(out_x + 3.80, yy, c, ha="right", va="center",
                fontsize=12, family="monospace", fontweight=fw, color=C_TEXT)

# ── Save ────────────────────────────────────────────────────────────────────
plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, hspace=0.10)
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "us_workflow.png")
plt.savefig(out_path, dpi=120, bbox_inches="tight",
            facecolor="white", edgecolor="none")
print(f"Saved {out_path}")
