#!/usr/bin/env python3
"""
tree_workflow.png: the tree index (uniqsketch --tree) and its query, two
columns, three stages each. Abstract references (R1..R6), one coarse level
shown (L0_1..L0_4; the near-identical references share L0_4). Regenerate with
    python3 docs/img/tree_workflow.py
(matplotlib only; writes tree_workflow.png beside this script).
"""
import os
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(ROOT, "tree_workflow.png")

fig, (ax_idx, ax_qry) = plt.subplots(1, 2, figsize=(16, 8.5))
for ax in (ax_idx, ax_qry):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")

# ---- palette ----------------------------------------------------------------
UQBLUE   = "#124678"
UQORANGE = "#DE8028"
UQGREY   = "#505050"
C_REF    = "#4ECDC4"
C_CLOSE  = "#FF6B6B"
C_TREE   = "#FFEAA7"
C_SKETCH = "#96CEB4"
C_HIT    = "#98D8C8"
C_PRUNE  = "#E5E5E5"
C_READ   = "#DDA0DD"
ARROW    = "#666666"

# ---- helpers ----------------------------------------------------------------
def box(ax, x, y, w, h, text, color, fontsize=11, bold=False, italic=False,
        edge="#333", lw=1.2):
    r = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                       facecolor=color, edgecolor=edge, linewidth=lw)
    ax.add_patch(r)
    weight = "bold" if bold else "normal"
    style  = "italic" if italic else "normal"
    ax.text(x + w/2, y + h/2, text, ha="center", va="center",
            fontsize=fontsize, weight=weight, style=style, color="#1A1A1A")

def arrow(ax, x1, y1, x2, y2, lw=2.0, color=ARROW):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw))

def stage_label(ax, y, text, n):
    ax.add_patch(Circle((0.6, y), 0.35, facecolor=UQORANGE,
                        edgecolor="none", zorder=3))
    ax.text(0.6, y, str(n), ha="center", va="center",
            fontsize=13, fontweight="bold", color="white", zorder=4)
    ax.text(1.2, y, text, ha="left", va="center",
            fontsize=14, fontweight="bold", color=UQBLUE)

# ============================================================
# LEFT: INDEXING
# ============================================================
ax_idx.set_title("Indexing", fontsize=18, fontweight="bold",
                 color=UQBLUE, pad=14, loc="left")

# ── Stage 1: reference universe ────────────────────────────
stage_label(ax_idx, 9.1, "Reference universe", 1)
labels_distant = [r"R$_1$", r"R$_2$", r"R$_3$"]
labels_close   = [r"R$_4$", r"R$_5$", r"R$_6$"]
y_refs = 8.0
for i, r in enumerate(labels_distant):
    box(ax_idx, 1.0 + i*1.30, y_refs, 1.10, 0.7, r, C_REF, fontsize=13)
for i, r in enumerate(labels_close):
    box(ax_idx, 5.5 + i*1.30, y_refs, 1.10, 0.7, r, C_CLOSE, fontsize=13)
ax_idx.text(2.65, 7.6, "distant", ha="center", fontsize=10,
            color=UQGREY, style="italic")
ax_idx.text(7.15, 7.6, "near-identical", ha="center", fontsize=10,
            color=UQGREY, style="italic")

arrow(ax_idx, 5.0, 7.4, 5.0, 6.55)

# ── Stage 2: hierarchical clustering tree ───────────────────
stage_label(ax_idx, 6.3, "Cluster by k-mer similarity", 2)

# Level-0 ribbons removed for simplicity. Instead, use 4 named L0 nodes.
tree_top_y = 5.05
l0_names = [r"L0$_1$", r"L0$_2$", r"L0$_3$"]
for i, lbl in enumerate(l0_names):
    cx = 1.55 + i*1.30
    box(ax_idx, cx - 0.55, tree_top_y, 1.10, 0.6, lbl, C_TREE, fontsize=11,
        bold=True)
    arrow(ax_idx, cx, tree_top_y, cx, 4.05, lw=1.2)

# Close-cluster Level-0 node = L0_4
cx_close = 6.80
box(ax_idx, cx_close - 1.20, tree_top_y, 2.40, 0.6, r"L0$_4$",
    C_TREE, fontsize=11, bold=True)
for i in range(3):
    arrow(ax_idx, cx_close, tree_top_y, 5.5 + 0.55 + i*1.30, 4.05, lw=1.2)

# Leaves
y_leaves = 3.55
for i, r in enumerate(labels_distant):
    box(ax_idx, 1.0 + i*1.30, y_leaves, 1.10, 0.5, r, C_REF, fontsize=11)
for i, r in enumerate(labels_close):
    box(ax_idx, 5.5 + i*1.30, y_leaves, 1.10, 0.5, r, C_CLOSE, fontsize=11)

ax_idx.text(0.02, tree_top_y + 0.3, "Level 0", ha="left",
            fontsize=9, color=UQGREY, fontweight="bold")
ax_idx.text(0.02, y_leaves + 0.25, "Leaves", ha="left",
            fontsize=9, color=UQGREY, fontweight="bold")

arrow(ax_idx, 5.0, 3.45, 5.0, 2.65)

# ── Stage 3: per-level sketches ────────────────────────────
stage_label(ax_idx, 2.4, "Per-level sketches", 3)
box(ax_idx, 0.6, 0.7, 4.3, 1.05,
    "Coarse sketch per level\n(node representatives)",
    C_SKETCH, fontsize=12)
box(ax_idx, 5.1, 0.7, 4.3, 1.05,
    "Leaf sketch (every reference)\n+ within-group signatures",
    C_SKETCH, fontsize=12)

# ============================================================
# RIGHT: QUERYING
# ============================================================
ax_qry.set_title("Querying", fontsize=18, fontweight="bold",
                 color=UQBLUE, pad=14, loc="left")

stage_label(ax_qry, 9.1, "Sequencing reads", 1)
box(ax_qry, 2.6, 7.85, 4.8, 0.85, "r1.fq.gz   +   r2.fq.gz",
    C_READ, fontsize=13, bold=True)

arrow(ax_qry, 5.0, 7.6, 5.0, 6.55)

stage_label(ax_qry, 6.3, "Top-down traversal", 2)

y_top = 5.05
group_w, group_h = 1.6, 0.65
group_xs = [0.5, 2.35, 4.20, 6.05, 7.90]
status   = ["✓",  "✗",   "✓",   "✗",  "✗"]
for x, st in zip(group_xs, status):
    color = C_HIT if st == "✓" else C_PRUNE
    edge  = UQBLUE if st == "✓" else "#999"
    label = f"group {st}"
    box(ax_qry, x, y_top, group_w, group_h, label,
        color, fontsize=11, bold=(st == "✓"), edge=edge,
        lw=1.5 if st == "✓" else 0.8)

ax_qry.text(9.5, 4.72, "match the coarse sketch,\nkeep the surviving groups",
            ha="right", va="top", fontsize=10, style="italic", color=UQGREY)

y_leaf = 3.30
leaf_w, leaf_h = 0.85, 0.55
sx1 = group_xs[0] + group_w/2
box(ax_qry, sx1 - leaf_w/2, y_leaf, leaf_w, leaf_h, r"R$_a$",
    C_HIT, fontsize=11)
arrow(ax_qry, sx1, y_top, sx1, y_leaf + leaf_h, lw=1.4, color=UQBLUE)

sx2_center = group_xs[2] + group_w/2
leaf_offsets = [-1.0, 0.0, 1.0]
for off in leaf_offsets:
    lx = sx2_center + off
    box(ax_qry, lx - leaf_w/2, y_leaf, leaf_w, leaf_h, "",
        C_HIT, fontsize=10)
    arrow(ax_qry, sx2_center, y_top, lx, y_leaf + leaf_h,
          lw=1.4, color=UQBLUE)
ax_qry.text(sx2_center, y_leaf - 0.38, "fine resolution within the close cluster",
            ha="center", fontsize=9, style="italic", color=UQGREY)

arrow(ax_qry, 5.0, 2.85, 5.0, 2.2)

stage_label(ax_qry, 1.95, "Per-reference abundance", 3)
box(ax_qry, 0.8, 0.45, 8.4, 1.10,
    r"abundance$_i$ = (count$_i$ / sigs$_i$) $\,/\,$ $\sum_j$ (count$_j$ / sigs$_j$)"
    "\nover the references and families reported",
    C_SKETCH, fontsize=12)

plt.tight_layout(w_pad=3)
plt.savefig(OUT, dpi=150, bbox_inches="tight",
            facecolor="white")
print("Saved:", OUT)
