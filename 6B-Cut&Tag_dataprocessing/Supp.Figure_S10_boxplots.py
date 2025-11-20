#!/usr/bin/env python3
"""
Boxplots of % 5mC and % 5hmC in genomic windows (IGV snapshots) .
- Axis labels and overall header are removed (clean layout).
- Boxes are colored by col6 using the provided hex palette.
"""

import argparse
import gzip
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import math
import re

# --------- Color palette (hex, without leading # accepted) ----------
PALETTE = {
    "H3K4Me1.mc":  "46982E",
    "H3K4Me1.hmc": "A3CC98",
    "H3K27ac.mc":  "2629CE",
    "H3K27ac.hmc": "9294E6",
    "H3K4Me3.mc":  "FF9300",
    "H3K4Me3.hmc": "F8FE67",
    "H3K27Me3.mc": "EA3323",
    "H3K27Me3.hmc":"F09890",
    "IgG.mc":      "9C9C9C",
    "IgG.hmc":     "CACACA",
}
# --------------------------------------------------------------------

def smart_open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")

def sanitize_filename(name: str) -> str:
    name = re.sub(r"[^\w\-.]+", "_", str(name))
    return name.strip("_") or "group"

def to_hex(c):
    c = str(c).lstrip("#")
    return f"#{c}" if len(c) in (3, 6) else "#808080"  # fallback gray

def color_for(cat6: str) -> str:
    return to_hex(PALETTE.get(cat6, "808080"))

def color_boxes(bp, labels):
    """Apply face/edge colors to each box according to its label."""
    for box, lab in zip(bp["boxes"], labels):
        col = color_for(lab)
        box.set_facecolor(col)
        box.set_edgecolor("black")
        box.set_linewidth(1.1)
    # Style medians/whiskers/caps consistently
    for k in ("medians", "whiskers", "caps"):
        for artist in bp[k]:
            artist.set_color("black")
            artist.set_linewidth(1.1)

def make_box(ax, data, labels, showfliers):
    bp = ax.boxplot(
        data,
        labels=labels,
        patch_artist=True,
        showfliers=showfliers,
        medianprops=dict(linewidth=1.5),
        boxprops=dict(linewidth=1.2),
        whiskerprops=dict(linewidth=1.2),
        capprops=dict(linewidth=1.2),
    )
    color_boxes(bp, labels)
    return bp

def main():
    ap = argparse.ArgumentParser(
        description=" Boxplots of % 5mC and % 5hmC in genomic windows (IGV snapshots) "
    )
    ap.add_argument("tsv", help="Input TSV (optionally .gz), no header.")
    ap.add_argument("--outdir", default="plots", help="Directory to write plots.")
    ap.add_argument("--logy", action="store_true", help="Log10 scale on Y axis.")
    ap.add_argument("--showfliers", action="store_true", help="Show fliers on boxplots.")
    ap.add_argument("--dpi", type=int, default=300, help="PNG resolution.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Read file
    with smart_open(args.tsv) as fh:
        df = pd.read_csv(
            fh,
            sep=r"\t",
            header=None,
            engine="python",
            names=["chr", "start", "end", "value", "cat5", "cat6"],
        )

    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["value", "cat5", "cat6"])
    if df.empty:
        raise SystemExit("No valid rows after parsing 'value', 'cat5', and 'cat6'.")

    unique_cat5 = list(pd.unique(df["cat5"]))

    # One plot per cat5 (split), boxes = cat6 (group)
    for c5 in unique_cat5:
        sub = df[df["cat5"] == c5].copy()
        if sub.empty:
            continue

        # Order cat6 by median value (descending)
        cat6_medians = sub.groupby("cat6")["value"].median().sort_values(ascending=False)
        order = list(cat6_medians.index)
        data = [sub.loc[sub["cat6"] == c6, "value"].dropna().values for c6 in order]
        if all(len(x) == 0 for x in data):
            continue

        fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(order)), 5))

        plot_data = data
        if args.logy:
            min_pos = sub["value"][sub["value"] > 0].min()
            eps = (min_pos / 10) if pd.notna(min_pos) else 1e-6
            plot_data = [x + eps for x in data]

        make_box(ax, plot_data, order, args.showfliers)

        if len(order) > 6:
            plt.setp(ax.get_xticklabels(), rotation=30, ha="right")

        # Remove axis labels and overall header; keep per-panel title
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title(c5)
        if args.logy:
            ax.set_yscale("log")

        fig.tight_layout()
        base = sanitize_filename(c5)
        fig.savefig(outdir / f"boxplot_{base}.png", dpi=args.dpi)
        fig.savefig(outdir / f"boxplot_{base}.pdf")
        plt.close(fig)

    # Multi-panel summary: one panel per cat5
    n = len(unique_cat5)
    if n > 1:
        cols = min(3, n)
        rows = math.ceil(n / cols)
        fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 4.5 * rows), squeeze=False)

        for i, c5 in enumerate(unique_cat5):
            r, c = divmod(i, cols)
            ax = axes[r][c]
            sub = df[df["cat5"] == c5]
            if sub.empty:
                ax.axis("off")
                continue

            order = list(sub.groupby("cat6")["value"].median().sort_values(ascending=False).index)
            data = [sub.loc[sub["cat6"] == c6, "value"].dropna().values for c6 in order]
            if all(len(x) == 0 for x in data):
                ax.axis("off")
                continue

            plot_data = data
            if args.logy:
                min_pos = sub["value"][sub["value"] > 0].min()
                eps = (min_pos / 10) if pd.notna(min_pos) else 1e-6
                plot_data = [x + eps for x in data]

            make_box(ax, plot_data, order, args.showfliers)
            if len(order) > 5:
                plt.setp(ax.get_xticklabels(), rotation=30, ha="right")

            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.set_title(str(c5))
            if args.logy:
                ax.set_yscale("log")

        # Hide any unused subplots
        for j in range(i + 1, rows * cols):
            r, c = divmod(j, cols)
            axes[r][c].axis("off")

        fig.tight_layout(rect=[0, 0, 1, 1])
        fig.savefig(outdir / "box.png", dpi=args.dpi)
        fig.savefig(outdir / "box.pdf")
        plt.close(fig)

if __name__ == "__main__":
    main()

