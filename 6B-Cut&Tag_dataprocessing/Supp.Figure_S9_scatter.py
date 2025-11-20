#!/usr/bin/env python3
"""
Scatter plot with regression line from a tab-delimited file.

Requirements:
  - Python 3.8+
  - pandas
  - matplotlib
  - numpy

Usage:
  python Supp.Figure_S9_scatter.py <INPUT_FILE>

Notes:
  - X = 4th column; Y = 5th column (0-based indices 3 and 4).
  - Axis labels are derived from those column names with 'bam' removed (case-insensitive).
  - Output figure is saved as '<INPUT_FILE>.pdf' ('.pdf' is appended, not replacing existing extensions).
"""

import sys
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


def clean_label(s: str) -> str:
    """Remove 'bam' or '.bam' (case-insensitive) anywhere in the string."""
    return re.sub(r'(?i)\.?bam', '', str(s)).strip()


def main():
    # --- Parse CLI args ---
    if len(sys.argv) < 2:
        sys.exit("Usage: python Supp.Figure_S9_scatter.py <INPUT_FILE>")

    input_path = Path(sys.argv[1])
    INPUT_FILE = str(input_path)

    # Output file: append .pdf to the provided filename (keep original extensions)
    OUTPUT_FILE = input_path.parent / (input_path.name + ".pdf")

    # --- Styling for publication-quality output ---
    mpl.rcParams.update({
        "pdf.fonttype": 42,           # Embed TrueType fonts (better for Illustrator/Inkscape)
        "ps.fonttype": 42,
        "font.size": 11,
        "font.family": "DejaVu Serif",
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "axes.linewidth": 0.8,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.minor.size": 2,
        "ytick.minor.size": 2,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })

    # --- Load data ---
    try:
        df = pd.read_csv(INPUT_FILE, sep="\t")
    except FileNotFoundError:
        sys.exit(f"Error: '{INPUT_FILE}' not found.")
    except Exception as e:
        sys.exit(f"Error reading '{INPUT_FILE}': {e}")

    # --- Ensure enough columns ---
    if df.shape[1] < 5:
        sys.exit("Error: Input file must have at least 5 columns (need 4th and 5th for X and Y).")

    # Select 4th and 5th columns (0-based indices 3 and 4)
    x_col_raw = df.columns[3]
    y_col_raw = df.columns[4]

    # Cleaned labels (remove 'bam')
    X_COL_LABEL = clean_label(x_col_raw)
    Y_COL_LABEL = clean_label(y_col_raw)

    # --- Clean numeric data ---
    data = df[[x_col_raw, y_col_raw]].copy()
    data[x_col_raw] = pd.to_numeric(data[x_col_raw], errors="coerce")
    data[y_col_raw] = pd.to_numeric(data[y_col_raw], errors="coerce")
    data = data.dropna()
    if data.empty:
        sys.exit("Error: No valid numeric rows after cleaning.")

    x = data[x_col_raw].to_numpy()
    y = data[y_col_raw].to_numpy()

    # --- Correlation (Pearson) ---
    r = np.corrcoef(x, y)[0, 1]
    print (r)
    r2 = r ** 2

    # --- Linear regression (least squares) ---
    m, b = np.polyfit(x, y, deg=1)

    # --- Plot ---
    fig = plt.figure(figsize=(3.3, 3.3))  # ~84mm square
    ax = plt.gca()

    # Scatter
    ax.scatter(x, y, s=6, alpha=0.6, linewidths=0, rasterized=(len(x) > 50000))

    # Regression line
    x_min, x_max = np.min(x), np.max(x)
    x_line = np.linspace(x_min, x_max, 200)
    y_line = m * x_line + b
    ax.plot(x_line, y_line, linewidth=1.2)

    # Labels & title
    ax.set_xlabel(X_COL_LABEL)
    ax.set_ylabel(Y_COL_LABEL)
    ax.set_title(f"{Y_COL_LABEL} vs {X_COL_LABEL}")

    # Ticks & spines
    ax.tick_params(direction="out", which="both")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # Limits with padding
    x_pad = 0.02 * (x_max - x_min if x_max > x_min else 1.0)
    y_min, y_max = np.min(y), np.max(y)
    y_pad = 0.02 * (y_max - y_min if y_max > y_min else 1.0)
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)

    # Annotation (R^2)
    ax.text(0.03, 0.97, f"$R^2$ = {r2:.4f}", transform=ax.transAxes, va="top", ha="left")

    # Save
    fig.tight_layout()
    fig.savefig(OUTPUT_FILE)  # PDF (vector)
    plt.close(fig)

    print(f"Using columns: X='{x_col_raw}' (label '{X_COL_LABEL}'), Y='{y_col_raw}' (label '{Y_COL_LABEL}')")
    print(f"Pearson r: {r:.6f} (R^2 = {r2:.6f})")
    print(f"Figure saved to: {Path(OUTPUT_FILE).resolve()}")

if __name__ == "__main__":
    main()

