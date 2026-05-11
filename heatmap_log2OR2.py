#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact

# =====================================================
# INPUT
# =====================================================

input_file = "transcript_parental_with_DE_flags.tsv"

# =====================================================
# LOAD
# =====================================================

df = pd.read_csv(input_file, sep="\t")

de_cols = df.columns[2:]

classes = ["Ss", "So", "Both"]

# =====================================================
# CONTRAST LABELS
# =====================================================

label_map = {

    "DE_CB49_INOC_L_vs_CB49_CTRL_L_filtered":
        "Susceptible to L. xyly",

    "DE_IAC66_INOC_S_vs_IAC66_CTRL_S_filtered":
        "Susceptible to S. scitamineum",

    "DE_SP78_INOC_G_vs_SP78_CTRL_X_filtered":
        "Susceptible to G. diazotrophicus",

    "DE_SP78_INOC_XG_vs_SP78_CTRL_X_filtered":
        "Susceptible to G. diazotrophicus + X. albilineans",

    "DE_SP78_INOC_X_vs_SP78_CTRL_X_filtered":
        "Susceptible to X. albilineans",

    "DE_SP80_INOC_G_vs_SP80_CTRL_X_filtered":
        "Tolerant to G. diazotrophicus",

    "DE_SP80_INOC_L_vs_SP80_CTRL_L_filtered":
        "Tolerant to L. xyly",

    "DE_SP80_INOC_S_vs_SP80_CTRL_S_filtered":
        "Tolerant to S. scitamineum",

    "DE_SP80_INOC_XG_vs_SP80_CTRL_X_filtered":
        "Tolerant to G. diazotrophicus + X. albilineans",

    "DE_SP80_INOC_X_vs_SP80_CTRL_X_filtered":
        "Tolerant to X. albilineans"
}

pretty_labels = [label_map[col] for col in de_cols]

# =====================================================
# BACKGROUND
# =====================================================

total_transcripts = len(df)

background = {}

for cls in classes:

    background[cls] = (
        df["Classification"] == cls
    ).sum()

# =====================================================
# CALCULATE LOG2(OR) + PVALUES
# =====================================================

heatmap_data = []
pval_data = []

for cls in classes:

    row = []
    prow = []

    for col in de_cols:

        de_df = df[df[col] == 1]

        total_de = len(de_df)

        # contingency table

        a = (
            de_df["Classification"] == cls
        ).sum()

        b = background[cls] - a

        c = total_de - a

        d = (
            total_transcripts
            - background[cls]
            - c
        )

        table = [
            [a, b],
            [c, d]
        ]

        oddsratio, pvalue = fisher_exact(table)

        # avoid infinite
        if oddsratio == 0:
            oddsratio = 1e-10

        log2_or = np.log2(oddsratio)

        row.append(log2_or)
        prow.append(pvalue)

    heatmap_data.append(row)
    pval_data.append(prow)

# =====================================================
# DATAFRAMES
# =====================================================

heatmap_df = pd.DataFrame(
    heatmap_data,
    index=classes,
    columns=pretty_labels
)

pval_df = pd.DataFrame(
    pval_data,
    index=classes,
    columns=pretty_labels
)

print("\nLOG2(OR):")
print(heatmap_df)

print("\nPVALUES:")
print(pval_df)

# =====================================================
# COLOR SCALE
# =====================================================

max_abs = np.max(np.abs(heatmap_df.values))

# round upward
max_abs = np.ceil(max_abs * 10) / 10

print(f"\nUsing color scale: {-max_abs} to {max_abs}")

# =====================================================
# PLOT
# =====================================================

fig, ax = plt.subplots(figsize=(18, 5))

im = ax.imshow(
    heatmap_df.values,
    aspect='auto',
    cmap='bwr',
    vmin=-max_abs,
    vmax=max_abs
)

# =====================================================
# TICKS
# =====================================================

ax.set_xticks(np.arange(len(pretty_labels)))
ax.set_yticks(np.arange(len(classes)))

ax.set_xticklabels(
    pretty_labels,
    rotation=45,
    ha='right',
    fontsize=10,
    fontstyle='italic'
)

ax.set_yticklabels(
    classes,
    fontsize=12,
    fontweight='bold'
)

# =====================================================
# VALUES + SIGNIFICANCE
# =====================================================

for i in range(len(classes)):
    for j in range(len(pretty_labels)):

        value = heatmap_df.iloc[i, j]
        pval = pval_df.iloc[i, j]

        # significance stars

        if pval < 0.001:
            stars = "***"

        elif pval < 0.01:
            stars = "**"

        elif pval < 0.05:
            stars = "*"

        else:
            stars = ""

        label = f"{value:.2f}{stars}"

        ax.text(
            j,
            i,
            label,
            ha='center',
            va='center',
            fontsize=10,
            fontweight='bold',
            color='black'
        )

# =====================================================
# COLORBAR
# =====================================================

cbar = plt.colorbar(im)

cbar.set_label(
    "log2(Odds Ratio)",
    fontsize=12,
    fontweight='bold'
)

# =====================================================
# TITLE
# =====================================================

plt.title(
    "Subgenomic enrichment of DE lncRNAs",
    fontsize=15,
    fontweight='bold'
)

# =====================================================
# LAYOUT
# =====================================================

plt.tight_layout()

# =====================================================
# SAVE
# =====================================================

plt.savefig(
    "heatmap_log2OR.png",
    dpi=300,
    bbox_inches='tight'
)

plt.savefig(
    "heatmap_log2OR.svg",
    bbox_inches='tight'
)

# =====================================================
# SHOW
# =====================================================

plt.show()
