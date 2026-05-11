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
# CALCULATE LOG2(OR)
# =====================================================

heatmap_data = []

for cls in classes:

    row = []

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

    heatmap_data.append(row)

# =====================================================
# DATAFRAME
# =====================================================

heatmap_df = pd.DataFrame(
    heatmap_data,
    index=classes,
    columns=pretty_labels
)

print(heatmap_df)

# =====================================================
# PLOT
# =====================================================

fig, ax = plt.subplots(figsize=(16, 5))

im = ax.imshow(
    heatmap_df.values,
    aspect='auto',
    cmap='bwr',
    vmin=-0.4,
    vmax=0.4
)

# ticks

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

# values inside cells

for i in range(len(classes)):
    for j in range(len(pretty_labels)):

        value = heatmap_df.iloc[i, j]

        ax.text(
            j,
            i,
            f"{value:.2f}",
            ha='center',
            va='center',
            fontsize=10,
            fontweight='bold',
            color='black'
        )

# colorbar

cbar = plt.colorbar(im)

cbar.set_label(
    "log2(Odds Ratio)",
    fontsize=12,
    fontweight='bold'
)

# title

plt.title(
    "Subgenomic enrichment of DE lncRNAs",
    fontsize=14,
    fontweight='bold'
)

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

plt.show()
