#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact

# =====================================================
# INPUT
# =====================================================

input_file = "transcript_parental_with_DE_direction.tsv"

# =====================================================
# BACKGROUND TOTALS
# =====================================================

background = {
    "Ss": 2833,
    "So": 10328,
    "Both": 9843
}

total_transcripts = 23084

classes = ["Ss", "So", "Both"]

# =====================================================
# LOAD
# =====================================================

df = pd.read_csv(
    input_file,
    sep="\t",
    index_col=0
)

# keep full table
full_df = df.copy()

# only biological classes
df_classes = df.loc[
    ~df.index.str.contains("Total")
]

# =====================================================
# LABELS
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

contrasts = list(df.columns)

pretty_labels = [
    label_map[x]
    for x in contrasts
]

# =====================================================
# FUNCTION
# =====================================================

def calculate_matrix(direction):

    heatmap = []
    pvals = []

    for cls in classes:

        row = []
        prow = []

        for contrast in contrasts:

            # observed DE transcripts

            a = df_classes.loc[
                f"{cls}_{direction}",
                contrast
            ]

            # total UP or DOWN transcripts

            total_direction = full_df.loc[
                f"Total_{direction}",
                contrast
            ]

            # contingency table

            b = background[cls] - a

            c = total_direction - a

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

            # avoid zero OR

            if oddsratio == 0:
                oddsratio = 1e-10

            log2or = np.log2(oddsratio)

            row.append(log2or)
            prow.append(pvalue)

        heatmap.append(row)
        pvals.append(prow)

    heatmap_df = pd.DataFrame(
        heatmap,
        index=classes,
        columns=pretty_labels
    )

    pval_df = pd.DataFrame(
        pvals,
        index=classes,
        columns=pretty_labels
    )

    return heatmap_df, pval_df

# =====================================================
# BUILD MATRICES
# =====================================================

up_df, up_p = calculate_matrix("UP")

down_df, down_p = calculate_matrix("DOWN")

# =====================================================
# PLOT FUNCTION
# =====================================================

def plot_heatmap(data, pvals, title, outfile):

    vmax = np.max(np.abs(data.values))

    vmax = np.ceil(vmax * 10) / 10

    fig, ax = plt.subplots(figsize=(18, 5))

    im = ax.imshow(
        data.values,
        cmap="bwr",
        aspect="auto",
        vmin=-vmax,
        vmax=vmax
    )

    # ==========================================
    # ticks
    # ==========================================

    ax.set_xticks(np.arange(len(data.columns)))
    ax.set_yticks(np.arange(len(data.index)))

    ax.set_xticklabels(
        data.columns,
        rotation=45,
        ha="right",
        fontsize=10,
        fontstyle="italic"
    )

    ax.set_yticklabels(
        data.index,
        fontsize=12,
        fontweight="bold"
    )

    # ==========================================
    # values + significance
    # ==========================================

    for i in range(len(data.index)):
        for j in range(len(data.columns)):

            value = data.iloc[i, j]

            pval = pvals.iloc[i, j]

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
                ha="center",
                va="center",
                fontsize=10,
                fontweight="bold",
                color="black"
            )

    # ==========================================
    # colorbar
    # ==========================================

    cbar = plt.colorbar(im)

    cbar.set_label(
        "log2(Odds Ratio)",
        fontsize=12,
        fontweight="bold"
    )

    # ==========================================
    # title
    # ==========================================

    plt.title(
        title,
        fontsize=15,
        fontweight="bold"
    )

    plt.tight_layout()

    # ==========================================
    # save
    # ==========================================

    plt.savefig(
        f"{outfile}.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.savefig(
        f"{outfile}.svg",
        bbox_inches="tight"
    )

# =====================================================
# GENERATE FIGURES
# =====================================================

plot_heatmap(
    up_df,
    up_p,
    "Subgenomic enrichment of UP-regulated lncRNAs",
    "heatmap_UP_log2OR"
)

plot_heatmap(
    down_df,
    down_p,
    "Subgenomic enrichment of DOWN-regulated lncRNAs",
    "heatmap_DOWN_log2OR"
)

# =====================================================
# DONE
# =====================================================

print("\nDone.")
print("Generated files:")
print("- heatmap_UP_log2OR.png")
print("- heatmap_UP_log2OR.svg")
print("- heatmap_DOWN_log2OR.png")
print("- heatmap_DOWN_log2OR.svg")
