#!/usr/bin/env python3

import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency

# =====================================================
# arquivo de entrada
# =====================================================

input_file = "transcript_parental_with_DE_flags.tsv"

# =====================================================
# carregar dados
# =====================================================

df = pd.read_csv(input_file, sep="\t")

# colunas de comparações DE
de_cols = df.columns[2:]

# classes consideradas
# removi NA porque não possui DE
classes = ["Ss", "So", "Both"]

# =====================================================
# background global
# =====================================================

total_transcripts = len(df)

background = {}

for cls in classes:

    background[cls] = (
        df["Classification"] == cls
    ).sum()

print("\nBACKGROUND")
print("="*80)

for cls in classes:
    print(f"{cls}: {background[cls]}")

print(f"\nTotal transcripts: {total_transcripts}")

# =====================================================
# análises
# =====================================================

for col in de_cols:

    print("\n")
    print("="*80)
    print(col)
    print("="*80)

    # -------------------------------------------------
    # transcripts DE
    # -------------------------------------------------

    de_df = df[df[col] == 1]

    total_de = len(de_df)

    print(f"\nTotal DE: {total_de}")

    # -------------------------------------------------
    # contagens observadas
    # -------------------------------------------------

    observed = []

    print("\nObserved counts:")

    for cls in classes:

        count = (
            de_df["Classification"] == cls
        ).sum()

        observed.append(count)

        print(f"{cls}: {count}")

    # =================================================
    # CHI-SQUARE TEST
    # =================================================

    print("\n--- Chi-square test ---")

    # esperado proporcional ao background
    expected = []

    for cls in classes:

        exp = total_de * (
            background[cls] / total_transcripts
        )

        expected.append(exp)

    # tabela 2 x N
    table = [observed, expected]

    chi2, pval, dof, exp = chi2_contingency(table)

    print("\nExpected counts:")

    for cls, expv in zip(classes, expected):
        print(f"{cls}: {expv:.2f}")

    print(f"\nChi2 statistic : {chi2:.4f}")
    print(f"p-value        : {pval:.6e}")
    print(f"Degrees freedom: {dof}")

    if pval < 0.05:
        print("Result         : Significant")
    else:
        print("Result         : Not significant")

    # =================================================
    # FISHER EXACT TEST
    # =================================================

    print("\n--- Fisher exact tests ---")

    for cls in classes:

        # DE da classe
        a = (
            de_df["Classification"] == cls
        ).sum()

        # não-DE da classe
        b = background[cls] - a

        # DE do resto
        c = total_de - a

        # não-DE do resto
        d = (
            total_transcripts
            - background[cls]
            - c
        )

        contingency = [
            [a, b],
            [c, d]
        ]

        oddsratio, pvalue = fisher_exact(contingency)

        print("\n")
        print(f"Class: {cls}")

        print("Contingency table:")
        print(contingency)

        print(f"Odds ratio : {oddsratio:.4f}")
        print(f"p-value    : {pvalue:.6e}")

        # interpretação
        if pvalue < 0.05:

            if oddsratio > 1:
                print("Interpretation: Enriched among DE")

            elif oddsratio < 1:
                print("Interpretation: Depleted among DE")

            else:
                print("Interpretation: Neutral")

        else:
            print("Interpretation: Not significant")
