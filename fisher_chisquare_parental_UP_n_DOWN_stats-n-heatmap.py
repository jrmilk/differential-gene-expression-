# =====================================================
# PRINT STATISTICS
# =====================================================

from scipy.stats import chi2_contingency

def print_statistics(direction):

    print("\n")
    print("=" * 80)
    print(f"{direction}-regulated lncRNAs")
    print("=" * 80)

    for contrast in contrasts:

        print("\n")
        print("=" * 80)
        print(contrast)
        print("=" * 80)

        total_de = full_df.loc[
            f"Total_{direction}",
            contrast
        ]

        print(f"\nTotal {direction}: {total_de}\n")

        observed = []

        print("Observed counts:")

        for cls in classes:

            val = df_classes.loc[
                f"{cls}_{direction}",
                contrast
            ]

            observed.append(val)

            print(f"{cls}: {val}")

        # ==========================================
        # Expected
        # ==========================================

        expected = []

        for cls in classes:

            exp = (
                background[cls]
                / total_transcripts
            ) * total_de

            expected.append(exp)

        print("\n--- Chi-square test ---\n")

        print("Expected counts:")

        for cls, exp in zip(classes, expected):

            print(f"{cls}: {exp:.2f}")

        # contingency table

        chi2_table = []

        for cls in classes:

            a = df_classes.loc[
                f"{cls}_{direction}",
                contrast
            ]

            b = total_de - a

            chi2_table.append([a, b])

        chi2, pval, dof, exp = chi2_contingency(
            chi2_table
        )

        print(f"\nChi2 statistic : {chi2:.4f}")
        print(f"p-value        : {pval:.6e}")
        print(f"Degrees freedom: {dof}")

        if pval < 0.05:
            print("Result         : Significant")
        else:
            print("Result         : Not significant")

        # ==========================================
        # Fisher tests
        # ==========================================

        print("\n--- Fisher exact tests ---\n")

        for cls in classes:

            a = df_classes.loc[
                f"{cls}_{direction}",
                contrast
            ]

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

            oddsratio, fisher_p = fisher_exact(table)

            print(f"\nClass: {cls}")

            print("Contingency table:")
            print(table)

            print(f"Odds ratio : {oddsratio:.4f}")
            print(f"p-value    : {fisher_p:.6e}")

            if fisher_p < 0.05:

                if oddsratio > 1:
                    interpretation = "Enriched among DE"

                else:
                    interpretation = "Depleted among DE"

            else:
                interpretation = "Not significant"

            print(f"Interpretation: {interpretation}")
