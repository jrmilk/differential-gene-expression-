#!/usr/bin/env python3

import pandas as pd

# -------------------------------------------------------
# entrada
# -------------------------------------------------------

input_file = "transcript_parental_with_DE_flags.tsv"

# saída
output_file = "parental_DE_summary_matrix.tsv"

# -------------------------------------------------------
# carregar tabela
# -------------------------------------------------------

df = pd.read_csv(input_file, sep="\t")

# colunas DE
de_cols = df.columns[2:]

# categorias
categories = ["Ss", "So", "Both", "NA"]

# -------------------------------------------------------
# construir matriz resumo
# -------------------------------------------------------

summary = {}

for col in de_cols:

    counts = {}

    # DE
    de_df = df[df[col] == 1]

    # not DE
    nde_df = df[df[col] == 0]

    for cat in categories:

        counts[f"{cat}_DE"] = (de_df["Classification"] == cat).sum()

        counts[f"{cat}_notDE"] = (nde_df["Classification"] == cat).sum()

    # totais
    counts["Total_DE"] = len(de_df)
    counts["Total_notDE"] = len(nde_df)
    counts["Total"] = len(df)

    summary[col] = counts

# -------------------------------------------------------
# converter para dataframe
# -------------------------------------------------------

summary_df = pd.DataFrame(summary)

# ordem desejada das linhas
ordered_rows = []

for cat in categories:
    ordered_rows.append(f"{cat}_DE")
    ordered_rows.append(f"{cat}_notDE")

ordered_rows += [
    "Total_DE",
    "Total_notDE",
    "Total"
]

summary_df = summary_df.loc[ordered_rows]

# -------------------------------------------------------
# salvar
# -------------------------------------------------------

summary_df.to_csv(output_file, sep="\t")

print(summary_df)

print(f"\nArquivo salvo: {output_file}")
