#!/usr/bin/env python3

import pandas as pd
import glob
import os

# -------------------------------------------------------
# arquivos
# -------------------------------------------------------

classification_file = "transcript_parental_classification.txt"

csv_folder = "/home/jnunes/Documents/epp/difexp/filtrados"

output_file = "parental_UP_DOWN_summary.tsv"

# -------------------------------------------------------
# carregar classificação parental
# -------------------------------------------------------

class_df = pd.read_csv(classification_file, sep="\t")

# garantir nomes
class_df.columns = ["Transcript", "Classification"]

# -------------------------------------------------------
# preparar resultados
# -------------------------------------------------------

results = {}

# categorias
classes = ["Ss", "So", "Both", "NA"]

# -------------------------------------------------------
# percorrer csvs
# -------------------------------------------------------

csv_files = sorted(glob.glob(os.path.join(csv_folder, "*.csv")))

for csv_file in csv_files:

    comp = os.path.basename(csv_file).replace(".csv", "")

    print(f"Processando: {comp}")

    df = pd.read_csv(csv_file)

    # ---------------------------------------------------
    # detectar coluna transcript
    # ---------------------------------------------------

    transcript_col = None

    for c in df.columns:
        cl = c.lower()

        if "transcript" in cl or "gene" in cl or "id" in cl:
            transcript_col = c
            break

    if transcript_col is None:
        print(f"ERRO transcript column: {csv_file}")
        continue

    # ---------------------------------------------------
    # detectar logFC
    # ---------------------------------------------------

    logfc_col = None

    for c in df.columns:
        cl = c.lower()

        if "log2foldchange" in cl \
        or "logfc" in cl \
        or "log2fc" in cl:
            logfc_col = c
            break

    if logfc_col is None:
        print(f"ERRO logFC column: {csv_file}")
        continue

    # ---------------------------------------------------
    # merge
    # ---------------------------------------------------

    merged = pd.merge(
        class_df,
        df[[transcript_col, logfc_col]],
        left_on="Transcript",
        right_on=transcript_col,
        how="left"
    )

    # transcripts DE
    de_df = merged[~merged[logfc_col].isna()].copy()

    # UP / DOWN
    up_df = de_df[de_df[logfc_col] > 0]
    down_df = de_df[de_df[logfc_col] < 0]

    # ---------------------------------------------------
    # armazenar contagens
    # ---------------------------------------------------

    counts = {}

    for cls in classes:

        counts[f"{cls}_UP"] = (
            up_df["Classification"] == cls
        ).sum()

        counts[f"{cls}_DOWN"] = (
            down_df["Classification"] == cls
        ).sum()

    counts["Total_UP"] = len(up_df)
    counts["Total_DOWN"] = len(down_df)
    counts["Total_DE"] = len(de_df)

    results[comp] = counts

# -------------------------------------------------------
# converter para matriz
# -------------------------------------------------------

summary = pd.DataFrame(results)

# ordem das linhas
ordered_rows = []

for cls in classes:
    ordered_rows.append(f"{cls}_UP")
    ordered_rows.append(f"{cls}_DOWN")

ordered_rows += [
    "Total_UP",
    "Total_DOWN",
    "Total_DE"
]

summary = summary.loc[ordered_rows]

# -------------------------------------------------------
# salvar
# -------------------------------------------------------

summary.to_csv(output_file, sep="\t")

print(summary)

print(f"\nArquivo salvo: {output_file}")
