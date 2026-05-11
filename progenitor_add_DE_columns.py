#!/usr/bin/env python3

import pandas as pd
import glob
import os

# -------------------------------------------------------
# arquivos de entrada
# -------------------------------------------------------

classification_file = "transcript_parental_classification.txt"

csv_folder = "/home/jnunes/Documents/epp/difexp/filtrados"

output_file = "transcript_parental_with_DE_flags.tsv"

# -------------------------------------------------------
# carregar classificação parental
# -------------------------------------------------------

df = pd.read_csv(classification_file, sep="\t")

# garantir nome da coluna
df.columns = ["Transcript", "Classification"]

# -------------------------------------------------------
# percorrer CSVs
# -------------------------------------------------------

csv_files = glob.glob(os.path.join(csv_folder, "*.csv"))

for csv_file in sorted(csv_files):

    # nome da coluna = nome do arquivo sem extensão
    col_name = os.path.basename(csv_file).replace(".csv", "")

    print(f"Processando: {col_name}")

    # carregar csv
    tmp = pd.read_csv(csv_file)

    # procurar coluna contendo transcript
    possible_cols = [
        c for c in tmp.columns
        if "transcript" in c.lower()
        or "gene" in c.lower()
        or "id" in c.lower()
    ]

    if len(possible_cols) == 0:
        print(f"ERRO: nenhuma coluna de transcript encontrada em {csv_file}")
        continue

    transcript_col = possible_cols[0]

    # conjunto de transcripts presentes no csv
    transcripts_in_csv = set(tmp[transcript_col].astype(str))

    # adicionar coluna binária
    df[col_name] = df["Transcript"].apply(
        lambda x: 1 if x in transcripts_in_csv else 0
    )

# -------------------------------------------------------
# salvar
# -------------------------------------------------------

df.to_csv(output_file, sep="\t", index=False)

print(f"\nArquivo salvo: {output_file}")
