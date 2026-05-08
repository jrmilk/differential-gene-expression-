import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import glob
import os

# ── Thresholds ────────────────────────────────────────────────────────────────
fdr_thresh   = 0.01   # FDR / padj máximo
logfc_thresh = 1.0    # |log2FoldChange| mínimo (equivale a 2× fold change)
top_n_labels = 10     # quantidade de genes anotados no plot

# ── Colunas obrigatórias ───────────────────────────────────────────────────────
REQUIRED_COLS = {"padj", "log2FoldChange"}


def classify_gene(lfc: float, padj: float,
                  logfc_thresh: float, fdr_thresh: float) -> str:
    """Classifica um gene em up-regulated, down-regulated ou não significativo."""
    if padj >= fdr_thresh:
        return "ns"
    if lfc >= logfc_thresh:
        return "up"
    if lfc <= -logfc_thresh:
        return "down"
    return "ns"


COLOR_MAP = {
    "up":   "#E05C5C",   # vermelho — up-regulated
    "down": "#4A90D9",   # azul     — down-regulated
    "ns":   "#CCCCCC",   # cinza    — não significativo
}

ZORDER_MAP = {"up": 3, "down": 3, "ns": 1}


# ── Loop principal ────────────────────────────────────────────────────────────
csv_files = glob.glob("DE_*.csv")
if not csv_files:
    print("Nenhum arquivo DE_*.csv encontrado.")

for csv in csv_files:
    print(f"\nProcessando: {csv}")

    df = pd.read_csv(csv)

    # Verifica colunas obrigatórias
    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        print(f"  [AVISO] Colunas ausentes {missing} — pulando arquivo.")
        continue

    # Remove linhas com NA nas colunas essenciais
    df = df.dropna(subset=["padj", "log2FoldChange"])

    # Nome da coluna de ID (primeira coluna ou 'gene_id' se existir)
    id_col = "gene_id" if "gene_id" in df.columns else df.columns[0]

    # ── Classifica cada gene ──────────────────────────────────────────────────
    df["_class"] = df.apply(
        lambda r: classify_gene(r["log2FoldChange"], r["padj"],
                                logfc_thresh, fdr_thresh),
        axis=1,
    )

    n_up   = (df["_class"] == "up").sum()
    n_down = (df["_class"] == "down").sum()
    n_ns   = (df["_class"] == "ns").sum()
    print(f"  Up: {n_up}  |  Down: {n_down}  |  NS: {n_ns}")

    df["_color"]  = df["_class"].map(COLOR_MAP)
    df["_zorder"] = df["_class"].map(ZORDER_MAP)
    df["_neglog10padj"] = -np.log10(df["padj"].clip(lower=1e-300))

    # ── Plot ─────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))

    for grp in ["ns", "up", "down"]:
        sub = df[df["_class"] == grp]
        ax.scatter(
            sub["log2FoldChange"],
            sub["_neglog10padj"],
            c=sub["_color"],
            alpha=0.75 if grp == "ns" else 0.90,
            edgecolors="none",
            s=14 if grp == "ns" else 22,
            zorder=ZORDER_MAP[grp],
            rasterized=True,
        )

    # Linhas de threshold
    ax.axhline(-np.log10(fdr_thresh), color="#888888",
               linestyle="--", linewidth=0.9, zorder=2)
    ax.axvline( logfc_thresh,  color="#888888",
               linestyle="--", linewidth=0.9, zorder=2)
    ax.axvline(-logfc_thresh,  color="#888888",
               linestyle="--", linewidth=0.9, zorder=2)

    # ── Anotação dos top N genes mais significativos ──────────────────────────
    sig = df[df["_class"].isin(["up", "down"])].copy()
    if not sig.empty:
        sig_sorted = sig.nsmallest(top_n_labels, "padj")
        for _, row in sig_sorted.iterrows():
            ax.annotate(
                str(row[id_col]),
                xy=(row["log2FoldChange"], row["_neglog10padj"]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=6.5,
                color="#333333",
                zorder=4,
            )

    # ── Legenda ───────────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(color=COLOR_MAP["up"],
                       label=f"Up-regulated (n={n_up})"),
        mpatches.Patch(color=COLOR_MAP["down"],
                       label=f"Down-regulated (n={n_down})"),
        mpatches.Patch(color=COLOR_MAP["ns"],
                       label=f"Not significant (n={n_ns})"),
    ]
    ax.legend(handles=legend_patches, fontsize=9, framealpha=0.7,
              loc="upper left")

    # ── Rótulos e título ──────────────────────────────────────────────────────
    ax.set_title(
        f"Volcano Plot — {os.path.basename(csv)}\n"
        f"FDR < {fdr_thresh}  |  |log₂FC| ≥ {logfc_thresh}",
        fontsize=12, pad=12,
    )
    ax.set_xlabel("log₂ Fold Change", fontsize=11)
    ax.set_ylabel("-log₁₀(FDR / padj)", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)

    # ── Salvar ────────────────────────────────────────────────────────────────
    outname = csv.replace(".csv", "_volcano_fdr.png")
    fig.savefig(outname, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Salvo: {outname}")

print("\nConcluído!")
