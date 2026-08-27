"""
Cross-dataset dysregulation consistency for the 17 creatine-related genes,
built directly from the primary DE tables in tsvs/ (limma outputs already
filtered to significant probes per comparison — see README for cutoffs used
per dataset), NOT from the narrative prose in README.md (which paraphrases
some of this via generative AI and is not treated as a data source here).

Datasets covered (array-based, GEO):
  - GSE1563 / GDS724 kidney transplant biopsies: TX, AR, NR vs Control, plus
    pairwise AR-vs-TX, AR-vs-NR, TX-vs-NR.
  - GSE6344 clear renal cell carcinoma: Stage 1 and Stage 2 tumor vs normal.

Not included (no raw DE table checked into tsvs/, only prose + PNG in README):
  - GDS3712 nephrosclerosis (GSE20602)
  - GSE137570 CKD/fibrosis RNA-seq + eGFR subgroups
These would need to be pulled fresh from GEO if this analysis is extended.

Output: geo_dysregulation_17genes.csv (long format) +
        geo_dysregulation_summary.csv (per-gene consistency summary)
"""
import glob
import re
import pandas as pd
import numpy as np

GENES = ["SLC2A4", "IGF1", "GATM", "SLC6A8", "AKT1", "AKT2", "AKT3",
         "PSMB5", "PSMD3", "SGK1", "PRPS1", "CKB", "CKM", "CKMT1A",
         "CKMT1B", "CKMT2", "SLC16A12"]

FILES = {
    "TXvsControl": "tsvs/TXxN.tsv",
    "TXvsNR": "tsvs/TXxNR.tsv",
    "ARvsControl": "tsvs/ARxN.tsv",
    "ARvsNR": "tsvs/ARxNR.tsv",
    "ARvsTX": "tsvs/ARxTX.tsv",
    "ST1_Tumor_vs_Normal": "tsvs/ST1-TxN.tsv",
    "ST2_Tumor_vs_Normal": "tsvs/ST2-TxN.tsv",
}


def load_comparison(path, comparison_label):
    df = pd.read_csv(path, sep="\t")
    fc_col = [c for c in df.columns if c.startswith("log2(fold change)")][0]
    p_col = [c for c in df.columns if c.startswith("-log10(Pvalue)")][0]
    df = df.rename(columns={fc_col: "log2FC", p_col: "neglog10P"})
    df["comparison"] = comparison_label
    return df[["Gene.symbol", "log2FC", "neglog10P", "comparison"]]


def main():
    frames = [load_comparison(path, label) for label, path in FILES.items()]
    all_de = pd.concat(frames, ignore_index=True)

    hits = all_de[all_de["Gene.symbol"].isin(GENES)].copy()
    hits["direction"] = np.where(hits["log2FC"] > 0, "up", "down")
    # multiple probes per gene/comparison: keep the one with highest -log10P
    hits = hits.sort_values("neglog10P", ascending=False)
    hits_dedup = hits.drop_duplicates(subset=["Gene.symbol", "comparison"], keep="first")
    hits_dedup.to_csv("geo_dysregulation_17genes.csv", index=False)

    rows = []
    for gene in GENES:
        sub = hits_dedup[hits_dedup["Gene.symbol"] == gene]
        n_sig = len(sub)
        n_up = int((sub["direction"] == "up").sum())
        n_down = int((sub["direction"] == "down").sum())
        consistency = max(n_up, n_down) / n_sig if n_sig > 0 else np.nan
        mean_abs_fc = sub["log2FC"].abs().mean() if n_sig > 0 else np.nan
        rows.append({
            "gene": gene,
            "n_comparisons_significant": n_sig,
            "n_comparisons_total": len(FILES),
            "n_up": n_up,
            "n_down": n_down,
            "direction_consistency": consistency,   # 1.0 = always same direction
            "mean_abs_log2FC": mean_abs_fc,
            "comparisons": ";".join(sub["comparison"] + ":" + sub["direction"]),
        })

    summary = pd.DataFrame(rows).sort_values(
        ["n_comparisons_significant", "direction_consistency"], ascending=[False, False]
    ).reset_index(drop=True)
    summary.to_csv("geo_dysregulation_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
