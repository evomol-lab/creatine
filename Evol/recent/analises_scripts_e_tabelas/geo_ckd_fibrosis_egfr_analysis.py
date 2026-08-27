"""
CKD fibrosis/eGFR continuous analysis (GSE137570, Doyle et al.) for the 17
creatine-related genes — pulled fresh from GEO supplementary normalized
counts (no raw DE table was checked into tsvs/ for this dataset).

IMPORTANT deviation from README: the README describes a binary split
("Fibrosis>50%", "EGFR>50% vs EGFR<50%"), but the script that produced those
groupings is not in this repo (RNAseq.R here is a GEO2R boilerplate for a
DIFFERENT dataset, GSE82291 — not reproducible for GSE137570). Rather than
guess the original threshold, this uses "Cohort 1" (n=24 samples with both
GFR and % tubulointerstitial fibrosis [TIF] recorded) and computes Spearman
correlation of each gene's normalized expression against GFR and against
TIF, continuously — no arbitrary cutoff invented.

Output: geo_ckd_fibrosis_egfr_17genes.csv
"""
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, false_discovery_control

GENES = ["SLC2A4", "IGF1", "GATM", "SLC6A8", "AKT1", "AKT2", "AKT3",
         "PSMB5", "PSMD3", "SGK1", "PRPS1", "CKB", "CKM", "CKMT1A",
         "CKMT1B", "CKMT2", "SLC16A12"]


def main():
    xl = pd.ExcelFile("geo_raw/GSE137570_counts.xlsx")
    df = xl.parse("COHORT 1_NORMALIZED READ COUNTS")
    df = df.set_index("SAMPLE ID")

    gfr = df.loc["GFR"].astype(float)
    tif = df.loc["TIF"].astype(float)

    rows = []
    for gene in GENES:
        if gene not in df.index:
            rows.append({"gene": gene, "present": False})
            continue
        expr = df.loc[gene].astype(float)
        rho_gfr, p_gfr = spearmanr(expr, gfr)
        rho_tif, p_tif = spearmanr(expr, tif)
        rows.append({
            "gene": gene, "present": True, "n_samples": len(expr),
            "spearman_rho_vs_GFR": rho_gfr, "pvalue_vs_GFR": p_gfr,
            "spearman_rho_vs_TIF_fibrosis": rho_tif, "pvalue_vs_TIF_fibrosis": p_tif,
        })

    out = pd.DataFrame(rows)
    present = out["present"]
    out.loc[present, "padj_vs_GFR_BH"] = false_discovery_control(out.loc[present, "pvalue_vs_GFR"])
    out.loc[present, "padj_vs_TIF_BH"] = false_discovery_control(out.loc[present, "pvalue_vs_TIF_fibrosis"])
    out = out.sort_values("pvalue_vs_TIF_fibrosis")
    out.to_csv("geo_ckd_fibrosis_egfr_17genes.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
