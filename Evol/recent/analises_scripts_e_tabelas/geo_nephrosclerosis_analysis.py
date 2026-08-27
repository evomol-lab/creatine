"""
Nephrosclerosis (GDS3712 / GSE20602) DE analysis for the 17 creatine-related
genes, pulled fresh from GEO (curated GDS SOFT file, which ships pre-mapped
gene symbols and transformed/log2 expression values) — this dataset had no
raw DE table checked into tsvs/, only prose+PNG in README.

Groups: TN (tumor-free control, n=4) vs NSC (nephrosclerosis, n=14).
Stats: Welch's t-test per gene (probe-level, log2 scale) + BH FDR correction
across the 17 genes tested. This is a lighter-weight substitute for the
original limma pipeline (not available here without R/Bioconductor), but
operates on the same underlying processed values GEO provides.
"""
import pandas as pd
import numpy as np
from scipy import stats

GENES = ["SLC2A4", "IGF1", "GATM", "SLC6A8", "AKT1", "AKT2", "AKT3",
         "PSMB5", "PSMD3", "SGK1", "PRPS1", "CKB", "CKM", "CKMT1A",
         "CKMT1B", "CKMT2", "SLC16A12"]

TN_SAMPLES = ["GSM517369", "GSM517370", "GSM517371", "GSM517372"]
NSC_SAMPLES = ["GSM517355", "GSM517356", "GSM517357", "GSM517358", "GSM517359",
               "GSM517360", "GSM517361", "GSM517362", "GSM517363", "GSM517364",
               "GSM517365", "GSM517366", "GSM517367", "GSM517368"]


def load_gds_table(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("!dataset_table_begin")) + 1
    end = next(i for i, l in enumerate(lines) if l.startswith("!dataset_table_end"))
    from io import StringIO
    table_txt = "".join(lines[start:end])
    return pd.read_csv(StringIO(table_txt), sep="\t")


def main():
    df = load_gds_table("geo_raw/GDS3712.soft")
    df = df[df["IDENTIFIER"].isin(GENES)].copy()

    rows = []
    for gene in GENES:
        sub = df[df["IDENTIFIER"] == gene]
        if sub.empty:
            rows.append({"gene": gene, "n_probes": 0})
            continue
        for _, probe in sub.iterrows():
            tn_vals = probe[TN_SAMPLES].astype(float).values
            nsc_vals = probe[NSC_SAMPLES].astype(float).values
            # log2FC: NSC vs TN (disease vs control), matching direction convention
            # used elsewhere in this analysis (positive = up in disease)
            log2fc = nsc_vals.mean() - tn_vals.mean()
            tstat, pval = stats.ttest_ind(nsc_vals, tn_vals, equal_var=False)
            rows.append({
                "gene": gene,
                "probe": probe["ID_REF"],
                "log2FC_NSCvsTN": log2fc,
                "pvalue": pval,
                "tn_mean": tn_vals.mean(),
                "nsc_mean": nsc_vals.mean(),
            })

    out = pd.DataFrame(rows)
    # BH FDR across the tested probes (17-gene-list scope, not genome-wide)
    valid = out["pvalue"].notna()
    from scipy.stats import false_discovery_control
    out.loc[valid, "padj_BH"] = false_discovery_control(out.loc[valid, "pvalue"])
    out = out.sort_values("pvalue")
    out.to_csv("geo_nephrosclerosis_17genes.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
