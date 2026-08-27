"""
Tissue-specificity (tau) for the 17 creatine-related genes, from GTEx v8 median TPM.

Pulls median gene expression per tissue directly from the public GTEx Portal API
(dataset gtex_v8, matching the release used for the GTex-*.svg figures already in
this repo), then computes:

  - tau (Yanai et al. 2005, Bioinformatics): tissue-specificity index, 0 = ubiquitous,
    1 = restricted to a single tissue.
  - kidney_mean_tpm: mean of Kidney_Cortex and Kidney_Medulla median TPM.
  - kidney_rank_pct: percentile rank of kidney expression among all profiled tissues
    for that gene (1.0 = kidney is the single highest-expressing tissue).

Output: gtex_tau_17genes.csv + gtex_tau_barplot.png
"""
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

GENES = ["SLC2A4", "IGF1", "GATM", "SLC6A8", "AKT1", "AKT2", "AKT3",
         "PSMB5", "PSMD3", "SGK1", "PRPS1", "CKB", "CKM", "CKMT1A",
         "CKMT1B", "CKMT2", "SLC16A12"]

DATASET = "gtex_v8"
BASE = "https://gtexportal.org/api/v2"


def fetch_gencode_ids(genes):
    params = [("geneId", g) for g in genes]
    params += [("itemsPerPage", 50), ("format", "json")]
    r = requests.get(f"{BASE}/reference/gene", params=params, timeout=30)
    r.raise_for_status()
    data = r.json()["data"]
    id_map = {d["geneSymbol"]: d["gencodeId"] for d in data}
    missing = set(genes) - set(id_map)
    if missing:
        raise RuntimeError(f"gencodeId not found for: {missing}")
    return id_map


def fetch_median_expression(gencode_ids, dataset=DATASET):
    params = [("gencodeId", g) for g in gencode_ids]
    params += [("datasetId", dataset), ("itemsPerPage", 2000), ("format", "json")]
    r = requests.get(f"{BASE}/expression/medianGeneExpression", params=params, timeout=60)
    r.raise_for_status()
    return pd.DataFrame(r.json()["data"])


def tau(values):
    """Yanai et al. 2005 tissue-specificity index. values: array of non-negative TPMs."""
    x = np.asarray(values, dtype=float)
    xmax = x.max()
    if xmax <= 0:
        return np.nan
    return float(np.sum(1 - x / xmax) / (len(x) - 1))


def main():
    id_map = fetch_gencode_ids(GENES)
    df = fetch_median_expression(list(id_map.values()))
    df["geneSymbol"] = df["geneSymbol"].astype(str)

    rows = []
    for gene in GENES:
        sub = df[df["geneSymbol"] == gene]
        if sub.empty:
            print(f"WARNING: no expression data returned for {gene}")
            continue
        vals = sub.set_index("tissueSiteDetailId")["median"]
        t = tau(vals.values)

        kidney_tissues = [c for c in vals.index if c.startswith("Kidney")]
        kidney_mean = vals[kidney_tissues].mean() if kidney_tissues else np.nan
        kidney_rank_pct = float((vals <= kidney_mean).sum()) / len(vals) if kidney_tissues else np.nan

        rows.append({
            "gene": gene,
            "n_tissues": len(vals),
            "tau": t,
            "kidney_mean_tpm": kidney_mean,
            "max_tissue": vals.idxmax(),
            "max_tpm": vals.max(),
            "kidney_rank_pct": kidney_rank_pct,
        })

    out = pd.DataFrame(rows).sort_values("tau", ascending=False).reset_index(drop=True)
    out.to_csv("gtex_tau_17genes.csv", index=False)
    print(out.to_string(index=False))

    fig, ax = plt.subplots(figsize=(9, 6))
    colors = ["#2ecc71" if k else "#3498db" for k in (out["max_tissue"].str.startswith("Kidney"))]
    ax.barh(out["gene"], out["tau"], color=colors, edgecolor="black")
    ax.set_xlabel("tau (tissue-specificity index)")
    ax.set_title("Tissue specificity (GTEx v8) — 17 creatine-related genes\n(green = kidney is the top-expressing tissue)")
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig("gtex_tau_barplot.png", dpi=200)
    print("\nSaved: gtex_tau_17genes.csv, gtex_tau_barplot.png")


if __name__ == "__main__":
    main()
