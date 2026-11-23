# 06_score_exhaustion.py: Freeze CD8 exhaustion signatures (discovery) and score cells via pyUCell.
# Usage: python scripts/06_score_exhaustion.py [metadata/signatures.yaml]
#
# 1. Load TEXP-42 genes/rules from signatures.yaml.
# 2. Flag genes with background noise (myeloid or hepatocyte expression >= CD8 in discovery).
# 3. Save 4 frozen signature sets (V1_full42, V2_noprolif, V3_noprolif_lineage, V4_compact6) before scoring.
# 4. Score CD8 cells across all 4 variants + 2 ambient-RNA scores and save results.

import datetime
import os
import sys
from importlib.metadata import version

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml
from pyucell import compute_ucell_scores

cfg_path = sys.argv[1] if len(sys.argv) > 1 else "metadata/signatures.yaml"
cfg = yaml.safe_load(open(cfg_path))
P = cfg["paths"]
os.makedirs(P["out_tables"], exist_ok=True)
os.makedirs(os.path.dirname(P["frozen_out"]) or ".", exist_ok=True)

texp = list(cfg["texp42"]["genes"])
assert len(texp) == len(set(texp)), "duplicate genes in texp42 list"
if "expected_n" in cfg["texp42"]:
    assert len(texp) == cfg["texp42"]["expected_n"]
prolif = set(cfg["texp42"]["proliferation_genes"])
assert prolif <= set(texp), f"proliferation genes not in the signature: {prolif - set(texp)}"
compact = list(cfg["compact6"]["genes"])
ratio = float(cfg["lineage_filter"]["max_nonT_to_CD8_ratio"])
print(f"TEXP signature genes: {len(texp)} | proliferation-annotated: {len(prolif)}")

# 1. lineage means in the discovery cohort 
cd8_names = ad.read_h5ad(P["cd8_clean"], backed="r").obs_names
mye_names = ad.read_h5ad(P["myeloid_clean"], backed="r").obs_names

raw = sc.read_h5ad(P["raw_counts"])
X = raw.X.tocsr() if sp.issparse(raw.X) else sp.csr_matrix(raw.X)
total = np.asarray(X.sum(axis=1)).ravel()
present = [g for g in texp if g in raw.var_names]
missing = [g for g in texp if g not in raw.var_names]
if missing:
    print("WARNING: signature genes absent from the data:", missing)
idx = [raw.var_names.get_loc(g) for g in present]
expr = np.log1p(X[:, idx].toarray() / total[:, None] * 1e4)      # same recipe as normalize_total(1e4) + log1p
obs_index = pd.Index(raw.obs_names)
hep_names = raw.obs_names[(raw.obs[P["coarse_col"]] == "Hepatocyte").to_numpy()]
groups = {"CD8_clean": cd8_names, "Myeloid_clean": mye_names, "Hepatocyte": hep_names}
means = {}
for name, bc in groups.items():
    rows = obs_index.get_indexer(bc)
    assert (rows >= 0).all(), f"{name}: barcodes missing from the raw object"
    means[name] = expr[rows].mean(axis=0)
lin = pd.DataFrame(means, index=present)
lin["nonT_max"] = lin[["Myeloid_clean", "Hepatocyte"]].max(axis=1)
lin["ratio_nonT_to_CD8"] = lin["nonT_max"] / lin["CD8_clean"].clip(lower=1e-3)
min_nonT = float(cfg["lineage_filter"].get("min_nonT_mean", 0.05))   # ignore noise-level differences
lin["lineage_ambiguous"] = (lin["nonT_max"] >= ratio * lin["CD8_clean"]) & (lin["nonT_max"] >= min_nonT)
lin["low_in_CD8"] = lin["CD8_clean"] < min_nonT
lin["proliferation_annotated"] = [g in prolif for g in lin.index]
lin = lin.round(3)
lin.to_csv(os.path.join(P["out_tables"], "TEXP42_gene_lineage_means_GSE149614.csv"))
print("Gene lineage means (log-normalised), discovery cohort:")
print(lin.sort_values("ratio_nonT_to_CD8", ascending=False).to_string())

# 2. variants, frozen before scoring 
ambiguous = set(lin.index[lin["lineage_ambiguous"]])
v2 = [g for g in texp if g not in prolif]
variants = {
    "V1_full42": texp,
    "V2_noprolif": v2,
    "V3_noprolif_lineage": [g for g in v2 if g not in ambiguous],
    "V4_compact6": compact,
}
frozen = {
    "frozen_on": datetime.date.today().isoformat(),
    "discovery_cohort": "GSE149614",
    "rule": f"lineage-ambiguous if max(mean_myeloid, mean_hepatocyte) >= {ratio} x mean_CD8 (log-normalised)",
    "primary_variant_proposed": cfg.get("primary_variant"),
    "dropped_proliferation": sorted(prolif),
    "dropped_lineage_ambiguous": sorted(ambiguous - prolif),
    "variants": {k: list(v) for k, v in variants.items()},
    "ambient": cfg["ambient"],
}
with open(P["frozen_out"], "w") as fh:
    yaml.safe_dump(frozen, fh, sort_keys=False)
print(f"Frozen signatures written to {P['frozen_out']}")
for k, v in variants.items():
    print(f"  {k}: {len(v)} genes")
print("  dropped as proliferation:", sorted(prolif))
print("  dropped as lineage-ambiguous (excluding proliferation genes):", sorted(ambiguous - prolif))
del raw, X, expr

# 3. scoring
cd8 = ad.read_h5ad(P["cd8_clean"])
for col in ("flag_cd8_expressed", "flag_cd3d_negative", "flag_cd4_only", "cluster_label", "patient", "site"):
    assert col in cd8.obs, f"{col} missing from the CD8 object"
sigs = {}
for name, genes in variants.items():
    keep = [g for g in genes if g in cd8.var_names]
    if len(keep) < len(genes):
        print(f"WARNING: {name}: {len(genes) - len(keep)} genes not in the CD8 object")
    sigs[name] = keep
for name, genes in cfg["ambient"].items():
    sigs[f"ambient_{name}"] = [g for g in genes if g in cd8.var_names]

compute_ucell_scores(cd8, sigs, max_rank=1500, missing_genes="skip", n_jobs=int(cfg.get("n_jobs", 4)), device="cpu")
score_cols = [f"{n}_UCell" for n in sigs]
cd8.obs["gate_strict"] = (cd8.obs["flag_cd8_expressed"] & ~cd8.obs["flag_cd3d_negative"] & ~cd8.obs["flag_cd4_only"])
cd8.uns["signatures_frozen"] = yaml.safe_dump(frozen)
cd8.uns["pyucell_version"] = version("pyucell")

# 4. summaries
pd.set_option("display.width", 200)
print(f"Cells: {cd8.n_obs:,} | strict gate: {int(cd8.obs['gate_strict'].sum()):,}")
print("Median UCell score by cluster (all cells):")
print(cd8.obs.groupby("cluster_label", observed=True)[score_cols].median().round(3).to_string())
print("Median UCell score by site (strict gate):")
print(cd8.obs[cd8.obs["gate_strict"]].groupby("site", observed=True)[score_cols].median().round(3).to_string())

tum = cd8.obs[cd8.obs["site"].isin(["Tumor", "PVTT"]) & cd8.obs["gate_strict"]]
print("Strict-gate CD8 cells in Tumor + PVTT per patient (power check):")
print(tum.groupby("patient", observed=True).size().to_string())
print("Spearman correlation, Tumor + PVTT strict gate (ambient scores vs variants):")
print(tum[score_cols].corr(method="spearman").loc[[c for c in score_cols if "ambient" in c],
                                                    [c for c in score_cols if "ambient" not in c]].round(2).to_string())

cd8.obs[["patient", "site", "cluster_label", "gate_strict"] + score_cols].to_csv(
    os.path.join(P["out_tables"], "exhaustion_scores_GSE149614.csv"))
cd8.write_h5ad(P["out_scored"], compression="gzip")
print(f"Saved scored object to {P['out_scored']} and per-cell scores to {P['out_tables']}/")