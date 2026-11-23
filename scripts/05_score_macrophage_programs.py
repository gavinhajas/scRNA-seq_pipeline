# 05_score_macrophage_programs.py: Measure macrophage program activity per cell and patient.
# Usage: python scripts/05_score_macrophage_programs.py [metadata/macrophage_programs.yaml]
#
# - Scores each cell across 3 overlapping programs (FOLR2_C1Q, SPP1_MMP, TREM2_GPNMB).
#   Creates a single "top program" label for charts only—never use it for stats.
# - Saves raw ("<program>") and cleaned ("<program>__filtered") scores. Cleaned scores drop
#   genes likely coming from hepatocyte or T/NK background noise (e.g., SEPP1/SELENOP), locked in before validation.
# - Outputs average macrophage scores and cell counts per patient/site.

import datetime
import os
import sys
from importlib.metadata import version

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from pyucell import compute_ucell_scores

cfg_path = sys.argv[1] if len(sys.argv) > 1 else "metadata/macrophage_programs.yaml"
cfg = yaml.safe_load(open(cfg_path))
P = cfg["paths"]
os.makedirs(P["out_tables"], exist_ok=True)
tag = P.get("dataset_tag", "dataset")

a = ad.read_h5ad(P["myeloid_clean"])
for col in ("patient", "site", "cluster_label"):
    assert col in a.obs, f"{col} missing from the myeloid object"
print(f"{a.n_obs:,} myeloid cells x {a.n_vars:,} genes")


def resolve(genes):
    used, missing = [], []
    for g in genes:
        found = next((x for x in [g] + list(cfg.get("aliases", {}).get(g, [])) if x in a.var_names), None)
        (used if found else missing).append(found or g)
        if found and found != g:
            print(f"  note: {g} -> {found}")
    return used, missing


programs = {}
for name, genes in cfg["programs"].items():
    used, missing = resolve(genes)
    print(f"{name}: using {used}" + (f" | MISSING {missing}" if missing else ""))
    assert len(used) >= 2, f"{name}: fewer than 2 genes available"
    programs[name] = used
amb, amb_missing = resolve(cfg["ambient_hepatocyte"])
if amb_missing:
    print("ambient genes missing:", amb_missing)

# macrophage-like cells (frozen discovery clusters)
a.obs["macrophage_like"] = a.obs["cluster_label"].astype(str).isin(cfg["macrophage_cluster_labels"])
mac = a.obs["macrophage_like"].to_numpy()
assert mac.sum() > 0, "no macrophage-like cells: check cluster_label values"
print(f"macrophage-like cells: {int(mac.sum()):,} of {a.n_obs:,}")

# lineage check in the discovery cohort
lc = cfg["lineage_filter"]
raw = ad.read_h5ad(P["raw_counts"])
X = raw.X.tocsr() if sp.issparse(raw.X) else sp.csr_matrix(raw.X)
total = np.asarray(X.sum(axis=1)).ravel()
all_genes = sorted({g for gs in programs.values() for g in gs})
expr = np.log1p(X[:, [raw.var_names.get_loc(g) for g in all_genes]].toarray() / total[:, None] * 1e4)
names = pd.Index(raw.obs_names)
groups = {"macrophage_like": a.obs_names[mac]}
for lab in lc["nonmac_coarse_labels"]:
    groups[lab] = raw.obs_names[(raw.obs[P["coarse_col"]] == lab).to_numpy()]
means = {}
for k, bc in groups.items():
    rows = names.get_indexer(bc)
    assert (rows >= 0).all(), f"{k}: barcodes missing from the raw object"
    means[k] = expr[rows].mean(axis=0)
lin = pd.DataFrame(means, index=all_genes)
nonmac = lin[lc["nonmac_coarse_labels"]].max(axis=1)
lin["nonmac_max"] = nonmac
lin["ratio_nonmac_to_mac"] = nonmac / lin["macrophage_like"].clip(lower=1e-3)
lin["lineage_ambiguous"] = (nonmac >= float(lc["max_nonmac_to_mac_ratio"]) * lin["macrophage_like"]) & \
                           (nonmac >= float(lc["min_nonmac_mean"]))
lin = lin.round(3)
lin.to_csv(os.path.join(P["out_tables"], f"macrophage_program_gene_lineage_means_{tag}.csv"))
pd.set_option("display.width", 200)
print("Program gene lineage means (log-normalised), discovery cohort:")
print(lin.sort_values("ratio_nonmac_to_mac", ascending=False).to_string())

filtered = {}
for name, genes in programs.items():
    kept = [g for g in genes if not lin.loc[g, "lineage_ambiguous"]]
    if len(kept) < 2:
        print(f"WARNING: {name}: fewer than 2 genes survive the lineage filter; the filtered set keeps the original genes")
        kept = list(genes)
    filtered[name] = kept
frozen = {"frozen_on": datetime.date.today().isoformat(), "discovery_cohort": tag,
          "rule": f"lineage-ambiguous if max(mean in {lc['nonmac_coarse_labels']}) >= {lc['max_nonmac_to_mac_ratio']} x mean in "
                  f"macrophage-like cells and >= {lc['min_nonmac_mean']}",
          "approved_sets": programs, "filtered_sets": filtered,
          "dropped": {n: sorted(set(programs[n]) - set(filtered[n])) for n in programs}}
with open(P["frozen_out"], "w") as fh:
    yaml.safe_dump(frozen, fh, sort_keys=False)
print(f"Frozen program sets written to {P['frozen_out']}")
for n in programs:
    print(f"  {n}: approved {programs[n]} -> filtered {filtered[n]}")
del raw, X, expr

# scoring 
sigs = dict(programs)
sigs.update({f"{n}__filtered": g for n, g in filtered.items()})
sigs["ambient_hepatocyte"] = amb
compute_ucell_scores(a, sigs, max_rank=1500, missing_genes="skip", n_jobs=int(cfg.get("n_jobs", 4)), device="cpu")
pcols = [f"{p}_UCell" for p in programs]
fcols = [f"{p}__filtered_UCell" for p in programs]
acol = "ambient_hepatocyte_UCell"

# display-only label from the approved sets
z = (a.obs.loc[mac, pcols] - a.obs.loc[mac, pcols].mean()) / a.obs.loc[mac, pcols].std(ddof=0)
top = z.idxmax(axis=1).str.replace("_UCell", "", regex=False)
top[z.max(axis=1) < float(cfg["display_min_z"])] = "mixed_or_low"
a.obs["display_label"] = "not_macrophage_like"
a.obs.loc[mac, "display_label"] = top

# summaries 
print("Median program score by discovery cluster (approved sets, then filtered sets):")
print(a.obs.groupby("cluster_label", observed=True)[pcols + [acol]].median().round(3).to_string())
print(a.obs.groupby("cluster_label", observed=True)[fcols].median().round(3).to_string())
print("Spearman correlation between programs (macrophage-like cells, approved sets):")
print(a.obs.loc[mac, pcols].corr(method="spearman").round(2).to_string())
print("Spearman correlation between programs (macrophage-like cells, filtered sets):")
print(a.obs.loc[mac, fcols].corr(method="spearman").round(2).to_string())
print("Display labels by site (macrophage-like cells, approved sets):")
print(pd.crosstab(a.obs.loc[mac, "display_label"], a.obs.loc[mac, "site"]).to_string())

m = a.obs[mac]
g = m.groupby(["patient", "site"], observed=True)
table = g[pcols + fcols + [acol]].mean().round(4)
table.insert(0, "n_macrophage_cells", g.size())
table.insert(1, "n_myeloid_cells", a.obs.groupby(["patient", "site"], observed=True).size())
table.to_csv(os.path.join(P["out_tables"], f"macrophage_program_patient_site_{tag}.csv"))
print("Tumor macrophage-like cells per patient (mean program scores; approved | filtered | ambient):")
print(table.xs("Tumor", level="site").to_string())

keep_cols = ["patient", "site", "cluster_label", "macrophage_like", "display_label"] + pcols + fcols + [acol]
a.obs[keep_cols].to_csv(os.path.join(P["out_tables"], f"macrophage_program_scores_{tag}.csv"))
a.uns["pyucell_version"] = version("pyucell")
a.uns["programs_frozen"] = yaml.safe_dump(frozen)
a.write_h5ad(P["out_scored"], compression="gzip")
print(f"Saved {P['out_scored']} and tables in {P['out_tables']}/")