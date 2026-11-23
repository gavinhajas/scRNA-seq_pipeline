# 08_score_validation_cohorts.py: Apply frozen discovery rules/signatures to GSE156625 and GSE140228.
# Usage: python scripts/08_score_validation_cohorts.py [GSE156625 GSE140228]
#
# 1. Load frozen gene lists + primary ('macrophage_like') and program-free ('macrophage_like_noprogram') gates.
# 2. Standardize patient/site labels in memory (GSE156625: drop HN1, map T/N, pull unscaled .raw; GSE140228: Donor/Tissue).
# 3. Score myeloid & CD8 cells (pyUCell, max_rank=1500), apply strict CD8 gate, and save CSVs, summaries, and *_scored.h5ad.

import os
import sys
from importlib.metadata import version

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from pyucell import compute_ucell_scores

COHORT_SPECS = {
    "GSE156625": {
        "myeloid_clean": "data/processed/GSE156625_myeloid_clean.h5ad",
        "cd8_clean": "data/processed/GSE156625_cd8_clean.h5ad",
        "patient_col": "patientno",
        "exclude_patients": ["HN1"],
        "site_col": "NormalvsTumor",
        "site_map": {"T": "Tumor", "N": "Normal"},
        "cluster_col": "legacy_label",
    },
    "GSE140228": {
        "myeloid_clean": "data/processed/GSE140228_myeloid_clean.h5ad",
        "cd8_clean": "data/processed/GSE140228_cd8_clean.h5ad",
        "patient_col": "Donor",
        "exclude_patients": [],
        "site_col": "Tissue",
        "site_map": {},
        "cluster_col": "celltype_sub",
    },
}


def prepare_expression_adata(adata_in: ad.AnnData, label: str) -> ad.AnnData:
    X_sample = adata_in.X[: min(500, adata_in.n_obs)]
    min_val = float(X_sample.min()) if sp.issparse(X_sample) else float(np.nanmin(X_sample))
    if min_val < -1e-6:
        assert adata_in.raw is not None, f"{label}: .X has negative (z-scored) values and .raw is None"
        print(f"  [{label}] Detected z-scored .X (min={min_val:.2f}); extracting unscaled log-normalised .raw")
        work = adata_in.raw.to_adata()
        work.obs = adata_in.obs.copy()
    else:
        work = adata_in.copy()

    work.var_names = work.var_names.astype(str)
    if not work.var_names.is_unique:
        work.var_names_make_unique()
    return work


def standardize_obs(adata_obj: ad.AnnData, spec: dict, label: str) -> ad.AnnData:
    pcol, scol, ccol = spec["patient_col"], spec["site_col"], spec["cluster_col"]
    assert pcol in adata_obj.obs, f"{label}: missing patient column '{pcol}'"
    assert scol in adata_obj.obs, f"{label}: missing site column '{scol}'"

    obs = adata_obj.obs.copy()
    obs["patient"] = obs[pcol].astype(str)
    if spec["exclude_patients"]:
        keep_mask = ~obs["patient"].isin(spec["exclude_patients"])
        dropped = int((~keep_mask).sum())
        if dropped > 0:
            print(f"  [{label}] Excluding {dropped} cells from non-HCC donors {spec['exclude_patients']}")
            adata_obj = adata_obj[keep_mask].copy()
            obs = adata_obj.obs.copy()
            obs["patient"] = obs[pcol].astype(str)

    raw_site = obs[scol].astype(str)
    obs["site"] = raw_site.replace(spec["site_map"]) if spec["site_map"] else raw_site
    obs["cluster_label"] = obs[ccol].astype(str) if ccol in obs else "unlabeled"
    adata_obj.obs = obs
    return adata_obj


def apply_gate(obs: pd.DataFrame, rule: dict) -> pd.Series:
    pan_key = rule.get("pan_panel_key", "gate_pan_macrophage")
    pan_col = f"{pan_key}_UCell"
    non_mac_cols = [f"{k}_UCell" for k in rule["panels"] if k != pan_key and f"{k}_UCell" in obs.columns]
    non_mac_max = obs[non_mac_cols].max(axis=1)
    return (obs[pan_col] >= float(rule["min_pan_macrophage"])) & (non_mac_max <= float(rule["max_non_macrophage"]))


def score_myeloid(gse: str, spec: dict, prog_frozen: dict, sig_frozen: dict, out_tables: str, out_h5ad_dir: str):
    print(f"{'=' * 70}Scoring Myeloid: {gse}{'=' * 70}")
    raw_ad = ad.read_h5ad(spec["myeloid_clean"])
    a = prepare_expression_adata(raw_ad, f"{gse}_myeloid")
    a = standardize_obs(a, spec, f"{gse}_myeloid")
    print(f"  Cells retained: {a.n_obs:,} across {a.obs['patient'].nunique()} patients")

    rule = prog_frozen.get("macrophage_like_rule")
    rule_np = prog_frozen.get("macrophage_like_rule_noprogram")
    assert rule is not None, "macrophage_like_rule missing; run 07b_calibrate_macrophage_rule.py first"

    sigs = {}
    for panel_name, genes in rule["panels"].items():
        sigs[panel_name] = [g for g in genes if g in a.var_names]
    if rule_np is not None:
        for panel_name, genes in rule_np["panels"].items():
            if panel_name not in sigs:
                sigs[panel_name] = [g for g in genes if g in a.var_names]
    for prog_name, genes in prog_frozen["approved_sets"].items():
        sigs[prog_name] = [g for g in genes if g in a.var_names]
    for prog_name, genes in prog_frozen["filtered_sets"].items():
        sigs[f"{prog_name}__filtered"] = [g for g in genes if g in a.var_names]
    sigs["ambient_hepatocyte"] = [g for g in sig_frozen["ambient"]["hepatocyte"] if g in a.var_names]

    compute_ucell_scores(a, sigs, max_rank=1500, missing_genes="skip", n_jobs=4, device="cpu")

    a.obs["macrophage_like"] = apply_gate(a.obs, rule)
    mac = a.obs["macrophage_like"].to_numpy()
    print(f"  Macrophage-like cells (primary gate): {int(mac.sum()):,} / {a.n_obs:,} ({100.0 * mac.mean():.1f}%)")
    if "cluster_label" in a.obs and a.obs["cluster_label"].nunique() > 1:
        print("  Primary gate breakdown by cluster/subcluster:")
        print(pd.crosstab(a.obs["cluster_label"], a.obs["macrophage_like"]).to_string())

    if rule_np is not None:
        a.obs["macrophage_like_noprogram"] = apply_gate(a.obs, rule_np)
        mac_np = a.obs["macrophage_like_noprogram"].to_numpy()
        print(f"  Macrophage-like cells (program-free gate): {int(mac_np.sum()):,} / {a.n_obs:,} ({100.0 * mac_np.mean():.1f}%)")
        if "cluster_label" in a.obs and a.obs["cluster_label"].nunique() > 1:
            print("  Program-free gate breakdown by cluster/subcluster:")
            print(pd.crosstab(a.obs["cluster_label"], a.obs["macrophage_like_noprogram"]).to_string())

    pcols = [f"{p}_UCell" for p in prog_frozen["approved_sets"]]
    fcols = [f"{p}__filtered_UCell" for p in prog_frozen["filtered_sets"]]
    acol = "ambient_hepatocyte_UCell"

    #for figures, not used in stats
    a.obs["display_label"] = "not_macrophage_like"
    if mac.sum() > 0:
        stds = a.obs.loc[mac, pcols].std(ddof=0).replace(0, np.nan)
        z = (a.obs.loc[mac, pcols] - a.obs.loc[mac, pcols].mean()) / stds
        top = z.idxmax(axis=1).fillna("mixed_or_low").str.replace("_UCell", "", regex=False)
        top[z.max(axis=1).fillna(-np.inf) < 0.5] = "mixed_or_low"
        a.obs.loc[mac, "display_label"] = top

    m = a.obs[mac]
    g = m.groupby(["patient", "site"], observed=True)
    table = g[pcols + fcols + [acol]].mean().round(4)
    table.insert(0, "n_macrophage_cells", g.size())
    table.insert(1, "n_myeloid_cells", a.obs.groupby(["patient", "site"], observed=True).size())
    table.to_csv(os.path.join(out_tables, f"macrophage_program_patient_site_{gse}.csv"))

    if "Tumor" in m["site"].unique():
        print("  Tumor macrophage-like cells per patient (mean program scores, primary gate):")
        print(table.xs("Tumor", level="site").to_string())

    keep_cols = ["patient", "site", "cluster_label", "macrophage_like"]
    if "macrophage_like_noprogram" in a.obs:
        keep_cols.append("macrophage_like_noprogram")
    keep_cols += ["display_label"] + pcols + fcols + [acol]
    a.obs[keep_cols].to_csv(os.path.join(out_tables, f"macrophage_program_scores_{gse}.csv"))

    a.uns["pyucell_version"] = version("pyucell")
    a.uns["programs_frozen"] = yaml.safe_dump(prog_frozen)
    out_h5ad = os.path.join(out_h5ad_dir, f"{gse}_myeloid_scored.h5ad")
    a.write_h5ad(out_h5ad, compression="gzip")
    print(f"  Saved {out_h5ad} and CSV tables in {out_tables}/")


def score_cd8(gse: str, spec: dict, sig_frozen: dict, out_tables: str, out_h5ad_dir: str):
    """Score CD8 compartment with frozen exhaustion variants V1-V4, ambient scores, and strict CD8 gate."""
    print(f"{'=' * 70}Scoring CD8: {gse}{'=' * 70}")
    raw_ad = ad.read_h5ad(spec["cd8_clean"])
    cd8 = prepare_expression_adata(raw_ad, f"{gse}_cd8")
    cd8 = standardize_obs(cd8, spec, f"{gse}_cd8")
    print(f"  Cells retained: {cd8.n_obs:,} across {cd8.obs['patient'].nunique()} patients")

    for flag in ("flag_cd8_expressed", "flag_cd3d_negative", "flag_cd4_only"):
        assert flag in cd8.obs, f"{gse}_cd8: missing expected flag column '{flag}'"
    cd8.obs["gate_strict"] = (
        cd8.obs["flag_cd8_expressed"].astype(bool)
        & ~cd8.obs["flag_cd3d_negative"].astype(bool)
        & ~cd8.obs["flag_cd4_only"].astype(bool)
    )
    print(f"  Strict-gate CD8 cells: {int(cd8.obs['gate_strict'].sum()):,} / {cd8.n_obs:,}")

    sigs = {}
    for vname, genes in sig_frozen["variants"].items():
        keep = [g for g in genes if g in cd8.var_names]
        if len(keep) < len(genes):
            missing = sorted(set(genes) - set(keep))
            print(f"  Note: {vname} missing {len(missing)} gene(s) in {gse}: {missing}")
        sigs[vname] = keep
    for amb_name, genes in sig_frozen["ambient"].items():
        sigs[f"ambient_{amb_name}"] = [g for g in genes if g in cd8.var_names]

    compute_ucell_scores(cd8, sigs, max_rank=1500, missing_genes="skip", n_jobs=4, device="cpu")
    score_cols = [f"{n}_UCell" for n in sigs]

    tum = cd8.obs[(cd8.obs["site"] == "Tumor") & cd8.obs["gate_strict"]]
    print("  Tumor strict-gate CD8 cells per patient:")
    print(tum.groupby("patient", observed=True).size().to_string())

    cd8.obs[["patient", "site", "cluster_label", "gate_strict"] + score_cols].to_csv(
        os.path.join(out_tables, f"exhaustion_scores_{gse}.csv")
    )
    cd8.uns["signatures_frozen"] = yaml.safe_dump(sig_frozen)
    cd8.uns["pyucell_version"] = version("pyucell")
    out_h5ad = os.path.join(out_h5ad_dir, f"{gse}_cd8_scored.h5ad")
    cd8.write_h5ad(out_h5ad, compression="gzip")
    print(f"  Saved {out_h5ad} and exhaustion_scores_{gse}.csv")


def main():
    with open("metadata/macrophage_programs_frozen_GSE149614.yaml") as fh:
        prog_frozen = yaml.safe_load(fh)
    with open("metadata/signatures_frozen_GSE149614.yaml") as fh:
        sig_frozen = yaml.safe_load(fh)

    targets = sys.argv[1:] if len(sys.argv) > 1 else ["GSE156625", "GSE140228"]
    out_tables = "results/tables"
    out_h5ad_dir = "data/processed"
    os.makedirs(out_tables, exist_ok=True)
    os.makedirs(out_h5ad_dir, exist_ok=True)

    for gse in targets:
        assert gse in COHORT_SPECS, f"Unknown cohort '{gse}'. Available: {list(COHORT_SPECS)}"
        spec = COHORT_SPECS[gse]
        score_myeloid(gse, spec, prog_frozen, sig_frozen, out_tables, out_h5ad_dir)
        score_cd8(gse, spec, sig_frozen, out_tables, out_h5ad_dir)


if __name__ == "__main__":
    main()