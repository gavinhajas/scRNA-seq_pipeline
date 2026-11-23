# 07b_calibrate_macrophage_rule.py: Tune and freeze UCell macrophage gates on GSE149614 before validation.
# Usage: python scripts/07b_calibrate_macrophage_rule.py [metadata/macrophage_programs_frozen_GSE149614.yaml]
#
# - Gate rule: pan_macrophage_UCell >= min_pan_mac AND max(monocyte_FCN1, cDC1, cDC2) <= max_non_mac,
#   optimized (max F1, Jaccard tie-break) against discovery macrophage clusters {2, 3, 5}.
# - Variants: Primary ('macrophage_like_rule') + circularity-free ('macrophage_like_rule_noprogram', drops program genes).
# - Outputs: Verifies rho concordance with script 07, freezes rules in YAML, and updates GSE149614 scores CSV.

import datetime
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from pyucell import compute_ucell_scores
from scipy.stats import spearmanr

frozen_path = sys.argv[1] if len(sys.argv) > 1 else "metadata/macrophage_programs_frozen_GSE149614.yaml"
stats_plan_path = sys.argv[2] if len(sys.argv) > 2 else "metadata/stats_plan.yaml"

with open(frozen_path) as fh:
    frozen = yaml.safe_load(fh)
with open(stats_plan_path) as fh:
    plan = yaml.safe_load(fh)

myeloid_path = (
    "data/processed/GSE149614_myeloid_scored.h5ad"
    if os.path.exists("data/processed/GSE149614_myeloid_scored.h5ad")
    else "data/processed/GSE149614_myeloid_clean.h5ad"
)
print(f"Loading discovery myeloid object: {myeloid_path}")
a = ad.read_h5ad(myeloid_path)
for col in ("patient", "site", "cluster_label"):
    assert col in a.obs, f"Required column '{col}' missing from {myeloid_path}"

true_mac_labels = ["TAM_SPP1_MMP9", "TAM_C1Q_APOE_FOLR2", "Kupffer_CD5L_MARCO"]
y_true = a.obs["cluster_label"].astype(str).isin(true_mac_labels).to_numpy()
assert y_true.sum() > 0, "No ground-truth macrophage-like cells found in GSE149614"
print(f"Discovery myeloid cells: {a.n_obs:,} | Ground-truth macrophage-like: {int(y_true.sum()):,}")

# Verify program genes so the program-free panel is guaranteed disjoint from all program genes
all_program_genes = {
    g for gs in list(frozen["approved_sets"].values()) + list(frozen["filtered_sets"].values()) for g in gs
}

# Primary pole-balanced panel (retains SPP1+ TAMs alongside C1Q/FOLR2 TAMs)
# Sensitivity program-free panel (zero overlap with any macrophage program gene)
gate_panels = {
    "gate_pan_macrophage": ["CD68", "APOE", "TYROBP", "FCER1G", "C1QA", "SPP1", "TREM2"],
    "gate_pan_mac_noprogram": ["CD68", "APOE", "TYROBP", "FCER1G", "CSF1R", "AIF1", "MS4A7"],
    "gate_monocyte_FCN1": ["FCN1", "S100A8", "S100A9", "VCAN"],
    "gate_cDC1": ["CLEC9A", "XCR1", "CADM1"],
    "gate_cDC2": ["CD1C", "FCER1A", "CLEC10A"],
    "gate_mregDC": ["LAMP3", "CCR7", "FSCN1"],
}
overlap = set(gate_panels["gate_pan_mac_noprogram"]) & all_program_genes
assert not overlap, f"Program-free panel overlaps program genes: {overlap}"

resolved_panels = {}
for k, genes in gate_panels.items():
    present = [g for g in genes if g in a.var_names]
    missing = [g for g in genes if g not in a.var_names]
    print(f"  {k}: {present}" + (f" (missing: {missing})" if missing else ""))
    assert len(present) >= 2, f"{k} has fewer than 2 genes in GSE149614"
    resolved_panels[k] = present

compute_ucell_scores(a, resolved_panels, max_rank=1500, missing_genes="skip", n_jobs=4, device="cpu")

non_mac_keys = ["gate_monocyte_FCN1", "gate_cDC1", "gate_cDC2", "gate_mregDC"]
non_mac_cols = [f"{k}_UCell" for k in non_mac_keys]
non_mac_score = a.obs[non_mac_cols].max(axis=1).to_numpy()

print("Median lineage gate scores by discovery cluster:")
summary_cols = ["gate_pan_macrophage_UCell", "gate_pan_mac_noprogram_UCell"] + non_mac_cols
print(a.obs.groupby("cluster_label", observed=True)[summary_cols].median().round(3).to_string())


def calibrate_gate(pan_key: str, label: str, out_csv: str) -> tuple[dict, np.ndarray]:
    """Grid-search (min_pan, max_non_mac) on GSE149614 maximizing F1 with Jaccard tie-break."""
    pan_score = a.obs[f"{pan_key}_UCell"].to_numpy()
    pan_grid = np.round(np.arange(0.15, 0.65, 0.02), 2)
    non_mac_grid = np.round(np.arange(0.15, 0.55, 0.02), 2)

    records = []
    for t_pan in pan_grid:
        for t_non in non_mac_grid:
            y_pred = (pan_score >= t_pan) & (non_mac_score <= t_non)
            tp = int(np.sum(y_true & y_pred))
            fp = int(np.sum(~y_true & y_pred))
            fn = int(np.sum(y_true & ~y_pred))
            if tp == 0:
                continue
            prec = tp / (tp + fp)
            rec = tp / (tp + fn)
            f1 = 2 * prec * rec / (prec + rec)
            jaccard = tp / (tp + fp + fn)
            records.append({
                "min_pan_macrophage": float(t_pan),
                "max_non_macrophage": float(t_non),
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "n_selected": tp + fp,
                "precision": round(prec, 4),
                "recall": round(rec, 4),
                "f1": round(f1, 4),
                "jaccard": round(jaccard, 4),
            })

    grid_df = pd.DataFrame(records).sort_values(["f1", "jaccard", "precision"], ascending=False).reset_index(drop=True)
    os.makedirs("results/tables", exist_ok=True)
    grid_df.to_csv(out_csv, index=False)

    best = grid_df.iloc[0]
    t_pan_best = float(best["min_pan_macrophage"])
    t_non_best = float(best["max_non_macrophage"])
    y_best = (pan_score >= t_pan_best) & (non_mac_score <= t_non_best)

    print(f"{'=' * 72}Calibrated Rule: {label}{'=' * 72}")
    print(grid_df.head(5).to_string(index=False))
    print("Confusion breakdown by discovery cluster:")
    ct = pd.crosstab(a.obs["cluster_label"], y_best, rownames=["cluster"], colnames=["selected"])
    print(ct.to_string())

    rule_panels = {pan_key: resolved_panels[pan_key]}
    for k in non_mac_keys:
        rule_panels[k] = resolved_panels[k]

    rule_dict = {
        "calibrated_on": "GSE149614",
        "calibrated_date": datetime.date.today().isoformat(),
        "ucell_max_rank": 1500,
        "pan_panel_key": pan_key,
        "panels": rule_panels,
        "min_pan_macrophage": t_pan_best,
        "max_non_macrophage": t_non_best,
        "rule_description": f"{pan_key}_UCell >= {t_pan_best} and max({', '.join(non_mac_cols)}) <= {t_non_best}",
        "discovery_metrics": {
            "precision": float(best["precision"]),
            "recall": float(best["recall"]),
            "f1": float(best["f1"]),
            "jaccard": float(best["jaccard"]),
        },
    }
    return rule_dict, y_best


rule_primary, y_primary = calibrate_gate(
    "gate_pan_macrophage",
    "PRIMARY (pole-balanced)",
    "results/tables/macrophage_rule_calibration_GSE149614.csv",
)
rule_noprogram, y_noprogram = calibrate_gate(
    "gate_pan_mac_noprogram",
    "SENSITIVITY (program-free, zero overlap with FOLR2/SPP1/TREM2 genes)",
    "results/tables/macrophage_rule_calibration_noprogram_GSE149614.csv",
)

# Verify discovery Tumor Spearman rho across cluster-based, primary agnostic, and program-free agnostic rules
prog_csv = plan["paths"]["program_scores"]
ex_csv = plan["paths"]["exhaustion_scores"]
if os.path.exists(prog_csv) and os.path.exists(ex_csv):
    pr = pd.read_csv(prog_csv, index_col=0)
    ex = pd.read_csv(ex_csv, index_col=0)
    pr["macrophage_like_agnostic"] = pd.Series(y_primary, index=a.obs_names).reindex(pr.index).astype(bool)
    pr["macrophage_like_noprogram"] = pd.Series(y_noprogram, index=a.obs_names).reindex(pr.index).astype(bool)
    pr.to_csv(prog_csv)
    print(f"Appended macrophage_like_noprogram to {prog_csv}")

    pm = plan["primary"]
    e_tum = ex[ex["site"].isin(pm["sites"]) & ex["gate_strict"]]
    cd8_med = e_tum.groupby("patient")[pm["exhaustion_variant"]].median()
    cd8_n = e_tum.groupby("patient").size()

    print("Discovery primary Tumor Spearman rho check (cluster vs. primary agnostic vs. program-free agnostic):")
    for prog in pm["programs"]:
        out_str = []
        for col_name, short in [
            ("macrophage_like", "cluster"),
            ("macrophage_like_agnostic", "primary_gate"),
            ("macrophage_like_noprogram", "noprogram_gate"),
        ]:
            m_sub = pr[pr["site"].isin(pm["sites"]) & pr[col_name]]
            df = pd.DataFrame({
                "y": cd8_med,
                "n_cd8": cd8_n,
                "x": m_sub.groupby("patient")[prog].mean(),
                "n_mac": m_sub.groupby("patient").size(),
            }).dropna()
            df = df[(df["n_cd8"] >= pm["min_cd8_cells"]) & (df["n_mac"] >= pm["min_macrophage_cells"])]
            r_val = spearmanr(df["x"], df["y"]).statistic
            out_str.append(f"{short}={r_val:+.3f} (n={len(df)})")
        print(f"  {prog:28s}: " + " | ".join(out_str))

frozen["macrophage_like_rule"] = rule_primary
frozen["macrophage_like_rule_noprogram"] = rule_noprogram

with open(frozen_path, "w") as fh:
    yaml.safe_dump(frozen, fh, sort_keys=False)

print(f"Frozen primary and program-free macrophage-like rules saved to {frozen_path}")