# 07_patient_level_stats.py: Patient-level correlation of macrophage programs vs. CD8 exhaustion (discovery).
# Usage: python scripts/07_patient_level_stats.py [metadata/stats_plan.yaml]
# - Design: Patient-level Spearman correlation (exact permutation p; Monte Carlo if >9 patients).
# - Variables: Primary tumor mean macrophage program score vs. median strict-gate CD8 exhaustion score.
# - Scope: Primary = V3 exhaustion x 3 filtered programs (Holm-corrected); others = sensitivity.
# - Slight problem: Underpowered at n ~ 8 (needs |rho| ~ 0.74 for p < 0.05) basically descriptive until multi-cohort pooling.

import itertools
import os
import sys

import numpy as np
import pandas as pd
import yaml
from scipy.stats import rankdata

plan_path = sys.argv[1] if len(sys.argv) > 1 else "metadata/stats_plan.yaml"
plan = yaml.safe_load(open(plan_path))
P = plan["paths"]
rng = np.random.default_rng(plan.get("random_seed", 0))
os.makedirs(P["out_dir"], exist_ok=True)
tag = P["tag"]


def rank_corr_perm(x, y, w=None, n_mc=200000):
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    w = np.ones(n) if w is None else np.asarray(w, float)
    w = w / w.sum()
    xr, yr = rankdata(x), rankdata(y)

    def wcorr(xmat):
        cx = xmat - (xmat * w).sum(axis=-1, keepdims=True)
        cy = yr - (yr * w).sum()
        cov = (cx * (w * cy)).sum(axis=-1)
        vx = (w * cx ** 2).sum(axis=-1)
        vy = (w * cy ** 2).sum()
        return cov / np.sqrt(vx * vy)

    r_obs = float(wcorr(xr[None, :])[0])
    if n <= 9:
        perms = np.array(list(itertools.permutations(range(n))))
    else:
        perms = np.array([rng.permutation(n) for _ in range(n_mc)])
    r_perm = wcorr(xr[perms])
    p = float(np.mean(np.abs(r_perm) >= abs(r_obs) - 1e-12))
    return r_obs, p, ("exact" if n <= 9 else "monte_carlo")


def partial_rank_corr(x, y, z):
    xr, yr, zr = rankdata(x), rankdata(y), rankdata(z)
    A = np.column_stack([np.ones(len(zr)), zr])
    rx = xr - A @ np.linalg.lstsq(A, xr, rcond=None)[0]
    ry = yr - A @ np.linalg.lstsq(A, yr, rcond=None)[0]
    return float(np.corrcoef(rx, ry)[0, 1])


# data
ex = pd.read_csv(P["exhaustion_scores"], index_col=0)
pr = pd.read_csv(P["program_scores"], index_col=0)
meta = pd.read_csv(P["clinical_metadata"], sep="\t", keep_default_na=False)
stage = meta.groupby("patient")["stage"].first().map(plan["stage_order"])
assert stage.notna().all(), "unmapped stage values"
pr_all_cols = [c for c in pr.columns if c.endswith("_UCell") and not c.startswith("ambient")]


def patient_table(sites, ex_col, prog_cols, min_cd8, min_mac):
    e = ex[ex["site"].isin(sites) & ex["gate_strict"]]
    m = pr[pr["site"].isin(sites) & pr["macrophage_like"]]
    out = pd.DataFrame({"n_cd8": e.groupby("patient").size(), "exhaustion": e.groupby("patient")[ex_col].median(),
                        "n_mac": m.groupby("patient").size(), "ambient_hepatocyte": m.groupby("patient")["ambient_hepatocyte_UCell"].mean()})
    for c in prog_cols:
        out[c] = m.groupby("patient")[c].mean()
    out["stage"] = stage
    out = out.dropna(subset=["exhaustion"] + prog_cols)
    return out[(out["n_cd8"] >= min_cd8) & (out["n_mac"] >= min_mac)]


def run_tests(tab, prog_cols, label, weighted=False):
    rows = []
    for c in prog_cols:
        r, p, how = rank_corr_perm(tab[c], tab["exhaustion"], tab["n_cd8"] if weighted else None, plan["n_monte_carlo"])
        rows.append({"analysis": label, "program": c, "n_patients": len(tab), "rho": round(r, 3), "p_perm": round(p, 4), "p_method": how})
    return rows


pm = plan["primary"]
progs = pm["programs"]
tab = patient_table(pm["sites"], pm["exhaustion_variant"], progs, pm["min_cd8_cells"], pm["min_macrophage_cells"])
pd.set_option("display.width", 220)
print("Patients included in the primary analysis (tumor only):")
print(tab.round(3).to_string())
excluded = sorted(set(ex["patient"].unique()) - set(tab.index))
print("Excluded by thresholds or missing data:", excluded)
tab.to_csv(os.path.join(P["out_dir"], f"patient_level_table_primary_{tag}.csv"))

# primary
prim = pd.DataFrame(run_tests(tab, progs, "PRIMARY"))
order = np.argsort(prim["p_perm"].to_numpy())
m_tests = len(prim)
holm = np.empty(m_tests)
running = 0.0
for rank, i in enumerate(order):
    running = max(running, (m_tests - rank) * prim.loc[i, "p_perm"])
    holm[i] = min(1.0, running)
prim["p_holm"] = holm.round(4)
print("PRIMARY (V3 exhaustion x filtered programs; Spearman, permutation p, Holm over 3 programs):")
print(prim.to_string(index=False))
prim.to_csv(os.path.join(P["out_dir"], f"stats_primary_{tag}.csv"), index=False)

# sensitivity
sens = plan["sensitivity"]
rows = []
for ex_col in sens["exhaustion_variants"]:
    for pset in sens["program_sets"]:
        pcols = [c.replace("__filtered", "") if pset == "approved" else c for c in progs]
        t = patient_table(pm["sites"], ex_col, pcols, pm["min_cd8_cells"], pm["min_macrophage_cells"])
        for r in run_tests(t, pcols, f"variant={ex_col}|programs={pset}"):
            rows.append(r)
for sites in sens.get("extra_site_sets", []):
    t = patient_table(sites, pm["exhaustion_variant"], progs, pm["min_cd8_cells"], pm["min_macrophage_cells"])
    rows += run_tests(t, progs, f"sites={'+'.join(sites)}")
if sens.get("no_thresholds"):
    t = patient_table(pm["sites"], pm["exhaustion_variant"], progs, 1, 1)
    rows += run_tests(t, progs, "no thresholds (all patients with data)")
if sens.get("weighted"):
    rows += run_tests(tab, progs, "weighted by strict-gate CD8 cell count", weighted=True)
for cov in sens.get("covariates_one_at_a_time", []):
    for c in progs:
        pr_rho = partial_rank_corr(tab[c], tab["exhaustion"], tab[cov])
        rows.append({"analysis": f"partial rho controlling {cov} (descriptive)", "program": c, "n_patients": len(tab),
                     "rho": round(pr_rho, 3), "p_perm": np.nan, "p_method": "none"})
if sens.get("leave_one_patient_out"):
    for c in progs:
        rhos = []
        for pat in tab.index:
            t = tab.drop(index=pat)
            rhos.append(rank_corr_perm(t[c], t["exhaustion"], None, 1000)[0])
        rows.append({"analysis": f"leave-one-patient-out range of rho (min to max {round(min(rhos), 2)} to {round(max(rhos), 2)})",
                     "program": c, "n_patients": len(tab) - 1, "rho": round(float(np.median(rhos)), 3), "p_perm": np.nan, "p_method": "none"})
sens_df = pd.DataFrame(rows)
print("SENSITIVITY (not corrected; descriptive):")
print(sens_df.to_string(index=False))
sens_df.to_csv(os.path.join(P["out_dir"], f"stats_sensitivity_{tag}.csv"), index=False)