# 09_validation_and_pooled_stats.py: Validation stats, cohort pooling, and circularity checks.
# Usage: python scripts/09_validation_and_pooled_stats.py [metadata/stats_plan.yaml]
#
# - Pooling: Excludes n < 5 cohorts (GSE140228) from Fisher-z; uses sqrt(n)-weighted signed Stouffer on exact
#   permutation p-values to include all cohorts safely. Sets random-effects estimate to NaN if k < 3 or I^2 > 75%.
# - Tiers: CONFIRMATORY_VALIDATION_ONLY (held-out cohorts, Holm-corrected) vs. SECONDARY_ALL_COHORTS (+ discovery).
# - Re-evaluates all cohorts with program-free gate ('macrophage_like_noprogram') to rule out circularity.

import itertools
import os
import sys

import numpy as np
import pandas as pd
import yaml
from scipy.stats import norm, rankdata

plan_path = sys.argv[1] if len(sys.argv) > 1 else "metadata/stats_plan.yaml"
with open(plan_path) as fh:
    plan = yaml.safe_load(fh)

out_dir = plan["paths"]["out_dir"]
os.makedirs(out_dir, exist_ok=True)
rng = np.random.default_rng(plan.get("random_seed", 0))

ALL_COHORTS = ["GSE149614", "GSE156625", "GSE140228"]
VALIDATION_COHORTS = ["GSE156625", "GSE140228"]

pool_cfg = plan.get("pooling", {})
MIN_FISHER_Z_PATIENTS = int(pool_cfg.get("min_fisher_z_patients", 5))
MIN_STOUFFER_PATIENTS = int(pool_cfg.get("min_stouffer_patients", 4))
MIN_COHORT_MEAN_EXPOSURE = float(pool_cfg.get("min_cohort_mean_exposure", 0.05))
MAX_I2_RANDOM_EFFECTS = float(pool_cfg.get("max_i2_for_random_effects", 75.0))

GSE140228_ID_MAP = {
    "DSN09": "D20171109",
    "DEJ08": "D20180108",
    "DEJ10": "D20180110",
    "DEJ16": "D20180116",
}


def load_validation_clinical(gse: str) -> pd.DataFrame | None:
    candidates = [
        plan["paths"].get("validation_clinical_metadata", "metadata/validation_clinical.csv"),
        "validation_clinical.csv",
    ]
    clin_path = next((p for p in candidates if p and os.path.exists(p)), None)
    if clin_path is None:
        return None
    clin = pd.read_csv(clin_path, keep_default_na=False)
    sub = clin[clin["cohort"] == gse].copy()
    if sub.empty:
        return None
    if gse == "GSE140228":
        sub["obs_patient"] = sub["patient_id"].replace(GSE140228_ID_MAP)
        sub = sub.rename(columns={"patient_id": "clinical_patient_id"}).set_index("obs_patient")
    else:
        sub = sub.set_index("patient_id")
    sub["stage"] = sub["tumor_stage"].map(plan["stage_order"])
    return sub


def rank_corr_perm(x, y, w=None, n_mc=200000):
    # Rank (Spearman-type) correlation of x and y, optionally weighted, with a permutation p-value. 
    # Exposure ranks x are permuted across patients while each patient's (outcome, weight) stays together.
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan, np.nan, "insufficient_n"
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
    # Spearman correlation of x and y after removing the linear effect of rank(z) from both ranks.
    xr, yr, zr = rankdata(x), rankdata(y), rankdata(z)
    A = np.column_stack([np.ones(len(zr)), zr])
    rx = xr - A @ np.linalg.lstsq(A, xr, rcond=None)[0]
    ry = yr - A @ np.linalg.lstsq(A, yr, rcond=None)[0]
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def holm_adjust(pvals: np.ndarray) -> np.ndarray:
    # Holm step-down adjustment handling NaN entries cleanly.
    pvals = np.asarray(pvals, float)
    out = np.full_like(pvals, np.nan)
    valid = np.where(~np.isnan(pvals))[0]
    if len(valid) == 0:
        return out
    sub = pvals[valid]
    order = np.argsort(sub)
    m = len(sub)
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * sub[idx])
        adj[idx] = min(1.0, running)
    out[valid] = adj
    return out


def build_patient_table(
    ex: pd.DataFrame,
    pr: pd.DataFrame,
    sites: list,
    ex_col: str,
    prog_cols: list,
    min_cd8: int,
    min_mac: int,
    mac_col: str = "macrophage_like",
    stage_series=None,
    clin_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    #Construct patient-level exposure/outcome table for a given site filter, macrophage gate, and cutoffs."""
    e = ex[ex["site"].isin(sites) & ex["gate_strict"].astype(bool)]
    m = pr[pr["site"].isin(sites) & pr[mac_col].astype(bool)]
    out = pd.DataFrame({
        "n_cd8": e.groupby("patient").size(),
        "exhaustion": e.groupby("patient")[ex_col].median(),
        "n_mac": m.groupby("patient").size(),
        "ambient_hepatocyte": m.groupby("patient")["ambient_hepatocyte_UCell"].mean(),
    })
    for c in prog_cols:
        out[c] = m.groupby("patient")[c].mean()
    if stage_series is not None:
        out["stage"] = stage_series
    if clin_df is not None:
        for col in ["clinical_patient_id", "sex", "age", "hbv_status", "tumor_stage", "stage"]:
            if col in clin_df.columns:
                out[col] = clin_df[col]
    out = out.dropna(subset=["exhaustion"] + prog_cols)
    return out[(out["n_cd8"] >= min_cd8) & (out["n_mac"] >= min_mac)]


def run_tests(tab: pd.DataFrame, prog_cols: list, label: str, weighted: bool = False) -> list:
    #Run rank correlation + permutation p-value across all requested programs.
    rows = []
    for c in prog_cols:
        r, p, how = rank_corr_perm(
            tab[c], tab["exhaustion"], tab["n_cd8"] if weighted else None, plan["n_monte_carlo"]
        )
        if len(tab) > 0 and tab["n_mac"].sum() > 0:
            w_mean_exp = float(np.average(tab[c], weights=tab["n_mac"]))
        else:
            w_mean_exp = np.nan
        rows.append({
            "analysis": label,
            "program": c,
            "n_patients": int(len(tab)),
            "mean_exposure": round(w_mean_exp, 4) if not np.isnan(w_mean_exp) else np.nan,
            "rho": round(r, 4) if not np.isnan(r) else np.nan,
            "p_perm": round(p, 4) if not np.isnan(p) else np.nan,
            "p_method": how,
        })
    return rows


def analyze_cohort(gse: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    #Run primary and sensitivity patient-level statistics for one validation cohort
    print(f"{'=' * 76}Patient-Level Analysis: {gse}{'=' * 76}")
    ex = pd.read_csv(os.path.join(out_dir, f"exhaustion_scores_{gse}.csv"), index_col=0)
    pr = pd.read_csv(os.path.join(out_dir, f"macrophage_program_scores_{gse}.csv"), index_col=0)
    ex["patient"] = ex["patient"].astype(str)
    pr["patient"] = pr["patient"].astype(str)

    stage_series = None
    clin_df = None
    if gse == "GSE149614" and os.path.exists(plan["paths"]["clinical_metadata"]):
        meta = pd.read_csv(plan["paths"]["clinical_metadata"], sep="\t", keep_default_na=False)
        stage_series = meta.groupby("patient")["stage"].first().map(plan["stage_order"])
    else:
        clin_df = load_validation_clinical(gse)

    pm = plan["primary"]
    progs = pm["programs"]
    min_mac = (
        int(pm.get("validation_min_macrophage_cells", pm["min_macrophage_cells"]))
        if gse != "GSE149614"
        else int(pm["min_macrophage_cells"])
    )
    tab = build_patient_table(
        ex, pr, pm["sites"], pm["exhaustion_variant"], progs,
        pm["min_cd8_cells"], min_mac, "macrophage_like", stage_series, clin_df,
    )
    print(f"{'=' * 76}Included patients (Tumor only; min_cd8={pm['min_cd8_cells']}, min_mac={min_mac}):{'=' * 76}")
    print(tab.round(3).to_string())
    excluded = sorted(set(ex["patient"].unique()) - set(tab.index))
    print("Excluded by cell thresholds or missing site data:", excluded)
    tab.to_csv(os.path.join(out_dir, f"patient_level_table_primary_{gse}.csv"))

    # Primary tests
    prim = pd.DataFrame(run_tests(tab, progs, "PRIMARY"))
    prim.insert(0, "cohort", gse)
    prim["p_holm"] = np.round(holm_adjust(prim["p_perm"].to_numpy()), 4)
    print(f"{'=' * 76}PRIMARY results:{'=' * 76}")
    print(prim.to_string(index=False))
    prim.to_csv(os.path.join(out_dir, f"stats_primary_{gse}.csv"), index=False)

    # Sensitivity tests
    sens = plan["sensitivity"]
    rows = []
    for ex_col in sens["exhaustion_variants"]:
        for pset in sens["program_sets"]:
            pcols = [c.replace("__filtered", "") if pset == "approved" else c for c in progs]
            t = build_patient_table(
                ex, pr, pm["sites"], ex_col, pcols, pm["min_cd8_cells"], min_mac,
                "macrophage_like", stage_series, clin_df,
            )
            rows += run_tests(t, pcols, f"variant={ex_col}|programs={pset}")

    if sens.get("program_free_macrophage_gate") and "macrophage_like_noprogram" in pr.columns:
        t_np = build_patient_table(
            ex, pr, pm["sites"], pm["exhaustion_variant"], progs, pm["min_cd8_cells"], min_mac,
            "macrophage_like_noprogram", stage_series, clin_df,
        )
        np_rows = run_tests(t_np, progs, "program_free_macrophage_gate")
        rows += np_rows
        print(f"{'=' * 76}SENSITIVITY (program-free macrophage gate — excludes SPP1, TREM2, and all program genes):{'=' * 76}")
        print(pd.DataFrame(np_rows).to_string(index=False))

    for sites in sens.get("extra_site_sets", []):
        if any(s in ex["site"].unique() for s in sites if s != "Tumor"):
            t = build_patient_table(
                ex, pr, sites, pm["exhaustion_variant"], progs, pm["min_cd8_cells"], min_mac,
                "macrophage_like", stage_series, clin_df,
            )
            rows += run_tests(t, progs, f"sites={'+'.join(sites)}")

    if sens.get("no_thresholds"):
        t = build_patient_table(
            ex, pr, pm["sites"], pm["exhaustion_variant"], progs, 1, 1,
            "macrophage_like", stage_series, clin_df,
        )
        rows += run_tests(t, progs, "no thresholds (all patients with data)")

    if sens.get("weighted"):
        rows += run_tests(tab, progs, "weighted by strict-gate CD8 cell count", weighted=True)

    for cov in sens.get("covariates_one_at_a_time", []):
        if cov not in tab.columns or tab[cov].isna().all():
            print(f"  Note: covariate '{cov}' not available for {gse}; skipping descriptive partial correlation.")
            continue
        if tab[cov].nunique(dropna=True) < 2:
            vals = tab[cov].dropna().unique().tolist()
            print(f"  Note: covariate '{cov}' is constant ({vals}) across included {gse} patients; skipping descriptive partial correlation.")
            continue
        for c in progs:
            pr_rho = partial_rank_corr(tab[c], tab["exhaustion"], tab[cov])
            rows.append({
                "analysis": f"partial rho controlling {cov} (descriptive)",
                "program": c,
                "n_patients": int(len(tab)),
                "mean_exposure": round(float(np.average(tab[c], weights=tab["n_mac"])), 4),
                "rho": round(pr_rho, 4) if not np.isnan(pr_rho) else np.nan,
                "p_perm": np.nan,
                "p_method": "none",
            })

    if sens.get("leave_one_patient_out") and len(tab) >= 4:
        for c in progs:
            rhos = []
            for pat in tab.index:
                t = tab.drop(index=pat)
                rhos.append(rank_corr_perm(t[c], t["exhaustion"], None, 1000)[0])
            rows.append({
                "analysis": f"leave-one-patient-out range of rho (min to max {round(min(rhos), 2)} to {round(max(rhos), 2)})",
                "program": c,
                "n_patients": int(len(tab) - 1),
                "mean_exposure": round(float(np.average(tab[c], weights=tab["n_mac"])), 4),
                "rho": round(float(np.median(rhos)), 4),
                "p_perm": np.nan,
                "p_method": "none",
            })

    sens_df = pd.DataFrame(rows)
    sens_df.insert(0, "cohort", gse)
    sens_df.to_csv(os.path.join(out_dir, f"stats_sensitivity_{gse}.csv"), index=False)
    return prim, sens_df


def stouffer_perm_pool(sub_df: pd.DataFrame) -> dict:
    #Non-parametric signed Stouffer combination of exact two-sided permutation p-values (weights = sqrt(n)). 
    # Avoids infinite Fisher-z when n = 4 and rho = +/-1.0.
    elig = sub_df[
        (sub_df["n_patients"] >= MIN_STOUFFER_PATIENTS)
        & (sub_df["mean_exposure"] >= MIN_COHORT_MEAN_EXPOSURE)
        & sub_df["rho"].notna()
        & sub_df["p_perm"].notna()
    ].copy()
    if elig.empty:
        return {
            "stouffer_k": 0,
            "stouffer_cohorts": "",
            "stouffer_n": 0,
            "stouffer_z": np.nan,
            "stouffer_p": np.nan,
        }
    rhos = elig["rho"].to_numpy(float)
    p_two = elig["p_perm"].to_numpy(float)
    ns = elig["n_patients"].to_numpy(float)
    weights = np.sqrt(ns)

    # Exact Spearman permutation distribution is symmetric: one-sided p = p_two / 2
    z_list = []
    for r, p2 in zip(rhos, p_two):
        if r == 0.0:
            z_list.append(0.0)
        else:
            p1 = float(np.clip(p2 / 2.0, 1e-12, 0.5))
            z_list.append(float(np.sign(r) * norm.isf(p1)))
    z_arr = np.asarray(z_list, float)
    z_comb = float(np.sum(weights * z_arr) / np.sqrt(np.sum(weights ** 2)))
    p_comb = float(2.0 * norm.sf(abs(z_comb)))
    return {
        "stouffer_k": int(len(elig)),
        "stouffer_cohorts": ",".join(elig["cohort"].tolist()),
        "stouffer_n": int(ns.sum()),
        "stouffer_z": round(z_comb, 4),
        "stouffer_p": round(p_comb, 5),
    }


def fisher_z_pool(sub_df: pd.DataFrame, tier: str, label: str, program: str) -> dict:
    #Inverse-variance Fisher-z meta-analysis (n >= 5 only; no +/-1 clipping) + signed Stouffer check.
    eligible = sub_df[
        (sub_df["n_patients"] >= MIN_FISHER_Z_PATIENTS)
        & (sub_df["mean_exposure"] >= MIN_COHORT_MEAN_EXPOSURE)
        & sub_df["rho"].notna()
        & (sub_df["rho"].abs() < 1.0)
    ].copy()

    excluded_cohorts = sorted(set(sub_df["cohort"]) - set(eligible["cohort"]))
    st = stouffer_perm_pool(sub_df)

    k = len(eligible)
    if k == 0:
        return {
            "tier": tier,
            "analysis": label,
            "program": program,
            "k_cohorts": 0,
            "cohorts_pooled": "",
            "cohorts_excluded_fisher_z": ",".join(excluded_cohorts),
            "total_patients": 0,
            "rho_fixed": np.nan,
            "ci95_low_fixed": np.nan,
            "ci95_high_fixed": np.nan,
            "p_fixed": np.nan,
            "rho_random": np.nan,
            "ci95_low_random": np.nan,
            "ci95_high_random": np.nan,
            "p_random": np.nan,
            "cochran_Q": np.nan,
            "I2_percent": np.nan,
            "random_effects_note": "no eligible cohorts (n >= 5)",
            **st,
        }

    r = eligible["rho"].to_numpy(float)
    n = eligible["n_patients"].to_numpy(float)
    z = np.arctanh(r)
    v = 1.0 / (n - 3.0)
    w = 1.0 / v

    z_fe = float(np.sum(w * z) / np.sum(w))
    se_fe = float(np.sqrt(1.0 / np.sum(w)))
    p_fe = float(2.0 * norm.sf(abs(z_fe / se_fe)))
    rho_fe = float(np.tanh(z_fe))
    ci_lo_fe = float(np.tanh(z_fe - 1.96 * se_fe))
    ci_hi_fe = float(np.tanh(z_fe + 1.96 * se_fe))

    q_stat, i2 = np.nan, np.nan
    rho_re, ci_lo_re, ci_hi_re, p_re = np.nan, np.nan, np.nan, np.nan
    re_note = "single cohort (k=1)"

    if k >= 2:
        q_stat = float(np.sum(w * (z - z_fe) ** 2))
        df_q = k - 1
        i2 = max(0.0, 100.0 * (q_stat - df_q) / q_stat) if q_stat > 0 else 0.0
        if k < 3:
            re_note = "suppressed (k < 3; report Q and I^2)"
        elif i2 > MAX_I2_RANDOM_EFFECTS:
            re_note = f"suppressed (I^2={i2:.1f}% > {MAX_I2_RANDOM_EFFECTS:.0f}%)"
        else:
            c_denom = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
            tau2 = max(0.0, (q_stat - df_q) / c_denom) if c_denom > 0 else 0.0
            w_re = 1.0 / (v + tau2)
            z_re = float(np.sum(w_re * z) / np.sum(w_re))
            se_re = float(np.sqrt(1.0 / np.sum(w_re)))
            p_re = float(2.0 * norm.sf(abs(z_re / se_re)))
            rho_re = float(np.tanh(z_re))
            ci_lo_re = float(np.tanh(z_re - 1.96 * se_re))
            ci_hi_re = float(np.tanh(z_re + 1.96 * se_re))
            re_note = "reported"

    return {
        "tier": tier,
        "analysis": label,
        "program": program,
        "k_cohorts": int(k),
        "cohorts_pooled": ",".join(eligible["cohort"].tolist()),
        "cohorts_excluded_fisher_z": ",".join(excluded_cohorts),
        "total_patients": int(n.sum()),
        "rho_fixed": round(rho_fe, 4),
        "ci95_low_fixed": round(ci_lo_fe, 4),
        "ci95_high_fixed": round(ci_hi_fe, 4),
        "p_fixed": round(p_fe, 5),
        "rho_random": round(rho_re, 4) if not np.isnan(rho_re) else np.nan,
        "ci95_low_random": round(ci_lo_re, 4) if not np.isnan(ci_lo_re) else np.nan,
        "ci95_high_random": round(ci_hi_re, 4) if not np.isnan(ci_hi_re) else np.nan,
        "p_random": round(p_re, 5) if not np.isnan(p_re) else np.nan,
        "cochran_Q": round(q_stat, 4) if not np.isnan(q_stat) else np.nan,
        "I2_percent": round(i2, 1) if not np.isnan(i2) else np.nan,
        "random_effects_note": re_note,
        **st,
    }


def add_discovery_noprogram_sensitivity(sens_all: list[pd.DataFrame]):
    #Compute program-free macrophage gate sensitivity for GSE149614 if column exists in discovery scores.
    pr_path = os.path.join(out_dir, "macrophage_program_scores_GSE149614.csv")
    ex_path = os.path.join(out_dir, "exhaustion_scores_GSE149614.csv")
    if not (os.path.exists(pr_path) and os.path.exists(ex_path)):
        return
    pr = pd.read_csv(pr_path, index_col=0)
    if "macrophage_like_noprogram" not in pr.columns:
        return
    ex = pd.read_csv(ex_path, index_col=0)
    ex["patient"] = ex["patient"].astype(str)
    pr["patient"] = pr["patient"].astype(str)
    pm = plan["primary"]
    t_np = build_patient_table(
        ex, pr, pm["sites"], pm["exhaustion_variant"], pm["programs"],
        pm["min_cd8_cells"], int(pm["min_macrophage_cells"]), "macrophage_like_noprogram",
    )
    np_df = pd.DataFrame(run_tests(t_np, pm["programs"], "program_free_macrophage_gate"))
    np_df.insert(0, "cohort", "GSE149614")
    sens_all.append(np_df)


def main():
    # 1. Run validation cohorts
    for gse in VALIDATION_COHORTS:
        analyze_cohort(gse)

    # 2. Load primary and sensitivity results across all 3 cohorts
    prim_all = []
    sens_all = []
    for gse in ALL_COHORTS:
        p_df = pd.read_csv(os.path.join(out_dir, f"stats_primary_{gse}.csv"))
        if "cohort" not in p_df.columns:
            p_df.insert(0, "cohort", gse)
        if "mean_exposure" not in p_df.columns:
            ptab = pd.read_csv(os.path.join(out_dir, f"patient_level_table_primary_{gse}.csv"), index_col=0)
            p_df["mean_exposure"] = p_df["program"].map(
                lambda c: float(np.average(ptab[c], weights=ptab["n_mac"]))
            )
        prim_all.append(p_df)

        s_df = pd.read_csv(os.path.join(out_dir, f"stats_sensitivity_{gse}.csv"))
        if "cohort" not in s_df.columns:
            s_df.insert(0, "cohort", gse)
        if "mean_exposure" not in s_df.columns:
            ptab = pd.read_csv(os.path.join(out_dir, f"patient_level_table_primary_{gse}.csv"), index_col=0)
            s_df["mean_exposure"] = s_df["program"].map(
                lambda c: float(np.average(ptab[c], weights=ptab["n_mac"])) if c in ptab.columns else MIN_COHORT_MEAN_EXPOSURE
            )
        sens_all.append(s_df)

    add_discovery_noprogram_sensitivity(sens_all)

    prim_cat = pd.concat(prim_all, ignore_index=True)
    sens_cat = pd.concat(sens_all, ignore_index=True)

    # 3. Confirmatory (Validation-Only) and Secondary (All Cohorts incl. Discovery) Pooling
    pm = plan["primary"]
    val_rows = []
    all_rows = []
    for prog in pm["programs"]:
        sub_val = prim_cat[(prim_cat["program"] == prog) & (prim_cat["cohort"].isin(VALIDATION_COHORTS))]
        sub_all = prim_cat[prim_cat["program"] == prog]
        val_rows.append(fisher_z_pool(sub_val, "CONFIRMATORY_VALIDATION_ONLY", "POOLED_VALIDATION_ONLY", prog))
        all_rows.append(fisher_z_pool(sub_all, "SECONDARY_ALL_COHORTS", "POOLED_ALL_INCL_DISCOVERY", prog))

    df_val = pd.DataFrame(val_rows)
    df_val["p_holm_fixed"] = np.round(holm_adjust(df_val["p_fixed"].to_numpy()), 5)
    df_val["p_holm_stouffer"] = np.round(holm_adjust(df_val["stouffer_p"].to_numpy()), 5)

    df_all = pd.DataFrame(all_rows)
    df_all["p_holm_fixed"] = np.round(holm_adjust(df_all["p_fixed"].to_numpy()), 5)
    df_all["p_holm_stouffer"] = np.round(holm_adjust(df_all["stouffer_p"].to_numpy()), 5)

    pooled_prim = pd.concat([df_val, df_all], ignore_index=True)

    show_cols = [
        "tier", "program", "k_cohorts", "cohorts_pooled", "total_patients",
        "rho_fixed", "ci95_low_fixed", "ci95_high_fixed", "p_fixed", "p_holm_fixed",
        "cochran_Q", "I2_percent", "stouffer_cohorts", "stouffer_n", "stouffer_z", "stouffer_p", "p_holm_stouffer",
    ]
    print(f"{'=' * 90}")
    print("CONFIRMATORY (VALIDATION ONLY) & SECONDARY (ALL COHORTS) META-ANALYSIS")
    print(f"Fisher-z requires n >= {MIN_FISHER_Z_PATIENTS} (excludes n=4 GSE140228); Stouffer pools exact perm p (n >= {MIN_STOUFFER_PATIENTS})")
    print(f"{'=' * 90}")
    print(pooled_prim[show_cols].to_string(index=False))
    pooled_prim.to_csv(os.path.join(out_dir, "stats_pooled_primary.csv"), index=False)

    # 4. Program-free macrophage gate sensitivity summary across all 3 cohorts
    np_sub = sens_cat[sens_cat["analysis"] == "program_free_macrophage_gate"]
    if not np_sub.empty:
        print(f"{'=' * 90}")
        print("SENSITIVITY: PROGRAM-FREE MACROPHAGE GATE (ZERO SPP1/TREM2/FOLR2 GENES)")
        print(f"{'=' * 90}")
        print(np_sub[["cohort", "program", "n_patients", "mean_exposure", "rho", "p_perm"]].to_string(index=False))

    # 5. Pooled sensitivity analyses (secondary all-cohort + validation-only for program-free gate)
    pooled_sens_rows = []
    shared_analyses = sens_cat[sens_cat["p_method"] != "none"][["analysis", "program"]].drop_duplicates()
    for _, row in shared_analyses.iterrows():
        ana, prog = row["analysis"], row["program"]
        sub = sens_cat[(sens_cat["analysis"] == ana) & (sens_cat["program"] == prog)]
        pooled_sens_rows.append(fisher_z_pool(sub, "SECONDARY_SENSITIVITY", ana, prog))
        if ana == "program_free_macrophage_gate":
            sub_v = sub[sub["cohort"].isin(VALIDATION_COHORTS)]
            pooled_sens_rows.append(fisher_z_pool(sub_v, "CONFIRMATORY_SENSITIVITY", ana, prog))

    pooled_sens = pd.DataFrame(pooled_sens_rows)
    pooled_sens.to_csv(os.path.join(out_dir, "stats_pooled_sensitivity.csv"), index=False)
    print(f"Saved pooled primary and sensitivity tables to {out_dir}/")


if __name__ == "__main__":
    main()