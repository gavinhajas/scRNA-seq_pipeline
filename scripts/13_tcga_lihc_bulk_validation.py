# 13_tcga_lihc_bulk_validation.py: Validate macrophage-exhaustion links & Overall Survival in TCGA-LIHC (n ~ 370).
# Usage: python scripts/13_tcga_lihc_bulk_validation.py
#
# - Data: Cached UCSC Xena TCGA-LIHC Primary Tumors ('-01', 1/patient) with expression, OS/PFI survival, and stage.
# - Scores: Frozen filtered macrophage programs + V3 exhaustion (V1/V2/V4 sensitivity), pan-mac/CD8 bulk controls,
#   and purity-adjusted indices (z-score polarization diffs + CD8-residualized exhaustion).
# - Stats: Raw & partial Spearman (+/- CD8, pan-mac, stage; Holm-corrected), multivariable OLS, and univariable/
#   stage-adjusted Cox PH (Breslow ties) + log-rank OS tests. Saves results/tables/tcga_lihc_*.csv & PDF figure.

import os
import re
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})

import numpy as np
import pandas as pd
import yaml
from scipy.optimize import minimize
from scipy.stats import chi2, norm, rankdata, spearmanr, t as student_t

CACHE_DIR = "data/raw/TCGA_LIHC"
OUT_TABLES = "results/tables"
OUT_FIGS = "figures"
os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs(OUT_TABLES, exist_ok=True)
os.makedirs(OUT_FIGS, exist_ok=True)

# UCSC Xena TCGA-LIHC Hub URLs (public, stable S3 endpoints)
XENA_FILES = {
    "expr": (
        "TCGA_LIHC_HiSeqV2.gz",
        "https://tcga.xenahubs.net/download/TCGA.LIHC.sampleMap%2FHiSeqV2.gz",
    ),
    "survival": (
        "TCGA_LIHC_survival.txt",
        "https://tcga.xenahubs.net/download/survival%2FLIHC_survival.txt",
    ),
    "clinical": (
        "TCGA_LIHC_clinicalMatrix.tsv",
        "https://tcga.xenahubs.net/download/TCGA.LIHC.sampleMap%2FLIHC_clinicalMatrix",
    ),
}

# Bulk lineage-abundance reference panels to control tumor purity / immune infiltration
PAN_MAC_BULK = ["CD68", "CSF1R", "C1QA", "C1QB", "C1QC", "TYROBP", "FCER1G"]
CD8_INFILT_BULK = ["CD8A", "CD8B", "CD3D", "CD3E"]
TREG_BULK = ["FOXP3", "IL2RA"]  # Sensitivity control for bulk Treg overlap with TIGIT/CTLA4 in V3

def fetch_file(local_name: str, url: str) -> str:
    #Download file once into CACHE_DIR if not already cached.
    dest = os.path.join(CACHE_DIR, local_name)
    if not os.path.exists(dest):
        print(f"  Downloading {local_name} from UCSC Xena...")
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as resp, open(dest, "wb") as fh:
            fh.write(resp.read())
    return dest


def holm_adjust(pvals: np.ndarray) -> np.ndarray:
    #Holm step-down adjustment handling NaN entries cleanly
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


def parse_ajcc_stage(stage_str: str) -> float:
    #Map AJCC pathologic stage strings ('Stage I', 'Stage IIIB', etc.) to ordinal 1..4.
    if not isinstance(stage_str, str):
        return np.nan
    s = stage_str.strip().upper()
    if not s or "NOT REPORTED" in s or "DISCREPANCY" in s or s == "NA":
        return np.nan
    m = re.search(r"STAGE\s*(IV|III|II|I)\b", s)
    if not m:
        m = re.search(r"STAGE\s*(IV|III|II|I)", s)
    if not m:
        return np.nan
    return {"I": 1.0, "II": 2.0, "III": 3.0, "IV": 4.0}[m.group(1)]


def score_mean_z(expr_df: pd.DataFrame, genes: list[str], label: str) -> pd.Series:
    #Compute per-sample signature score as mean gene-wise z-score across primary tumor samples.
    present = [g for g in genes if g in expr_df.index]
    missing = sorted(set(genes) - set(present))
    if missing:
        print(f"  Note [{label}]: missing {len(missing)} gene(s) in TCGA-LIHC: {missing}")
    assert len(present) >= 2, f"{label}: fewer than 2 genes found in TCGA-LIHC"
    sub = expr_df.loc[present].astype(float)
    std = sub.std(axis=1, ddof=0).replace(0, np.nan)
    z = sub.sub(sub.mean(axis=1), axis=0).div(std, axis=0)
    return z.mean(axis=0)


def partial_spearman(x: np.ndarray, y: np.ndarray, covs: np.ndarray) -> tuple[float, float, int]:
    #Partial Spearman correlation of x and y controlling for one or more covariates (via rank residuals).
    mat = np.column_stack([x, y, covs])
    valid = ~np.isnan(mat).any(axis=1)
    mat = mat[valid]
    n = mat.shape[0]
    p_cov = mat.shape[1] - 2
    if n <= p_cov + 3:
        return np.nan, np.nan, n
    xr = rankdata(mat[:, 0])
    yr = rankdata(mat[:, 1])
    Z = np.column_stack([np.ones(n)] + [rankdata(mat[:, 2 + j]) for j in range(p_cov)])
    rx = xr - Z @ np.linalg.lstsq(Z, xr, rcond=None)[0]
    ry = yr - Z @ np.linalg.lstsq(Z, yr, rcond=None)[0]
    r = float(np.corrcoef(rx, ry)[0, 1])
    r_clip = np.clip(r, -0.999999, 0.999999)
    df = n - 2 - p_cov
    t_stat = r_clip * np.sqrt(df / (1.0 - r_clip ** 2))
    pval = float(2.0 * student_t.sf(abs(t_stat), df))
    return r, pval, n


def ols_regression(y: np.ndarray, X: np.ndarray, col_names: list[str]) -> dict:
    #Standard OLS regression returning standardized beta, SE, and two-sided t p-value for X[:, 0].
    mat = np.column_stack([y, X])
    valid = ~np.isnan(mat).any(axis=1)
    mat = mat[valid]
    n = mat.shape[0]
    yv = mat[:, 0]
    Xv = mat[:, 1:]
    # Standardize y and predictors so beta is directly comparable across programs
    yv = (yv - yv.mean()) / yv.std(ddof=0)
    Xv = (Xv - Xv.mean(axis=0)) / Xv.std(axis=0, ddof=0)
    design = np.column_stack([np.ones(n), Xv])
    beta, _, _, _ = np.linalg.lstsq(design, yv, rcond=None)
    resid = yv - design @ beta
    df = n - design.shape[1]
    s2 = float(np.sum(resid ** 2) / df)
    cov = s2 * np.linalg.inv(design.T @ design)
    se = np.sqrt(np.diag(cov))
    t_val = beta[1] / se[1]
    pval = float(2.0 * student_t.sf(abs(t_val), df))
    return {
        "predictor": col_names[0],
        "n": int(n),
        "std_beta": float(beta[1]),
        "se": float(se[1]),
        "p_ols": pval,
    }


def fit_cox_ph(time: np.ndarray, event: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    #Exact Breslow Cox Proportional Hazards regression via Newton-Raphson (BFGS + analytical grad/Hess).
    #Returns (beta, se, wald_p) for each column of X.
    t = np.asarray(time, float)
    e = np.asarray(event, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]

    # Sort descending by time so cumulative sums form risk sets cleanly
    order = np.argsort(-t, kind="mergesort")
    t, e, X = t[order], e[order], X[order]
    n, p = X.shape

    # Identify unique event times for Breslow tie handling
    def neg_log_lik_and_grad(beta):
        eta = np.clip(X @ beta, -30.0, 30.0)
        exp_eta = np.exp(eta)
        # Risk set sums at each index (ties share the risk set of the last observation with the same time)
        cum_exp = np.cumsum(exp_eta)
        cum_x_exp = np.cumsum(X * exp_eta[:, None], axis=0)

        # For tied times, use the maximum index among ties (largest risk set >= t_i)
        _, inv_idx = np.unique(-t, return_inverse=True)
        # Last occurrence of each unique time in sorted-descending order
        last_pos = np.zeros(len(np.unique(t)), dtype=int)
        for idx_pos, grp in enumerate(inv_idx):
            last_pos[grp] = idx_pos
        risk_idx = last_pos[inv_idx]

        denom = cum_exp[risk_idx]
        num_x = cum_x_exp[risk_idx]

        ev_mask = e == 1.0
        ll = np.sum(eta[ev_mask] - np.log(denom[ev_mask]))
        grad = np.sum(X[ev_mask] - num_x[ev_mask] / denom[ev_mask, None], axis=0)
        return -float(ll), -grad

    def hessian(beta):
        eta = np.clip(X @ beta, -30.0, 30.0)
        exp_eta = np.exp(eta)
        cum_exp = np.cumsum(exp_eta)
        cum_x_exp = np.cumsum(X * exp_eta[:, None], axis=0)
        # Outer products X_i X_i^T * exp(eta_i)
        outer = (X[:, :, None] * X[:, None, :]) * exp_eta[:, None, None]
        cum_xx_exp = np.cumsum(outer, axis=0)

        _, inv_idx = np.unique(-t, return_inverse=True)
        last_pos = np.zeros(len(np.unique(t)), dtype=int)
        for idx_pos, grp in enumerate(inv_idx):
            last_pos[grp] = idx_pos
        risk_idx = last_pos[inv_idx]

        ev_mask = e == 1.0
        r_idx = risk_idx[ev_mask]
        denom = cum_exp[r_idx]
        mean_x = cum_x_exp[r_idx] / denom[:, None]
        mean_xx = cum_xx_exp[r_idx] / denom[:, None, None]
        info = np.sum(mean_xx - mean_x[:, :, None] * mean_x[:, None, :], axis=0)
        return info

    res = minimize(neg_log_lik_and_grad, np.zeros(p), jac=True, method="BFGS")
    beta = res.x
    info = hessian(beta)
    cov = np.linalg.inv(info)
    se = np.sqrt(np.diag(cov))
    z = beta / se
    pvals = 2.0 * norm.sf(np.abs(z))
    return beta, se, pvals


def logrank_test(time: np.ndarray, event: np.ndarray, group: np.ndarray) -> tuple[float, float]:
    #Mantel-Haenszel two-group log-rank test (group is boolean/0-1). Returns (chi2_stat, p_value).
    t = np.asarray(time, float)
    e = np.asarray(event, float)
    g = np.asarray(group, int)
    event_times = np.unique(t[e == 1])
    O1 = 0.0
    E1 = 0.0
    V = 0.0
    for ti in event_times:
        at_risk = t >= ti
        n_total = np.sum(at_risk)
        n1 = np.sum(at_risk & (g == 1))
        d_total = np.sum((t == ti) & (e == 1))
        d1 = np.sum((t == ti) & (e == 1) & (g == 1))
        O1 += d1
        E1 += d_total * (n1 / n_total)
        if n_total > 1:
            V += (n1 * (n_total - n1) * d_total * (n_total - d_total)) / (n_total ** 2 * (n_total - 1))
    if V <= 0:
        return np.nan, np.nan
    stat = float(((O1 - E1) ** 2) / V)
    pval = float(chi2.sf(stat, df=1))
    return stat, pval


def kaplan_meier_curve(time: np.ndarray, event: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    #Compute step-function Kaplan-Meier survival curve.
    t = np.asarray(time, float)
    e = np.asarray(event, float)
    uniq = np.sort(np.unique(t[e == 1]))
    times = [0.0]
    surv = [1.0]
    s = 1.0
    for ti in uniq:
        n_risk = np.sum(t >= ti)
        d = np.sum((t == ti) & (e == 1))
        if n_risk > 0:
            s *= (1.0 - d / n_risk)
        times.append(float(ti))
        surv.append(float(s))
    return np.array(times), np.array(surv)


def main():
    print("=" * 76)
    print("TCGA-LIHC Bulk Validation & Survival Analysis (Frozen Discovery Signatures)")
    print("=" * 76)

    # 1. Load frozen discovery signatures
    with open("metadata/macrophage_programs_frozen_GSE149614.yaml") as fh:
        prog_frozen = yaml.safe_load(fh)
    with open("metadata/signatures_frozen_GSE149614.yaml") as fh:
        sig_frozen = yaml.safe_load(fh)

    # 2. Fetch / load UCSC Xena TCGA-LIHC files
    expr_path = fetch_file(*XENA_FILES["expr"])
    surv_path = fetch_file(*XENA_FILES["survival"])
    clin_path = fetch_file(*XENA_FILES["clinical"])

    print("Loading TCGA-LIHC matrices")
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    surv = pd.read_csv(surv_path, sep="\t", index_col=0)
    clin = pd.read_csv(clin_path, sep="\t", index_col=0)
    print(f"  Raw expression matrix: {expr.shape[0]:,} genes x {expr.shape[1]:,} samples")

    # 3. Restrict strictly to Primary Solid Tumors ('-01') and 1 sample per patient
    prim_samples = [s for s in expr.columns if re.search(r"-01[A-Z]?$", str(s))]
    expr_prim = expr[prim_samples].copy()
    # Deduplicate by 12-char TCGA patient barcode if any patient has multiple '-01' vials
    patient_ids = pd.Series([s[:12] for s in expr_prim.columns], index=expr_prim.columns)
    keep_samples = patient_ids.drop_duplicates(keep="first").index.tolist()
    expr_prim = expr_prim[keep_samples]
    print(f"  Primary Tumor samples (unique patients): {expr_prim.shape[1]:,}")

    # 4. Score frozen macrophage programs, exhaustion variants, and bulk infiltration controls
    scores = pd.DataFrame(index=expr_prim.columns)
    scores["patient"] = [s[:12] for s in expr_prim.columns]

    # Macrophage programs (filtered = primary; approved = sensitivity)
    for name, genes in prog_frozen["filtered_sets"].items():
        # Handle SEPP1 -> SELENOP alias if needed
        resolved = [g if g in expr_prim.index else ("SELENOP" if g == "SEPP1" and "SELENOP" in expr_prim.index else g) for g in genes]
        scores[f"{name}__filtered"] = score_mean_z(expr_prim, resolved, f"{name}__filtered")
    for name, genes in prog_frozen["approved_sets"].items():
        resolved = [g if g in expr_prim.index else ("SELENOP" if g == "SEPP1" and "SELENOP" in expr_prim.index else g) for g in genes]
        scores[f"{name}__approved"] = score_mean_z(expr_prim, resolved, f"{name}__approved")

    # Exhaustion signatures (V3 = primary)
    for vname, genes in sig_frozen["variants"].items():
        scores[vname] = score_mean_z(expr_prim, genes, vname)

    # Bulk lineage-abundance controls
    scores["pan_macrophage_bulk"] = score_mean_z(expr_prim, PAN_MAC_BULK, "pan_macrophage_bulk")
    scores["cd8_infiltration_bulk"] = score_mean_z(expr_prim, CD8_INFILT_BULK, "cd8_infiltration_bulk")
    scores["treg_bulk"] = score_mean_z(expr_prim, TREG_BULK, "treg_bulk")

    # Relative polarization indices & CD8-normalized exhaustion
    scores["pol_SPP1_vs_FOLR2"] = scores["SPP1_MMP__filtered"] - scores["FOLR2_C1Q__filtered"]
    scores["pol_TREM2_vs_FOLR2"] = scores["TREM2_GPNMB__filtered"] - scores["FOLR2_C1Q__filtered"]
    # Residualize V3 exhaustion on CD8 infiltration to capture exhaustion-per-CD8-cell
    x_cd8 = scores["cd8_infiltration_bulk"].to_numpy()
    y_v3 = scores["V3_noprolif_lineage"].to_numpy()
    slope, intercept = np.polyfit(x_cd8, y_v3, 1)
    scores["V3_exhaustion_per_CD8"] = y_v3 - (slope * x_cd8 + intercept)

    # 5. Merge survival and clinical stage
    scores = scores.join(surv[["OS", "OS.time", "PFI", "PFI.time"]], how="left")
    stage_col = next(
        (c for c in ["pathologic_stage", "ajcc_pathologic_tumor_stage", "tumor_stage"] if c in clin.columns),
        None,
    )
    if stage_col is not None:
        scores["raw_stage"] = clin.reindex(scores.index)[stage_col].astype(str)
        scores["stage_num"] = scores["raw_stage"].map(parse_ajcc_stage)
    else:
        scores["raw_stage"] = ""
        scores["stage_num"] = np.nan

    if "age_at_initial_pathologic_diagnosis" in clin.columns:
        scores["age"] = pd.to_numeric(clin.reindex(scores.index)["age_at_initial_pathologic_diagnosis"], errors="coerce")
    if "gender" in clin.columns:
        scores["sex"] = clin.reindex(scores.index)["gender"].astype(str)

    n_stage = int(scores["stage_num"].notna().sum())
    n_os = int(scores[["OS", "OS.time"]].dropna().query("`OS.time` > 0").shape[0])
    print(f"  Patients with AJCC pathologic stage (I-IV): {n_stage:,} / {len(scores):,}")
    print(f"  Patients with valid Overall Survival (OS.time > 0): {n_os:,} / {len(scores):,}")

    scores.to_csv(os.path.join(OUT_TABLES, "tcga_lihc_patient_scores.csv"))

    # 6. Program vs. V3 Exhaustion Associations (Raw, Infiltration-Adjusted, Infiltration+Stage-Adjusted)
    primary_progs = ["FOLR2_C1Q__filtered", "SPP1_MMP__filtered", "TREM2_GPNMB__filtered"]
    extra_progs = ["pol_SPP1_vs_FOLR2", "pol_TREM2_vs_FOLR2"]
    assoc_rows = []

    for prog in primary_progs + extra_progs:
        is_primary = prog in primary_progs
        x = scores[prog].to_numpy()
        y = scores["V3_noprolif_lineage"].to_numpy()

        # Raw Spearman
        r_raw, p_raw = spearmanr(x, y)

        # Partial Spearman controlling CD8 + pan-macrophage infiltration
        cov_inf = scores[["cd8_infiltration_bulk", "pan_macrophage_bulk"]].to_numpy()
        r_inf, p_inf, n_inf = partial_spearman(x, y, cov_inf)

        # Partial Spearman controlling CD8 + pan-macrophage + AJCC stage
        cov_full = scores[["cd8_infiltration_bulk", "pan_macrophage_bulk", "stage_num"]].to_numpy()
        r_full, p_full, n_full = partial_spearman(x, y, cov_full)

        # Sensitivity: Partial Spearman additionally controlling bulk Treg score (FOXP3, IL2RA)
        cov_treg = scores[["cd8_infiltration_bulk", "pan_macrophage_bulk", "stage_num", "treg_bulk"]].to_numpy()
        r_treg, p_treg, n_treg = partial_spearman(x, y, cov_treg)

        # Multivariable OLS: V3 ~ prog + CD8 + pan-macrophage + stage
        ols_res = ols_regression(
            y,
            np.column_stack([x, cov_full]),
            [prog, "cd8_infiltration_bulk", "pan_macrophage_bulk", "stage_num"],
        )

        # Correlation against CD8-normalized exhaustion residual
        r_per_cd8, p_per_cd8 = spearmanr(x, scores["V3_exhaustion_per_CD8"].to_numpy())

        assoc_rows.append({
            "tier": "PRIMARY" if is_primary else "POLARIZATION_INDEX",
            "program": prog,
            "n_all": len(scores),
            "rho_raw": round(float(r_raw), 4),
            "p_raw": float(p_raw),
            "rho_adj_infiltration": round(float(r_inf), 4),
            "p_adj_infiltration": float(p_inf),
            "n_adj_infiltration_stage": n_full,
            "rho_adj_infiltration_stage": round(float(r_full), 4),
            "p_adj_infiltration_stage": float(p_full),
            "rho_adj_inf_stage_treg": round(float(r_treg), 4),
            "p_adj_inf_stage_treg": float(p_treg),
            "ols_std_beta_full": round(ols_res["std_beta"], 4),
            "ols_p_full": ols_res["p_ols"],
            "rho_vs_V3_per_CD8": round(float(r_per_cd8), 4),
            "p_vs_V3_per_CD8": float(p_per_cd8),
        })

    assoc_df = pd.DataFrame(assoc_rows)
    prim_mask = assoc_df["tier"] == "PRIMARY"
    assoc_df["p_holm_adj_infiltration_stage"] = np.nan
    assoc_df.loc[prim_mask, "p_holm_adj_infiltration_stage"] = holm_adjust(
        assoc_df.loc[prim_mask, "p_adj_infiltration_stage"].to_numpy()
    )
    assoc_df.to_csv(os.path.join(OUT_TABLES, "tcga_lihc_program_exhaustion_associations.csv"), index=False)

    print("TCGA-LIHC: Macrophage Programs vs. Frozen V3 CD8 Exhaustion")
    show_cols = [
        "tier", "program", "rho_raw", "rho_adj_infiltration",
        "rho_adj_infiltration_stage", "p_adj_infiltration_stage", "p_holm_adj_infiltration_stage",
        "rho_adj_inf_stage_treg", "p_adj_inf_stage_treg",
        "ols_std_beta_full", "rho_vs_V3_per_CD8",
    ]
    print(assoc_df[show_cols].to_string(index=False))

    # 7. Overall Survival Analysis (Univariable & Stage-Adjusted Cox PH + Log-Rank)
    os_df = scores.dropna(subset=["OS", "OS.time"]).copy()
    os_df = os_df[os_df["OS.time"] > 0].copy()
    os_df["OS_years"] = os_df["OS.time"] / 365.25

    surv_features = primary_progs + extra_progs + ["V3_noprolif_lineage", "V3_exhaustion_per_CD8"]
    surv_rows = []
    for feat in surv_features:
        is_primary = feat in primary_progs
        x_raw = os_df[feat].to_numpy(float)
        x_z = (x_raw - x_raw.mean()) / x_raw.std(ddof=0)  # HR per 1-SD increase
        t_arr = os_df["OS.time"].to_numpy(float)
        e_arr = os_df["OS"].to_numpy(float)

        # Univariable Cox PH
        b_uni, se_uni, p_uni = fit_cox_ph(t_arr, e_arr, x_z[:, None])
        hr_uni = float(np.exp(b_uni[0]))
        ci_lo_uni = float(np.exp(b_uni[0] - 1.96 * se_uni[0]))
        ci_hi_uni = float(np.exp(b_uni[0] + 1.96 * se_uni[0]))

        # Median-split Log-Rank test
        hi_grp = (x_raw >= np.median(x_raw)).astype(int)
        lr_stat, lr_p = logrank_test(t_arr, e_arr, hi_grp)

        # Stage-adjusted Multivariable Cox PH
        sub_st = os_df.dropna(subset=["stage_num"]).copy()
        xs = sub_st[feat].to_numpy(float)
        xs_z = (xs - xs.mean()) / xs.std(ddof=0)
        st_arr = sub_st["stage_num"].to_numpy(float)
        b_adj, se_adj, p_adj = fit_cox_ph(
            sub_st["OS.time"].to_numpy(float),
            sub_st["OS"].to_numpy(float),
            np.column_stack([xs_z, st_arr]),
        )
        hr_adj = float(np.exp(b_adj[0]))
        ci_lo_adj = float(np.exp(b_adj[0] - 1.96 * se_adj[0]))
        ci_hi_adj = float(np.exp(b_adj[0] + 1.96 * se_adj[0]))

        surv_rows.append({
            "tier": "PRIMARY" if is_primary else "SECONDARY",
            "feature": feat,
            "n_os": int(len(os_df)),
            "events_os": int(e_arr.sum()),
            "HR_univar_per_SD": round(hr_uni, 3),
            "CI95_univar": f"{ci_lo_uni:.3f}-{ci_hi_uni:.3f}",
            "p_cox_univar": float(p_uni[0]),
            "p_logrank_median": float(lr_p),
            "n_stage_adj": int(len(sub_st)),
            "HR_stage_adj_per_SD": round(hr_adj, 3),
            "CI95_stage_adj": f"{ci_lo_adj:.3f}-{ci_hi_adj:.3f}",
            "p_cox_stage_adj": float(p_adj[0]),
        })

    surv_res = pd.DataFrame(surv_rows)
    pm_s = surv_res["tier"] == "PRIMARY"
    surv_res["p_holm_stage_adj"] = np.nan
    surv_res.loc[pm_s, "p_holm_stage_adj"] = holm_adjust(surv_res.loc[pm_s, "p_cox_stage_adj"].to_numpy())
    surv_res.to_csv(os.path.join(OUT_TABLES, "tcga_lihc_survival_cox.csv"), index=False)

    print("TCGA-LIHC Overall Survival (Cox PH per 1-SD & Median-Split Log-Rank)")
    print(surv_res.to_string(index=False))
    print(f"Saved TCGA-LIHC tables to {OUT_TABLES}/")

# 8. Generate 4-panel validation figure (figures/tcga_lihc_validation.pdf)
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 9.6))

    # Panel A: Partial Spearman rho (adj. CD8 + Pan-Mac + Stage)
    ax = axes[0, 0]
    label_map = {
        "FOLR2_C1Q__filtered": "FOLR2_C1Q",
        "SPP1_MMP__filtered": "SPP1_MMP",
        "TREM2_GPNMB__filtered": "TREM2_GPNMB",
        "pol_SPP1_vs_FOLR2": "SPP1 - FOLR2",
        "pol_TREM2_vs_FOLR2": "TREM2 - FOLR2",
    }
    colors_a = ["#2B6CB0", "#C53030", "#D69E2E", "#8B1E3F", "#B7791F"]
    y_pos = np.arange(len(assoc_df))
    rhos = assoc_df["rho_adj_infiltration_stage"].to_numpy(float)
    pvals = assoc_df["p_adj_infiltration_stage"].to_numpy(float)
    ax.barh(y_pos, rhos, color=colors_a[:len(y_pos)], alpha=0.88, height=0.6)
    ax.axvline(0.0, color="#4A5568", lw=1.0, ls="--")
    for yi, r_val, p_val in zip(y_pos, rhos, pvals):
        offset = 0.015 if r_val >= 0 else -0.015
        ha = "left" if r_val >= 0 else "right"
        ax.text(r_val + offset, yi, f"r={r_val:+.2f} (p={p_val:.1e})", va="center", ha=ha, fontsize=8.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([label_map.get(p, p) for p in assoc_df["program"]])
    ax.invert_yaxis()
    ax.set_xlim(-0.42, 0.48)
    ax.set_xlabel("Partial Spearman rho (adj. CD8 + Pan-Mac + Stage)")
    ax.set_title("A. Bulk Program vs. V3 Exhaustion (TCGA-LIHC)", fontweight="bold")

    # Panel B: Polarization vs. CD8-Normalized V3 Exhaustion Residual
    ax = axes[0, 1]
    x_pol = scores["pol_SPP1_vs_FOLR2"].to_numpy(float)
    y_res = scores["V3_exhaustion_per_CD8"].to_numpy(float)
    r_b, p_b = spearmanr(x_pol, y_res)
    ax.scatter(x_pol, y_res, s=22, color="#4A5568", alpha=0.45, edgecolors="none")
    m_b, b_b = np.polyfit(x_pol, y_res, 1)
    x_grid = np.linspace(x_pol.min(), x_pol.max(), 80)
    ax.plot(x_grid, m_b * x_grid + b_b, color="#C53030", lw=2.2)
    ax.set_xlabel("Macrophage Polarization (z[SPP1_MMP] - z[FOLR2_C1Q])")
    ax.set_ylabel("CD8-Normalized V3 Exhaustion (Residual)")
    ax.set_title(f"B. Polarization vs. Exhaustion/CD8 (rho={r_b:+.3f}, p={p_b:.2e})", fontweight="bold")

    # Panels C & D: Kaplan-Meier Overall Survival curves (median split)
    km_specs = [
        (axes[1, 0], "SPP1_MMP__filtered", "SPP1_MMP", "C. TCGA-LIHC OS: SPP1_MMP"),
        (axes[1, 1], "pol_SPP1_vs_FOLR2", "SPP1 > FOLR2", "D. TCGA-LIHC OS: SPP1 - FOLR2 Axis"),
    ]
    for ax, col, leg_prefix, title_prefix in km_specs:
        vals = os_df[col].to_numpy(float)
        hi = vals >= np.median(vals)
        t_yrs = os_df["OS_years"].to_numpy(float)
        ev = os_df["OS"].to_numpy(float)
        _, lr_p = logrank_test(os_df["OS.time"].to_numpy(float), ev, hi.astype(int))

        t_hi, s_hi = kaplan_meier_curve(t_yrs[hi], ev[hi])
        t_lo, s_lo = kaplan_meier_curve(t_yrs[~hi], ev[~hi])
        ax.step(t_hi, s_hi, where="post", color="#B22222", lw=2.0, label=f"{leg_prefix} High (n={int(hi.sum())})")
        ax.step(t_lo, s_lo, where="post", color="#2B6CB0", lw=2.0, label=f"{leg_prefix} Low (n={int((~hi).sum())})")
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel("Overall Survival (Years)")
        ax.set_ylabel("Survival Probability")
        ax.set_title(f"{title_prefix} (log-rank p={lr_p:.3g})", fontweight="bold")
        ax.legend(frameon=False, loc="upper right")

    plt.tight_layout()
    pdf_path = os.path.join(OUT_FIGS, "tcga_lihc_validation.pdf")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved TCGA-LIHC validation figure to {pdf_path}")


if __name__ == "__main__":
    main()