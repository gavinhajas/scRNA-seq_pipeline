# 14_make_summary_figures.py: Generate publication summary figures (PDF & 300-dpi PNG in figures/).
# Usage: python scripts/14_make_summary_figures.py [--tables-dir results/tables] [--figs-dir figures]
#
# Outputs:
# 1. fig1_scrna_patient_scatter_grid.{pdf,png}: 3x3 grid of patient macrophage program means vs. CD8 V3 exhaustion
#    across GSE149614 (n=8), GSE156625 (n=9), and GSE140228 (n=4), point-sized by n_mac with size legend.
# 2. fig2_pooled_fisher_z_forest.{pdf,png}: 3-panel Spearman correlation forest plot (excludes n<5 from Fisher-z),
#    separating Validation-only (Fisher-z & Stouffer exact-perm p) and Secondary All-Cohorts (k=2, reporting I^2).
# 3. fig3_cross_cohort_communication.{pdf,png}: NicheNet v2 scatter + 2-cohort faceted CellChat paired dot plot.

from __future__ import annotations

import argparse
import os

import matplotlib
matplotlib.use("Agg")
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from adjustText import adjust_text  # type: ignore
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False

# Configuration & Constants
DEFAULT_TABLES_DIR = "results/tables"
DEFAULT_FIGS_DIR = "figures"

MIN_MEAN_THRESHOLD = 0.05
MIN_FISHER_N = 5
OLS_P_THRESHOLD = 0.10      # Only draw OLS visual guide when n >= MIN_FISHER_N and p_perm < OLS_P_THRESHOLD
NICHENET_TOP_N = 6          # Number of top shared NicheNet ligands to label cleanly
SIZE_LEGEND_VALUES = (50, 250, 1000)

PROG_SPECS = [
    ("FOLR2_C1Q__filtered_UCell", "FOLR2 / C1Q Program", "#2B6CB0"),
    ("SPP1_MMP__filtered_UCell", "SPP1 / MMP Program", "#C53030"),
    ("TREM2_GPNMB__filtered_UCell", "TREM2 / GPNMB Program", "#D69E2E"),
]

COHORT_SPECS = [
    ("GSE149614", "GSE149614 (Discovery, unsorted)"),
    ("GSE156625", "GSE156625 (Validation, 1:1 CD45)"),
    ("GSE140228", "GSE140228 (Supportive, CD45+ Stage II)"),
]

MAC_PROG_COLORS = {
    "Mac_FOLR2_C1Q": "#2B6CB0",
    "Mac_SPP1_MMP": "#C53030",
    "Mac_TREM2_GPNMB": "#D69E2E",
    "Mac_mixed_or_low": "#718096",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
})


# Shared Utilities
def mac_size_scale(n_mac: np.ndarray | float) -> np.ndarray:
    """Map macrophage cell counts to scatter marker area (shared by plots and legend)."""
    arr = np.asarray(n_mac, dtype=float)
    return 35.0 + 210.0 * np.sqrt(np.clip(arr, 20.0, 1000.0) / 1000.0)


def save_fig(fig: plt.Figure, figs_dir: str, stem: str) -> None:
    """Save editable vector PDF (TrueType fonttype 42) and 300-dpi PNG with safe padding."""
    pdf_path = os.path.join(figs_dir, f"{stem}.pdf")
    png_path = os.path.join(figs_dir, f"{stem}.png")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight", pad_inches=0.18)
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)
    print(f"  Saved {pdf_path} and {png_path}")


def fisher_ci(rho: float, n: int) -> tuple[float, float]:
    """Compute 95% Fisher-z CI only when n >= MIN_FISHER_N and |rho| < 1."""
    if np.isnan(rho) or n < MIN_FISHER_N or abs(rho) >= 1.0:
        return np.nan, np.nan
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - 3.0)
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def derive_exclusion_note(n_pat: int, prog_mean: float, stats_row: pd.Series | None = None) -> str:
    """Derive cohort/program status note dynamically from table metrics."""
    if stats_row is not None:
        for col in ("mean_exposure", "mean_program", "program_mean", "mean_ucell"):
            if col in stats_row.index and pd.notna(stats_row[col]):
                prog_mean = float(stats_row[col])
                break
    if np.isfinite(prog_mean) and prog_mean < MIN_MEAN_THRESHOLD:
        return f"excl.: mean<{MIN_MEAN_THRESHOLD:.2f}"
    if n_pat < MIN_FISHER_N:
        return f"n<{MIN_FISHER_N}: Stouffer only"
    return ""


def annotate_non_overlapping(
    ax: plt.Axes,
    xs: np.ndarray,
    ys: np.ndarray,
    labels: list[str],
    fontsize: float = 7.8,
    fontweight: str = "normal",
    color: str = "#2D3748",
) -> None:
    """Place point labels cleanly using adjustText when available, or a greedy offset fallback."""
    valid_pts = [
        (float(x), float(y), str(lbl))
        for x, y, lbl in zip(xs, ys, labels)
        if np.isfinite(x) and np.isfinite(y)
    ]
    if not valid_pts:
        return

    if HAS_ADJUST_TEXT:
        texts = [
            ax.text(x, y, lbl, fontsize=fontsize, fontweight=fontweight, color=color, zorder=6)
            for x, y, lbl in valid_pts
        ]
        adjust_text(
            texts,
            x=[p[0] for p in valid_pts],
            y=[p[1] for p in valid_pts],
            ax=ax,
            expand=(1.25, 1.35),
            arrowprops=dict(arrowstyle="-", color="#A0AEC0", lw=0.6, shrinkA=2, shrinkB=2),
        )
        return

    # Deterministic pure-matplotlib fallback: choose the candidate offset farthest from
    # nearby points and previously placed labels in normalized axis coordinates.
    x_lo, x_hi = ax.get_xlim()
    y_lo, y_hi = ax.get_ylim()
    x_span = max(x_hi - x_lo, 1e-9)
    y_span = max(y_hi - y_lo, 1e-9)

    pts_norm = np.array([[(x - x_lo) / x_span, (y - y_lo) / y_span] for x, y, _ in valid_pts])
    # (dx_pts, dy_pts, ha, va, approx_norm_dx, approx_norm_dy)
    candidates = [
        (6, 5, "left", "bottom", 0.045, 0.045),
        (-6, 5, "right", "bottom", -0.045, 0.045),
        (6, -6, "left", "top", 0.045, -0.045),
        (-6, -6, "right", "top", -0.045, -0.045),
        (0, 8, "center", "bottom", 0.0, 0.065),
        (0, -9, "center", "top", 0.0, -0.065),
        (9, 0, "left", "center", 0.065, 0.0),
        (-9, 0, "right", "center", -0.065, 0.0),
    ]
    placed_norm: list[np.ndarray] = []

    for i, (x, y, lbl) in enumerate(valid_pts):
        p = pts_norm[i]
        best_cand = candidates[0]
        best_score = -1e9

        for cand in candidates:
            _, _, _, _, ndx, ndy = cand
            pos = p + np.array([ndx, ndy])
            # Penalize placing a label outside the plot area, near the left spine, or inside the stat box
            edge_pen = 0.0
            if pos[0] < 0.14 or pos[0] > 0.90 or pos[1] < 0.05 or pos[1] > 0.93:
                edge_pen -= 0.65
            if pos[0] < 0.42 and pos[1] > 0.68:
                edge_pen -= 0.35

            # Minimum distance to all other points
            other_pts = np.delete(pts_norm, i, axis=0)
            min_pt_dist = float(np.min(np.linalg.norm(other_pts - pos, axis=1))) if len(other_pts) else 1.0

            # Minimum distance to already-placed labels (weighted slightly higher)
            if placed_norm:
                min_lbl_dist = float(np.min(np.linalg.norm(np.vstack(placed_norm) - pos, axis=1)))
            else:
                min_lbl_dist = 1.0

            score = min(min_pt_dist, 1.15 * min_lbl_dist) + edge_pen
            if score > best_score:
                best_score = score
                best_cand = cand

        dx, dy, ha, va, ndx, ndy = best_cand
        placed_norm.append(p + np.array([ndx, ndy]))
        ax.annotate(
            lbl,
            (x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            va=va,
            fontsize=fontsize,
            fontweight=fontweight,
            color=color,
            zorder=6,
        )


# Figure 1: 3 x 3 Single-Cell Patient-Level Scatter Grid
def make_fig1_scatter_grid(tables_dir: str, figs_dir: str) -> None:
    print("\nBuilding Figure 1: scRNA-seq patient-level scatter grid...")
    fig, axes = plt.subplots(3, 3, figsize=(13.8, 11.4), sharey="row")

    for r_idx, (gse, cohort_title) in enumerate(COHORT_SPECS):
        ptab = pd.read_csv(os.path.join(tables_dir, f"patient_level_table_primary_{gse}.csv"), index_col=0)
        stats = pd.read_csv(os.path.join(tables_dir, f"stats_primary_{gse}.csv")).set_index("program")

        n_mac = ptab["n_mac"].to_numpy(float)
        sizes = mac_size_scale(n_mac)
        print(f"  {gse}: {len(ptab)} patients (n_mac range: {np.nanmin(n_mac):.0f}–{np.nanmax(n_mac):.0f})")

        for c_idx, (prog_col, prog_label, color) in enumerate(PROG_SPECS):
            ax = axes[r_idx, c_idx]
            x = ptab[prog_col].to_numpy(float)
            y = ptab["exhaustion"].to_numpy(float)
            valid = np.isfinite(x) & np.isfinite(y)
            xv, yv = x[valid], y[valid]
            if not np.all(valid):
                print(f"    [WARN] {gse} / {prog_col}: dropped {int((~valid).sum())} non-finite patient(s)")

            row = stats.loc[prog_col]
            rho = float(row["rho"])
            p_perm = float(row["p_perm"])
            n_pat = int(row["n_patients"])
            prog_mean = float(np.nanmean(xv)) if len(xv) > 0 else np.nan
            note = derive_exclusion_note(n_pat, prog_mean, row)

            # Draw OLS visual guide only when n >= MIN_FISHER_N and p_perm < OLS_P_THRESHOLD
            if len(xv) >= MIN_FISHER_N and np.std(xv) > 0 and p_perm < OLS_P_THRESHOLD:
                slope, intercept = np.polyfit(xv, yv, 1)
                x_line = np.linspace(xv.min(), xv.max(), 60)
                ax.plot(x_line, slope * x_line + intercept, color=color, lw=1.6, ls="--", alpha=0.70, zorder=2)

            ax.scatter(x, y, s=sizes, color=color, alpha=0.84, edgecolors="white", linewidth=0.9, zorder=3)

            # Pad both axes before placing labels so edge points (e.g., D18-0116) never clip
            x_min, x_max = float(np.nanmin(xv)), float(np.nanmax(xv))
            x_span = max(0.02, x_max - x_min)
            ax.set_xlim(x_min - 0.14 * x_span, x_max + 0.16 * x_span)

            y_min, y_max = float(np.nanmin(yv)), float(np.nanmax(yv))
            y_span = max(0.02, y_max - y_min)
            ax.set_ylim(y_min - 0.18 * y_span, y_max + 0.38 * y_span)

            short_labels = [
                str(pat).replace("D2017", "D17-").replace("D2018", "D18-")
                for pat in ptab.index
            ]
            annotate_non_overlapping(ax, x, y, short_labels, fontsize=7.6, color="#2D3748")

            box_txt = f"rho = {rho:+.2f}\np = {p_perm:.3f} (n={n_pat})"
            if note:
                box_txt += f"\n[{note}]"
            # Keep stat box clear of points in the upper-left corner (e.g., GSE140228 FOLR2 D18-0116)
            x_lo, x_hi = ax.get_xlim()
            y_lo, y_hi = ax.get_ylim()
            ul_occupied = np.any(
                ((xv - x_lo) / max(x_hi - x_lo, 1e-9) < 0.38)
                & ((yv - y_lo) / max(y_hi - y_lo, 1e-9) > 0.68)
            )
            box_x, box_ha = (0.96, "right") if ul_occupied else (0.04, "left")
            ax.text(
                box_x, 0.95, box_txt, transform=ax.transAxes, va="top", ha=box_ha, fontsize=8.4,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#CBD5E0", alpha=0.92),
                zorder=7,
            )

            ax.grid(True, ls=":", alpha=0.35)
            if r_idx == 0:
                ax.set_title(prog_label, fontsize=11.5, fontweight="bold", color=color, pad=8)
            if r_idx == 2:
                ax.set_xlabel("Patient Mean Macrophage Program (UCell)", fontsize=9.5)
            if c_idx == 0:
                ax.set_ylabel(f"{cohort_title}\nMedian CD8 V3 Exhaustion", fontsize=9.2, fontweight="bold")

    # Shared bubble-size legend for n_mac
    size_handles = [
        ax_leg_marker(mac_size_scale(val), f"{val:,} cells")
        for val in SIZE_LEGEND_VALUES
    ]
    fig.legend(
        handles=size_handles,
        title="Tumor Macrophage Count (n_mac)",
        loc="lower center",
        bbox_to_anchor=(0.5, 0.028),
        ncol=len(SIZE_LEGEND_VALUES),
        frameon=True,
        facecolor="white",
        edgecolor="#CBD5E0",
        fontsize=8.6,
        title_fontsize=8.8,
        handletextpad=0.7,
        columnspacing=1.6,
    )

    fig.suptitle(
        "Patient-Level Macrophage Programs vs. Strict-Gate CD8+ T Cell Exhaustion (Tumor Only)",
        fontsize=13, fontweight="bold", y=0.992,
    )
    fig.text(
        0.5, 0.005,
        f"Dashed lines: OLS visual guides (shown only for n >= {MIN_FISHER_N} and p < {OLS_P_THRESHOLD:.2f}). "
        "Tests: exact permutation Spearman correlation.",
        ha="center", fontsize=8.8, color="#4A5568",
    )
    plt.tight_layout(rect=[0, 0.085, 1, 0.968])
    save_fig(fig, figs_dir, "fig1_scrna_patient_scatter_grid")


def ax_leg_marker(area: float, label: str) -> Line2D:
    """Create a legend handle matching matplotlib scatter marker area (in pt^2)."""
    return Line2D(
        [0], [0],
        marker="o",
        color="none",
        markerfacecolor="#718096",
        markeredgecolor="white",
        markeredgewidth=0.9,
        markersize=float(np.sqrt(area)),
        alpha=0.84,
        label=label,
    )


# Figure 2: Validation-Only (Replication) & Secondary All-Cohort Forest Plot

def make_fig2_forest_plot(tables_dir: str, figs_dir: str) -> None:
    print("\nBuilding Figure 2: Validation-only (replication) & secondary Fisher-z forest plot...")
    pooled = pd.read_csv(os.path.join(tables_dir, "stats_pooled_primary.csv"))
    val_pool = pooled[pooled["tier"] == "CONFIRMATORY_VALIDATION_ONLY"].set_index("program")
    all_pool = pooled[pooled["tier"] == "SECONDARY_ALL_COHORTS"].set_index("program")

    cohort_stats = {}
    cohort_means = {}
    for gse, _ in COHORT_SPECS:
        cohort_stats[gse] = pd.read_csv(os.path.join(tables_dir, f"stats_primary_{gse}.csv")).set_index("program")
        ptab = pd.read_csv(os.path.join(tables_dir, f"patient_level_table_primary_{gse}.csv"), index_col=0)
        cohort_means[gse] = {
            prog_col: float(np.nanmean(ptab[prog_col].to_numpy(float)))
            for prog_col, _, _ in PROG_SPECS
        }

    fig, axes = plt.subplots(3, 1, figsize=(12.2, 10.8), sharex=True)

    for idx, (prog_col, prog_label, color) in enumerate(PROG_SPECS):
        ax = axes[idx]
        vrow = val_pool.loc[prog_col]
        arow = all_pool.loc[prog_col]

        y_labels: list[str] = []
        y_positions: list[int] = []
        y_cursor = 0

        # Individual cohorts
        val_plus_supp_n = 0
        for gse, cohort_title in COHORT_SPECS:
            crow = cohort_stats[gse].loc[prog_col]
            rho = float(crow["rho"])
            p_perm = float(crow["p_perm"])
            n_pat = int(crow["n_patients"])
            mean_exp = (
                float(crow["mean_exposure"])
                if "mean_exposure" in crow.index and pd.notna(crow["mean_exposure"])
                else cohort_means[gse][prog_col]
            )
            if gse in ("GSE156625", "GSE140228") and mean_exp >= MIN_MEAN_THRESHOLD:
                val_plus_supp_n += n_pat
            lo, hi = fisher_ci(rho, n_pat)

            y_positions.append(y_cursor)
            y_labels.append(f"{cohort_title} (n={n_pat})")

            if not np.isnan(lo):
                ax.plot([lo, hi], [y_cursor, y_cursor], color=color, lw=1.8, alpha=0.85)
                ax.scatter([rho], [y_cursor], color=color, marker="o", s=55 + 4 * n_pat, zorder=4)
                txt = f"{rho:+.2f} [{lo:+.2f}, {hi:+.2f}]   p={p_perm:.3f}"
            else:
                ax.scatter([rho], [y_cursor], color="#718096", marker="s", s=55, zorder=4)
                note = derive_exclusion_note(n_pat, cohort_means[gse][prog_col], crow)
                txt = f"{rho:+.2f} (no CI; {note})   p={p_perm:.3f}"
            ax.text(1.06, y_cursor, txt, va="center", ha="left", fontsize=8.4, color="#1A202C")
            y_cursor += 1

        ax.axhline(y_cursor - 0.5, color="#CBD5E0", lw=1.0)

        # Validation-only (replication): Fisher-z
        rho_v = float(vrow["rho_fixed"])
        lo_v, hi_v = float(vrow["ci95_low_fixed"]), float(vrow["ci95_high_fixed"])
        n_v = int(vrow["total_patients"]) if "total_patients" in vrow.index and pd.notna(vrow["total_patients"]) else 9
        k_v = (
            int(vrow["n_cohorts"])
            if "n_cohorts" in vrow.index and pd.notna(vrow["n_cohorts"])
            else (int(vrow["k_cohorts"]) if "k_cohorts" in vrow.index and pd.notna(vrow["k_cohorts"]) else 1)
        )
        p_v = float(vrow["p_fixed"]) if "p_fixed" in vrow.index and pd.notna(vrow["p_fixed"]) else np.nan

        y_positions.append(y_cursor)
        y_labels.append(f"Validation: Fisher-z (k={k_v}, N={n_v})")
        ax.plot([lo_v, hi_v], [y_cursor, y_cursor], color="#1A202C", lw=2.2)
        ax.scatter([rho_v], [y_cursor], color=color, edgecolors="#1A202C", marker="D", s=95, zorder=5)
        v_txt = f"{rho_v:+.2f} [{lo_v:+.2f}, {hi_v:+.2f}]"
        if np.isfinite(p_v):
            v_txt += f"   p={p_v:.3f}"
        ax.text(1.06, y_cursor, v_txt, va="center", ha="left", fontsize=8.4, fontweight="bold")
        y_cursor += 1

        # Validation-only (replication): Stouffer exact-perm
        p_st_v, ph_st_v = float(vrow["stouffer_p"]), float(vrow["p_holm_stouffer"])
        st_n_v = int(vrow["stouffer_n"])
        if st_n_v != val_plus_supp_n:
            print(
                f"  [WARN] {prog_label}: stats_pooled_primary.csv has stouffer_n={st_n_v}, "
                f"but GSE156625 + GSE140228 has N={val_plus_supp_n}. "
                "Verify whether GSE140228 was excluded in script 09 (e.g., mean threshold)."
            )

        y_positions.append(y_cursor)
        y_labels.append(f"Validation: Stouffer exact-perm (N={st_n_v})")
        ax.text(
            1.06, y_cursor,
            f"Stouffer p={p_st_v:.3f}  (Holm {ph_st_v:.3f})",
            va="center", ha="left", fontsize=8.4, fontweight="bold", color="#2D3748",
        )
        y_cursor += 1

        # Secondary: All Cohorts incl. Discovery
        rho_a = float(arow["rho_fixed"])
        lo_a, hi_a = float(arow["ci95_low_fixed"]), float(arow["ci95_high_fixed"])
        p_a, ph_a = float(arow["p_fixed"]), float(arow["p_holm_fixed"])
        i2_a = float(arow["I2_percent"]) if pd.notna(arow["I2_percent"]) else np.nan
        n_a = int(arow["total_patients"])
        i2_str = f"{i2_a:.1f}%" if np.isfinite(i2_a) else "N/A"

        y_positions.append(y_cursor)
        y_labels.append(f"All cohorts: Fisher-z (N={n_a}, $I^2$={i2_str}*)")
        ax.plot([lo_a, hi_a], [y_cursor, y_cursor], color="#4A5568", lw=1.8, ls="--")
        ax.scatter([rho_a], [y_cursor], color="white", edgecolors="#1A202C", marker="D", s=85, lw=1.6, zorder=5)
        ax.text(
            1.06, y_cursor,
            f"{rho_a:+.2f} [{lo_a:+.2f}, {hi_a:+.2f}]   p={p_a:.3f} (Holm {ph_a:.3f})",
            va="center", ha="left", fontsize=8.4, color="#2D3748",
        )

        ax.axvline(0.0, color="#718096", lw=1.0, ls=":")
        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=8.9)
        ax.invert_yaxis()
        ax.set_xlim(-1.05, 1.05)
        ax.set_title(prog_label, loc="left", fontsize=11.2, fontweight="bold", color=color, pad=6)
        ax.grid(axis="x", ls=":", alpha=0.35)

        if idx == 0:
            ax.text(
                1.06, -0.78, "Estimate [95% CI] / p-value",
                va="bottom", ha="left", fontsize=8.6, fontweight="bold", color="#4A5568",
            )

    axes[-1].set_xlabel(f"Spearman Correlation (rho) with 95% Fisher-z CI (n >= {MIN_FISHER_N})", fontsize=10.2)
    fig.suptitle(
        "Single-Cell Cross-Cohort Summary: Validation-Only (Replication) vs. Secondary All-Cohort Pool",
        fontsize=12.4, fontweight="bold", y=0.985,
    )
    fig.text(
        0.5, 0.012,
        "*All-cohort $I^2$ is descriptive (k=2 cohorts with n >= 5); n < 5 cohorts with mean >= 0.05 contribute to Stouffer exact-perm p-values.",
        ha="center", fontsize=8.6, color="#4A5568",
    )
    plt.subplots_adjust(left=0.30, right=0.64, top=0.92, bottom=0.085, hspace=0.40)
    save_fig(fig, figs_dir, "fig2_pooled_fisher_z_forest")



# Figure 3: Cross-Cohort CellChat & NicheNet Mechanistic Replication
def make_fig3_communication(tables_dir: str, figs_dir: str) -> None:
    print("\nBuilding Figure 3: Cross-cohort NicheNet & CellChat mechanistic replication...")
    nn1_path = os.path.join(tables_dir, "nichenet_ligand_activities_V3_GSE149614.csv")
    nn2_path = os.path.join(tables_dir, "nichenet_ligand_activities_V3_GSE156625.csv")
    cc1_path = os.path.join(tables_dir, "cellchat_pathways_mac_to_cd8_GSE149614.csv")
    cc2_path = os.path.join(tables_dir, "cellchat_pathways_mac_to_cd8_GSE156625.csv")

    if not all(os.path.exists(p) for p in [nn1_path, nn2_path, cc1_path, cc2_path]):
        print("  Skipping Figure 3 (missing CellChat or NicheNet CSV tables).")
        return

    nn1 = pd.read_csv(nn1_path)
    nn2 = pd.read_csv(nn2_path)
    cc1 = pd.read_csv(cc1_path)
    cc2 = pd.read_csv(cc2_path)

    fig = plt.figure(figsize=(15.4, 6.4))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.35], wspace=0.26)
    ax1 = fig.add_subplot(outer[0])
    inner_b = outer[1].subgridspec(1, 2, wspace=0.12)
    ax2_a = fig.add_subplot(inner_b[0])
    ax2_b = fig.add_subplot(inner_b[1], sharey=ax2_a)

    # --- Panel A: NicheNet (Spearman correlation for skewed AUPR values) ---
    m_nn = nn1[["test_ligand", "aupr_corrected", "dominant_mac_program"]].merge(
        nn2[["test_ligand", "aupr_corrected", "dominant_mac_program"]],
        on="test_ligand",
        suffixes=("_149614", "_156625"),
    )
    print(f"  NicheNet shared ligands: {len(m_nn)} (from {len(nn1)} GSE149614, {len(nn2)} GSE156625)")

    for prog_name, col in MAC_PROG_COLORS.items():
        sub = m_nn[m_nn["dominant_mac_program_149614"] == prog_name]
        if sub.empty:
            continue
        ax1.scatter(
            sub["aupr_corrected_149614"], sub["aupr_corrected_156625"],
            c=col, s=50, alpha=0.86, edgecolors="white", lw=0.6, label=prog_name, zorder=3,
        )

    # Pad limits before annotating top ligands so HLA-E never clips
    x_nn = m_nn["aupr_corrected_149614"].to_numpy(float)
    y_nn = m_nn["aupr_corrected_156625"].to_numpy(float)
    x_pad = 0.12 * max(float(np.nanmax(x_nn) - np.nanmin(x_nn)), 1e-4)
    y_pad = 0.10 * max(float(np.nanmax(y_nn) - np.nanmin(y_nn)), 1e-4)
    ax1.set_xlim(float(np.nanmin(x_nn)) - 0.5 * x_pad, float(np.nanmax(x_nn)) + x_pad)
    ax1.set_ylim(float(np.nanmin(y_nn)) - 0.5 * y_pad, float(np.nanmax(y_nn)) + y_pad)

    top_both = (
        m_nn.assign(mean_aupr=0.5 * (m_nn["aupr_corrected_149614"] + m_nn["aupr_corrected_156625"]))
        .sort_values("mean_aupr", ascending=False)
        .head(NICHENET_TOP_N)
    )
    annotate_non_overlapping(
        ax1,
        top_both["aupr_corrected_149614"].to_numpy(float),
        top_both["aupr_corrected_156625"].to_numpy(float),
        top_both["test_ligand"].astype(str).tolist(),
        fontsize=8.6,
        fontweight="bold",
        color="#1A202C",
    )

    rho_nn = m_nn["aupr_corrected_149614"].corr(m_nn["aupr_corrected_156625"], method="spearman")
    ax1.set_xlabel("GSE149614 Ligand Activity (aupr_corrected)")
    ax1.set_ylabel("GSE156625 Ligand Activity (aupr_corrected)")
    ax1.set_title(
        f"A. NicheNet V3 Exhaustion Ligands\n(Spearman rho={rho_nn:.2f}, n={len(m_nn)} shared; top {NICHENET_TOP_N} labeled)",
        fontweight="bold", fontsize=10.4, pad=8,
    )
    ax1.legend(title="Discovery Dominant Program", loc="upper left", frameon=False, fontsize=8.4, title_fontsize=8.4)
    ax1.grid(True, ls=":", alpha=0.35)

    # --- Panel B: 2-Cohort Faceted CellChat Paired Dot / Dumbbell Plot ---
    mac_order = ["Mac_FOLR2_C1Q", "Mac_SPP1_MMP", "Mac_TREM2_GPNMB"]

    def prep_cc(df: pd.DataFrame, cohort_label: str) -> tuple[pd.DataFrame, float]:
        sub = df[(df["target"].isin(["CD8_Exh_High", "CD8_Exh_Low"])) & (df["pathway_name"] != "MHC-I")].copy()
        sub = sub[sub["source"].isin(mac_order)]
        sub["cohort"] = cohort_label
        high_sub = sub[sub["target"] == "CD8_Exh_High"]
        max_p = float(high_sub["prob"].max()) if len(high_sub) > 0 else 1.0
        sub["rel_prob"] = sub["prob"] / max_p if max_p > 0 else 0.0
        return sub, max_p

    cc1_sub, max_p1 = prep_cc(cc1, "GSE149614")
    cc2_sub, max_p2 = prep_cc(cc2, "GSE156625")
    cc_all = pd.concat([cc1_sub, cc2_sub], ignore_index=True)

    cc_high = cc_all[cc_all["target"] == "CD8_Exh_High"]
    cc_low = cc_all[cc_all["target"] == "CD8_Exh_Low"]
    has_low = not cc_low.empty

    top_paths = (
        cc_high.groupby("pathway_name")["rel_prob"].sum().sort_values(ascending=False).head(8).index.tolist()
    )
    piv_high = (
        cc_high[cc_high["pathway_name"].isin(top_paths)]
        .pivot_table(index="pathway_name", columns=["cohort", "source"], values="rel_prob", fill_value=0.0)
        .reindex(top_paths)
    )
    piv_low = (
        cc_low[cc_low["pathway_name"].isin(top_paths)]
        .pivot_table(index="pathway_name", columns=["cohort", "source"], values="rel_prob", fill_value=0.0)
        .reindex(top_paths)
        if has_low else pd.DataFrame(index=top_paths)
    )

    y_base = np.arange(len(top_paths))
    offsets = [-0.23, 0.0, 0.23]
    cohort_panels = [
        (ax2_a, "GSE149614", max_p1, True),
        (ax2_b, "GSE156625", max_p2, False),
    ]

    for ax_c, cohort_name, max_p, show_ylabels in cohort_panels:
        # Subtle alternating horizontal bands to separate pathways
        for i in range(len(top_paths)):
            if i % 2 == 1:
                ax_c.axhspan(i - 0.45, i + 0.45, color="#F7FAFC", zorder=0)

        for m_idx, mac in enumerate(mac_order):
            col = MAC_PROG_COLORS[mac]
            ys = y_base + offsets[m_idx]
            vh = (
                piv_high[(cohort_name, mac)].to_numpy()
                if (cohort_name, mac) in piv_high.columns
                else np.zeros(len(top_paths))
            )
            vl = (
                piv_low[(cohort_name, mac)].to_numpy()
                if has_low and (cohort_name, mac) in piv_low.columns
                else np.zeros(len(top_paths))
            )

            for yi, h_val, l_val in zip(ys, vh, vl):
                if h_val <= 0 and l_val <= 0:
                    continue
                # Baseline stem from 0 to max signaling so magnitude is immediately readable
                ax_c.plot([0.0, max(h_val, l_val)], [yi, yi], color="#CBD5E0", lw=1.0, zorder=2)
                # Colored dumbbell segment between Exh_High and Exh_Low when both are detected
                if h_val > 0 and l_val > 0:
                    ax_c.plot([l_val, h_val], [yi, yi], color=col, lw=2.2, alpha=0.85, zorder=3)

            m_h = vh > 0
            if np.any(m_h):
                ax_c.scatter(
                    vh[m_h], ys[m_h],
                    marker="o", s=46, color=col, edgecolors="white", lw=0.6, zorder=5,
                )
            m_l = vl > 0
            if has_low and np.any(m_l):
                ax_c.scatter(
                    vl[m_l], ys[m_l],
                    marker="o", s=42, facecolors="white", edgecolors=col, lw=1.4, zorder=4,
                )

        ax_c.set_xlim(-0.04, 1.08)
        ax_c.set_xticks([0.0, 0.25, 0.50, 0.75, 1.0])
        ax_c.set_xlabel("Relative Signaling Prob.", fontsize=9.2)
        ax_c.set_title(f"{cohort_name}\n(max prob = {max_p:.3f})", fontsize=9.8, fontweight="bold", pad=6)
        ax_c.grid(axis="x", ls=":", alpha=0.35)

        if show_ylabels:
            ax_c.set_yticks(y_base)
            ax_c.set_yticklabels(top_paths, fontsize=9.4, fontweight="bold")
            ax_c.invert_yaxis()
        else:
            ax_c.tick_params(axis="y", left=False, labelleft=False)

    # Panel B super-title centered over both cohort facets
    pos_a = ax2_a.get_position()
    pos_b = ax2_b.get_position()
    b_center_x = 0.5 * (pos_a.x0 + pos_b.x1)
    fig.text(
        b_center_x, 0.965,
        "B. CellChat Non-MHC-I Pathways -> CD8 (Not Exhaustion-Specific)",
        ha="center", va="top", fontsize=10.6, fontweight="bold",
    )

    legend_handles = [
        Line2D([0], [0], marker="o", color=MAC_PROG_COLORS[m], lw=2.0, markersize=6.5, label=m)
        for m in mac_order
    ]
    if has_low:
        legend_handles.extend([
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#2D3748",
                   markeredgecolor="white", markersize=6.8, label="Target: CD8_Exh_High (filled)"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
                   markeredgecolor="#2D3748", markeredgewidth=1.4, markersize=6.5, label="Target: CD8_Exh_Low (open)"),
        ])
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(b_center_x, 0.042),
        ncol=3,
        frameon=True,
        facecolor="white",
        edgecolor="#CBD5E0",
        fontsize=8.0,
        columnspacing=1.2,
    )

    fig.text(
        0.5, 0.006,
        "Note: Panel A rank concordance reflects shared NicheNet prior and fixed V3 target genes; "
        "Panel B scales each cohort to its CD8_Exh_High max prob (Exh_High and Exh_Low receive comparable signaling).",
        ha="center", fontsize=8.3, color="#4A5568",
    )

    fig.subplots_adjust(left=0.06, right=0.985, top=0.84, bottom=0.19)
    save_fig(fig, figs_dir, "fig3_cross_cohort_communication")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate publication summary figures (Figures 1-3).")
    parser.add_argument("--tables-dir", default=DEFAULT_TABLES_DIR, help="Input directory containing CSV tables.")
    parser.add_argument("--figs-dir", default=DEFAULT_FIGS_DIR, help="Output directory for PDF/PNG figures.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.figs_dir, exist_ok=True)
    print("Generating Updated Summary Figures")
    print(f"  tables_dir : {args.tables_dir}")
    print(f"  figs_dir   : {args.figs_dir}")
    print(f"  adjustText : {'enabled' if HAS_ADJUST_TEXT else 'not installed (using greedy fallback)'}")
    make_fig1_scatter_grid(args.tables_dir, args.figs_dir)
    make_fig2_forest_plot(args.tables_dir, args.figs_dir)
    make_fig3_communication(args.tables_dir, args.figs_dir)
    print("\nAll summary figures generated successfully.")


if __name__ == "__main__":
    main()