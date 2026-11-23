# 10_export_for_r_communication.py: Export sparse 10x-style matrices & metadata for CellChat/NicheNet in R (<1 GB RAM).
# Usage: python scripts/10_export_for_r_communication.py [GSE149614 GSE156625]
#
# - Outputs: matrix.mtx.gz, genes.tsv, barcodes.tsv, metadata.csv (avoids brittle R/h5ad readers).
# - Tumor macrophages (macrophage_like): Mac_FOLR2_C1Q, Mac_SPP1_MMP, Mac_TREM2_GPNMB, Mac_mixed_or_low (+ cluster_label).
# - Strict CD8 cells: Tumor split at median V3_noprolif_lineage_UCell (CD8_Exh_High / CD8_Exh_Low); Normal as CD8_Normal.

import gzip
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


def ensure_log1p_sparse(adata_obj: ad.AnnData, label: str) -> sp.csc_matrix:
    #Verify expression matrix is non-negative log1p-normalized (1e4 scale) and return genes x cells sparse."""
    X = adata_obj.X
    if sp.issparse(X):
        X = X.tocsr().astype(np.float32)
        min_val = float(X.min()) if X.nnz > 0 else 0.0
        max_val = float(X.max()) if X.nnz > 0 else 0.0
    else:
        X = np.asarray(X, dtype=np.float32)
        min_val = float(np.nanmin(X))
        max_val = float(np.nanmax(X))
        X = sp.csr_matrix(X)

    assert min_val >= -1e-5, f"{label}: negative values detected (min={min_val:.3f}); expected unscaled log1p."

    # Check if raw integer counts were stored instead of log1p(1e4)
    if max_val > 30.0:
        print(f"  [{label}] Max value {max_val:.1f} > 30; applying library-size normalization (1e4) + log1p.")
        row_sums = np.asarray(X.sum(axis=1)).ravel()
        row_sums[row_sums == 0] = 1.0
        inv = sp.diags(1e4 / row_sums)
        X = inv @ X
        X.data = np.log1p(X.data)

    return X.tocsr()


def export_cohort(gse: str, out_root: str = "data/processed/r_export"):
    #Subset, align genes, and export Tumor macrophages + Tumor/Normal strict-gate CD8 cells.
    mye_path = f"data/processed/{gse}_myeloid_scored.h5ad"
    cd8_path = f"data/processed/{gse}_cd8_scored.h5ad"
    assert os.path.exists(mye_path), f"Missing {mye_path}"
    assert os.path.exists(cd8_path), f"Missing {cd8_path}"

    print(f"{'=' * 70}Exporting {gse} for R (CellChat / NicheNet){'=' * 70}")
    mye = ad.read_h5ad(mye_path)
    cd8 = ad.read_h5ad(cd8_path)

    # Subset macrophages: Tumor primary site + macrophage_like
    mye_mask = (mye.obs["site"].astype(str) == "Tumor") & (mye.obs["macrophage_like"].astype(bool))
    mye_sub = mye[mye_mask].copy()

    # Subset CD8: Tumor and Normal sites + strict CD8 gate
    cd8_mask = cd8.obs["site"].astype(str).isin(["Tumor", "Normal"]) & (cd8.obs["gate_strict"].astype(bool))
    cd8_sub = cd8[cd8_mask].copy()

    print(f"  Tumor macrophage-like cells: {mye_sub.n_obs:,}")
    print(f"  Strict-gate CD8 cells (Tumor + Normal): {cd8_sub.n_obs:,} "
          f"(Tumor={int((cd8_sub.obs['site'] == 'Tumor').sum()):,}, "
          f"Normal={int((cd8_sub.obs['site'] == 'Normal').sum()):,})")

    # Align on shared gene symbols
    shared_genes = mye_sub.var_names.intersection(cd8_sub.var_names)
    assert len(shared_genes) > 5000, f"{gse}: unexpectedly few shared genes ({len(shared_genes)})"
    mye_sub = mye_sub[:, shared_genes].copy()
    cd8_sub = cd8_sub[:, shared_genes].copy()
    print(f"  Shared genes exported: {len(shared_genes):,}")

    X_mye = ensure_log1p_sparse(mye_sub, f"{gse}_myeloid")
    X_cd8 = ensure_log1p_sparse(cd8_sub, f"{gse}_cd8")

    # Build unified metadata table
    mye_meta = pd.DataFrame(index=[f"mac__{b}" for b in mye_sub.obs_names])
    mye_meta["orig_barcode"] = mye_sub.obs_names.astype(str).values
    mye_meta["compartment"] = "Macrophage"
    mye_meta["patient"] = mye_sub.obs["patient"].astype(str).values
    mye_meta["site"] = mye_sub.obs["site"].astype(str).values
    mye_meta["cluster_label"] = mye_sub.obs["cluster_label"].astype(str).values
    mye_meta["display_label"] = mye_sub.obs["display_label"].astype(str).values
    mye_meta["cell_group"] = ["Mac_" + x for x in mye_meta["display_label"]]
    for col in ("FOLR2_C1Q__filtered_UCell", "SPP1_MMP__filtered_UCell", "TREM2_GPNMB__filtered_UCell"):
        mye_meta[col] = mye_sub.obs[col].to_numpy(dtype=float) if col in mye_sub.obs else np.nan
    mye_meta["V3_noprolif_lineage_UCell"] = np.nan

    cd8_meta = pd.DataFrame(index=[f"cd8__{b}" for b in cd8_sub.obs_names])
    cd8_meta["orig_barcode"] = cd8_sub.obs_names.astype(str).values
    cd8_meta["compartment"] = "CD8"
    cd8_meta["patient"] = cd8_sub.obs["patient"].astype(str).values
    cd8_meta["site"] = cd8_sub.obs["site"].astype(str).values
    cd8_meta["cluster_label"] = cd8_sub.obs["cluster_label"].astype(str).values
    cd8_meta["display_label"] = "CD8_strict"
    for col in ("FOLR2_C1Q__filtered_UCell", "SPP1_MMP__filtered_UCell", "TREM2_GPNMB__filtered_UCell"):
        cd8_meta[col] = np.nan
    v3_scores = cd8_sub.obs["V3_noprolif_lineage_UCell"].to_numpy(dtype=float)
    cd8_meta["V3_noprolif_lineage_UCell"] = v3_scores

    # Split Tumor CD8 cells at cohort Tumor median V3 exhaustion score
    tum_cd8_mask = cd8_meta["site"] == "Tumor"
    v3_med = float(np.median(cd8_meta.loc[tum_cd8_mask, "V3_noprolif_lineage_UCell"]))
    cd8_groups = np.where(
        cd8_meta["site"] == "Normal",
        "CD8_Normal",
        np.where(cd8_meta["V3_noprolif_lineage_UCell"] >= v3_med, "CD8_Exh_High", "CD8_Exh_Low"),
    )
    cd8_meta["cell_group"] = cd8_groups
    print(f"  Tumor CD8 V3 median threshold: {v3_med:.4f}")

    meta = pd.concat([mye_meta, cd8_meta], axis=0)
    print("  Exported cell_group counts:")
    print(meta["cell_group"].value_counts().to_string())

    # Combine matrices (cells x genes -> transpose to genes x cells for R Matrix::readMM)
    X_all = sp.vstack([X_mye, X_cd8], format="csr").T.tocsc()
    assert X_all.shape == (len(shared_genes), len(meta)), "Matrix shape mismatch"

    out_dir = os.path.join(out_root, gse)
    os.makedirs(out_dir, exist_ok=True)

    with gzip.open(os.path.join(out_dir, "matrix.mtx.gz"), "wb") as fh:
        sio.mmwrite(fh, X_all)
    pd.Series(shared_genes.astype(str)).to_csv(os.path.join(out_dir, "genes.tsv"), index=False, header=False)
    pd.Series(meta.index.astype(str)).to_csv(os.path.join(out_dir, "barcodes.tsv"), index=False, header=False)
    meta.to_csv(os.path.join(out_dir, "metadata.csv"))
    print(f"  Saved R export bundle to {out_dir}/")


def main():
    targets = sys.argv[1:] if len(sys.argv) > 1 else ["GSE149614", "GSE156625"]
    for gse in targets:
        export_cohort(gse)


if __name__ == "__main__":
    main()