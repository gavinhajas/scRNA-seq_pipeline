#engine for 03.py and 04.py
import json
import os
import platform
from importlib.metadata import version

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml
from scipy.io import mmread

CFG_PATH = "metadata/datasets.yaml"
OUT_DIR = "data/processed"
LOG_DIR = "results/logs"


def load_cfg(dataset):
    with open(CFG_PATH) as fh:
        cfg = yaml.safe_load(fh)
    if dataset not in cfg:
        raise KeyError(f"{dataset} not in {CFG_PATH}; available: {[k for k in cfg if k != 'defaults']}")
    return cfg.get("defaults") or {}, cfg[dataset]


def load_dataset(c):
    """Returns (AnnData with log-normalised .X, bool: still needs normalisation)."""
    d, kind = c["raw_dir"], c["loader"]
    if kind == "mtx":
        mtx = mmread(os.path.join(d, c["counts"])).T.tocsr().astype(np.float32)
        genes = pd.read_csv(os.path.join(d, c["genes"]), sep="\t")
        cellinfo = pd.read_csv(os.path.join(d, c["cellinfo"]), sep="\t")
        assert mtx.shape == (len(cellinfo), len(genes)), (
            f"matrix {mtx.shape} vs cellinfo {len(cellinfo)} / genes {len(genes)}: "
            "check whether the genes/cellinfo files have a header row")
        gene_col = genes.columns[1] if genes.shape[1] > 1 else genes.columns[0]
        a = ad.AnnData(X=mtx, obs=cellinfo.set_index(c["barcode_col"]),
                       var=pd.DataFrame(index=genes[gene_col].astype(str).values))
        a.var_names_make_unique()
        a.layers["counts"] = a.X.copy()
        return a, True
    a = sc.read_h5ad(os.path.join(d, c["counts"]))
    if kind == "h5ad_raw":  # .X is z-scored; unscaled log-norm values live in .raw
        u = a.raw.to_adata()
        u.obs = a.obs.copy()
        if not sp.issparse(u.X):
            u.X = sp.csr_matrix(u.X)
        return u, False
    if kind == "h5ad_counts":  # raw counts in .X
        a.layers["counts"] = a.X.copy()
        return a, True
    raise ValueError(f"unknown loader '{kind}'")


def _vals(v):
    return [str(x) for x in (v if isinstance(v, list) else [v])]


def subset(a, spec, name):
    """spec keys: include / exclude (column -> values), include_prefix / exclude_prefix (column -> prefix)."""
    mask = np.ones(a.n_obs, bool)
    for col, vals in (spec.get("include") or {}).items():
        s = a.obs[col].astype(str)
        missing = sorted(set(_vals(vals)) - set(s.unique()))
        if missing:
            raise ValueError(f"[{name}] {missing} not found in obs['{col}']. Available: {sorted(s.unique())}")
        mask &= s.isin(_vals(vals)).to_numpy()
    for col, vals in (spec.get("exclude") or {}).items():
        mask &= ~a.obs[col].astype(str).isin(_vals(vals)).to_numpy()
    for col, prefix in (spec.get("include_prefix") or {}).items():
        mask &= a.obs[col].astype(str).str.startswith(str(prefix)).to_numpy()
    for col, prefix in (spec.get("exclude_prefix") or {}).items():
        mask &= ~a.obs[col].astype(str).str.startswith(str(prefix)).to_numpy()
    if not mask.any():
        raise ValueError(f"[{name}] selection is empty")
    return a[mask].copy()


def recluster(a, p, seed):
    sc.pp.highly_variable_genes(a, n_top_genes=p["n_hvg"])
    rx = p.get("hvg_exclude_regex")  # such as heat-shock / immediate-early / ribosomal genes that drive artefact clusters
    if rx:
        drop = np.asarray(a.var_names.str.match(rx))
        n_before = int(a.var["highly_variable"].sum())
        a.var["highly_variable"] = a.var["highly_variable"].to_numpy() & ~drop
        print(f"HVG: removed {n_before - int(a.var['highly_variable'].sum())} of {n_before} genes matching {rx}")
    sc.pp.pca(a, mask_var="highly_variable", random_state=seed)
    a.obsm["X_pca"] = np.ascontiguousarray(a.obsm["X_pca"], dtype=np.float32)
    kw = dict(p.get("neighbors_kwargs") or {})
    if p.get("harmony_key"):  # batch-correct the PCA by patient before building the graph
        import harmonypy  # called directly: scanpy's wrapper breaks with some harmonypy versions
        ho = harmonypy.run_harmony(a.obsm["X_pca"], a.obs, [p["harmony_key"]], random_state=seed)
        z = np.asarray(ho.Z_corr)
        if z.shape[0] != a.n_obs:  # older harmonypy returns PCs x cells
            z = z.T
        assert z.shape == a.obsm["X_pca"].shape, f"unexpected Harmony output shape {z.shape}"
        a.obsm["X_pca_harmony"] = np.ascontiguousarray(z, dtype=np.float32)
        kw["use_rep"] = "X_pca_harmony"
    sc.pp.neighbors(a, n_pcs=p["n_pcs"], random_state=seed, **kw)
    sc.tl.leiden(a, resolution=p["resolution"], key_added=p["key"], flavor="igraph",
                 n_iterations=2, directed=False, random_state=seed)


def check_markers(a, key, checks, on_fail):
    """Each entry says: cluster X must have the highest mean of gene G (catches renumbered clusters)."""
    genes = sorted({g for g in checks.values() if g in a.var_names})
    means = sc.get.obs_df(a, keys=[key] + genes).groupby(key, observed=False).mean()
    problems = []
    for cl, g in checks.items():
        if g not in means.columns:
            problems.append(f"{g} not in data")
            continue
        top = str(means[g].idxmax())
        if top != str(cl):
            problems.append(f"cluster {cl}: expected highest mean {g}, found cluster {top}")
    if problems:
        msg = "marker check failed (cluster numbering may have changed): " + "; ".join(problems)
        if on_fail == "error":
            raise AssertionError(msg)
        print("WARNING:", msg)
    else:
        print(f"marker checks passed for '{key}': {checks}")


def print_marker_table(a, key, genes):
    """Mean log-normalised expression of chosen genes per cluster (for choosing clusters to keep)."""
    present = [g for g in genes if g in a.var_names]
    missing = sorted(set(genes) - set(present))
    df = sc.get.obs_df(a, keys=[key] + present).groupby(key, observed=False).mean().round(2)
    df.insert(0, "n_cells", a.obs[key].value_counts().reindex(df.index).values)
    print(f"\nMarker means per {key} (log-normalised):")
    print(df.to_string())
    if missing:
        print("genes not in data:", missing)


def print_top_degs(a, key, n):
    sc.tl.rank_genes_groups(a, groupby=key, method="wilcoxon", n_genes=n)
    names = pd.DataFrame(a.uns["rank_genes_groups"]["names"])
    print(f"\nTop {n} DEGs per {key} (Wilcoxon):")
    for group, row in names.T.iterrows():
        print(f"{group}: " + ", ".join(map(str, row)))


def _expressed(a, gene):
    if gene not in a.var_names:
        return np.zeros(a.n_obs, bool)
    x = a[:, gene].X
    x = x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel()
    return x > 0


def add_flags(a):
    """Annotation only: dropout makes hard gating on single genes lose true CD8 cells."""
    cd8 = _expressed(a, "CD8A") | _expressed(a, "CD8B")
    a.obs["flag_cd8_expressed"] = cd8
    a.obs["flag_cd4_only"] = _expressed(a, "CD4") & ~cd8
    a.obs["flag_cd3d_negative"] = ~_expressed(a, "CD3D")
    a.obs["flag_mki67"] = _expressed(a, "MKI67")
    print("\nCD8 flags (cells):", {k: int(a.obs[k].sum()) for k in
          ["flag_cd8_expressed", "flag_cd4_only", "flag_cd3d_negative", "flag_mki67"]})


def report(a, key, c):
    for col in (c.get("patient_col"), c.get("tissue_col")):
        if key in a.obs and col in a.obs:
            print(f"\n{key} x {col}:")
            print(pd.crosstab(a.obs[key], a.obs[col]).to_string())


def run(dataset, compartment):
    defaults, c = load_cfg(dataset)
    seed = defaults.get("seed", 0)
    np.random.seed(seed)
    comp = c["compartments"][compartment]

    a, needs_norm = load_dataset(c)
    print(f"[{dataset}] loaded {a.n_obs:,} cells x {a.n_vars:,} genes")
    if c.get("extra_obs"):  # per-cell columns computed by another script; annotation only, no cells removed
        extra = pd.read_csv(c["extra_obs"], index_col=0)
        extra = extra.drop(columns=[x for x in extra.columns if x in a.obs.columns])
        missing = a.obs_names.difference(extra.index)
        assert len(missing) == 0, f"{len(missing)} cells missing from {c['extra_obs']}"
        a.obs = a.obs.join(extra.loc[a.obs_names])
        print(f"merged columns from {c['extra_obs']}: {list(extra.columns)}")
    a = subset(a, c.get("prefilter") or {}, f"{dataset} prefilter")
    print(f"[{dataset}] {a.n_obs:,} cells after prefilter")
    if needs_norm:
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)

    last_key = None
    reported = None
    for i, step in enumerate(comp["steps"], 1):
        if step.get("select"):
            a = subset(a, step["select"], f"{dataset} {compartment} step {i}")
        print(f"\nstep {i}: {a.n_obs:,} cells after selection")
        rc = step.get("recluster")
        if rc:
            recluster(a, rc, seed)
            last_key = rc["key"]
            print(a.obs[last_key].value_counts().sort_index().to_string())
            if step.get("marker_checks"):
                check_markers(a, last_key, step["marker_checks"], step.get("on_check_fail", "error"))
            report(a, last_key, c)
            reported = (last_key, a.n_obs)
            if step.get("print_markers"):
                print_marker_table(a, last_key, step["print_markers"])
            if step.get("print_degs"):
                print_top_degs(a, last_key, int(step["print_degs"]))
            if step.get("keep"):
                a = a[a.obs[last_key].astype(str).isin(_vals(step["keep"]))].copy()
                print(f"kept {_vals(step['keep'])}: {a.n_obs:,} cells")
            if step.get("cluster_labels"):
                m = {str(k): v for k, v in step["cluster_labels"].items()}
                a.obs["cluster_label"] = a.obs[last_key].astype(str).map(m)  # descriptive, reference only
            if step.get("legacy_labels"):
                m = {str(k): v for k, v in step["legacy_labels"].items()}
                a.obs["legacy_label"] = a.obs[last_key].astype(str).map(m)  # reference only

    if comp.get("flags"):
        add_flags(a)
        for f in comp.get("drop_flags") or []:
            a = a[~a.obs[f]].copy()
            print(f"dropped cells with {f}: {a.n_obs:,} remain")
    final_key = comp.get("report_key") or last_key
    if (final_key, a.n_obs) != reported:
        report(a, final_key, c)

    exp = comp.get("expected_n")
    if exp is not None:
        status = "OK" if a.n_obs == exp else "MISMATCH"
        print(f"\nexpected_n check: {status} (got {a.n_obs:,}, expected {exp:,})")

    prov = {"dataset": dataset, "compartment": compartment, "seed": seed,
            "python": platform.python_version(), "scanpy": version("scanpy"),
            "anndata": version("anndata"), "numpy": version("numpy"), "pandas": version("pandas")}
    a.uns["provenance"] = json.dumps(prov)
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    tag = "explore" if comp.get("exploration") else "clean"  # exploration outputs are never mistaken for clean ones
    out = os.path.join(OUT_DIR, f"{dataset}_{compartment}_{tag}.h5ad")
    a.write_h5ad(out, compression="gzip")
    with open(os.path.join(LOG_DIR, f"{dataset}_{compartment}{'' if tag == 'clean' else '_explore'}.json"), "w") as fh:
        json.dump({**prov, "n_cells": int(a.n_obs), "n_genes": int(a.n_vars)}, fh, indent=2)
    print(f"\nSaved {a.n_obs:,} cells to {out}")