import scanpy as sc

SRC = "data/processed/GSE149614_all_raw.h5ad"   # raw integer counts in .X
OUT = "metadata/GSE149614_doublets.csv"

a = sc.read_h5ad(SRC)
print(f"loaded {a.n_obs:,} cells x {a.n_vars:,} genes")
sc.pp.scrublet(a, batch_key="sample", random_state=0)

df = a.obs[["doublet_score", "predicted_doublet"]].copy()
df.index.name = "Cell"
df.to_csv(OUT)
print(f"saved {OUT}\n")

summary = a.obs.groupby("sample", observed=True)["predicted_doublet"].agg(cells="size", doublets="sum", rate="mean")
print("Predicted doublets per sample (0% or >15% suggests threshold failure):")  
print(summary.round(3).to_string())
print("Predicted doublet rate by cell type:")
print(a.obs.groupby("celltype", observed=True)["predicted_doublet"].mean().round(3).to_string())