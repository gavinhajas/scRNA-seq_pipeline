import gzip, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

counts = "data/raw/GSE149614/GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
meta_f = "data/raw/GSE149614/GSE149614_HCC.metadata.updated.txt.gz"

meta = pd.read_csv(meta_f, sep="\t", index_col="Cell")

with gzip.open(counts, "rt") as fh:
    cells = fh.readline().rstrip("\n").split("\t")
assert len(cells) == len(set(cells)) == 71915
assert set(cells) == set(meta.index), "barcodes differ between counts and metadata"

blocks, genes = [], []
reader = pd.read_csv(counts, sep="\t", header=None, skiprows=1, index_col=0,
                     names=["gene"] + cells, dtype={c: "float32" for c in cells},
                     chunksize=500)
for chunk in reader:
    genes += chunk.index.tolist()
    blocks.append(sp.csr_matrix(chunk.to_numpy()))

X = sp.vstack(blocks, format="csr").T.tocsr()
adata = ad.AnnData(X=X, obs=meta.loc[cells].copy(), var=pd.DataFrame(index=genes))
n_dup = adata.var_names.duplicated().sum()
adata.var_names_make_unique()

vals = adata.X.data[:100000]
print(adata)
print("patients:", adata.obs["patient"].nunique(), "| duplicated gene names:", n_dup)
print("integer-valued counts:", bool(np.allclose(vals, np.round(vals))))
adata.write_h5ad("data/processed/GSE149614_all_raw.h5ad", compression="gzip")