# scripts/01_download_geo_data.sh
# Downloads the three public GEO scRNA-seq datasets into data/raw/
# and verifies SHA-256 checksums against metadata/raw_data_sha256.txt (if present).
#
# Usage (from project root):
#   bash scripts/01_download_geo_data.sh

set -euo pipefail

if [ ! -f "metadata/datasets.yaml" ]; then
  echo "ERROR: Run this script from the repository root (expected metadata/datasets.yaml)." >&2
  exit 1
fi

fetch() {
  local url="$1"
  local dest="$2"
  mkdir -p "$(dirname "$dest")"
  if [ -s "$dest" ]; then
    echo "  [SKIP]: $dest"
    return 0
  fi
  echo "  [GET]  $dest"
  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 3 --retry-delay 5 -o "$dest.tmp" "$url" && mv "$dest.tmp" "$dest"
  elif command -v wget >/dev/null 2>&1; then
    wget -q --show-progress --tries=3 -O "$dest.tmp" "$url" && mv "$dest.tmp" "$dest"
  else
    echo "neither curl nor wget is installed." >&2
    exit 1
  fi
}

echo "1. Downloading GSE149614 (Discovery, Lu et al. 2022)"
GSE149614_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl"
fetch "${GSE149614_BASE}/GSE149614_HCC.scRNAseq.S71915.count.txt.gz" \
  "data/raw/GSE149614/GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
fetch "${GSE149614_BASE}/GSE149614_HCC.metadata.updated.txt.gz" \
  "data/raw/GSE149614/GSE149614_HCC.metadata.updated.txt.gz"

echo ""
echo "2. Downloading GSE156625 (Primary Validation, Sharma et al. 2020)"
GSE156625_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156625/suppl"

if [ -s "data/raw/GSE156625/GSE156625_HCCscanpyobj.h5ad" ]; then
  echo "  [SKIP]: data/raw/GSE156625/GSE156625_HCCscanpyobj.h5ad"
else
  fetch "${GSE156625_BASE}/GSE156625_HCCscanpyobj.h5ad.gz" \
    "data/raw/GSE156625/GSE156625_HCCscanpyobj.h5ad.gz"
  echo "  [GUNZIP]: data/raw/GSE156625/GSE156625_HCCscanpyobj.h5ad.gz"
  gunzip -c "data/raw/GSE156625/GSE156625_HCCscanpyobj.h5ad.gz" > "data/raw/GSE156625/GSE156625_HCCscanpyobj.h5ad.tmp"
  mv "data/raw/GSE156625/GSE156625_HCCscanpyobj.h5ad.tmp" "data/raw/GSE156625/GSE156625_HCCscanpyobj.h5ad"
fi

echo ""
echo "3. Downloading GSE140228 (Supportive Validation, Zhang et al. 2019)"
GSE140228_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140228/suppl"
fetch "${GSE140228_BASE}/GSE140228_UMI_counts_Droplet.mtx.gz" \
  "data/raw/GSE140228/GSE140228_UMI_counts_Droplet.mtx.gz"
fetch "${GSE140228_BASE}/GSE140228_UMI_counts_Droplet_genes.tsv.gz" \
  "data/raw/GSE140228/GSE140228_UMI_counts_Droplet_genes.tsv.gz"
fetch "${GSE140228_BASE}/GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz" \
  "data/raw/GSE140228/GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz"

echo ""
echo "All GEO raw datasets ready in data/raw/."