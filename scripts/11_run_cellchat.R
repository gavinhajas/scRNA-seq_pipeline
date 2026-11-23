# 11_run_cellchat.R: CellChat ligand-receptor inference from macrophage subsets to strict-gate CD8+ T cells.
# Usage: Rscript scripts/11_run_cellchat.R <GSE149614 | GSE156625>
#
# - Input: Unscaled log1p-normalised Tumor cells from script 10 (keeps groups with >= 10 cells).
# - Settings: population.size = FALSE in computeCommunProb() to measure per-cell signaling rather than
#   cohort cell-type proportions (which vary between sorted and unsorted datasets).
# - Outputs:
#   * results/tables/cellchat_mac_to_cd8_<GSE>.csv          (significant Macrophage -> CD8 L-R pairs)
#   * results/tables/cellchat_pathways_mac_to_cd8_<GSE>.csv (aggregated signaling pathways)
#   * data/processed/r_export/<GSE>/cellchat_<GSE>.rds      (saved CellChat object)

args <- commandArgs(trailingOnly = TRUE)
gse <- if (length(args) >= 1) args[1] else "GSE149614"

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(CellChat)
})

set.seed(0)

export_dir <- file.path("data", "processed", "r_export", gse)
out_tables <- file.path("results", "tables")
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)

# 1. Load exported sparse log1p matrix and metadata
mtx_path <- file.path(export_dir, "matrix.mtx.gz")
if (!file.exists(mtx_path)) {
  stop(sprintf("Missing %s. Run: python scripts/10_export_for_r_communication.py %s", mtx_path, gse))
}

mat <- as(readMM(gzfile(mtx_path)), "CsparseMatrix")
rownames(mat) <- readLines(file.path(export_dir, "genes.tsv"))
colnames(mat) <- readLines(file.path(export_dir, "barcodes.tsv"))

meta <- read.csv(file.path(export_dir, "metadata.csv"), row.names = 1, stringsAsFactors = FALSE)
stopifnot(identical(colnames(mat), rownames(meta)))

# Restrict to primary Tumor cells and groups with >= 10 cells
tum_idx <- which(meta$site == "Tumor")
mat <- mat[, tum_idx, drop = FALSE]
meta <- meta[tum_idx, , drop = FALSE]

keep_groups <- names(which(table(meta$cell_group) >= 10))
keep_idx <- which(meta$cell_group %in% keep_groups)
mat <- mat[, keep_idx, drop = FALSE]
meta <- meta[keep_idx, , drop = FALSE]
meta$cell_group <- factor(meta$cell_group)

# 2. Initialize CellChat & run inference
cellchat <- createCellChat(object = mat, meta = meta, group.by = "cell_group")
ann <- unique(CellChatDB.human$interaction$annotation)
cellchat@DB <- subsetDB(CellChatDB.human, search = setdiff(ann, "Non-protein Signaling"), key = "annotation")

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 3. Save Macrophage -> CD8 tables & CellChat object
mac_groups <- grep("^Mac_", levels(meta$cell_group), value = TRUE)
cd8_groups <- grep("^CD8_", levels(meta$cell_group), value = TRUE)

df_lr <- subsetCommunication(cellchat, sources.use = mac_groups, targets.use = cd8_groups)
df_path <- subsetCommunication(cellchat, sources.use = mac_groups, targets.use = cd8_groups, slot.name = "netP")

df_lr <- if (!is.null(df_lr) && nrow(df_lr) > 0) arrange(df_lr, desc(prob), pval) else data.frame()
df_path <- if (!is.null(df_path) && nrow(df_path) > 0) arrange(df_path, desc(prob), pval) else data.frame()

write.csv(df_lr, file.path(out_tables, sprintf("cellchat_mac_to_cd8_%s.csv", gse)), row.names = FALSE)
write.csv(df_path, file.path(out_tables, sprintf("cellchat_pathways_mac_to_cd8_%s.csv", gse)), row.names = FALSE)
saveRDS(cellchat, file = file.path(export_dir, sprintf("cellchat_%s.rds", gse)))