# 12_run_nichenet.R: NicheNet ligand activity inference from Tumor macrophage programs to strict-gate CD8+ T cells.
# Usage: Rscript scripts/12_run_nichenet.R <GSE149614 | GSE156625>
#
# 1. Loads script 10 log1p matrix/metadata and caches NicheNet v2 human prior models in data/raw/nichenet/.
# 2. Selects candidate ligands expressed (>= 10% cells) in any Tumor macrophage subset (Mac_FOLR2_C1Q,
#    Mac_SPP1_MMP, Mac_TREM2_GPNMB) with cognate receptors expressed (>= 10% cells) in Tumor strict CD8 cells.
# 3. Evaluates ligand activity against two target gene sets:
#    - Primary: frozen 19-gene V3_noprolif_lineage exhaustion signature (metadata/signatures_frozen_GSE149614.yaml).
#    - Sensitivity: top 150 unbiased DEGs in Tumor CD8_Exh_High vs Normal CD8_Normal (Wilcoxon, BH p < 0.05,
#      logFC >= 0.25, pct >= 10%).
#  4. Annotates top ligands (mean expression, % detection, dominant macrophage program) and saves CSVs.


args <- commandArgs(trailingOnly = TRUE)
gse <- if (length(args) >= 1) args[1] else "GSE149614"

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(yaml)
  library(nichenetr)
})

set.seed(0)
options(timeout = 600)

export_dir <- file.path("data", "processed", "r_export", gse)
cache_dir <- file.path("data", "raw", "nichenet")
out_tables <- file.path("results", "tables")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)

# 1. Download / load NicheNet v2 human prior models (cached locally)
fetch_rds <- function(filename, url) {
  dest <- file.path(cache_dir, filename)
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    message("Downloading ", filename, " ...")
    tryCatch(
      download.file(url, destfile = dest, mode = "wb", quiet = TRUE),
      error = function(e) {
        if (file.exists(dest)) unlink(dest)
        stop("Failed to download ", filename, " from ", url, ": ", conditionMessage(e), call. = FALSE)
      }
    )
  }
  readRDS(dest)
}

zenodo_base <- "https://zenodo.org/records/7074291/files"
lr_network <- fetch_rds(
  "lr_network_human_21122021.rds",
  paste0(zenodo_base, "/lr_network_human_21122021.rds?download=1")
) %>%
  distinct(from, to)

ligand_target_matrix <- fetch_rds(
  "ligand_target_matrix_nsga2r_final.rds",
  paste0(zenodo_base, "/ligand_target_matrix_nsga2r_final.rds?download=1")
)

# 2. Load exported sparse log1p matrix, metadata, and frozen V3 signature
mtx_path <- file.path(export_dir, "matrix.mtx.gz")
if (!file.exists(mtx_path)) {
  stop(sprintf("Missing %s. Run: python scripts/10_export_for_r_communication.py %s", mtx_path, gse))
}

mat <- as(readMM(gzfile(mtx_path)), "CsparseMatrix")
rownames(mat) <- readLines(file.path(export_dir, "genes.tsv"))
colnames(mat) <- readLines(file.path(export_dir, "barcodes.tsv"))
meta <- read.csv(file.path(export_dir, "metadata.csv"), row.names = 1, stringsAsFactors = FALSE)
stopifnot(identical(colnames(mat), rownames(meta)))

sig_frozen <- yaml::read_yaml("metadata/signatures_frozen_GSE149614.yaml")
v3_genes <- as.character(sig_frozen$variants$V3_noprolif_lineage)

# 3. Define sender (Tumor Macrophages) and receiver (Tumor strict-gate CD8) expressed genes
pct_expressed <- function(m, idx) {
  if (length(idx) == 0) return(setNames(numeric(nrow(m)), rownames(m)))
  Matrix::rowMeans(m[, idx, drop = FALSE] > 0)
}

mean_expr <- function(m, idx) {
  if (length(idx) == 0) return(setNames(numeric(nrow(m)), rownames(m)))
  Matrix::rowMeans(m[, idx, drop = FALSE])
}

mac_states <- intersect(c("Mac_FOLR2_C1Q", "Mac_SPP1_MMP", "Mac_TREM2_GPNMB"), unique(meta$cell_group))
cd8_tum_idx <- which(meta$compartment == "CD8" & meta$site == "Tumor")
cd8_high_idx <- which(meta$cell_group == "CD8_Exh_High")
cd8_norm_idx <- which(meta$cell_group == "CD8_Normal")

min_pct <- 0.10
sender_pct_mat <- sapply(mac_states, function(st) pct_expressed(mat, which(meta$cell_group == st)))
sender_mean_mat <- sapply(mac_states, function(st) mean_expr(mat, which(meta$cell_group == st)))
expressed_genes_sender <- rownames(sender_pct_mat)[apply(sender_pct_mat, 1, max) >= min_pct]

cd8_tum_pct <- pct_expressed(mat, cd8_tum_idx)
expressed_genes_receiver <- names(cd8_tum_pct)[cd8_tum_pct >= min_pct]
background_expressed_genes <- intersect(expressed_genes_receiver, rownames(ligand_target_matrix))

lr_network_expressed <- lr_network %>%
  filter(from %in% expressed_genes_sender & to %in% expressed_genes_receiver)
potential_ligands <- intersect(unique(lr_network_expressed$from), colnames(ligand_target_matrix))

annotate_ligands <- function(act_df) {
  rank_col <- if ("aupr_corrected" %in% colnames(act_df)) "aupr_corrected" else "pearson"
  act_df <- act_df %>% arrange(desc(.data[[rank_col]]))
  ligs <- act_df$test_ligand

  expr_df <- as.data.frame(sender_mean_mat[ligs, , drop = FALSE])
  colnames(expr_df) <- paste0("mean_", colnames(expr_df))
  pct_df <- as.data.frame(sender_pct_mat[ligs, , drop = FALSE])
  colnames(pct_df) <- paste0("pct_", colnames(pct_df))

  dom_prog <- colnames(sender_mean_mat)[max.col(sender_mean_mat[ligs, , drop = FALSE], ties.method = "first")]

  out <- bind_cols(
    act_df,
    tibble(dominant_mac_program = dom_prog),
    as_tibble(round(expr_df, 4)),
    as_tibble(round(pct_df, 4))
  )
  if (all(c("mean_Mac_SPP1_MMP", "mean_Mac_FOLR2_C1Q") %in% colnames(out))) {
    out <- mutate(out, logFC_SPP1_vs_FOLR2 = round(mean_Mac_SPP1_MMP - mean_Mac_FOLR2_C1Q, 4))
  }
  if (all(c("mean_Mac_TREM2_GPNMB", "mean_Mac_FOLR2_C1Q") %in% colnames(out))) {
    out <- mutate(out, logFC_TREM2_vs_FOLR2 = round(mean_Mac_TREM2_GPNMB - mean_Mac_FOLR2_C1Q, 4))
  }
  out
}

# 4. PRIMARY NicheNet run: Frozen V3_noprolif_lineage exhaustion gene set
geneset_v3 <- intersect(v3_genes, rownames(ligand_target_matrix))
lig_act_v3_ann <- predict_ligand_activities(
  geneset = geneset_v3,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
) %>% annotate_ligands()

write.csv(lig_act_v3_ann, file.path(out_tables, sprintf("nichenet_ligand_activities_V3_%s.csv", gse)), row.names = FALSE)

top20_ligands <- head(lig_act_v3_ann$test_ligand, 20)
active_targets_v3 <- top20_ligands %>%
  lapply(get_weighted_ligand_target_links, geneset = geneset_v3, ligand_target_matrix = ligand_target_matrix, n = 250) %>%
  bind_rows() %>%
  filter(!is.na(weight)) %>%
  arrange(desc(weight))
write.csv(active_targets_v3, file.path(out_tables, sprintf("nichenet_ligand_targets_V3_%s.csv", gse)), row.names = FALSE)

weighted_networks <- fetch_rds("weighted_networks_nsga2r_final.rds", paste0(zenodo_base, "/weighted_networks_nsga2r_final.rds"))
lr_links_v3 <- lr_network_expressed %>%
  filter(from %in% top20_ligands) %>%
  left_join(weighted_networks$lr_sig, by = c("from", "to")) %>%
  arrange(from, desc(weight))
write.csv(lr_links_v3, file.path(out_tables, sprintf("nichenet_ligand_receptors_V3_%s.csv", gse)), row.names = FALSE)
rm(weighted_networks); gc(verbose = FALSE)

# 5. SENSITIVITY NicheNet run: Unbiased Tumor CD8_Exh_High vs Normal CD8_Normal DEGs
  m_hi <- mean_expr(mat, cd8_high_idx)
  m_no <- mean_expr(mat, cd8_norm_idx)
  p_hi <- pct_expressed(mat, cd8_high_idx)
  logfc <- m_hi - m_no

  cand_genes <- intersect(names(logfc)[logfc >= 0.25 & p_hi >= 0.10], rownames(ligand_target_matrix))
  if (length(cand_genes) >= 10) {
    hi_sub <- as.matrix(mat[cand_genes, cd8_high_idx, drop = FALSE])
    no_sub <- as.matrix(mat[cand_genes, cd8_norm_idx, drop = FALSE])
    pvals <- vapply(seq_along(cand_genes), function(i) {
      suppressWarnings(wilcox.test(hi_sub[i, ], no_sub[i, ])$p.value)
    }, numeric(1))

    deg_df <- tibble(
      gene = cand_genes,
      mean_Exh_High = round(m_hi[cand_genes], 4),
      mean_Normal = round(m_no[cand_genes], 4),
      pct_Exh_High = round(p_hi[cand_genes], 4),
      logFC_diff = round(logfc[cand_genes], 4),
      pval = pvals,
      padj = p.adjust(pvals, method = "BH")
    ) %>%
      filter(padj < 0.05) %>%
      arrange(desc(logFC_diff))

    write.csv(deg_df, file.path(out_tables, sprintf("nichenet_cd8_tumor_vs_normal_DEGs_%s.csv", gse)), row.names = FALSE)
    geneset_deg <- head(deg_df$gene, 150)

    if (length(geneset_deg) >= 10) {
      lig_act_deg_ann <- predict_ligand_activities(
        geneset = geneset_deg,
        background_expressed_genes = background_expressed_genes,
        ligand_target_matrix = ligand_target_matrix,
        potential_ligands = potential_ligands
      ) %>% annotate_ligands()
      write.csv(lig_act_deg_ann, file.path(out_tables, sprintf("nichenet_ligand_activities_DEG_%s.csv", gse)), row.names = FALSE)
    }
  }
