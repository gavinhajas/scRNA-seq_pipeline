files <- list.files(
  pattern = "^cellchat_mac_to_cd8_(GSE149614|GSE156625)\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

for (f in files) {
  cat("===", f, "===")
  df <- read.csv(f)
  sub <- subset(df,
    (source == "Mac_FOLR2_C1Q" & pathway_name == "CXCL") |
    (source == "Mac_TREM2_GPNMB" & pathway_name == "GALECTIN")
  )
  sub <- sub[order(sub$pathway_name, -sub$prob), ]
  print(sub[, c("source", "target", "pathway_name", "interaction_name", "prob", "pval")])
}