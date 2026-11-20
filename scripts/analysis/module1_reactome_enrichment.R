#!/usr/bin/env Rscript
#===========================================================
# Module 1: Reactome 通路富集分析
#
# 对 module1 差异表达（共识显著蛋白）进行 Reactome 通路
# 富集分析，输出气泡图 + CSV 结果表。
# 使用 ReactomePA::enrichPathway（返回 enrichResult，兼容
# clusterProfiler 的 dotplot 等可视化）。
#
# 用法:
#   Rscript module1_reactome_enrichment.R                   # 默认 module1
#   Rscript module1_reactome_enrichment.R --tag main
#===========================================================
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))), "plot_functions.R"))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
})

# ========== 路径与参数 ==========
cli_tag <- parse_cli_tag()
cat(sprintf("  [CLI] tag: %s\n", ifelse(is.null(cli_tag), "default", cli_tag)))

BASE_DIR <- get_project_root()
OUTPUT_DIR <- file.path(BASE_DIR, "results",
                        if (!is.null(cli_tag)) paste0("module1_", cli_tag) else "module1")
DIR_CSV <- file.path(OUTPUT_DIR, "csv")
DIR_ENRICH <- file.path(OUTPUT_DIR, "enrichment")
DIR_FIG <- file.path(DIR_ENRICH, "fig")
DIR_CSV_OUT <- file.path(DIR_ENRICH, "csv")
invisible(lapply(c(DIR_ENRICH, DIR_FIG, DIR_CSV_OUT), dir.create, recursive = TRUE, showWarnings = FALSE))

# ========== DE 文件列表 ==========
de_files <- list.files(DIR_CSV, pattern = "^DE_.*\\.csv$", full.names = TRUE)
de_files <- grep("DE_all_comparisons_summary", de_files, value = TRUE, invert = TRUE)
if (length(de_files) == 0) stop("未找到 DE 结果文件")

# ========== 主循环 ==========
for (f in de_files) {
  comp_label <- gsub("__[^_]+$", "", gsub("\\.csv$", "", basename(f)))
  cat(sprintf("\n========== %s ==========\n", comp_label))

  de <- read.csv(f, stringsAsFactors = FALSE)
  sig_genes <- de %>% filter(consensus == TRUE) %>% pull(protein)
  sig_genes <- unique(gsub("/.*", "", sig_genes))
  cat(sprintf("  共识显著基因: %d 个\n", length(sig_genes)))

  if (length(sig_genes) < 5) {
    cat("  基因数太少，跳过\n")
    next
  }

  eg <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  cat(sprintf("  映射到 Entrez ID: %d 个\n", nrow(eg)))
  if (nrow(eg) < 5) {
    cat("  映射后基因数太少，跳过\n")
    next
  }

  react <- tryCatch(
    enrichPathway(gene = eg$ENTREZID, organism = "human",
                  pvalueCutoff = 0.05, readable = TRUE),
    error = function(e) { cat(sprintf("    Reactome 富集出错: %s\n", conditionMessage(e))); NULL }
  )

  if (!is.null(react) && nrow(react) > 0) {
    write.csv(as.data.frame(react), file.path(DIR_CSV_OUT, sprintf("Reactome_%s.csv", comp_label)), row.names = FALSE)
    cat(sprintf("    %d 条通路显著\n", nrow(react)))

    n_show <- min(20, nrow(react))
    p <- plot_go_dotplot(react, comp_label = comp_label, n_show = n_show,
                         title = paste0("Reactome - ", comp_label))
    sz <- enrichment_panel_size(react, n_show)
    ggsave(file.path(DIR_FIG, sprintf("dotplot_Reactome_%s.png", comp_label)),
           p, width = sz$w, height = sz$h, dpi = 150)
    cat(sprintf("    -> dotplot_Reactome_%s.png\n", comp_label))
  } else {
    cat("    无显著富集\n")
  }
}

cat("\n========================================\n")
cat("  Reactome 富集分析完成\n")
cat(sprintf("  输出: %s\n", DIR_ENRICH))
cat("========================================\n")
