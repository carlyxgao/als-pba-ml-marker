#!/usr/bin/env Rscript
#===========================================================
# Module 1: GO/KEGG 通路富集分析
#
# 对差异表达分析的显著蛋白进行 GO (BP/CC/MF) 和 KEGG
# 通路富集分析，输出气泡图和 CSV 结果表。
# 绘图复用共享函数库 plot_functions.R
#
# 用法:
#   Rscript module1_pathway_enrichment.R                           # 默认 module1
#   Rscript module1_pathway_enrichment.R --tag main
#===========================================================
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))), "plot_functions.R"))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# KEGG REST 偶发连接失败：重试而不是静默返回 0 条通路
options(timeout = max(120, getOption("timeout")))
enrich_kegg_retry <- function(gene, attempts = 3) {
  for (i in seq_len(attempts)) {
    kk <- tryCatch(
      enrichKEGG(gene = gene, organism = "hsa", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                 minGSSize = 3, maxGSSize = 500),
      error = function(e) NULL)
    if (!is.null(kk)) return(kk)
    Sys.sleep(2 * i)
  }
  NULL
}

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

# 定义 GO 本体
go_ontologies <- c("BP", "CC", "MF")

# ========== 主循环 ==========
for (f in de_files) {
  # 提取比较名称
  basename_f <- gsub("\\.csv$", "", basename(f))
  # 去掉 "__Jinan" 后缀用于标签
  comp_label <- gsub("__[^_]+$", "", basename_f)
  cat(sprintf("\n========== %s ==========\n", comp_label))

  # 读取 DE 结果
  de <- read.csv(f, stringsAsFactors = FALSE)

  # 筛选共识显著基因 (limma + 常规检验一致)
  sig_genes <- de %>% filter(consensus == TRUE) %>% pull(protein)
  # 处理多基因别名 (如 "ULBP2/RAET1L/RAET1G" → 只取第一个)
  sig_genes <- gsub("/.*", "", sig_genes)
  sig_genes <- unique(sig_genes)
  cat(sprintf("  共识显著基因: %d 个\n", length(sig_genes)))

  if (length(sig_genes) < 5) {
    cat("  基因数太少，跳过富集分析\n")
    next
  }

  # 基因名 → Entrez ID
  eg <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  cat(sprintf("  映射到 Entrez ID: %d 个\n", nrow(eg)))

  if (nrow(eg) < 5) {
    cat("  映射后基因数太少，跳过\n")
    next
  }

  gene_list <- eg$ENTREZID

  # ---- GO 富集 ----
  for (ont in go_ontologies) {
    cat(sprintf("  GO-%s ...\n", ont))
    go_res <- tryCatch(
      enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
               ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05,
               qvalueCutoff = 0.2, minGSSize = 3, maxGSSize = 500,
               readable = TRUE),
      error = function(e) { cat(sprintf("    GO-%s 出错: %s\n", ont, conditionMessage(e))); NULL }
    )

    if (!is.null(go_res) && nrow(go_res) > 0) {
      # 保存 CSV
      write.csv(as.data.frame(go_res), file.path(DIR_CSV_OUT, sprintf("GO_%s_%s.csv", ont, comp_label)), row.names = FALSE)
      cat(sprintf("    %d 条通路显著\n", nrow(go_res)))

      # 气泡图（Top 20）
      n_show <- min(20, nrow(go_res))
      p <- plot_go_dotplot(go_res, ont = ont, comp_label = comp_label, n_show = n_show)
      sz <- enrichment_panel_size(go_res, n_show)
      ggsave(file.path(DIR_FIG, sprintf("dotplot_GO_%s_%s.png", ont, comp_label)),
             p, width = sz$w, height = sz$h, dpi = 150)
      cat(sprintf("    -> dotplot_GO_%s_%s.png\n", ont, comp_label))
    } else {
      cat("    无显著富集\n")
    }
  }

  # ---- KEGG 富集 ----
  cat("  KEGG ...\n")
  kk <- enrich_kegg_retry(gene_list)
  if (is.null(kk)) cat("    KEGG 重试后仍失败\n")

  if (!is.null(kk) && nrow(kk) > 0) {
    # 添加基因符号
    kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(kk), file.path(DIR_CSV_OUT, sprintf("KEGG_%s.csv", comp_label)), row.names = FALSE)
    cat(sprintf("    %d 条通路显著\n", nrow(kk)))

    n_show <- min(20, nrow(kk))
    p <- plot_go_dotplot(kk, ont = NULL, comp_label = comp_label, n_show = n_show)
    sz <- enrichment_panel_size(kk, n_show)
    ggsave(file.path(DIR_FIG, sprintf("dotplot_KEGG_%s.png", comp_label)),
           p, width = sz$w, height = sz$h, dpi = 150)
    cat(sprintf("    -> dotplot_KEGG_%s.png\n", comp_label))
  } else {
    cat("    无显著富集\n")
  }
}

cat("\n========================================\n")
cat("  通路富集分析完成\n")
cat(sprintf("  输出: %s\n", DIR_ENRICH))
cat("========================================\n")
