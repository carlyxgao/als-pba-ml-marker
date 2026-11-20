#!/usr/bin/env Rscript
#===========================================================
# Module 1: 敏感性分析 — 多阈值扫描 + 稳定通路筛选
#
# 通过变化 log2FC 和 p.adjust 阈值，生成多组 DE 基因集，
# 分别进行 GO/KEGG 富集，筛选跨阈值稳定富集的通路。
# 绘图复用共享函数库 plot_functions.R
#
# 用法:
#   Rscript module1_sensitivity_analysis.R                           # 默认 module1
#   Rscript module1_sensitivity_analysis.R --tag main
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
DIR_SENS <- file.path(OUTPUT_DIR, "sensitivity")
DIR_FIG <- file.path(DIR_SENS, "fig")
DIR_CSV_OUT <- file.path(DIR_SENS, "csv")
invisible(lapply(c(DIR_SENS, DIR_FIG, DIR_CSV_OUT), dir.create, recursive = TRUE, showWarnings = FALSE))

# ========== 阈值网格 ==========
# 注意: 这些是 **log2FC** 阈值（筛选条件为 abs(limma_logFC) > fc），不是倍数；
#       0.585 = log2(1.5)。键名沿用 "FC=..." 以保持既有 CSV 列名不变。
fc_thresholds <- c(0.3, 0.5, 0.585, 1.0, 1.5, 2.0)
p_thresholds <- c(0.01, 0.05, 0.1, 0.2)

cat(sprintf("阈值网格: %d log2FC × %d p.adj = %d 组合\n",
            length(fc_thresholds), length(p_thresholds),
            length(fc_thresholds) * length(p_thresholds)))

# ========== DE 文件 ==========
de_files <- list.files(DIR_CSV, pattern = "^DE_.*\\.csv$", full.names = TRUE)
de_files <- grep("DE_all_comparisons_summary", de_files, value = TRUE, invert = TRUE)
if (length(de_files) == 0) stop("未找到 DE 结果文件")

# 结果收集
all_results <- list()
stable_list <- list()

for (f in de_files) {
  comp_label <- gsub("__[^_]+$", "", gsub("\\.csv$", "", basename(f)))
  cat(sprintf("\n========== %s ==========\n", comp_label))

  de <- read.csv(f, stringsAsFactors = FALSE)

  # 对每个阈值组合筛选基因 → 富集
  grid_results <- list()

  for (fc in fc_thresholds) {
    for (pv in p_thresholds) {
      # 筛选 limma 显著基因
      genes <- de %>% filter(limma_adj.P.Val < pv & abs(limma_logFC) > fc) %>% pull(protein)
      genes <- unique(gsub("/.*", "", genes))  # 处理多基因别名

      key <- sprintf("FC=%.2f_p=%.2f", fc, pv)
      cat(sprintf("  %s: %d 个差异基因\n", key, length(genes)))

      if (length(genes) < 5) {
        grid_results[[key]] <- list(n_genes = length(genes), GO_BP = 0, KEGG = 0,
                                    go_terms = "", kegg_terms = "")
        next
      }

      eg <- tryCatch(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db, drop = TRUE),
                     error = function(e) NULL)
      if (is.null(eg) || nrow(eg) < 5) {
        grid_results[[key]] <- list(n_genes = length(genes), GO_BP = 0, KEGG = 0,
                                    go_terms = "", kegg_terms = "")
        next
      }

      # GO BP 富集
      go <- tryCatch(enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db,
                              ont = "BP", pAdjustMethod = "BH",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                              minGSSize = 3, maxGSSize = 500, readable = TRUE),
                     error = function(e) NULL)
      n_go <- if (!is.null(go)) nrow(go) else 0
      go_terms <- if (n_go > 0) paste(head(go$Description, 10), collapse = "; ") else ""

      # KEGG 富集
      kk <- enrich_kegg_retry(eg$ENTREZID)
      n_kk <- if (!is.null(kk)) nrow(kk) else 0
      kegg_terms <- if (n_kk > 0) paste(head(kk$Description, 10), collapse = "; ") else ""

      grid_results[[key]] <- list(n_genes = length(genes), GO_BP = n_go, KEGG = n_kk,
                                  go_terms = go_terms, kegg_terms = kegg_terms)

      # 记录每个通路的出现（用于稳定性分析）
      if (n_go > 0) {
        for (desc in go$Description) {
          stable_list[[length(stable_list) + 1]] <- data.frame(
            Comparison = comp_label, Threshold = key, Type = "GO_BP",
            Pathway = desc, p.adjust = go$p.adjust[go$Description == desc][1],
            stringsAsFactors = FALSE)
        }
      }
      if (n_kk > 0) {
        for (desc in kk$Description) {
          stable_list[[length(stable_list) + 1]] <- data.frame(
            Comparison = comp_label, Threshold = key, Type = "KEGG",
            Pathway = desc, p.adjust = kk$p.adjust[kk$Description == desc][1],
            stringsAsFactors = FALSE)
        }
      }
    }
  }

  # 输出网格结果
  grid_df <- do.call(rbind, lapply(names(grid_results), function(k) {
    data.frame(Threshold = k, Comparison = comp_label,
               n_DE_genes = grid_results[[k]]$n_genes,
               n_GO_BP = grid_results[[k]]$GO_BP,
               n_KEGG = grid_results[[k]]$KEGG,
               GO_terms = grid_results[[k]]$go_terms,
               KEGG_terms = grid_results[[k]]$kegg_terms,
               stringsAsFactors = FALSE)
  }))
  write.csv(grid_df, file.path(DIR_CSV_OUT, sprintf("threshold_grid_%s.csv", comp_label)), row.names = FALSE)
  cat(sprintf("  网格结果已保存: threshold_grid_%s.csv\n", comp_label))

  # 阈值组合热图
  if (nrow(grid_df) > 0) {
    p_go <- plot_sensitivity_heatmap(grid_df, fc_thresholds, p_thresholds, "n_GO_BP",
                                     title = paste("GO BP -", comp_label))
    ggsave(file.path(DIR_FIG, sprintf("sensitivity_heatmap_GO_BP_%s.png", comp_label)),
           p_go, width = 8, height = 6, dpi = 150)
    cat(sprintf("  -> sensitivity_heatmap_GO_BP_%s.png\n", comp_label))

    p_kegg <- plot_sensitivity_heatmap(grid_df, fc_thresholds, p_thresholds, "n_KEGG",
                                       title = paste("KEGG -", comp_label))
    ggsave(file.path(DIR_FIG, sprintf("sensitivity_heatmap_KEGG_%s.png", comp_label)),
           p_kegg, width = 8, height = 6, dpi = 150)
    cat(sprintf("  -> sensitivity_heatmap_KEGG_%s.png\n", comp_label))
  }

  all_results[[comp_label]] <- grid_df
}

# ========== 跨阈值稳定性分析 ==========
if (length(stable_list) > 0) {
  cat("\n========== 稳定性分析 ==========\n")
  stable_df <- bind_rows(stable_list)

  # 每个通路在多少阈值组合中显著
  pathway_summary <- stable_df %>%
    group_by(Comparison, Type, Pathway) %>%
    summarise(
      n_significant = n(),
      total_thresholds = length(fc_thresholds) * length(p_thresholds),
      stability_ratio = n_significant / total_thresholds,
      min_padj = min(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(stability_ratio))

  write.csv(pathway_summary, file.path(DIR_CSV_OUT, "sensitivity_summary.csv"), row.names = FALSE)
  cat(sprintf("  稳定性汇总已保存: %d 条通路记录\n", nrow(pathway_summary)))

  # 稳定富集的通路（≥30% 阈值组合中显著）
  stable_pathways <- pathway_summary %>% filter(stability_ratio >= 0.3)
  write.csv(stable_pathways, file.path(DIR_CSV_OUT, "stable_pathways.csv"), row.names = FALSE)
  cat(sprintf("  稳定富集的通路 (≥%.0f%%): %d 条\n", 30, nrow(stable_pathways)))

  # 稳定性柱状图（每比较仅绘制 Top30，高度按实际绘制的条形数计算，避免超大图）
  if (nrow(stable_pathways) > 0) {
    n_comp <- length(unique(stable_pathways$Comparison))
    p <- plot_stability_barplot(stable_pathways)
    plot_h <- max(8, 30 * 0.3 * ceiling(n_comp / 2) + 2)
    ggsave(file.path(DIR_FIG, "stability_barplot.png"), p,
           width = 14, height = plot_h, dpi = 150, limitsize = FALSE)
    cat("  -> stability_barplot.png\n")
  }
}

cat("\n========================================\n")
cat("  敏感性分析完成\n")
cat(sprintf("  输出: %s\n", DIR_SENS))
cat("========================================\n")
