#!/usr/bin/env Rscript
#===========================================================
# Module 1 绘图脚本
# 读取 module1 分析结果 CSV + RDS，生成所有图形
#   火山图、Venn 图、核心标志物热图（发现中心 & 全样本）、核心标志物箱线图
# 可独立运行: Rscript module1_plots.R
# 依赖共享绘图函数库 plot_functions.R
#===========================================================
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))), "plot_functions.R"))

# ========== 1. 路径 ==========
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  BASE_DIR <- normalizePath(file.path(script_dir, "..", ".."))
} else {
  BASE_DIR <- getwd()
}
cli_tag <- parse_cli_tag()
OUTPUT_DIR <- file.path(BASE_DIR, "results",
                        if (!is.null(cli_tag)) paste0("module1_", cli_tag)
                        else "module1")
DIR_CSV <- file.path(OUTPUT_DIR, "csv")
DIR_FIG <- file.path(OUTPUT_DIR, "fig")
invisible(dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE))

# ========== 2. 读取 DE 结果 ==========
de_files <- list.files(DIR_CSV, pattern = "^DE_.*\\.csv$", full.names = TRUE)
de_files <- grep("DE_all_comparisons_summary", de_files, value = TRUE, invert = TRUE)
if (length(de_files) == 0) stop("未找到 DE 结果文件")

# ========== 3. 火山图 ==========
cat(">>> 生成火山图...\n")
for (f in de_files) {
  res <- read.csv(f, stringsAsFactors = FALSE)
  key <- gsub("\\.csv$", "", basename(f))
  p <- plot_volcano(res, title = key)
  ggsave(file.path(DIR_FIG, paste0("volcano_", key, ".png")),
         p, width = 8, height = 7, dpi = 150)
  cat(sprintf("  -> volcano_%s.png\n", key))
}

# ========== 4. Venn 图 ==========
cat("\n>>> 生成 Venn 图...\n")
venn_file <- file.path(OUTPUT_DIR, "venn_data.rds")
if (file.exists(venn_file)) {
  vd <- readRDS(venn_file)
  p <- plot_venn(vd)
  ggsave(file.path(DIR_FIG, "venn_core_markers.png"), p,
         width = 9, height = 8, dpi = 300)
  cat("  -> venn_core_markers.png\n")
} else {
  cat("  未找到 venn_data.rds，跳过 Venn 图\n")
}

# ========== 5. 核心标志物热图（发现中心 per-comparison + 全样本综合）==========
cat("\n>>> 生成核心标志物热图...\n")
heatmap_rds <- file.path(OUTPUT_DIR, "heatmap_data.rds")
if (file.exists(heatmap_rds)) {
  hm <- readRDS(heatmap_rds)

  # ---- 5a. 各比较热图（发现中心） ----
  for (name in names(hm$data)) {
    comp_data <- hm$data[[name]]
    mat <- as.matrix(comp_data$expr)
    mat <- mat[apply(mat, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]
    if (nrow(mat) < 2) next

    anno_col <- comp_data$annotation
    available_groups <- intersect(c("ALS","AD","DC","HC"), unique(anno_col$Group))
    anno_colors <- list(Group = ALL_GROUP_COLORS[available_groups])
    p <- plot_core_heatmap(mat, anno_col, anno_colors,
                           title = paste0("Core Markers: ", gsub("__.*", "", name)),
                           scale = "row", fontsize = 14, show_colnames = FALSE)
    ggsave(file.path(DIR_FIG,
                     paste0("heatmap_core_markers_", gsub("__", "_", name), ".png")),
           p, width = 10, height = max(6, nrow(mat) * 0.4), dpi = 150)
  }
  cat("  heatmap_core_markers_*.png (per-comparison) 已保存\n")

  # ---- 5b. 全样本综合热图 ----
  cat("  生成全样本综合热图...\n")
  all_expr <- hm$all_expr
  all_anno <- hm$all_anno

  if (!is.null(all_expr) && nrow(all_expr) >= 2) {
    mat_all <- as.matrix(all_expr)
    mat_all <- mat_all[apply(mat_all, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]

    # 分组和中心配色（动态：只使用实际存在的分组）
    avail_groups <- intersect(ALL_GROUP_LEVELS, unique(all_anno$Group))
    anno_colors <- list(Group = ALL_GROUP_COLORS[avail_groups], Center = CENTER_COLORS)

    # 按 Group 排序列
    all_anno$Group <- factor(all_anno$Group, levels = avail_groups)
    ord <- order(all_anno$Group, all_anno$Center)
    mat_all <- mat_all[, ord, drop = FALSE]
    all_anno <- all_anno[ord, , drop = FALSE]

    p <- plot_core_heatmap(mat_all, all_anno, anno_colors,
                           title = "Core Markers Across All Samples",
                           scale = "row", fontsize = 14, show_colnames = FALSE,
                           clustering_method = "ward.D2")
    ggsave(file.path(DIR_FIG, "heatmap_core_markers_all_samples.png"),
           p, width = 14, height = max(6, nrow(mat_all) * 0.5), dpi = 150)
    cat("  -> heatmap_core_markers_all_samples.png\n")
  } else {
    cat("  全样本热图跳过（core_set < 2 或无表达数据）\n")
  }
} else {
  cat("  未找到 heatmap_data.rds，跳过热图\n")
}

# ========== 6. 核心标志物分组箱线图（新增；Figure 2C 面板源）==========
cat("\n>>> 生成核心标志物箱线图...\n")
core_file <- file.path(DIR_CSV, "core_markers__Jinan.csv")
tmm_file <- file.path(BASE_DIR, "data", "csv", "03-2.Protein expression TMM nomalized.csv")
si_file <- file.path(BASE_DIR, "data", "csv", "Sample info merged.csv")
if (file.exists(core_file) && file.exists(tmm_file) && file.exists(si_file)) {
  core <- read.csv(core_file, stringsAsFactors = FALSE)
  tmm <- read.csv(tmm_file, row.names = 1, check.names = FALSE)
  si <- read.csv(si_file, check.names = FALSE)
  p <- plot_marker_boxplots(tmm, si, core$protein,
                            title = "12 Core Markers: Expression by Group")
  ggsave(file.path(DIR_FIG, "boxplot_core_markers.png"),
         p, width = 14, height = 12, dpi = 150)
  cat("  -> boxplot_core_markers.png\n")
} else {
  cat("  核心标志物/表达数据缺失，跳过箱线图\n")
}

cat("\n>>> Module 1 绘图完成\n")
