#!/usr/bin/env Rscript
#===========================================================
# Module 2 绘图脚本（独立版本，不依赖 plots.R / common.R）
# 读取 module2 分析结果 CSV，生成所有图形
# 依赖共享绘图函数库 plot_functions.R
#===========================================================
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))), "plot_functions.R"))

# ========== 1. 路径 ==========
BASE_DIR <- get_project_root()
cli_tag <- parse_cli_tag()
OUTPUT_DIR <- file.path(BASE_DIR, "results",
                        if (!is.null(cli_tag)) paste0("module2_", cli_tag)
                        else "module2")
DIR_CSV <- file.path(OUTPUT_DIR, "csv")
DIR_FIG <- file.path(OUTPUT_DIR, "fig")
invisible(dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE))

# 最佳多分类模型名称（SHAP 标题用）
best_mc_model <- "RF"
if (file.exists(file.path(DIR_CSV, "multiclass_performance_table.csv"))) {
  mt <- read.csv(file.path(DIR_CSV, "multiclass_performance_table.csv"), stringsAsFactors = FALSE)
  if ("Model" %in% names(mt) && "Accuracy" %in% names(mt)) {
    best_mc_model <- mt$Model[which.max(as.numeric(mt$Accuracy))]
  }
}

# ========== 1. ROC / PRC 曲线 ==========
cat(">>> 生成 ROC/PRC 曲线(来自多分类模型ALS概率)...\n")
pred_file <- file.path(DIR_CSV, "multiclass_predictions.csv")
if (file.exists(pred_file)) {
  multi_pred <- read.csv(pred_file, stringsAsFactors = FALSE)
  for (tn in unique(multi_pred$TestSet)) {
    p_roc <- plot_roc(multi_pred, test_set = tn, ci = TRUE)
    ggsave(file.path(DIR_FIG, paste0("ROC_", tn, ".png")), p_roc, width = 6, height = 5.5, dpi = 150)
    cat(sprintf("  ROC_%s.png 已保存\n", tn))

    p_prc <- plot_prc(multi_pred, test_set = tn, ci = TRUE)
    ggsave(file.path(DIR_FIG, paste0("PRC_", tn, ".png")), p_prc, width = 6, height = 5.5, dpi = 150)
    cat(sprintf("  PRC_%s.png 已保存\n", tn))
  }
}

# ========== 2. 混淆矩阵热图 ==========
cat(">>> 生成混淆矩阵...\n")
if (file.exists(pred_file)) {
  multi_pred <- read.csv(pred_file, stringsAsFactors = FALSE)
  for (tn in unique(multi_pred$TestSet)) {
    p_cm <- plot_confusion_matrix(multi_pred, test_set = tn)
    ggsave(file.path(DIR_FIG, paste0("confusion_matrix_", tn, ".png")), p_cm, width = 6.5, height = 5.5, dpi = 150)
    cat(sprintf("  confusion_matrix_%s.png 已保存\n", tn))
  }
}

# ========== 3. SHAP 分析图 ==========
cat(">>> 生成 SHAP 分析图...\n")
shap_imp_file <- file.path(DIR_CSV, "shap_importance_multi.csv")
shap_val_file <- file.path(DIR_CSV, "shap_values_multi_ALS.csv")
shap_feat_file <- file.path(DIR_CSV, "shap_feature_values_multi.csv")

if (file.exists(shap_imp_file) && file.exists(shap_val_file)) {
  # ---- 3a. 特征重要性柱状图 ----
  shap_imp <- read.csv(shap_imp_file, stringsAsFactors = FALSE) %>%
    filter(Class == "ALS") %>% select(Feature, MeanAbsSHAP)
  p_imp <- plot_shap_importance(shap_imp,
                                title = paste("SHAP Feature Importance (Multiclass -", best_mc_model, ")"))
  ggsave(file.path(DIR_FIG, "shap_importance.png"), p_imp,
         width = 9, height = max(5, nrow(shap_imp) * 0.4), dpi = 150)
  cat("  shap_importance.png 已保存\n")

  # ---- 3b. SHAP Beeswarm 图（Top 15 特征）----
  if (file.exists(shap_feat_file)) {
    shap_vals <- read.csv(shap_val_file, stringsAsFactors = FALSE)
    feat_vals <- read.csv(shap_feat_file, stringsAsFactors = FALSE)
    p_beeswarm <- plot_shap_beeswarm(shap_imp, shap_vals, feat_vals, top_n = 15,
                                     title = paste("SHAP Beeswarm Plot -", best_mc_model),
                                     subtitle = "Multiclass model: ALS class SHAP values")
    ggsave(file.path(DIR_FIG, "shap_beeswarm.png"), p_beeswarm,
           width = 12, height = max(6, min(15, nrow(shap_imp)) * 0.45), dpi = 150)
    cat("  shap_beeswarm.png 已保存\n")
  } else {
    cat("  警告: shap_feature_values.csv 不存在，跳过 beeswarm 图\n")
  }
}

# ========== 4. 多分类 SHAP 图 ==========
cat(">>> 生成多分类 SHAP 分析图...\n")
shap_m_imp_file <- file.path(DIR_CSV, "shap_importance_multi.csv")
if (file.exists(shap_m_imp_file)) {
  shap_m_imp <- read.csv(shap_m_imp_file, stringsAsFactors = FALSE)

  # ---- 4a. 多分类 SHAP 特征重要性（每个类别） ----
  p_multi_imp <- plot_shap_importance_multi(shap_m_imp,
                                            title = paste("Multiclass SHAP Feature Importance (", best_mc_model, ")"))
  n_feat <- length(unique(shap_m_imp$Feature))
  ggsave(file.path(DIR_FIG, "shap_importance_multi.png"), p_multi_imp,
         width = 12, height = max(5, n_feat * 0.35), dpi = 150)
  cat("  shap_importance_multi.png 已保存\n")

  # ---- 4b. 每个类别的特征重要性图 ----
  for (cl in unique(shap_m_imp$Class)) {
    p_cl <- plot_shap_importance_class(shap_m_imp, cl)
    ggsave(file.path(DIR_FIG, sprintf("shap_importance_multi_%s.png", cl)), p_cl,
           width = 9, height = max(5, n_feat * 0.35), dpi = 150)
    cat(sprintf("  shap_importance_multi_%s.png 已保存\n", cl))
  }
} else {
  cat("  多分类 SHAP 文件不存在，跳过\n")
}

# ========== 5. Benchmark 柱状图 + ROC ==========
cat(">>> 生成 Benchmark 柱状图...\n")
bench_file <- file.path(DIR_CSV, "benchmark_results.csv")
if (file.exists(bench_file)) {
  bm <- read.csv(bench_file, stringsAsFactors = FALSE)
  p_bench <- plot_benchmark_bar(bm)
  ggsave(file.path(DIR_FIG, "benchmark_auc_barplot.png"), p_bench, width = 10, height = 7, dpi = 150)
  cat("  -> benchmark_auc_barplot.png\n")

  # ---- 5b. ROC 曲线对比（最优模型的 3 个特征集） ----
  bm_data_file <- file.path(OUTPUT_DIR, "benchmark_data.rds")
  if (file.exists(bm_data_file)) {
    bmd <- readRDS(bm_data_file)
    bm$Avg_AUC <- (bm$Test1_AUC + bm$Test2_AUC) / 2
    best_model <- bm$Model[which.max(bm$Avg_AUC)]
    p_roc <- plot_benchmark_roc(bmd, best_model = best_model)
    ggsave(file.path(DIR_FIG, "benchmark_roc.png"), p_roc, width = 12, height = 6, dpi = 150)
    cat("  -> benchmark_roc.png\n")
  }
} else {
  cat("  benchmark_results.csv 不存在，跳过\n")
}

# ========== 6. ALS 亚组分析图 ==========
cat(">>> 生成 ALS 亚组分析图...\n")
subgroup_file <- file.path(DIR_CSV, "subgroup_analysis.csv")
if (file.exists(subgroup_file)) {
  sub_df <- read.csv(subgroup_file, stringsAsFactors = FALSE)
  p_sub <- plot_subgroup(sub_df,
                         title = paste("ALS Subgroup Analysis (", best_mc_model, " multiclass)"))
  ggsave(file.path(DIR_FIG, "subgroup_analysis.png"), p_sub, width = 14, height = 10, dpi = 150)
  cat("  -> subgroup_analysis.png\n")
} else {
  cat("  subgroup_analysis.csv 不存在，跳过\n")
}

# ========== 7. NFL vs 12-Marker 对比（新增；Figure 3F 面板源） ==========
cat(">>> 生成 NFL 对比图...\n")
nfl_auc_file <- file.path(DIR_CSV, "nfl_vs_marker_auc.csv")
multi_pred_file <- file.path(DIR_CSV, "multiclass_predictions.csv")
if (file.exists(nfl_auc_file)) {
  nfl_auc <- read.csv(nfl_auc_file, stringsAsFactors = FALSE)

  # ---- 7a. NFL vs 12-Marker AUC 柱状图 ----
  p_nfl <- plot_nfl_barplot(nfl_auc)
  ggsave(file.path(DIR_FIG, "nfl_vs_marker_barplot.png"), p_nfl, width = 8, height = 6, dpi = 150)
  cat("  -> nfl_vs_marker_barplot.png\n")

  # ---- 7b. NFL ROC vs 12-Marker RF ROC ----
  if (file.exists(multi_pred_file)) {
    multi_pred <- read.csv(multi_pred_file, stringsAsFactors = FALSE)
    bin_pred <- multi_pred %>%
      filter(Model == "RF") %>%
      transmute(Sample, Model, TestSet, TrueLabel, PredClass, Prob_ALS = ALS)
    nfl_vals <- tryCatch(
      assemble_nfl_values(file.path(BASE_DIR, "data", "raw", "nfl.xlsx"), multi_pred),
      error = function(e) { cat("  (NFL 值装配失败，使用 AUC 标注降级)\n"); NULL })

    roc_plots <- list()
    for (tn in unique(bin_pred$TestSet)) {
      nfl_row <- nfl_auc %>% filter(TestSet == tn, Model == "NFL") %>% head(1)
      p_roc <- plot_nfl_roc(bin_pred, nfl_vals, test_set = tn, nfl_auc_row = nfl_row)
      ggsave(file.path(DIR_FIG, paste0("ROC_NFL_vs_Marker_", tn, ".png")),
             p_roc, width = 6, height = 5.5, dpi = 150)
      cat(sprintf("  -> ROC_NFL_vs_Marker_%s.png\n", tn))
      roc_plots[[tn]] <- p_roc
    }
    if (length(roc_plots) >= 2) {
      combined <- ggpubr::ggarrange(plotlist = roc_plots, ncol = 2, labels = names(roc_plots),
                                    font.label = list(size = 12, face = "bold"))
      ggsave(file.path(DIR_FIG, "nfl_vs_marker_roc.png"), combined, width = 12, height = 5.5, dpi = 150)
      cat("  -> nfl_vs_marker_roc.png\n")
    }
  }
} else {
  cat("  nfl_vs_marker_auc.csv 不存在，跳过 NFL 对比图\n")
}

cat("\n>>> Module 2 绘图完成\n")
