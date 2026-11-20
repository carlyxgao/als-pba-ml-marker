#!/usr/bin/env Rscript
#===========================================================
# Module 1: 样本 Level 差异表达分析与特征筛选
#
# A. 差异表达分析
#    - limma: 经验贝叶斯线性模型
#    - 常规统计检验: 正态性检验→方差齐性检验→
#      Student's t-test / Welch's t-test / Mann-Whitney U test
# B. 特征筛选
#    - LASSO 回归 (L1 正则化, glmnet): 三组比较取并集
#    - 随机森林 (特征重要性, randomForest): 三组比较取并集
#    - 核心标志物 = LASSO 并集 ∩ RF 并集
# C. 核心标志物输出 + 分析报告
#===========================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
})

suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
})

# ========== 1. 路径与队列设置 ==========
# 自动获取项目根目录（脚本位于 scripts/analysis/）
get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    return(normalizePath(file.path(script_dir, "..", "..")))
  }
  getwd()
}
BASE_DIR <- get_project_root()
DATA_DIR <- file.path(BASE_DIR, "data", "csv")

# ========== 命令行参数 ==========
parse_cli_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  res <- list(exclude = c("AD"), tag = NULL)
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--exclude-groups" && i < length(args)) {
      res$exclude <- trimws(strsplit(args[i + 1], ",")[[1]]); i <- i + 2
    } else if (args[i] == "--tag" && i < length(args)) {
      res$tag <- args[i + 1]; i <- i + 2
    } else { i <- i + 1 }
  }
  res
}
cli <- parse_cli_args()
if (!is.null(cli$tag)) cat(sprintf("  [CLI] 输出标签: %s\n", cli$tag))

OUTPUT_DIR <- file.path(BASE_DIR, "results",
                        if (!is.null(cli$tag)) paste0("module1_", cli$tag)
                        else "module1")
DIR_CSV <- file.path(OUTPUT_DIR, "csv")
DIR_FIG <- file.path(OUTPUT_DIR, "fig")
invisible(lapply(c(OUTPUT_DIR, DIR_CSV, DIR_FIG), dir.create, recursive = TRUE, showWarnings = FALSE))

# 队列（中心/批次）配置
# center_col:  样本信息表中代表中心/队列的列名
# discovery_center: 用于发现标志物的中心名称
# 当样本信息中不存在 center_col 列时，所有样本自动归为同一中心 "All"
center_col <- "center"
discovery_center <- "Jinan"

# ========== 2. 读取数据与队列检测 ==========
cat(">>> 读取样本信息...\n")
sample_info <- read.csv(file.path(DATA_DIR, "Sample info merged.csv"),
                        stringsAsFactors = FALSE)
rownames(sample_info) <- sample_info$sample_name

# --- 中心检测 ---
if (!center_col %in% colnames(sample_info)) {
  cat(sprintf(">>> 未检测到中心信息列 '%s'，将所有样本归为同一中心\n", center_col))
  sample_info[[center_col]] <- "All"
}
cohorts_available <- unique(sample_info[[center_col]])
cat("  中心分布:\n"); print(table(sample_info[[center_col]]))

if (!discovery_center %in% cohorts_available) {
  discovery_center <- cohorts_available[1]
}
cat(sprintf("  发现中心: %s\n", discovery_center))

# 差异分析中心：当前仅在发现中心上进行
# 后续可扩展为 c(discovery_center, other_centers) 以进行多中心比较
cohorts_for_de <- discovery_center

cat(">>> 读取 TMM 标准化表达矩阵...\n")
expr_raw <- read.csv(file.path(DATA_DIR, "03-2.Protein expression TMM nomalized.csv"),
                     stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

common_samples <- intersect(colnames(expr_raw), sample_info$sample_name)
cat(sprintf("  共同样本数: %d\n", length(common_samples)))

expr <- expr_raw[, common_samples, drop = FALSE]
sample_info <- sample_info[common_samples, , drop = FALSE]

cat("  分组分布 (全样本):\n"); print(table(sample_info$group))
for (cohort in cohorts_available) {
  idx <- sample_info[[center_col]] == cohort
  cat(sprintf("  %s 分组分布:\n", cohort)); print(table(sample_info$group[idx]))
}

# --- 样本分组处理（内置）---
if (!is.null(cli$exclude)) {
  for (grp in cli$exclude) {
    n_before <- nrow(sample_info)
    sample_info <- sample_info[!sample_info$group %in% grp, , drop = FALSE]
    n_removed <- n_before - nrow(sample_info)
  }
  # 同步更新表达矩阵的样本
  common_samples <- intersect(colnames(expr), sample_info$sample_name)
  expr <- expr[, common_samples, drop = FALSE]
  expr_raw <- expr_raw[, common_samples, drop = FALSE]
  sample_info <- sample_info[common_samples, , drop = FALSE]
  cat("  分组分布:\n"); print(table(sample_info$group))
}

# ========== 3. 过滤低表达蛋白 ==========
cat("\n>>> 过滤低表达蛋白 (至少在 10 个样本中表达 > 0)...\n")
keep <- rowSums(expr > 0) >= 10
expr_filtered <- expr[keep, ]
cat(sprintf("  过滤前: %d 蛋白, 过滤后: %d 蛋白\n", nrow(expr), nrow(expr_filtered)))

expr_log <- log2(expr_filtered + 1)

# ========== 4. 定义比较（动态：ALS vs 其他保留的分组） ==========
available_groups <- setdiff(unique(sample_info$group), "ALS")
comparisons <- lapply(available_groups, function(g) c("ALS", g))
comp_names <- sapply(comparisons, paste, collapse = "_vs_")
cat(sprintf(">>> 比较方案: %s\n", paste(comp_names, collapse = ", ")))

# ========== 5. 差异表达分析 ==========

run_limma_de <- function(expr_mat, sample_info, group1, group2) {
  idx <- which(sample_info$group %in% c(group1, group2))
  sub_expr <- expr_mat[, idx, drop = FALSE]
  sub_info <- sample_info[idx, , drop = FALSE]
  sub_info$group <- factor(sub_info$group, levels = c(group1, group2))

  design <- model.matrix(~ group, data = sub_info)
  fit <- lmFit(sub_expr, design)
  fit <- eBayes(fit, trend = TRUE)
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  colnames(res) <- paste0("limma_", colnames(res))
  res$protein <- rownames(res)
  return(res)
}

select_and_run_test <- function(x, y) {
  n1 <- sum(!is.na(x)); n2 <- sum(!is.na(y))

  shapiro_pval1 <- NA; shapiro_pval2 <- NA
  normal1 <- TRUE; normal2 <- TRUE
  if (n1 >= 3 && n1 <= 5000) { shapiro_pval1 <- shapiro.test(x)$p.value; normal1 <- shapiro_pval1 > 0.05 }
  if (n2 >= 3 && n2 <= 5000) { shapiro_pval2 <- shapiro.test(y)$p.value; normal2 <- shapiro_pval2 > 0.05 }

  var_pval <- NA; var_equal <- NA

  if (normal1 && normal2) {
    vt <- var.test(x, y)
    var_pval <- vt$p.value
    var_equal <- var_pval > 0.05
    if (var_equal) {
      test <- t.test(x, y, var.equal = TRUE)
      method <- "Student_t"
    } else {
      test <- t.test(x, y, var.equal = FALSE)
      method <- "Welch_t"
    }
  } else {
    test <- wilcox.test(x, y, exact = FALSE)
    method <- "Mann_Whitney"
  }

  return(data.frame(
    pvalue    = test$p.value,
    statistic = as.numeric(test$statistic),
    method    = method,
    mean_g1   = mean(x, na.rm = TRUE),
    mean_g2   = mean(y, na.rm = TRUE),
    log2FC    = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
    normal1   = normal1,
    normal2   = normal2,
    shapiro_pval1 = shapiro_pval1,
    shapiro_pval2 = shapiro_pval2,
    var_pval  = var_pval,
    var_equal = var_equal
  ))
}

run_conventional_de <- function(expr_mat, sample_info, group1, group2) {
  idx <- which(sample_info$group %in% c(group1, group2))
  sub_expr <- expr_mat[, idx, drop = FALSE]
  sub_info <- sample_info[idx, , drop = FALSE]

  g1 <- rownames(sub_info)[sub_info$group == group1]
  g2 <- rownames(sub_info)[sub_info$group == group2]

  results_list <- apply(sub_expr, 1, function(row) {
    select_and_run_test(row[g1], row[g2])
  })
  results <- do.call(rbind, results_list)
  results$protein <- rownames(sub_expr)
  rownames(results) <- NULL
  return(results)
}

# ========== 6. 执行差异分析（按队列） ==========
all_results <- list()

for (cohort in cohorts_for_de) {
  cat(sprintf("\n########################################\n"))
  cat(sprintf("  队列: %s\n", cohort))
  cat(sprintf("########################################\n"))

  cohort_idx <- which(sample_info[[center_col]] == cohort)
  cohort_samples <- sample_info$sample_name[cohort_idx]
  cohort_samples <- intersect(cohort_samples, colnames(expr_log))
  cohort_expr <- expr_log[, cohort_samples, drop = FALSE]
  cohort_info <- sample_info[cohort_samples, , drop = FALSE]

  for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    name <- comp_names[i]
    g1 <- comp[1]; g2 <- comp[2]
    cat(sprintf("\n========== 差异分析 [%s]: %s vs %s ==========\n", cohort, g1, g2))

    n1 <- sum(cohort_info$group == g1)
    n2 <- sum(cohort_info$group == g2)
    cat(sprintf("  样本量: %s=%d, %s=%d\n", g1, n1, g2, n2))

    if (n1 < 2 || n2 < 2) {
      cat(sprintf("  跳过: 至少一组样本量不足\n"))
      next
    }

    cat("  [limma] 运行线性模型 + 经验贝叶斯...\n")
    limma_res <- run_limma_de(cohort_expr, cohort_info, g1, g2)

    cat("  [常规检验] 正态性→方差齐性→t/Welch/Mann-Whitney...\n")
    conv_res <- run_conventional_de(cohort_expr, cohort_info, g1, g2)

    merged <- merge(limma_res, conv_res, by = "protein", all = TRUE)

    merged$normal1 <- conv_res$normal1[match(merged$protein, conv_res$protein)]
    merged$normal2 <- conv_res$normal2[match(merged$protein, conv_res$protein)]
    merged$shapiro_pval1 <- conv_res$shapiro_pval1[match(merged$protein, conv_res$protein)]
    merged$shapiro_pval2 <- conv_res$shapiro_pval2[match(merged$protein, conv_res$protein)]
    merged$var_pval <- conv_res$var_pval[match(merged$protein, conv_res$protein)]
    merged$var_equal <- conv_res$var_equal[match(merged$protein, conv_res$protein)]

    merged$conv_padj <- p.adjust(merged$pvalue, method = "BH")
    merged$conv_sig <- merged$conv_padj < 0.05 & abs(merged$log2FC) > log2(1.5)
    merged$limma_adj.P.Val <- p.adjust(merged$limma_P.Value, method = "BH")
    merged$limma_sig <- merged$limma_adj.P.Val < 0.05 & abs(merged$limma_logFC) > log2(1.5)
    merged$consensus <- merged$conv_sig & merged$limma_sig

    merged <- merged[order(merged$pvalue), ]

    cohort_tag <- gsub("[^A-Za-z0-9]", "_", cohort)
    out_file <- file.path(DIR_CSV, sprintf("DE_%s_vs_%s__%s.csv", g1, g2, cohort_tag))
    write.csv(merged, out_file, row.names = FALSE)
    cat(sprintf("  -> %s\n", out_file))
    cat(sprintf("  常规显著: %d | limma显著: %d | 共识: %d\n",
                sum(merged$conv_sig, na.rm = TRUE),
                sum(merged$limma_sig, na.rm = TRUE),
                sum(merged$consensus, na.rm = TRUE)))

    key <- paste0(name, "__", cohort_tag)
    all_results[[key]] <- merged
  }
}

# ========== 7. 差异分析结果汇总 ==========
cat("\n========== 汇总差异分析结果 ==========\n")
all_proteins <- unique(unlist(lapply(all_results, `[[`, "protein")))
master <- data.frame(protein = all_proteins, stringsAsFactors = FALSE)

for (key in names(all_results)) {
  res <- all_results[[key]]
  m <- match(master$protein, res$protein)
  master[[paste0(key, "_conv_pval")]]    <- res$pvalue[m]
  master[[paste0(key, "_conv_padj")]]    <- res$conv_padj[m]
  master[[paste0(key, "_conv_log2FC")]]  <- res$log2FC[m]
  master[[paste0(key, "_conv_method")]]  <- res$method[m]
  master[[paste0(key, "_limma_PValue")]] <- res$limma_P.Value[m]
  master[[paste0(key, "_limma_adjP")]]   <- res$limma_adj.P.Val[m]
  master[[paste0(key, "_limma_logFC")]]  <- res$limma_logFC[m]
  master[[paste0(key, "_consensus")]]    <- res$consensus[m]
}

consensus_cols <- grep("_consensus$", colnames(master), value = TRUE)
master$n_consensus <- rowSums(master[, consensus_cols, drop = FALSE], na.rm = TRUE)
master <- master[order(master$n_consensus, decreasing = TRUE), ]

write.csv(master, file.path(DIR_CSV, "DE_all_comparisons_summary.csv"), row.names = FALSE)
cat(sprintf("  汇总表已保存: %d 个蛋白\n", nrow(master)))

# ========== 9. 特征筛选：LASSO（发现队列） ==========
cohort_tag_disc <- gsub("[^A-Za-z0-9]", "_", discovery_center)
disc_idx <- which(sample_info[[center_col]] == discovery_center)
disc_samples <- sample_info$sample_name[disc_idx]
disc_samples <- intersect(disc_samples, colnames(expr_log))
disc_expr <- expr_log[, disc_samples, drop = FALSE]
disc_info <- sample_info[disc_samples, , drop = FALSE]

cat(sprintf("========== LASSO 特征筛选（发现中心: %s） ==========\n", discovery_center))
cat(sprintf("  发现中心样本数: %d\n", ncol(disc_expr)))

run_lasso <- function(expr_mat, sample_info, group1, group2) {
  idx <- which(sample_info$group %in% c(group1, group2))
  sub_expr <- t(expr_mat[, idx, drop = FALSE])
  sub_info <- sample_info[idx, , drop = FALSE]
  y <- ifelse(sub_info$group == group1, 1, 0)

  set.seed(42)
  cv_fit <- cv.glmnet(sub_expr, y, family = "binomial", alpha = 1,
                      nfolds = min(10, min(table(y))), type.measure = "auc")
  lasso_min <- glmnet(sub_expr, y, family = "binomial", alpha = 1,
                      lambda = cv_fit$lambda.min)
  lasso_1se <- glmnet(sub_expr, y, family = "binomial", alpha = 1,
                      lambda = cv_fit$lambda.1se)

  coef_min <- as.matrix(coef(lasso_min))[-1, , drop = FALSE]
  coef_1se <- as.matrix(coef(lasso_1se))[-1, , drop = FALSE]

  proteins_min <- rownames(coef_min)[coef_min[, 1] != 0]
  proteins_1se <- rownames(coef_1se)[coef_1se[, 1] != 0]

  return(list(
    proteins_min = proteins_min,
    proteins_1se = proteins_1se,
    cv_fit       = cv_fit,
    coef_min     = coef_min[coef_min[, 1] != 0, , drop = FALSE],
    coef_1se     = coef_1se[coef_1se[, 1] != 0, , drop = FALSE]
  ))
}

lasso_results <- list()
for (name in comp_names) {
  parts <- strsplit(name, "_vs_")[[1]]
  g1 <- parts[1]; g2 <- parts[2]
  if (sum(disc_info$group == g1) < 2 || sum(disc_info$group == g2) < 2) {
      cat(sprintf("  LASSO 跳过 %s: 样本量不足\n", name))
      next
    }

  cat(sprintf("  LASSO [%s]: %s\n", discovery_center, name))
  lasso_res <- run_lasso(disc_expr, disc_info, g1, g2)
  lasso_results[[name]] <- lasso_res
  cat(sprintf("    lambda.min: %d 个蛋白\n", length(lasso_res$proteins_min)))
  cat(sprintf("    lambda.1se: %d 个蛋白\n", length(lasso_res$proteins_1se)))
  if (length(lasso_res$proteins_min) > 0) {
    cat("    ", paste(head(lasso_res$proteins_min, 15), collapse = ", "),
        ifelse(length(lasso_res$proteins_min) > 15, "...", ""), "\n")
  }

  write.csv(lasso_res$coef_min,
            file.path(DIR_CSV, sprintf("LASSO_coef_min_%s__%s.csv", name, cohort_tag_disc)))
  write.csv(lasso_res$coef_1se,
            file.path(DIR_CSV, sprintf("LASSO_coef_1se_%s__%s.csv", name, cohort_tag_disc)))
}

if (length(lasso_results) > 0) {
  lasso_union <- Reduce(union, lapply(lasso_results, `[[`, "proteins_min"))
  cat(sprintf("\n  LASSO 三组比较并集 [%s]: %d 个蛋白\n", discovery_center, length(lasso_union)))
} else {
  lasso_union <- character(0)
  cat("  LASSO: 无结果\n")
}

write.csv(data.frame(protein = lasso_union, center = discovery_center),
          file.path(DIR_CSV, paste0("LASSO_union_proteins__", cohort_tag_disc, ".csv")), row.names = FALSE)

# ========== 10. 特征筛选：随机森林（发现中心） ==========
cat(sprintf("\n========== 随机森林特征筛选（发现中心: %s） ==========\n", discovery_center))

run_rf <- function(expr_mat, sample_info, group1, group2, n_top = 30) {
  idx <- which(sample_info$group %in% c(group1, group2))
  sub_expr <- t(expr_mat[, idx, drop = FALSE])
  sub_info <- sample_info[idx, , drop = FALSE]
  y <- factor(ifelse(sub_info$group == group1, group1, group2),
              levels = c(group1, group2))

  set.seed(42)
  rf <- randomForest(sub_expr, y, ntree = 1000, importance = TRUE,
                     sampsize = rep(min(table(y)), 2))

  imp <- as.data.frame(rf$importance)
  imp$protein <- rownames(imp)
  imp <- imp[order(imp$MeanDecreaseAccuracy, decreasing = TRUE), ]
  top_acc <- head(imp$protein, n_top)

  imp <- imp[order(imp$MeanDecreaseGini, decreasing = TRUE), ]
  top_gini <- head(imp$protein, n_top)

  return(list(
    rf          = rf,
    importance  = imp,
    top_acc     = top_acc,
    top_gini    = top_gini
  ))
}

rf_results <- list()
for (name in comp_names) {
  parts <- strsplit(name, "_vs_")[[1]]
  g1 <- parts[1]; g2 <- parts[2]
  if (sum(disc_info$group == g1) < 2 || sum(disc_info$group == g2) < 2) {
    cat(sprintf("  RF 跳过 %s: 样本量不足\n", name))
    next
  }

  cat(sprintf("  Random Forest [%s]: %s\n", discovery_center, name))
  rf_res <- run_rf(disc_expr, disc_info, g1, g2, n_top = 30)
  rf_results[[name]] <- rf_res

  top10 <- head(rf_res$importance[order(rf_res$importance$MeanDecreaseAccuracy, decreasing = TRUE), ], 10)
  cat("    Top10 (MeanDecreaseAccuracy):\n")
  for (k in seq_len(nrow(top10))) {
    cat(sprintf("      %s: %.3f\n", top10$protein[k], top10$MeanDecreaseAccuracy[k]))
  }

  write.csv(rf_res$importance,
            file.path(DIR_CSV, sprintf("RF_importance_%s__%s.csv", name, cohort_tag_disc)), row.names = FALSE)
}

if (length(rf_results) > 0) {
  rf_acc_lists <- lapply(rf_results, `[[`, "top_acc")
  rf_union_acc <- Reduce(union, rf_acc_lists)
  cat(sprintf("\n  RF Top30(Acc) 三组比较并集 [%s]: %d 个蛋白\n", discovery_center, length(rf_union_acc)))
} else {
  rf_union_acc <- character(0)
  cat("  RF: 无结果\n")
}

write.csv(data.frame(protein = rf_union_acc, center = discovery_center),
          file.path(DIR_CSV, paste0("RF_union_proteins__", cohort_tag_disc, ".csv")), row.names = FALSE)

# ========== 11. 核心标志物 ==========
cat(sprintf("\n========== 核心标志物（发现中心: %s） ==========\n", discovery_center))

# 三重交集：共识差异蛋白 (≥1个比较) ∩ LASSO 并集 ∩ RF 并集
de_consensus_proteins <- master$protein[which(master$n_consensus >= 1)]
core_set <- Reduce(intersect, list(de_consensus_proteins, lasso_union, rf_union_acc))

cat(sprintf("  DE共识蛋白 (≥1个比较) [%s]: %d 个蛋白\n", discovery_center, length(de_consensus_proteins)))
cat(sprintf("  LASSO 三组并集 [%s]: %d 个蛋白\n", discovery_center, length(lasso_union)))
cat(sprintf("  RF Top30 三组并集 [%s]: %d 个蛋白\n", discovery_center, length(rf_union_acc)))
cat(sprintf("  三重交集 [%s]: %d 个蛋白\n", discovery_center, length(core_set)))
if (length(core_set) > 0) {
  cat("  ", paste(sort(core_set), collapse = ", "), "\n")
}

write.csv(data.frame(protein = core_set, center = discovery_center),
          file.path(DIR_CSV, paste0("core_markers__", cohort_tag_disc, ".csv")), row.names = FALSE)

# 保存 Venn 图数据（供绘图脚本使用）
venn_data <- list(
  de_consensus = de_consensus_proteins,
  lasso_union  = lasso_union,
  rf_union     = rf_union_acc,
  core_set     = core_set,
  center       = discovery_center
)
saveRDS(venn_data, file.path(OUTPUT_DIR, "venn_data.rds"))
cat("  Venn 数据已保存 -> venn_data.rds\n")

# ========== 12. 保存核心标志物表达数据（供绘图脚本使用）==========
cat(sprintf(">>> 保存核心标志物表达数据（发现中心: %s）...\n", discovery_center))
if (length(core_set) > 1) {
  heatmap_data_list <- list()
  for (name in comp_names) {
    parts <- strsplit(name, "_vs_")[[1]]
    idx <- which(disc_info$group %in% parts)
    sub_expr <- disc_expr[rownames(disc_expr) %in% core_set, idx, drop = FALSE]
    sub_info <- disc_info[idx, , drop = FALSE]
    heatmap_data_list[[name]] <- list(
      expr = sub_expr,
      annotation = data.frame(Group = sub_info$group, row.names = rownames(sub_info))
    )
  }
  # 新增：全样本热图数据（所有中心，core markers 表达 + 分组 + 中心信息）
  all_idx <- which(sample_info$sample_name %in% colnames(expr_log))
  all_expr <- expr_log[rownames(expr_log) %in% core_set, all_idx, drop = FALSE]
  all_anno <- data.frame(
    Group  = sample_info[colnames(all_expr), "group"],
    Center = sample_info[colnames(all_expr), center_col],
    row.names = colnames(all_expr)
  )
  saveRDS(list(
    data         = heatmap_data_list,
    all_expr     = all_expr,
    all_anno     = all_anno,
    center       = discovery_center,
    core_set     = core_set
  ), file.path(OUTPUT_DIR, "heatmap_data.rds"))
  cat("  表达数据已保存 -> heatmap_data.rds\n")
}

# ========== 13. 分析报告 ==========
cat("\n>>> 生成分析报告...\n")

report_lines <- c(
  "=========================================================================",
  "  EV-PBA 模块一：差异表达分析与特征筛选报告",
  "=========================================================================",
  "",
  paste("  生成时间:", Sys.time()),
  paste("  发现中心:", discovery_center),
  "",
  "-----------------------------------------------------------------------",
  "  数据概览",
  "-----------------------------------------------------------------------",
  "",
  paste("  总样本数:", length(common_samples)),
  paste("  总蛋白数:", nrow(expr_filtered)),
  paste("  中心分布:", paste(paste0(names(table(sample_info[[center_col]])), "=", table(sample_info[[center_col]])), collapse = ", ")),
  paste("  分组: ALS =", sum(sample_info$group == "ALS"),
        ", AD =", sum(sample_info$group == "AD"),
        ", HC =", sum(sample_info$group == "HC"),
        ", DC =", sum(sample_info$group == "DC")),
  "",
  "-----------------------------------------------------------------------",
  "  差异表达分析结果 (共识差异蛋白数)",
  "-----------------------------------------------------------------------",
  ""
)

for (key in names(all_results)) {
  res <- all_results[[key]]
  n_conv <- sum(res$conv_sig, na.rm = TRUE)
  n_limma <- sum(res$limma_sig, na.rm = TRUE)
  n_cons <- sum(res$consensus, na.rm = TRUE)
  report_lines <- c(report_lines,
    paste0("  ", key, ":"),
    paste0("    常规检验 (正态→t/Welch/MW): ", n_conv),
    paste0("    limma:                    ", n_limma),
    paste0("    共识差异蛋白:              ", n_cons),
    ""
  )
}

report_lines <- c(report_lines,
  "",
  paste("  共识差异蛋白总数 (≥1个比较):", sum(master$n_consensus >= 1)),
  paste("  共识差异蛋白总数 (≥2个比较):", sum(master$n_consensus >= 2)),
  paste("  共识差异蛋白总数 (全部3个比较):", sum(master$n_consensus >= 3)),
  "",
  "-----------------------------------------------------------------------",
  "  统计检验方法使用统计",
  "-----------------------------------------------------------------------",
  ""
)

for (key in names(all_results)) {
  res <- all_results[[key]]
  method_counts <- table(res$method, useNA = "always")
  report_lines <- c(report_lines,
    paste0("  ", key, ":"),
    paste0("    ", names(method_counts), ": ", method_counts, collapse = "\n"),
    ""
  )
}

report_lines <- c(report_lines,
  "",
  "-----------------------------------------------------------------------",
  "  正态性与方差齐性检验汇总",
  "-----------------------------------------------------------------------",
  ""
)

for (key in names(all_results)) {
  res <- all_results[[key]]
  n_total <- nrow(res)
  n_both_normal <- sum(res$normal1 & res$normal2, na.rm = TRUE)
  n_any_non_normal <- sum(!(res$normal1 & res$normal2), na.rm = TRUE)
  n_var_equal <- sum(res$var_equal, na.rm = TRUE)
  n_var_unequal <- sum(!res$var_equal, na.rm = TRUE)
  n_var_na <- sum(is.na(res$var_equal))

  report_lines <- c(report_lines,
    paste0("  ", key, " (总蛋白: ", n_total, "):"),
    paste0("    两组均正态:            ", n_both_normal, " (", sprintf("%.1f%%", n_both_normal/n_total*100), ")"),
    paste0("    任一非正态:            ", n_any_non_normal, " (", sprintf("%.1f%%", n_any_non_normal/n_total*100), ")"),
    paste0("    方差齐 (F-test p>0.05): ", n_var_equal, " (", sprintf("%.1f%%", n_var_equal/n_total*100), ")"),
    paste0("    方差不齐:              ", n_var_unequal, " (", sprintf("%.1f%%", n_var_unequal/n_total*100), ")"),
    paste0("    未做方差检验 (非正态):  ", n_var_na, " (", sprintf("%.1f%%", n_var_na/n_total*100), ")"),
    ""
  )
}

report_lines <- c(report_lines,
  "-----------------------------------------------------------------------",
  "  LASSO 特征筛选（发现中心）",
  "-----------------------------------------------------------------------",
  ""
)

for (name in names(lasso_results)) {
  n_min <- length(lasso_results[[name]]$proteins_min)
  n_1se <- length(lasso_results[[name]]$proteins_1se)
  report_lines <- c(report_lines,
    paste0("  [", discovery_center, "] ", name, ": lambda.min=", n_min, " 蛋白, lambda.1se=", n_1se, " 蛋白")
  )
}
report_lines <- c(report_lines,
  "",
  paste0("  [", discovery_center, "] LASSO 三组并集: ", length(lasso_union), " 个蛋白"),
  "",
  "-----------------------------------------------------------------------",
  ""
)

report_lines <- c(report_lines,
  "",
  "-----------------------------------------------------------------------",
  "  随机森林特征筛选（发现中心）",
  "-----------------------------------------------------------------------",
  ""
)

for (name in names(rf_results)) {
  imp <- rf_results[[name]]$importance
  top5 <- head(imp[order(imp$MeanDecreaseAccuracy, decreasing = TRUE), ]$protein, 5)
  report_lines <- c(report_lines,
    paste0("  [", discovery_center, "] ", name, " Top5 (MeanDecreaseAccuracy):"),
    paste0("    ", paste(top5, collapse = ", "))
  )
}

report_lines <- c(report_lines,
  paste0("  [", discovery_center, "] RF Top30(Acc) 三组并集: ", length(rf_union_acc), " 个蛋白"),
  "",
  "-----------------------------------------------------------------------",
  "  核心标志物",
  "-----------------------------------------------------------------------",
  "",
  paste("  发现中心:", discovery_center),
  paste("  核心标志物数量:", length(core_set)),
  "",
  "  筛选策略:",
  "    三重交集: DE共识蛋白 (≥1个比较) ∩ LASSO 并集 ∩ RF Top30 并集",
  "",
  "  核心标志物列表:"
)

if (length(core_set) > 0) {
  for (i in seq_along(core_set)) {
    report_lines <- c(report_lines, paste0("    ", i, ". ", core_set[i]))
  }
}

report_lines <- c(report_lines,
  "",
  "-----------------------------------------------------------------------",
  "  输出文件清单",
  "-----------------------------------------------------------------------",
  "",
  paste0("  发现队列标识: ", cohort_tag_disc),
  "  差异分析:",
  paste0("    DE_ALS_vs_HC__", cohort_tag_disc, ".csv, DE_ALS_vs_AD__", cohort_tag_disc, ".csv, DE_ALS_vs_DC__", cohort_tag_disc, ".csv"),
  "    DE_all_comparisons_summary.csv",
  "",
  "  可视化:",
  paste0("    volcano_ALS_vs_HC__", cohort_tag_disc, ".png, volcano_ALS_vs_AD__", cohort_tag_disc, ".png, volcano_ALS_vs_DC__", cohort_tag_disc, ".png"),
  paste0("    heatmap_core_markers_ALS_vs_HC__", cohort_tag_disc, ".png, heatmap_core_markers_ALS_vs_AD__", cohort_tag_disc, ".png, heatmap_core_markers_ALS_vs_DC__", cohort_tag_disc, ".png"),
  "    venn_core_markers.png: 核心标志物筛选 Venn 图",
  "    heatmap_core_markers_all_samples.png: 全样本核心标志物热图",
  "",
  "  LASSO:",
  paste0("    LASSO_coef_min_*__", cohort_tag_disc, ".csv, LASSO_coef_1se_*__", cohort_tag_disc, ".csv, LASSO_union_proteins__", cohort_tag_disc, ".csv"),
  "",
  "  随机森林:",
  paste0("    RF_importance_*__", cohort_tag_disc, ".csv, RF_union_proteins__", cohort_tag_disc, ".csv"),
  "",
  "  核心标志物:",
  paste0("    core_markers__", cohort_tag_disc, ".csv"),
  "",
  "  本报告:",
  "    module1_report.txt",
  "",
  "========================================================================="
)

writeLines(report_lines, file.path(OUTPUT_DIR, "module1_report.txt"))
cat("  报告已保存 ->", file.path(OUTPUT_DIR, "module1_report.txt"), "\n")

# ========== 14. 绘图 ==========
cat("\n>>> 调用绘图脚本生成图形...\n")
source(file.path(BASE_DIR, "scripts", "analysis", "module1_plots.R"))

# ========== 15. 完成 ==========
cat("\n========================================\n")
cat("  模块一分析完成\n")
cat(sprintf("  结果目录: %s\n", OUTPUT_DIR))
cat(sprintf("  发现中心: %s\n", discovery_center))
cat("  输出文件:\n")
  cat(sprintf("    - DE_{ALSvsXX}__%s.csv: 各比较差异分析结果\n", cohort_tag_disc))
  cat("    - DE_all_comparisons_summary.csv: 汇总表\n")
  cat(sprintf("    - volcano_{ALSvsXX}__%s.png: 火山图 (fig/)\n", cohort_tag_disc))
  cat(sprintf("    - LASSO_coef_*__%s.csv: LASSO 系数\n", cohort_tag_disc))
  cat(sprintf("    - LASSO_union_proteins__%s.csv: LASSO 并集特征集\n", cohort_tag_disc))
  cat(sprintf("    - RF_importance_*__%s.csv: RF 特征重要性\n", cohort_tag_disc))
  cat(sprintf("    - RF_union_proteins__%s.csv: RF 并集特征集\n", cohort_tag_disc))
  cat(sprintf("    - core_markers__%s.csv: 最终核心标志物\n", cohort_tag_disc))
  cat(sprintf("    - heatmap_core_markers_*.png: 核心标志物热图 (fig/)\n", cohort_tag_disc))
  cat("    - module1_report.txt: 详细分析报告\n")
cat("========================================\n")