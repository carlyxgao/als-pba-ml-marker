#!/usr/bin/env Rscript

# ============================================================
# 模块二：机器学习建模与验证
# - 二分类: ALS vs Others
# - 多分类: ALS / AD / HC / DC
# - 模型: 逻辑回归(LR) / 随机森林(RF) / SVM(RBF) / XGBoost
# ============================================================

# ========== 1. 加载包 ==========
suppressPackageStartupMessages({
  library(caret)
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ROCR)
  library(ggthemes)
  library(showtext)
  library(fastshap)
  library(xgboost)
})

# ========== CLI 参数 ==========
parse_cli_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  res <- list(exclude = c("AD"), tag = NULL, markers = NULL)
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--exclude-groups" && i < length(args)) {
      res$exclude <- trimws(strsplit(args[i + 1], ",")[[1]]); i <- i + 2
    } else if (args[i] == "--tag" && i < length(args)) {
      res$tag <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--markers" && i < length(args)) {
      res$markers <- args[i + 1]; i <- i + 2
    } else { i <- i + 1 }
  }
  res
}
cli <- parse_cli_args()
if (!is.null(cli$tag)) cat(sprintf("  [CLI] 输出标签: %s\n", cli$tag))
if (!is.null(cli$markers)) cat(sprintf("  [CLI] 标志物文件: %s\n", cli$markers))

# ========== 2. 路径与队列设置 ==========
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

showtext_auto(enable = TRUE)
theme_set(theme_few() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  text = element_text(size = 14)))

DATA_DIR <- file.path(BASE_DIR, "data", "csv")

# 标志物文件路径（优先 CLI 参数，否则默认从 module1 读取）
if (!is.null(cli$markers)) {
  marker_csv <- cli$markers
  if (!file.exists(marker_csv)) stop("标志物文件不存在: ", marker_csv)
} else {
  # 优先读取与 --tag 对应的模块一目录（如 results/module1_main/csv），
  # 避免在 module2_main 中静默复用  results/module1（全量标签）的旧 marker 集。
  m1_dir <- if (!is.null(cli$tag)) {
    file.path(BASE_DIR, "results", paste0("module1_", cli$tag), "csv")
  } else {
    file.path(BASE_DIR, "results", "module1", "csv")
  }
  marker_csv <- list.files(m1_dir, pattern = "core_markers.*\\.csv", full.names = TRUE)[1]
  if (is.na(marker_csv)) {
    marker_csv <- list.files(file.path(BASE_DIR, "results", "module1", "csv"),
                             pattern = "core_markers.*\\.csv", full.names = TRUE)[1]
  }
  if (is.na(marker_csv)) stop("未找到核心标志物文件")
}

OUTPUT_DIR <- file.path(BASE_DIR, "results",
                        if (!is.null(cli$tag)) paste0("module2_", cli$tag)
                        else "module2")
DIR_CSV <- file.path(OUTPUT_DIR, "csv")
DIR_FIG <- file.path(OUTPUT_DIR, "fig")
DIR_MODEL <- file.path(OUTPUT_DIR, "model")
invisible(lapply(c(OUTPUT_DIR, DIR_CSV, DIR_FIG, DIR_MODEL), dir.create, recursive = TRUE, showWarnings = FALSE))

# 队列（中心/批次）配置（与模块一保持一致）
center_col <- "center"
discovery_center <- "Jinan"

# ========== 3. 读取数据 ==========
cat(">>> 读取样本信息...\n")
sample_info <- read.csv(file.path(DATA_DIR, "Sample info merged.csv"),
                        stringsAsFactors = FALSE)
rownames(sample_info) <- sample_info$sample_name

# 中心检测逻辑
has_cohort <- center_col %in% colnames(sample_info)
if (has_cohort) {
  centers <- unique(sample_info[[center_col]])
  if (!discovery_center %in% centers) {
    discovery_center <- centers[1]
  }
  cat(sprintf(">>> 检测到中心信息列 '%s'\n", center_col))
  cat(sprintf("  中心分布:\n"))
  print(table(sample_info[[center_col]]))
  cat(sprintf("  发现中心: %s\n", discovery_center))
} else {
  cat(sprintf(">>> 未检测到中心信息列 '%s'，将所有样本归为同一中心\n", center_col))
  sample_info[[center_col]] <- "All"
  discovery_center <- "All"
  cat("  中心分布:\n")
  print(table(sample_info[[center_col]]))
  cat(sprintf("  发现中心: %s\n", discovery_center))
}

cat(sprintf(">>> 读取核心标志物列表: %s\n", marker_csv))
core_markers <- read.csv(marker_csv, stringsAsFactors = FALSE)
core_proteins <- core_markers$protein
cat(sprintf("  核心标志物数: %d\n", length(core_proteins)))
cat(sprintf("  标志物列表: %s\n", paste(core_proteins, collapse = ", ")))

cat(">>> 读取 TMM 标准化表达矩阵...\n")
expr_raw <- read.csv(file.path(DATA_DIR, "03-2.Protein expression TMM nomalized.csv"),
                     row.names = 1, check.names = FALSE)

# 对齐样本
common_samples <- intersect(colnames(expr_raw), sample_info$sample_name)
sample_info <- sample_info[common_samples, , drop = FALSE]
expr_raw <- expr_raw[, common_samples, drop = FALSE]

# 提取核心标志物表达 + log2 转换
expr_core <- expr_raw[intersect(core_proteins, rownames(expr_raw)), , drop = FALSE]
expr_core <- as.matrix(expr_core)
expr_log <- log2(expr_core + 1)

cat(sprintf("  共同样本数: %d\n", ncol(expr_log)))
cat(sprintf("  核心蛋白数: %d\n", nrow(expr_log)))

# 构建完整数据框（行=样本, 列=蛋白+组别+临床）
df_all <- as.data.frame(t(expr_log))
df_all$group <- sample_info[rownames(df_all), "group"]
df_all[[center_col]] <- sample_info[rownames(df_all), center_col]
df_all$gender <- sample_info[rownames(df_all), "gender"]
df_all$age <- as.numeric(sample_info[rownames(df_all), "age"])
# ALS 特有临床变量（用于亚组分析）
df_all$onset_site <- sample_info[rownames(df_all), "onset_site"]
df_all$kings_stage <- sample_info[rownames(df_all), "kings_stage"]
df_all$disease_duration <- as.numeric(sample_info[rownames(df_all), "disease_duration"])
df_all$functional_score <- as.numeric(sample_info[rownames(df_all), "functional_score"])

# 样本分组处理（内置）
if (!is.null(cli$exclude)) {
  for (grp in cli$exclude) {
    n_before <- nrow(df_all)
    df_all <- df_all[!df_all$group %in% grp, , drop = FALSE]
  }
}

cat(sprintf("  分组分布:\n"))
print(table(df_all$group))

# 定义 metadata 特征
metadata_cols <- c("age", "gender")
# 确保 metadata 无缺失（用中位数填补年龄，众数填补性别）
df_all$age[is.na(df_all$age)] <- median(df_all$age, na.rm = TRUE)
df_all$gender[is.na(df_all$gender)] <- names(sort(table(df_all$gender), decreasing = TRUE))[1]

# 三个特征集
feature_sets <- list(
  Markers = core_proteins,
  Metadata = metadata_cols,
  MarkersPlusMetadata = c(core_proteins, metadata_cols)
)

# ========== 4. 数据集划分 ==========
cat("\n>>> 数据集划分...\n")

if (has_cohort) {
  # 中心感知划分
  disc_samples <- rownames(df_all[df_all[[center_col]] == discovery_center, ])
  other_centers <- setdiff(unique(df_all[[center_col]]), discovery_center)
  if (length(other_centers) >= 2) {
    test1_center <- other_centers[1]
    test2_center <- other_centers[2]
  } else if (length(other_centers) == 1) {
    test1_center <- other_centers[1]
    # 从发现中心中分一部分作为 test2
    set.seed(42)
    holdout_idx <- sample(1:length(disc_samples), size = floor(length(disc_samples) * 0.2))
    test2_samples <- disc_samples[holdout_idx]
    disc_samples <- disc_samples[-holdout_idx]
    test2_center <- "Holdout"
  } else {
    stop("至少需要 1 个发现中心 + 1 个测试中心")
  }
  train_samples <- disc_samples
  test1_samples <- rownames(df_all[df_all[[center_col]] == test1_center, ])
  if (exists("test2_samples")) {
    # already defined
  } else {
    test2_samples <- rownames(df_all[df_all[[center_col]] == test2_center, ])
  }
} else {
  # 无队列信息：随机划分 60/20/20
  set.seed(42)
  n <- nrow(df_all)
  idx <- sample(1:n, n)
  n_train <- floor(n * 0.6)
  n_test1 <- floor(n * 0.2)
  train_samples <- rownames(df_all)[idx[1:n_train]]
  test1_samples <- rownames(df_all)[idx[(n_train + 1):(n_train + n_test1)]]
  test2_samples <- rownames(df_all)[idx[(n_train + n_test1 + 1):n]]
}

df_train <- df_all[train_samples, ]
df_test1 <- df_all[test1_samples, ]
df_test2 <- df_all[test2_samples, ]

cat(sprintf("  训练集: %d 样本\n", nrow(df_train)))
cat(sprintf("  测试集 (%s): %d 样本\n", test1_center, nrow(df_test1)))
cat(sprintf("  测试集 (%s): %d 样本\n", test2_center, nrow(df_test2)))

# 测试集名称标签（用于输出文件）
TEST1_NAME <- test1_center
TEST2_NAME <- test2_center

# ========== 5. 二分类 Benchmark: 3特征集 × 4模型 ==========
cat("\n========================================================\n")
cat("  二分类 Benchmark: 3特征集 × 4模型\n")
cat("========================================================\n")

# 构建特征矩阵并确保数据完整
prepare_features <- function(df, features) {
  avail <- intersect(features, colnames(df))
  # gender 转为 0/1
  if ("gender" %in% avail) {
    df$gender <- as.numeric(factor(df$gender, levels = c("Male", "Female"))) - 1
  }
  as.matrix(df[, avail, drop = FALSE])
}

# 标签
df_train_bin <- df_train; df_train_bin$label <- factor(ifelse(df_train_bin$group=="ALS","ALS","Other"), levels=c("ALS","Other"))
df_test1_bin <- df_test1; df_test1_bin$label <- factor(ifelse(df_test1_bin$group=="ALS","ALS","Other"), levels=c("ALS","Other"))
df_test2_bin <- df_test2; df_test2_bin$label <- factor(ifelse(df_test2_bin$group=="ALS","ALS","Other"), levels=c("ALS","Other"))

train_ctrl_bin <- trainControl(method="cv", number=5, classProbs=TRUE,
                                summaryFunction=twoClassSummary, savePredictions="final",
                                verboseIter=FALSE)

# 模型工厂
run_model <- function(method_name, train_x, train_y, test1_x, test2_x) {
  set.seed(42)
  if (method_name == "LR") {
    m <- train(x=train_x, y=train_y, method="glm", family="binomial", trControl=train_ctrl_bin, metric="ROC")
  } else if (method_name == "RF") {
    p <- ncol(train_x)
    grid <- expand.grid(mtry=unique(floor(seq(1, p, length.out=min(5, p)))))
    m <- train(x=train_x, y=train_y, method="rf", trControl=train_ctrl_bin, tuneGrid=grid, metric="ROC", ntree=500)
  } else if (method_name == "SVM") {
    grid <- expand.grid(C=c(0.1,1,10,100), sigma=c(0.001,0.01,0.1,1))
    m <- train(x=train_x, y=train_y, method="svmRadial", trControl=train_ctrl_bin, tuneGrid=grid, metric="ROC")
  }
  p1 <- predict(m, test1_x, type="prob")$ALS
  p2 <- predict(m, test2_x, type="prob")$ALS
  list(model=m, pred1=p1, pred2=p2)
}

# Benchmark 循环
benchmark_list <- list()
all_preds <- list()

for (fs_name in names(feature_sets)) {
  fs <- feature_sets[[fs_name]]
  cat(sprintf("\n--- 特征集: %s (%d 个特征) ---\n", fs_name, length(fs)))

  tx <- prepare_features(df_train_bin, fs)
  t1x <- prepare_features(df_test1_bin, fs)
  t2x <- prepare_features(df_test2_bin, fs)
  ty <- df_train_bin$label
  t1y <- df_test1_bin$label
  t2y <- df_test2_bin$label

  for (method in c("LR", "RF", "SVM")) {
    cat(sprintf("  [%s] ... ", method))
    result <- run_model(method, tx, ty, t1x, t2x)
    a1 <- roc(t1y, result$pred1)$auc; a2 <- roc(t2y, result$pred2)$auc
    cat(sprintf("AUC=%s=%.4f, %s=%.4f\n", TEST1_NAME, a1, TEST2_NAME, a2))
    benchmark_list[[paste(fs_name, method)]] <- data.frame(FeatureSet=fs_name, Model=method,
      Test1_AUC=as.numeric(a1), Test2_AUC=as.numeric(a2), stringsAsFactors=FALSE)
    all_preds[[paste(fs_name, method)]] <- list(pred1=result$pred1, pred2=result$pred2, true1=t1y, true2=t2y)
  }
}
benchmark_df <- bind_rows(benchmark_list)
write.csv(benchmark_df, file.path(DIR_CSV, "benchmark_results.csv"), row.names=FALSE)
cat(sprintf("\n>>> Benchmark 结果已保存: %s\n", file.path(DIR_CSV, "benchmark_results.csv")))

# ========== 6. 二分类 XGBoost（使用 markers 特征集） ==========
cat("\n========== XGBoost (Markers 特征集) ==========\n")
tx <- prepare_features(df_train_bin, core_proteins)
ty <- df_train_bin$label
t1x <- prepare_features(df_test1_bin, core_proteins)
t2x <- prepare_features(df_test2_bin, core_proteins)

dtrain <- xgb.DMatrix(tx, label=ifelse(ty=="ALS",1,0))
dtest1 <- xgb.DMatrix(t1x)
dtest2 <- xgb.DMatrix(t2x)

xgb_grid <- expand.grid(max_depth=c(3,5,7), eta=c(0.01,0.1), gamma=c(0,0.5), min_child_weight=c(1,5), nrounds=c(50,100,200))
set.seed(42)
xgb_best <- NULL; best_score <- -1
for (i in 1:nrow(xgb_grid)) {
  p <- xgb_grid[i,]
  cv <- xgb.cv(params=list(objective="binary:logistic", eval_metric="auc", max_depth=p$max_depth, eta=p$eta, gamma=p$gamma, min_child_weight=p$min_child_weight),
               data=dtrain, nrounds=p$nrounds, nfold=5, verbose=0, prediction=FALSE)
  sc <- max(cv$evaluation_log$test_auc_mean)
  if (sc > best_score) { best_score <- sc; xgb_best <- p }
}
cat(sprintf("  Best params: max_depth=%d, eta=%.2f, gamma=%.1f, min_child_weight=%d, nrounds=%d\n",
            xgb_best$max_depth, xgb_best$eta, xgb_best$gamma, xgb_best$min_child_weight, xgb_best$nrounds))
xgb_model <- xgb.train(params=list(objective="binary:logistic", eval_metric="auc", max_depth=xgb_best$max_depth, eta=xgb_best$eta, gamma=xgb_best$gamma, min_child_weight=xgb_best$min_child_weight),
                       data=dtrain, nrounds=xgb_best$nrounds, verbose=0)
xgb_p1 <- predict(xgb_model, dtest1)
xgb_p2 <- predict(xgb_model, dtest2)

benchmark_df <- rbind(benchmark_df, data.frame(FeatureSet="Markers", Model="XGB",
  Test1_AUC=as.numeric(roc(df_test1_bin$label, xgb_p1)$auc),
  Test2_AUC=as.numeric(roc(df_test2_bin$label, xgb_p2)$auc), stringsAsFactors=FALSE))
write.csv(benchmark_df, file.path(DIR_CSV, "benchmark_results.csv"), row.names=FALSE)
all_preds[["Markers XGB"]] <- list(pred1=xgb_p1, pred2=xgb_p2, true1=df_test1_bin$label, true2=df_test2_bin$label)
saveRDS(list(model=xgb_model, params=xgb_best), file.path(DIR_MODEL, "xgb_binary.rds"))
# ---- 旧二分类代码：当未使用新 benchmark 系统时运行 ----
if (!file.exists(file.path(DIR_CSV, "benchmark_results.csv"))) {

# ---- 6.1 逻辑回归 (glm, L2正则化) ----
cat("\n  [LR] 逻辑回归...\n")
set.seed(42)
lr_bin <- train(
  x = train_x, y = train_y,
  method = "glm",
  family = "binomial",
  trControl = train_ctrl_bin,
  metric = "ROC"
)
cat(sprintf("    CV ROC=%g\n", max(lr_bin$results$ROC)))

# ---- 6.2 随机森林 ----
cat("\n  [RF] 随机森林...\n")
p <- length(core_proteins)
rf_grid <- expand.grid(mtry = unique(floor(seq(1, p, length.out = 5))))
set.seed(42)
rf_bin <- train(
  x = train_x, y = train_y,
  method = "rf",
  trControl = train_ctrl_bin,
  tuneGrid = rf_grid,
  metric = "ROC",
  ntree = 500
)
cat(sprintf("    Best mtry=%d, CV ROC=%g\n",
            rf_bin$bestTune$mtry, max(rf_bin$results$ROC)))

# ---- 6.3 SVM (RBF核) ----
cat("\n  [SVM] 支持向量机 (RBF)...\n")
svm_grid <- expand.grid(C = c(0.1, 1, 10, 100), sigma = c(0.001, 0.01, 0.1, 1))
set.seed(42)
svm_bin <- train(
  x = train_x, y = train_y,
  method = "svmRadial",
  trControl = train_ctrl_bin,
  tuneGrid = svm_grid,
  metric = "ROC"
)
cat(sprintf("    Best sigma=%g, C=%g, CV ROC=%g\n",
            svm_bin$bestTune$sigma, svm_bin$bestTune$C, max(svm_bin$results$ROC)))

# ---- 6.4 XGBoost ----
cat("\n  [XGB] XGBoost (direct xgboost, caret兼容接口)...\n")

# caret 与 xgboost >= 3.0 存在 XGBAltrepPointerClass 兼容性问题
# 改用直接 xgboost 训练并封装为 caret 兼容的预测接口
xgb_bin <- NULL
tryCatch({
  set.seed(42)

  # 将数据转换为 xgb.DMatrix
  train_label <- ifelse(train_y == "ALS", 1, 0)
  test1_label <- ifelse(test1_y == "ALS", 1, 0)
  test2_label <- ifelse(test2_y == "ALS", 1, 0)

  dtrain <- xgb.DMatrix(train_x, label = train_label)
  dtest1 <- xgb.DMatrix(test1_x, label = test1_label)
  dtest2 <- xgb.DMatrix(test2_x, label = test2_label)

  # 参数网格搜索（简化版，使用 xgb.cv）
  best_score <- -Inf
  best_params <- NULL
  best_nrounds <- 50

  param_grid <- expand.grid(
    max_depth = c(3, 5),
    eta = c(0.05, 0.1),
    gamma = c(0, 0.1),
    subsample = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = c(1, 3),
    nrounds = c(50, 100)
  )

  for (i in seq_len(nrow(param_grid))) {
    pg <- param_grid[i, ]
    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = pg$max_depth,
      eta = pg$eta,
      gamma = pg$gamma,
      subsample = pg$subsample,
      colsample_bytree = pg$colsample_bytree,
      min_child_weight = pg$min_child_weight
    )
    cv_result <- xgb.cv(
      params = params,
      data = dtrain,
      nrounds = pg$nrounds,
      nfold = 5,
      early_stopping_rounds = 10,
      verbose = 0,
      prediction = FALSE
    )
    score <- max(cv_result$evaluation_log$test_auc_mean)
    if (score > best_score) {
      best_score <- score
      best_params <- params
      best_nrounds <- which.max(cv_result$evaluation_log$test_auc_mean)
    }
  }

  cat(sprintf("    最佳参数: max_depth=%d, eta=%g, gamma=%g, min_child_weight=%g, nrounds=%d\n",
              best_params$max_depth, best_params$eta, best_params$gamma,
              best_params$min_child_weight, best_nrounds))
  cat(sprintf("    CV AUC=%.6f\n", best_score))

  # 用最佳参数训练最终模型
  bst <- xgb.train(params = best_params, data = dtrain, nrounds = best_nrounds)

  # 封装为 caret 兼容的预测接口
  xgb_bin <- list(
    final_bst = bst,
    bestTune = data.frame(nrounds = best_nrounds,
                          max_depth = best_params$max_depth,
                          eta = best_params$eta,
                          gamma = best_params$gamma,
                          subsample = best_params$subsample,
                          colsample_bytree = best_params$colsample_bytree,
                          min_child_weight = best_params$min_child_weight),
    results = data.frame(ROC = best_score),
    method = "xgboost"
  )
  class(xgb_bin) <- "xgb_caret_compat"

  # 定义预测方法
  predict.xgb_caret_compat <- function(object, newdata, type = "raw") {
    if (is.matrix(newdata)) {
      dnew <- xgb.DMatrix(newdata)
    } else {
      dnew <- xgb.DMatrix(as.matrix(newdata))
    }
    prob <- predict(object$final_bst, dnew)
    if (type == "prob") {
      data.frame(Other = 1 - prob, ALS = prob)
    } else {
      factor(ifelse(prob > 0.5, "ALS", "Other"), levels = c("ALS", "Other"))
    }
  }

}, error = function(e) {
  cat(sprintf("    XGBoost 训练失败: %s\n", conditionMessage(e)))
})

# ========== 7. 二分类：评估 ==========

eval_binary <- function(model, test_x, test_y, test_name, model_name, test_samples) {
  pred_prob <- predict(model, test_x, type = "prob")$ALS
  pred_class <- predict(model, test_x)
  
  # ROC
  roc_obj <- roc(test_y, pred_prob, levels = c("Other", "ALS"), direction = "<", quiet = TRUE)
  auc_val <- auc(roc_obj)
  
  # PRC 曲线数据
  pred_obj <- prediction(pred_prob, test_y == "ALS")
  perf_prec <- performance(pred_obj, "prec", "rec")
  prec_vals <- perf_prec@y.values[[1]]
  rec_vals <- perf_prec@x.values[[1]]
  prec_vals[is.na(prec_vals)] <- 0
  
  # AUC-PRC（梯形积分）
  ord <- order(rec_vals)
  auc_prc <- sum(diff(rec_vals[ord]) * (prec_vals[ord][-1] + prec_vals[ord][-length(prec_vals[ord])]) / 2)
  
  # 混淆矩阵指标
  cm <- confusionMatrix(pred_class, test_y, positive = "ALS")
  
  # ROC 坐标（用于平滑 ROC 曲线）
  roc_coords <- coords(roc_obj, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
  
  list(
    model_name = model_name, test_name = test_name,
    test_samples = test_samples, true_label = as.character(test_y),
    pred_prob = pred_prob, pred_class = pred_class,
    roc_obj = roc_obj, auc_roc = auc_val,
    roc_fpr = 1 - roc_coords$specificity,
    roc_tpr = roc_coords$sensitivity,
    prc_data = data.frame(Recall = rec_vals, Precision = prec_vals),
    auc_prc = auc_prc,
    sens = cm$byClass["Sensitivity"],
    spec = cm$byClass["Specificity"],
    ppv = cm$byClass["Pos Pred Value"],
    npv = cm$byClass["Neg Pred Value"],
    f1 = cm$byClass["F1"],
    accuracy = cm$overall["Accuracy"]
  )
}

all_models_bin <- list(
  LR = lr_bin, RF = rf_bin, SVM = svm_bin, XGB = xgb_bin
)
all_models_bin <- all_models_bin[!sapply(all_models_bin, is.null)]

binary_results <- list()
for (nm in names(all_models_bin)) {
  binary_results[[paste0(nm, "_test1")]] <- eval_binary(
    all_models_bin[[nm]], test1_x, test1_y, TEST1_NAME, nm, rownames(test1_x))
  binary_results[[paste0(nm, "_test2")]] <- eval_binary(
    all_models_bin[[nm]], test2_x, test2_y, TEST2_NAME, nm, rownames(test2_x))
}



# ---- 7.2 二分类性能表格 ----
bin_perf <- do.call(rbind, lapply(binary_results, function(r) {
  data.frame(
    Model = r$model_name,
    TestSet = r$test_name,
    AUC_ROC = round(r$auc_roc, 4),
    AUC_PRC = round(r$auc_prc, 4),
    Sensitivity = round(r$sens, 4),
    Specificity = round(r$spec, 4),
    PPV = round(r$ppv, 4),
    NPV = round(r$npv, 4),
    F1 = round(r$f1, 4),
    Accuracy = round(r$accuracy, 4),
    stringsAsFactors = FALSE
  )
}))
write.csv(bin_perf, file.path(DIR_CSV, "binary_performance_table.csv"), row.names = FALSE)
cat(sprintf("  二分类性能表已保存\n"))

# ---- 7.3 输出预测值 ----
bin_pred_df <- do.call(rbind, lapply(binary_results, function(r) {
  data.frame(
    Sample = r$test_samples,
    Model = r$model_name,
    TestSet = r$test_name,
    TrueLabel = r$true_label,
    PredClass = as.character(r$pred_class),
    Prob_ALS = round(r$pred_prob, 6),
    stringsAsFactors = FALSE
  )
}))
write.csv(bin_pred_df, file.path(DIR_CSV, "binary_predictions.csv"), row.names = FALSE)
cat(sprintf("  二分类预测值已保存\n"))

# ---- 7.4 保存最佳模型 (按 Test1 AUC-ROC) ----
cat("\n>>> 保存最佳二分类模型...\n")
test1_perf <- bin_perf[bin_perf$TestSet == TEST1_NAME, ]
best_bin_name <- as.character(test1_perf$Model[which.max(test1_perf$AUC_ROC)])
best_bin_model <- all_models_bin[[best_bin_name]]
best_bin_model_save <- list(
  model = best_bin_model,
  metadata = list(
    model_name = best_bin_name,
    algorithm = best_bin_name,
    training_samples = nrow(df_train),
    features = length(core_proteins),
    feature_list = core_proteins,
    discovery_center = discovery_center,
    performance_Test1 = test1_perf[test1_perf$Model == best_bin_name, ],
    created_at = as.character(Sys.time())
  )
)
saveRDS(best_bin_model_save, file.path(DIR_MODEL, "best_binary_model.rds"))
cat(sprintf("  最佳模型: %s (AUC-ROC = %.4f)\n",
            best_bin_name, max(test1_perf$AUC_ROC)))
cat(sprintf("  模型已保存: model/best_binary_model.rds\n"))

} # ---- 旧二分类代码结束 ----

# ========== 8. 多分类准备 ==========
multiclass_levels <- unique(df_all$group)
cat(sprintf("\n========== 多分类: %s ==========\n", paste(multiclass_levels, collapse = " / ")))

df_train_multi <- df_train
df_train_multi$label <- factor(df_train_multi$group, levels = multiclass_levels)

df_test1_multi <- df_test1
df_test1_multi$label <- factor(df_test1_multi$group, levels = multiclass_levels)

df_test2_multi <- df_test2
df_test2_multi$label <- factor(df_test2_multi$group, levels = multiclass_levels)

train_x_m <- as.matrix(df_train_multi[, core_proteins, drop = FALSE])
test1_x_m <- as.matrix(df_test1_multi[, core_proteins, drop = FALSE])
test2_x_m <- as.matrix(df_test2_multi[, core_proteins, drop = FALSE])

train_y_m <- df_train_multi$label
test1_y_m <- df_test1_multi$label
test2_y_m <- df_test2_multi$label

cat(sprintf("  训练集分布:\n"))
print(table(train_y_m))

# ========== 9. 多分类：训练控制与网格搜索 ==========

train_ctrl_multi <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  savePredictions = "final",
  verboseIter = FALSE
)

# ---- 9.1 逻辑回归 (multinom, 多分类) ----
cat("\n  [LR] 多分类逻辑回归...\n")
set.seed(42)
lr_multi_ctrl <- trainControl(
  method = "cv", number = 5,
  classProbs = TRUE,
  savePredictions = "final", verboseIter = FALSE
)
lr_multi <- train(
  x = train_x_m, y = train_y_m,
  method = "multinom",
  trControl = lr_multi_ctrl,
  metric = "Accuracy",
  trace = FALSE
)
cat(sprintf("    CV Accuracy=%g\n", max(lr_multi$results$Accuracy)))

# ---- 9.2 随机森林 ----
cat("\n  [RF] 多分类随机森林...\n")
p <- ncol(train_x_m)
rf_grid <- expand.grid(mtry = unique(floor(seq(1, p, length.out = min(5, p)))))
set.seed(42)
rf_multi <- train(
  x = train_x_m, y = train_y_m,
  method = "rf",
  trControl = train_ctrl_multi,
  tuneGrid = rf_grid,
  metric = "Accuracy",
  ntree = 500
)
cat(sprintf("    Best mtry=%d, CV Accuracy=%g\n",
            rf_multi$bestTune$mtry, max(rf_multi$results$Accuracy)))

# ---- 9.3 SVM (RBF核) ----
cat("\n  [SVM] 多分类支持向量机...\n")
svm_grid <- expand.grid(C = c(0.1, 1, 10, 100), sigma = c(0.001, 0.01, 0.1, 1))
set.seed(42)
svm_multi <- train(
  x = train_x_m, y = train_y_m,
  method = "svmRadial",
  trControl = train_ctrl_multi,
  tuneGrid = svm_grid,
  metric = "Accuracy"
)
cat(sprintf("    Best sigma=%g, C=%g, CV Accuracy=%g\n",
            svm_multi$bestTune$sigma, svm_multi$bestTune$C, max(svm_multi$results$Accuracy)))

# ---- 9.4 XGBoost ----
cat("\n  [XGB] 多分类 XGBoost (direct xgboost)...\n")
xgb_multi <- NULL
tryCatch({
  set.seed(42)

  # 将数据转换为 xgb.DMatrix
  class_levels <- levels(train_y_m)
  n_class <- length(class_levels)
  train_label_m <- as.numeric(train_y_m) - 1
  test1_label_m <- as.numeric(test1_y_m) - 1
  test2_label_m <- as.numeric(test2_y_m) - 1

  dtrain_m <- xgb.DMatrix(train_x_m, label = train_label_m)
  dtest1_m <- xgb.DMatrix(test1_x_m, label = test1_label_m)
  dtest2_m <- xgb.DMatrix(test2_x_m, label = test2_label_m)

  # 参数搜索
  best_score_m <- -Inf
  best_params_m <- NULL
  best_nrounds_m <- 50

  param_grid_m <- expand.grid(
    max_depth = c(3, 5),
    eta = c(0.05, 0.1),
    gamma = c(0, 0.1),
    subsample = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = c(1, 3),
    nrounds = c(50, 100)
  )

  for (i in seq_len(nrow(param_grid_m))) {
    pg <- param_grid_m[i, ]
    params <- list(
      objective = "multi:softprob",
      eval_metric = "mlogloss",
      num_class = n_class,
      max_depth = pg$max_depth,
      eta = pg$eta,
      gamma = pg$gamma,
      subsample = pg$subsample,
      colsample_bytree = pg$colsample_bytree,
      min_child_weight = pg$min_child_weight
    )
    cv_result <- xgb.cv(
      params = params,
      data = dtrain_m,
      nrounds = pg$nrounds,
      nfold = 5,
      early_stopping_rounds = 10,
      verbose = 0,
      prediction = FALSE
    )
    score <- -min(cv_result$evaluation_log$test_mlogloss_mean)
    if (score > best_score_m) {
      best_score_m <- score
      best_params_m <- params
      best_nrounds_m <- which.min(cv_result$evaluation_log$test_mlogloss_mean)
    }
  }

  cat(sprintf("    最佳参数: max_depth=%d, eta=%g, gamma=%g, min_child_weight=%g, nrounds=%d\n",
              best_params_m$max_depth, best_params_m$eta, best_params_m$gamma,
              best_params_m$min_child_weight, best_nrounds_m))

  # 用最佳参数训练最终模型
  bst_m <- xgb.train(params = best_params_m, data = dtrain_m, nrounds = best_nrounds_m)

  # 封装为 caret 兼容的预测接口
  xgb_multi <- list(
    final_bst = bst_m,
    bestTune = data.frame(nrounds = best_nrounds_m,
                          max_depth = best_params_m$max_depth,
                          eta = best_params_m$eta),
    results = data.frame(Accuracy = exp(best_score_m) / (1 + exp(best_score_m))),
    method = "xgboost",
    class_levels = class_levels
  )
  class(xgb_multi) <- "xgb_caret_compat_multi"

  # 定义预测方法
  predict.xgb_caret_compat_multi <- function(object, newdata, type = "raw") {
    if (is.matrix(newdata)) {
      dnew <- xgb.DMatrix(newdata)
    } else {
      dnew <- xgb.DMatrix(as.matrix(newdata))
    }
    prob_raw <- predict(object$final_bst, dnew)
    # multi:softprob 的 predict() 已返回 (n × n_class) 矩阵；仅当返回向量时才需重排
    prob_matrix <- if (is.matrix(prob_raw)) prob_raw
                   else matrix(prob_raw, ncol = length(object$class_levels), byrow = TRUE)
    colnames(prob_matrix) <- object$class_levels
    if (type == "prob") {
      as.data.frame(prob_matrix)
    } else {
      pred_idx <- max.col(prob_matrix)
      factor(object$class_levels[pred_idx], levels = object$class_levels)
    }
  }

  cat(sprintf("    训练完成，类别: %s\n", paste(class_levels, collapse = ", ")))
}, error = function(e) {
  cat(sprintf("    多分类 XGBoost 训练失败: %s\n", conditionMessage(e)))
})

# ========== 10. 多分类：评估 ==========

eval_multiclass <- function(model, test_x, test_y, test_name, model_name, test_samples) {
  pred_class <- predict(model, test_x)
  pred_prob <- predict(model, test_x, type = "prob")
  
  cm <- confusionMatrix(pred_class, test_y)
  
  # Per-class metrics
  per_class <- as.data.frame(cm$byClass)
  per_class$Class <- rownames(per_class)
  
  list(
    model_name = model_name, test_name = test_name,
    test_samples = test_samples, true_label = as.character(test_y),
    pred_class = pred_class, pred_prob = pred_prob,
    confusion = cm$table,
    overall = cm$overall,
    per_class = per_class,
    accuracy = cm$overall["Accuracy"],
    kappa = cm$overall["Kappa"]
  )
}

all_models_multi <- list(
  LR = lr_multi, RF = rf_multi, SVM = svm_multi, XGB = xgb_multi
)
all_models_multi <- all_models_multi[!sapply(all_models_multi, is.null)]

multi_results <- list()
for (nm in names(all_models_multi)) {
  multi_results[[paste0(nm, "_test1")]] <- eval_multiclass(
    all_models_multi[[nm]], test1_x_m, test1_y_m, TEST1_NAME, nm, rownames(test1_x_m))
  multi_results[[paste0(nm, "_test2")]] <- eval_multiclass(
    all_models_multi[[nm]], test2_x_m, test2_y_m, TEST2_NAME, nm, rownames(test2_x_m))
}

# ---- 10.1 混淆矩阵热图（由绘图脚本生成）----

# ---- 10.2 多分类性能表格 ----
multi_perf <- do.call(rbind, lapply(multi_results, function(r) {
  data.frame(
    Model = r$model_name,
    TestSet = r$test_name,
    Accuracy = round(r$accuracy, 4),
    Kappa = round(r$kappa, 4),
    stringsAsFactors = FALSE
  )
}))

# 添加 per-class 指标
multiclass_class_names <- paste0("Class: ", multiclass_levels)
for (cls in multiclass_class_names) {
  col_sens <- paste0(gsub("Class: ", "", cls), "_Sensitivity")
  col_spec <- paste0(gsub("Class: ", "", cls), "_Specificity")
  multi_perf[[col_sens]] <- sapply(multi_results, function(r) {
    idx <- grep(cls, rownames(r$per_class))
    if (length(idx)) round(r$per_class[idx, "Sensitivity"], 4) else NA
  })
  multi_perf[[col_spec]] <- sapply(multi_results, function(r) {
    idx <- grep(cls, rownames(r$per_class))
    if (length(idx)) round(r$per_class[idx, "Specificity"], 4) else NA
  })
}

write.csv(multi_perf, file.path(DIR_CSV, "multiclass_performance_table.csv"), row.names = FALSE)
cat(sprintf("  多分类性能表已保存\n"))

# ---- 10.3 输出多分类预测值 ----
multi_pred_df <- do.call(rbind, lapply(multi_results, function(r) {
  df <- data.frame(
    Sample = r$test_samples,
    Model = r$model_name,
    TestSet = r$test_name,
    TrueLabel = r$true_label,
    PredClass = as.character(r$pred_class),
    stringsAsFactors = FALSE
  )
  # 添加各类的概率
  prob_cols <- as.data.frame(r$pred_prob)
  cbind(df, prob_cols)
}))
write.csv(multi_pred_df, file.path(DIR_CSV, "multiclass_predictions.csv"), row.names = FALSE)
cat(sprintf("  多分类预测值已保存\n"))

# ---- 10.3b 完整多分类性能表（Macro AUC / F1 等；供正文与图表核对）----
full_rows <- lapply(multi_results, function(r) {
  cm <- r$confusion
  probs <- as.data.frame(r$pred_prob)
  per <- lapply(multiclass_levels, function(cl) {
    tp <- cm[cl, cl]
    fp <- sum(cm[cl, ]) - tp
    fn <- sum(cm[, cl]) - tp
    tn <- sum(cm) - tp - fp - fn
    sens <- ifelse(tp + fn > 0, tp / (tp + fn), NA_real_)
    spec <- ifelse(tn + fp > 0, tn / (tn + fp), NA_real_)
    prec <- ifelse(tp + fp > 0, tp / (tp + fp), NA_real_)
    f1   <- ifelse(2 * tp + fp + fn > 0, 2 * tp / (2 * tp + fp + fn), NA_real_)
    auc  <- as.numeric(pROC::auc(pROC::roc(
      ifelse(r$true_label == cl, 1, 0), probs[[cl]], quiet = TRUE, direction = "<")))
    c(Sensitivity = sens, Specificity = spec, Precision = prec, F1 = f1, AUC = auc)
  })
  per_df <- as.data.frame(do.call(rbind, per))
  per_df$Class <- multiclass_levels
  macro <- colMeans(per_df[, c("Sensitivity", "Specificity", "Precision", "F1", "AUC")],
                    na.rm = TRUE)
  out <- data.frame(
    Model = r$model_name, TestSet = r$test_name,
    Accuracy = round(r$accuracy, 4), Kappa = round(r$kappa, 4),
    Macro_AUC = round(macro["AUC"], 4),
    Macro_Sensitivity = round(macro["Sensitivity"], 4),
    Macro_Specificity = round(macro["Specificity"], 4),
    Macro_Precision = round(macro["Precision"], 4),
    Macro_F1 = round(macro["F1"], 4),
    stringsAsFactors = FALSE)
  for (cl in multiclass_levels) {
    row <- per_df[per_df$Class == cl, ]
    out[[paste0(cl, "_Sensitivity")]] <- round(row$Sensitivity, 4)
    out[[paste0(cl, "_Specificity")]] <- round(row$Specificity, 4)
    out[[paste0(cl, "_F1")]] <- round(row$F1, 4)
  }
  out
})
full_perf <- bind_rows(full_rows)
write.csv(full_perf, file.path(DIR_CSV, "multiclass_full_performance.csv"), row.names = FALSE)
cat(sprintf("  完整多分类性能表已保存\n"))

# ---- 10.3c NfL vs 12-Marker AUC（按 Sample 配对，避免旧二分类文件）----
nfl_file <- file.path(BASE_DIR, "data", "raw", "nfl.xlsx")
if (file.exists(nfl_file) && requireNamespace("readxl", quietly = TRUE)) {
  nfl_raw <- as.data.frame(readxl::read_excel(nfl_file, sheet = 1))
  samp_col <- names(nfl_raw)[grepl("编号", names(nfl_raw))][1]
  if (is.na(samp_col)) {
    samp_col <- names(nfl_raw)[vapply(nfl_raw, function(x) {
      v <- as.character(x[!is.na(x)])
      length(v) > 0 && grepl("^[A-Z]", v[1])
    }, logical(1))][1]
  }
  nfl_col <- names(nfl_raw)[grepl("NFL", names(nfl_raw), ignore.case = TRUE)][1]
  if (!is.na(samp_col) && !is.na(nfl_col)) {
    nfl_df <- data.frame(
      Sample = as.character(nfl_raw[[samp_col]]),
      NFL = suppressWarnings(as.numeric(nfl_raw[[nfl_col]])),
      stringsAsFactors = FALSE) %>%
      filter(!is.na(Sample), !is.na(NFL)) %>%
      group_by(Sample) %>% summarise(NFL = mean(NFL), .groups = "drop")
    rf_pred <- multi_pred_df %>%
      filter(Model == "RF", TrueLabel %in% c("ALS", "DC")) %>%
      select(Sample, TestSet, TrueLabel, ALS)
    nfl_match <- inner_join(rf_pred, nfl_df, by = "Sample")
    nfl_rows <- list()
    for (ts in unique(nfl_match$TestSet)) {
      d <- nfl_match %>% filter(TestSet == ts)
      if (nrow(d) < 5 || length(unique(d$TrueLabel)) < 2) next
      roc_nfl <- pROC::roc(ifelse(d$TrueLabel == "ALS", 1, 0), d$NFL,
                           quiet = TRUE, direction = "<")
      roc_rf  <- pROC::roc(ifelse(d$TrueLabel == "ALS", 1, 0), d$ALS,
                           quiet = TRUE, direction = "<")
      ci_n <- as.numeric(pROC::ci.auc(roc_nfl, method = "delong"))
      ci_r <- as.numeric(pROC::ci.auc(roc_rf, method = "delong"))
      pv <- tryCatch(pROC::roc.test(roc_nfl, roc_rf, paired = TRUE, method = "delong")$p.value,
                     error = function(e) NA_real_)
      nfl_rows[[length(nfl_rows) + 1]] <- data.frame(
        Model = c("NFL", "12-Marker RF"), TestSet = ts,
        AUC = round(c(as.numeric(pROC::auc(roc_nfl)), as.numeric(pROC::auc(roc_rf))), 4),
        CI_low = round(c(ci_n[1], ci_r[1]), 4),
        CI_high = round(c(ci_n[3], ci_r[3]), 4),
        P_value = c(pv, pv), stringsAsFactors = FALSE)
    }
    if (length(nfl_rows) > 0) {
      write.csv(bind_rows(nfl_rows), file.path(DIR_CSV, "nfl_vs_marker_auc.csv"), row.names = FALSE)
      cat(sprintf("  NfL vs 12-Marker AUC 已保存（按 Sample 配对）\n"))
    }
  }
}

# ---- 10.4 保存最佳多分类模型 (按 Test1 Accuracy) ----
cat("\n>>> 保存最佳多分类模型...\n")
test1_m_perf <- multi_perf[multi_perf$TestSet == TEST1_NAME, ]
best_m_name <- as.character(test1_m_perf$Model[which.max(test1_m_perf$Accuracy)])
best_m_model <- all_models_multi[[best_m_name]]
best_m_model_save <- list(
  model = best_m_model,
  metadata = list(
    model_name = best_m_name,
    algorithm = best_m_name,
    training_samples = nrow(df_train_multi),
    features = length(core_proteins),
    feature_list = core_proteins,
    classes = multiclass_levels,
    discovery_center = discovery_center,
    performance_Test1 = test1_m_perf[test1_m_perf$Model == best_m_name, ],
    created_at = as.character(Sys.time())
  )
)
saveRDS(best_m_model_save, file.path(DIR_MODEL, "best_multiclass_model.rds"))
cat(sprintf("  最佳模型: %s (Accuracy = %.4f)\n",
            best_m_name, max(test1_m_perf$Accuracy)))
cat(sprintf("  模型已保存: model/best_multiclass_model.rds\n"))

# ========== 11. SHAP 分析 ==========
cat("\n========== SHAP 分析 ==========\n")

# ---- 11a. 二分类 SHAP（仅旧代码路径） ----
if (!file.exists(file.path(DIR_CSV, "benchmark_results.csv"))) {
cat(">>> 使用 fastshap::explain() 加速计算 SHAP 值...\n")

compute_shap <- function(model, test_x, bg_x, nsim = 50) {
  # 使用 fastshap::explain() 进行快速的 Monte Carlo SHAP 计算
  # fastshap 使用向量化操作，远快于手动循环
  pred_wrapper <- function(object, newdata) {
    predict(object, newdata = newdata, type = "prob")[, 2]
  }
  
  shap_vals <- fastshap::explain(
    object = model,
    X = as.data.frame(test_x),
    nsim = nsim,
    pred_wrapper = pred_wrapper,
    adjust = TRUE
  )
  
  feat_names <- colnames(test_x)
  baseline <- mean(pred_wrapper(model, as.data.frame(bg_x)))
  
  list(shap_values = as.matrix(shap_vals), baseline = baseline, feat_names = feat_names)
}

# 使用训练集作为背景数据
bg_x <- train_x
shap_result <- compute_shap(best_bin_model, rbind(test1_x, test2_x), bg_x, nsim = 80)

# 保存 SHAP 值到 CSV
shap_df <- as.data.frame(shap_result$shap_values)
shap_df$Sample <- c(rownames(test1_x), rownames(test2_x))
shap_df$TestSet <- c(rep(TEST1_NAME, nrow(test1_x)), rep(TEST2_NAME, nrow(test2_x)))
test_all_y <- c(as.character(test1_y), as.character(test2_y))
shap_df$TrueLabel <- test_all_y
write.csv(shap_df, file.path(DIR_CSV, "shap_values.csv"), row.names = FALSE)
cat(sprintf("  SHAP 值已保存\n"))

# 保存原始特征值（用于 beeswarm 图的特征值颜色映射）
test_combined <- rbind(test1_x, test2_x)
feat_df <- as.data.frame(test_combined)
feat_df$Sample <- c(rownames(test1_x), rownames(test2_x))
feat_df$TestSet <- c(rep(TEST1_NAME, nrow(test1_x)), rep(TEST2_NAME, nrow(test2_x)))
feat_df$TrueLabel <- test_all_y
write.csv(feat_df, file.path(DIR_CSV, "shap_feature_values.csv"), row.names = FALSE)
cat(sprintf("  原始特征值已保存（用于 SHAP beeswarm 图）\n"))

# 全局特征重要性（平均 |SHAP|）
shap_imp <- data.frame(
  Feature = shap_result$feat_names,
  MeanAbsSHAP = colMeans(abs(shap_result$shap_values)),
  stringsAsFactors = FALSE
) %>% arrange(desc(MeanAbsSHAP))

write.csv(shap_imp, file.path(DIR_CSV, "shap_importance.csv"), row.names = FALSE)
cat(sprintf("  SHAP 特征重要性已保存\n"))

} # ---- 二分类 SHAP 结束 ----

# ========== 11b. 多分类 SHAP 分析 ==========
cat("\n>>> 多分类 SHAP 分析（基于最佳多分类模型）...\n")

compute_shap_multi <- function(model, test_x, bg_x, classes, nsim = 50) {
  # 对每个类别分别计算 SHAP（fastshap 要求输出为数值向量）
  shap_by_class <- list()
  for (cl in classes) {
    pred_wrapper_cl <- function(object, newdata) {
      pred <- predict(object, newdata = newdata, type = "prob")
      pred[[cl]]  # 返回该类别的概率向量
    }
    shap_by_class[[cl]] <- fastshap::explain(
      object = model,
      X = as.data.frame(test_x),
      nsim = nsim,
      pred_wrapper = pred_wrapper_cl,
      adjust = TRUE
    )
  }
  list(shap_values = shap_by_class, feat_names = colnames(test_x))
}

if (!is.null(best_m_model)) {
  cat(sprintf("  分析模型: %s\n", best_m_name))
  classes <- levels(train_y_m)
  shap_m_result <- compute_shap_multi(best_m_model, rbind(test1_x_m, test2_x_m),
                                       train_x_m, classes, nsim = 50)
  
  # 保存多分类 SHAP 值
  for (cl in names(shap_m_result$shap_values)) {
    sv <- shap_m_result$shap_values[[cl]]
    df_sv <- as.data.frame(sv)
    df_sv$Sample <- c(rownames(test1_x_m), rownames(test2_x_m))
    df_sv$TestSet <- c(rep(TEST1_NAME, nrow(test1_x_m)),
                        rep(TEST2_NAME, nrow(test2_x_m)))
    write.csv(df_sv,
              file.path(DIR_CSV, paste0("shap_values_multi_", cl, ".csv")),
              row.names = FALSE)
  }

  # 保存多分类 SHAP 对应的原始特征值（beeswarm 着色需要；与 SHAP 行顺序一致）
  feat_m <- as.data.frame(rbind(test1_x_m, test2_x_m))
  feat_m$Sample <- c(rownames(test1_x_m), rownames(test2_x_m))
  feat_m$TestSet <- c(rep(TEST1_NAME, nrow(test1_x_m)),
                      rep(TEST2_NAME, nrow(test2_x_m)))
  feat_m$TrueLabel <- c(as.character(test1_y_m), as.character(test2_y_m))
  write.csv(feat_m, file.path(DIR_CSV, "shap_feature_values_multi.csv"), row.names = FALSE)

  # 全局特征重要性（多分类平均）
  shap_m_imp_list <- list()
  for (cl in names(shap_m_result$shap_values)) {
    sv <- shap_m_result$shap_values[[cl]]
    shap_m_imp_list[[cl]] <- data.frame(
      Feature = shap_m_result$feat_names,
      Class = cl,
      MeanAbsSHAP = colMeans(abs(sv)),
      stringsAsFactors = FALSE
    )
  }
  shap_m_imp <- bind_rows(shap_m_imp_list)
  write.csv(shap_m_imp, file.path(DIR_CSV, "shap_importance_multi.csv"),
            row.names = FALSE)
  cat(sprintf("  多分类 SHAP 已保存 (shap_values_multi_*.csv, shap_importance_multi.csv)\n"))
} else {
  cat("  多分类最佳模型为空，跳过 SHAP 分析\n")
}

# ========== 12. 从多分类模型派生二分类 ROC ==========
cat("\n========== 多分类 → 二分类 (ALS概率) ==========\n")
if (!is.null(best_m_model) && length(multi_results) > 0) {
  # 找到最佳模型的多分类结果
  mc_key1 <- paste0(best_m_name, "_test1")
  mc_key2 <- paste0(best_m_name, "_test2")

  for (kv in list(c(mc_key1, TEST1_NAME), c(mc_key2, TEST2_NAME))) {
    key <- kv[1]; tname <- kv[2]
    if (!key %in% names(multi_results)) next
    res <- multi_results[[key]]
    pp <- res$pred_prob
    if (is.null(pp) || !"ALS" %in% colnames(pp)) next
    als_prob <- pp[, "ALS"]
    df_orig <- if (tname == TEST1_NAME) df_test1 else df_test2
    true_bin <- factor(ifelse(df_orig$group == "ALS", "ALS", "Other"), levels = c("ALS", "Other"))
    roc_obj <- roc(true_bin, als_prob)
    cat(sprintf("  %s: AUC = %.4f (从%s多分类ALS概率派生)\n", tname, auc(roc_obj), best_m_name))
  }
}

# ========== 13. ALS 亚组分析 ==========
cat("\n========== ALS 亚组分析 ==========\n")
if (!is.null(best_m_model) && length(multi_results) > 0) {
  mc_key1 <- paste0(best_m_name, "_test1")
  mc_key2 <- paste0(best_m_name, "_test2")

  subgroup_results <- list()
  subgroup_vars <- list(
    age = list(name = "Age", breaks = c(0, 60, 120), labels = c("<60", ">=60")),
    gender = list(name = "Gender", breaks = NULL, labels = NULL),
    onset_site = list(name = "Onset", breaks = NULL, labels = NULL),
    kings_stage = list(name = "King's Stage", breaks = c(0, 2, 5), labels = c("Stage1-2", "Stage3-4")),
    disease_duration = list(name = "Duration", breaks = NULL, labels = NULL)
  )

  for (sv_name in names(subgroup_vars)) {
    sv <- subgroup_vars[[sv_name]]
    for (kv in list(c(mc_key1, TEST1_NAME), c(mc_key2, TEST2_NAME))) {
      key <- kv[1]; tname <- kv[2]
      if (!key %in% names(multi_results)) next
      res <- multi_results[[key]]
      pp <- res$pred_prob
      df_orig <- if (tname == TEST1_NAME) df_test1 else df_test2
      if (is.null(pp) || !"ALS" %in% colnames(pp) || sum(df_orig$group == "ALS") < 10) next

      als_only <- df_orig$group == "ALS"
      als_prob <- pp[als_only, "ALS"]

      if (sv_name == "gender") {
        groups <- df_orig$gender[als_only]
      } else if (sv_name == "onset_site") {
        groups <- df_orig$onset_site[als_only]
      } else if (sv_name == "kings_stage") {
        groups <- cut(as.numeric(df_orig$kings_stage[als_only]), breaks = sv$breaks, labels = sv$labels)
      } else if (sv_name %in% c("age", "disease_duration")) {
        val <- as.numeric(df_orig[[sv_name]][als_only])
        med <- median(val, na.rm = TRUE)
        groups <- ifelse(val >= med, paste0(">=", round(med)), paste0("<", round(med)))
      } else { next }

      valid <- !is.na(groups) & !is.na(als_prob)
      if (length(unique(groups[valid])) < 2 || sum(valid) < 10) next

      for (g in unique(groups[valid])) {
        idx <- which(groups == g & valid)
        subgroup_results[[length(subgroup_results) + 1]] <- data.frame(
          Variable = sv$name, Subgroup = g, TestSet = tname,
          N = length(idx), MeanProb = mean(als_prob[idx]),
          SD_Prob = sd(als_prob[idx]), stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(subgroup_results) > 0) {
    sub_df <- bind_rows(subgroup_results)
    write.csv(sub_df, file.path(DIR_CSV, "subgroup_analysis.csv"), row.names = FALSE)
    cat(sprintf("  亚组分析已保存: %s\n", file.path(DIR_CSV, "subgroup_analysis.csv")))
    print(sub_df, row.names = FALSE)
  }
}

# ========== 14. 保存 benchmark 预测值（供绘图脚本用） ==========
cat(">>> 保存 benchmark 预测值...\n")
bench_preds <- list()
for (nm in names(all_preds)) {
  p <- all_preds[[nm]]
  len1 <- length(p$pred1); len2 <- length(p$pred2)
  bench_preds[[nm]] <- data.frame(ModelSet = rep(nm, len1 + len2),
             Test = c(rep(TEST1_NAME, len1), rep(TEST2_NAME, len2)),
             Pred = c(p$pred1, p$pred2),
             True = c(as.character(p$true1), as.character(p$true2)),
             stringsAsFactors = FALSE)
}
saveRDS(list(
  benchmark = benchmark_df,
  predictions = do.call(rbind, bench_preds),
  feature_sets = feature_sets,
  core_proteins = core_proteins,
  test1_name = TEST1_NAME, test2_name = TEST2_NAME
), file = file.path(DIR_CSV, "..", "benchmark_data.rds"))
cat("  已保存: benchmark_data.rds\n")



# ========== 12. 输出模型使用说明 ==========
cat("\n>>> 生成模型使用说明...\n")
if (file.exists(file.path(DIR_CSV, "benchmark_results.csv"))) {
  # 新 pipeline: 使用多分类模型
  readme_lines <- c(
    "# EV-PBA 诊断模型使用说明（新版本 - 多分类统一模型）",
    "",
    "## 概述",
    "",
    paste("本模型基于", length(core_proteins), "个外泌体膜蛋白标志物，仅使用多分类模型进行预测。"),
    "二分类 ROC 从多分类模型中 ALS 类别的概率派生。",
    paste("训练集:", nrow(df_train), sprintf("样本, 测试集 (%s):", TEST1_NAME), nrow(df_test1),
          sprintf("样本, 测试集 (%s):", TEST2_NAME), nrow(df_test2), "样本"),
    "",
    paste("最佳多分类模型: **", best_m_name, "**"),
    paste(TEST1_NAME, "Accuracy:", round(max(multi_perf$Accuracy[multi_perf$TestSet == TEST1_NAME]), 4)),
    "",
    "## 特征集 Benchmark",
    "支持三种特征集: Markers, Metadata, Markers+Metadata",
    paste("Benchmark 结果:", file.path(DIR_CSV, "benchmark_results.csv")),
    "",
    "## 模型文件",
    paste0('  - best_multiclass_model.rds: "', file.path(DIR_MODEL, "best_multiclass_model.rds"), '"'),
    "",
    "## ALS 亚组分析",
    "支持按 Age, Gender, Onset Site, King's Stage, Disease Duration 分层分析",
    paste("亚组分析结果:", file.path(DIR_CSV, "subgroup_analysis.csv"))
  )
} else {
  # 旧 pipeline (有 AD 的完整数据)
  readme_lines <- c(
  "# EV-PBA 诊断模型使用说明",
  "",
  "## 概述",
  "",
  paste("本模型基于", length(core_proteins), "个外泌体膜蛋白标志物，用于 ALS 诊断分类。"),
  paste("训练集:", nrow(df_train), sprintf("样本, 测试集 (%s):", TEST1_NAME), nrow(df_test1), sprintf("样本, 测试集 (%s):", TEST2_NAME), nrow(df_test2), "样本"),
  "",
  "## 二分类模型（ALS vs Others）",
  "",
  paste("最佳模型: **", best_bin_name, "**"),
  paste(TEST1_NAME, "AUC-ROC:", round(max(test1_perf$AUC_ROC), 4)),
  "",
  "### 使用方法",
  "",
  "```r",
  "# 加载模型",
  paste0('  model <- readRDS("', file.path(DIR_MODEL, "best_binary_model.rds"), '")'),
  "",
  "# 准备表达矩阵 (行=蛋白, 列=样本)",
  paste0('  core_proteins <- c("', paste(core_proteins, collapse = '", "'), '")'),
  "expr_matrix <- ...  # log2(x+1) 标准化后的表达矩阵",
  "test_x <- as.matrix(expr_matrix[core_proteins, ])",
  "",
  "# 预测概率",
  "pred_prob <- predict(model, t(test_x), type = \"prob\")$ALS",
  "",
  "# 预测类别",
  "pred_class <- predict(model, t(test_x))",
  "```",
  "",
  "### 标志物列表",
  "",
  paste(sprintf("1. **%s**", core_proteins), collapse = "\n"),
  "",
  "## 多分类模型（ALS / AD / HC / DC）",
  "",
  paste("最佳模型: **", best_m_name, "**"),
  paste(TEST1_NAME, "Accuracy:", round(max(test1_m_perf$Accuracy), 4)),
  "",
  "## 文件说明",
  "",
  "| 文件 | 说明 |",
  "|------|------|",
  paste0("| `model/best_binary_model.rds` | 最佳二分类模型 (", best_bin_name, ") |"),
  paste0("| `model/best_multiclass_model.rds` | 最佳多分类模型 (", best_m_name, ") |"),
  "| `csv/binary_predictions.csv` | 二分类所有模型预测值 |",
  "| `csv/multiclass_predictions.csv` | 多分类所有模型预测值 |",
  "| `csv/shap_values.csv` | 二分类最佳模型 SHAP 值 |",
  "| `csv/shap_importance.csv` | SHAP 全局特征重要性 |",
  "",
  "## 依赖",
  "",
  "- R >= 4.0",
  "- caret, randomForest, e1071, xgboost, nnet"
  )
} # else branch 结束

writeLines(readme_lines, file.path(DIR_MODEL, "README.md"))
cat(sprintf("  模型使用说明已保存\n"))

# ========== 14. 网格搜索参数表 ==========
if (!file.exists(file.path(DIR_CSV, "benchmark_results.csv"))) {
# 二分类网格搜索参数表（LR 为无参模型，单独说明）
lr_grid_note <- data.frame(Model = "LR", Parameter = "N/A", Value = "Standard GLM (no tuning)",
                             ROC = round(max(lr_bin$results$ROC), 4), stringsAsFactors = FALSE)
rf_note <- data.frame(Model = "RF", Parameter = paste0("mtry=", rf_bin$results$mtry),
                      Value = "", ROC = round(rf_bin$results$ROC, 4), stringsAsFactors = FALSE)
svm_note <- data.frame(Model = "SVM", Parameter = paste0("sigma=", svm_bin$results$sigma, ",C=", svm_bin$results$C),
                       Value = "", ROC = round(svm_bin$results$ROC, 4), stringsAsFactors = FALSE)
xgb_note <- NULL
if (!is.null(xgb_bin)) {
  xgb_note <- data.frame(Model = "XGB", Parameter = paste0("nrounds=", xgb_bin$bestTune$nrounds, ",max_depth=", xgb_bin$bestTune$max_depth, ",eta=", xgb_bin$bestTune$eta),
                         Value = "", ROC = round(xgb_bin$results$ROC, 4), stringsAsFactors = FALSE)
}

all_grid <- rbind(lr_grid_note, rf_note, svm_note, xgb_note)
write.csv(all_grid, file.path(DIR_CSV, "binary_grid_search_params.csv"), row.names = FALSE)

# 多分类网格搜索参数表
lr_grid_m_note <- data.frame(Model = "LR", Parameter = "N/A",
                             Value = "Multinomial GLM (no tuning)",
                             Accuracy = round(max(lr_multi$results$Accuracy), 4),
                             stringsAsFactors = FALSE)
rf_m_note <- data.frame(Model = "RF", Parameter = paste0("mtry=", rf_multi$results$mtry),
                        Value = "", Accuracy = round(rf_multi$results$Accuracy, 4), stringsAsFactors = FALSE)
svm_m_note <- data.frame(Model = "SVM", Parameter = paste0("sigma=", svm_multi$results$sigma, ",C=", svm_multi$results$C),
                         Value = "", Accuracy = round(svm_multi$results$Accuracy, 4), stringsAsFactors = FALSE)
xgb_m_note <- NULL
if (!is.null(xgb_multi)) {
  xgb_m_note <- data.frame(Model = "XGB", Parameter = paste0("nrounds=", xgb_multi$bestTune$nrounds, ",max_depth=", xgb_multi$bestTune$max_depth, ",eta=", xgb_multi$bestTune$eta),
                           Value = "", Accuracy = round(xgb_multi$results$Accuracy, 4), stringsAsFactors = FALSE)
}

all_grid_m <- rbind(lr_grid_m_note, rf_m_note, svm_m_note, xgb_m_note)
write.csv(all_grid_m, file.path(DIR_CSV, "multiclass_grid_search_params.csv"), row.names = FALSE)
} # ---- 旧网格搜索参数表结束 ----

# ========== 15. 分析报告 ==========
cat("\n>>> 生成分析报告...\n")

report_lines <- c(
  "=========================================================================",
  "  EV-PBA 模块二：机器学习建模与验证报告",
  "=========================================================================",
  "",
  paste("  生成时间:", Sys.time()),
  paste("  发现中心:", discovery_center),
  "",
  "-----------------------------------------------------------------------",
  "  数据集概览",
  "-----------------------------------------------------------------------",
  "",
  paste("  总样本数:", nrow(df_all)),
  paste("  核心标志物数:", length(core_proteins)),
  paste("  标志物列表:", paste(core_proteins, collapse = ", ")),
  paste("  分组:", paste(names(table(df_all$group)), "=", table(df_all$group), collapse = ", ")),
  "",
  paste("  训练集:", nrow(df_train), "样本"),
  paste("  测试集(", TEST1_NAME, "):", nrow(df_test1), "样本"),
  paste("  测试集(", TEST2_NAME, "):", nrow(df_test2), "样本"),
  "",
  "-----------------------------------------------------------------------",
  "  多分类结果 (ALS / AD / HC / DC)",
  "-----------------------------------------------------------------------",
  "",
  "  性能总表:",
  "  Model  TestSet  Accuracy  Kappa"
)

for (i in seq_len(nrow(multi_perf))) {
  r <- multi_perf[i, ]
  report_lines <- c(report_lines,
    sprintf("  %-4s  %-7s  %-8.4f  %-5.4f",
            r$Model, r$TestSet, r$Accuracy, r$Kappa))
}

report_lines <- c(report_lines,
  "",
  "  每类 Sensitivity / Specificity:",
  "  Model  TestSet  ALS_Sens  ALS_Spec  AD_Sens  AD_Spec  HC_Sens  HC_Spec  DC_Sens  DC_Spec"
)

for (i in seq_len(nrow(multi_perf))) {
  r <- multi_perf[i, ]
  report_lines <- c(report_lines,
    sprintf("  %-4s  %-7s  %-8.4f  %-8.4f  %-7.4f  %-7.4f  %-7.4f  %-7.4f  %-7.4f  %-7.4f",
            r$Model, r$TestSet,
            r$ALS_Sensitivity, r$ALS_Specificity,
            r$AD_Sensitivity, r$AD_Specificity,
            r$HC_Sensitivity, r$HC_Specificity,
            r$DC_Sensitivity, r$DC_Specificity))
}

report_lines <- c(report_lines,
  "",
  "-----------------------------------------------------------------------",
  "  多分类 SHAP 分析（基于最佳多分类模型）",
  "-----------------------------------------------------------------------",
  "",
  paste("  分析模型:", best_m_name),
  "  已保存每个类别的 SHAP 值和全局特征重要性",
  "",
  "-----------------------------------------------------------------------",
  "  输出文件",
  "-----------------------------------------------------------------------",
  "",
  paste("  结果目录:", OUTPUT_DIR),
  "",
  "  图片 (fig/):",
  paste0("    - ROC_", TEST1_NAME, ".png / ROC_", TEST2_NAME, ".png: ROC 曲线"),
  paste0("    - PRC_", TEST1_NAME, ".png / PRC_", TEST2_NAME, ".png: PRC 曲线"),
  paste0("    - confusion_matrix_", TEST1_NAME, ".png / confusion_matrix_", TEST2_NAME, ".png: 混淆矩阵热图"),
  "    - shap_importance.png: SHAP 重要性柱状图",
  "    - shap_beeswarm.png: SHAP Beeswarm 图",
  "",
  "  数据 (csv/):",
  "    - binary_performance_table.csv: 二分类性能表",
  "    - binary_predictions.csv: 二分类预测值",
  "    - binary_grid_search_params.csv: 二分类网格搜索参数表",
  "    - multiclass_performance_table.csv: 多分类性能表",
  "    - multiclass_predictions.csv: 多分类预测值",
  "    - multiclass_grid_search_params.csv: 多分类网格搜索参数表",
  "    - shap_values.csv: 二分类 SHAP 值",
  "    - shap_feature_values.csv: 二分类原始特征值",
  "    - shap_importance.csv: 二分类 SHAP 重要性",
  "    - shap_values_multi_*.csv: 多分类 SHAP 值（每个类别一个文件）",
  "    - shap_importance_multi.csv: 多分类 SHAP 重要性",
  "",
  "  模型:",
  paste0("    - model/best_multiclass_model.rds: 最佳多分类模型 (", best_m_name, ")"),
  "",
  "  报告:",
  "    - module2_report.txt: 详细分析报告",
  "",
  "  其他:",
  "    - model/README.md: 模型使用说明",
  "",
  "========================================================================="
)

writeLines(report_lines, file.path(OUTPUT_DIR, "module2_report.txt"))
cat(sprintf("  报告已保存 -> %s\n", file.path(OUTPUT_DIR, "module2_report.txt")))

# 调用绘图脚本
cat("\n>>> 调用绘图脚本生成图形...\n")
source(file.path(BASE_DIR, "scripts", "analysis", "module2_plots.R"))

# ========== 15. 完成 ==========
cat(sprintf("\n========================================\n"))
cat(sprintf("  模块二分析完成\n"))
cat(sprintf("  结果目录: %s\n", OUTPUT_DIR))
cat(sprintf("  发现中心: %s\n", discovery_center))
cat("  输出文件:\n")
cat("    图片 (fig/):\n")
cat("    - ROC_*.png / PRC_*.png: 曲线图\n")
cat("    - confusion_matrix_*.png: 混淆矩阵\n")
cat("    - shap_importance.png: 重要性柱状图\n")
cat("    数据 (csv/):\n")
cat("    - binary_performance_table.csv: 二分类性能表\n")
cat("    - binary_predictions.csv: 二分类预测值\n")
cat("    - binary_grid_search_params.csv: 二分类网格搜索参数表\n")
cat("    - multiclass_performance_table.csv: 多分类性能表\n")
cat("    - multiclass_predictions.csv: 多分类预测值\n")
cat("    - multiclass_grid_search_params.csv: 多分类网格搜索参数表\n")
cat("    - shap_values.csv: 各样本 SHAP 值\n")
cat("    - shap_feature_values.csv: 原始特征值（用于 SHAP beeswarm 图）\n")
cat("    - shap_importance.csv: 全局特征重要性\n")
cat("    模型:\n")
cat(sprintf("    - model/best_multiclass_model.rds: 最佳多分类模型 (%s)\n", best_m_name))
cat("    报告:\n")
cat("    - module2_report.txt: 详细分析报告\n")
cat("    其他:\n")
cat("    - model/README.md: 模型使用说明\n")
cat(sprintf("========================================\n"))