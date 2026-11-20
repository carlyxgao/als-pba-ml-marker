#!/usr/bin/env Rscript
#===========================================================
# ALS-PBA-ML-Marker 主分析脚本（完整流程入口）
#
# 按顺序串联 bulk 蛋白组分析（论文 Fig 1-3 及相关补充图）：
#   数据整理 → 模块1(DE/GO/Reactome/敏感性/作图) → 模块2(ML/作图)
#   → 手稿专用 PDF（多分类 ROC/PRC、12-Marker RF vs NfL）
#
# 用法:
#   Rscript scripts/run_all.R                  # 跑完整流程
#   Rscript scripts/run_all.R --list           # 列出所有步骤
#   Rscript scripts/run_all.R --only module2_ml  # 只跑指定步骤
#   Rscript scripts/run_all.R --from module1_de  # 从指定步骤开始
#
# 说明:
#   * 每步用独立 Rscript 子进程运行（脚本用 --tag 参数区分结果目录）
#   * 任一步失败即中止（exit != 0），并打印失败的步骤名
#   * 数据准备: 将原始 PBA 报告 zip 与基线 xlsx、nfl.xlsx 放入 data/raw/
#     （详见 README），curation 步骤会解压并合并样本信息
#===========================================================

# ========== 配置（按需修改） ==========
TAG            <- "main"    # 结果目录后缀：results/module*_<TAG>

# ========== 路径 ==========
BASE_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) {
    normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", f[1]))), ".."))
  } else normalizePath(getwd())
})
RSCRIPT_BIN <- file.path(R.home("bin"), "Rscript")
MARKER_FILE <- file.path(BASE_DIR, "results", paste0("module1_", TAG), "csv", "core_markers__Jinan.csv")

# ========== 步骤定义（相对 scripts/ 的路径 + 参数） ==========
steps <- list(
  curation      = list("curation/merge_sample_info.R",
                       character()),
  module1_de     = list("analysis/module1_de_analysis.R",
                        c("--tag", TAG)),
  module1_go     = list("analysis/module1_pathway_enrichment.R",
                        c("--tag", TAG)),
  module1_react  = list("analysis/module1_reactome_enrichment.R",
                        c("--tag", TAG)),
  module1_sens   = list("analysis/module1_sensitivity_analysis.R",
                        c("--tag", TAG)),
  module1_plots  = list("analysis/module1_plots.R",
                        c("--tag", TAG)),
  module2_ml     = list("analysis/module2_ml_modeling.R",
                        c("--tag", TAG,
                          "--markers", MARKER_FILE)),
  module2_plots  = list("analysis/module2_plots.R",
                        c("--tag", TAG)),
  figure_roc_prc = list("analysis/figure3_abcd_pdf.R",
                        character()),
  figure_nfl_roc = list("analysis/nfl_roc_pdf.R",
                        character())
)

# ========== 运行助手 ==========
run_step <- function(name) {
  s <- steps[[name]]
  script <- file.path(BASE_DIR, "scripts", s[[1]])
  args <- s[[2]]
  cat(sprintf("\n%s\n>>> [%s] %s %s\n%s\n", strrep("=", 72), name,
              basename(script), paste(args, collapse = " "), strrep("=", 72)))
  if (!file.exists(script)) stop(sprintf("脚本不存在: %s", script))
  t0 <- Sys.time()
  # 模块二的 XGBoost 在多线程 BLAS 下可能长时间无输出；单线程稳定且更快
  env_vars <- if (name == "module2_ml")
    c("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1") else character(0)
  status <- system2(RSCRIPT_BIN, c(shQuote(script), args), env = env_vars)
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  if (status != 0) stop(sprintf("✗ [%s] 失败 (exit %d, %.1f min)", name, status, mins))
  cat(sprintf("✓ [%s] 完成 (%.1f min)\n", name, mins))
}

# ========== CLI 参数 ==========
cli <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag) {
  i <- match(flag, cli)
  if (!is.na(i) && i < length(cli)) cli[i + 1] else NULL
}
if ("--list" %in% cli) {
  cat("可用步骤:\n"); for (n in names(steps)) cat(sprintf("  %-18s %s\n", n, steps[[n]][[1]]))
  quit(status = 0)
}

todo <- names(steps)
only <- get_opt("--only")
from <- get_opt("--from")
if (!is.null(only)) {
  if (!only %in% names(steps)) stop(sprintf("未知步骤: %s（用 --list 查看）", only))
  todo <- only
} else if (!is.null(from)) {
  if (!from %in% names(steps)) stop(sprintf("未知步骤: %s（用 --list 查看）", from))
  todo <- names(steps)[seq(match(from, names(steps)), length(steps))]
}

cat(sprintf(">>> ALS-PBA-ML-Marker 主分析流程：%d 个步骤 (TAG=%s)\n",
            length(todo), TAG))
t_start <- Sys.time()
for (nm in todo) run_step(nm)
cat(sprintf("\n%s\n>>> 全部完成（总耗时 %.1f 分钟）\n%s\n",
            strrep("=", 72),
            as.numeric(difftime(Sys.time(), t_start, units = "mins")),
            strrep("=", 72)))
