#!/usr/bin/env Rscript
#===========================================================
# NfL vs 12-Marker RF ROC 对比 PDF
#
# 仅使用同时具有 NfL 和 PBA 数据的 ALS / DC 样本，按 Sample 配对；
# 输出 Qingdao、Xiamen 两个单独 PDF 及合并 PDF 到
# results/module2_main/pdf/
#===========================================================
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))), "plot_functions.R"))

suppressPackageStartupMessages({
  library(ggpubr)
  library(showtext)
})

.font_files <- list(
  regular = path.expand("~/.fonts/arial.ttf"),
  bold    = path.expand("~/.fonts/arialbd.ttf"),
  italic  = path.expand("~/.fonts/ariali.ttf"),
  bolditalic = path.expand("~/.fonts/arialbi.ttf"))
if (all(file.exists(unlist(.font_files)))) {
  do.call(sysfonts::font_add, c(list("Arial"), .font_files))
}
showtext_auto(enable = TRUE)

BASE_DIR <- get_project_root()
M2C <- file.path(BASE_DIR, "results", "module2_main", "csv")
OUT_DIR <- file.path(BASE_DIR, "results", "module2_main", "pdf")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

multi_pred <- read.csv(file.path(M2C, "multiclass_predictions.csv"), stringsAsFactors = FALSE)
nfl_auc <- read.csv(file.path(M2C, "nfl_vs_marker_auc.csv"), stringsAsFactors = FALSE)
nfl_vals <- assemble_nfl_values(file.path(BASE_DIR, "data", "raw", "nfl.xlsx"), multi_pred)
if (is.null(nfl_vals)) stop("NfL 数据装配失败")

rf <- multi_pred %>%
  filter(Model == "RF") %>%
  select(Sample, TestSet, TrueLabel, ALS)
paired <- nfl_vals %>%
  select(Sample, TestSet, TrueLabel, NFL) %>%
  inner_join(rf, by = c("Sample", "TestSet", "TrueLabel")) %>%
  filter(TrueLabel %in% c("ALS", "DC"))

make_panel <- function(ts) {
  d <- paired %>% filter(TestSet == ts)
  bin_d <- d %>%
    transmute(Sample, Model = "RF", TestSet, TrueLabel, Prob_ALS = ALS)
  nfl_d <- d %>% select(Sample, TestSet, TrueLabel, NFL)
  nfl_row <- nfl_auc %>% filter(TestSet == ts, Model == "NFL") %>% slice(1)
  p <- plot_nfl_roc(bin_d, nfl_d, test_set = ts, nfl_auc_row = nfl_row,
                    base_size = 11)
  p + labs(title = sprintf("%s (n = %d)", ts, nrow(d))) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          legend.position = c(0.72, 0.28))
}

p_q <- make_panel("Qingdao")
p_x <- make_panel("Xiamen")

ggsave(file.path(OUT_DIR, "ROC_NFL_vs_Marker_Qingdao.pdf"), p_q,
       width = 90, height = 85, units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(OUT_DIR, "ROC_NFL_vs_Marker_Xiamen.pdf"), p_x,
       width = 90, height = 85, units = "mm", device = cairo_pdf, bg = "white")

combined <- ggarrange(p_q, p_x, labels = c("A", "B"), ncol = 2, nrow = 1,
                      font.label = list(size = 11, face = "bold"))
ggsave(file.path(OUT_DIR, "ROC_NFL_vs_Marker_combined.pdf"), combined,
       width = 175, height = 85, units = "mm", device = cairo_pdf, bg = "white")

cat("  -> ROC_NFL_vs_Marker_Qingdao.pdf\n")
cat("  -> ROC_NFL_vs_Marker_Xiamen.pdf\n")
cat("  -> ROC_NFL_vs_Marker_combined.pdf\n")
