#!/usr/bin/env Rscript
#===========================================================
# Figure 3 A-D 单独导出（无参考虚线）
#
# 仅生成 Qingdao/Xiamen 的 ROC 与 PRC 四个面板：
#   A ROC-Qingdao, B PRC-Qingdao, C ROC-Xiamen, D PRC-Xiamen
# 输出到 results/module2_main/pdf/
#===========================================================
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))), "plot_functions.R"))

suppressPackageStartupMessages({
  library(ggpubr)
  library(showtext)
})

# 与 figure_publication.R 一致的 Arial / showtext 设置
.font_files <- list(
  regular = path.expand("~/.fonts/arial.ttf"),
  bold    = path.expand("~/.fonts/arialbd.ttf"),
  italic  = path.expand("~/.fonts/ariali.ttf"),
  bolditalic = path.expand("~/.fonts/arialbi.ttf"))
if (all(file.exists(unlist(.font_files)))) {
  do.call(sysfonts::font_add, c(list("Arial"), .font_files))
}
showtext_auto(enable = TRUE)
theme_set(theme_publication(base_size = PUB_TEXT_PT))

BASE_DIR <- get_project_root()
M2C <- file.path(BASE_DIR, "results", "module2_main", "csv")
OUT_DIR <- file.path(BASE_DIR, "results", "module2_main", "pdf")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

pred <- read.csv(file.path(M2C, "multiclass_predictions.csv"), stringsAsFactors = FALSE)

# 出版图中 Figure3 A-D 的基础字号；单独导出时可读性优先
BASE_PT <- 8
PANEL_W <- 90
PANEL_H <- 82

pA <- plot_roc(pred, test_set = "Qingdao", ci = TRUE, n_boot = 200,
               show_reference = FALSE, base_size = BASE_PT,
               legend.position = c(0.62, 0.18))
pB <- plot_prc(pred, test_set = "Qingdao", ci = TRUE, n_boot = 200,
               show_reference = FALSE, base_size = BASE_PT,
               legend.position = c(0.40, 0.18))
pC <- plot_roc(pred, test_set = "Xiamen", ci = TRUE, n_boot = 200,
               show_reference = FALSE, base_size = BASE_PT,
               legend.position = c(0.62, 0.18))
pD <- plot_prc(pred, test_set = "Xiamen", ci = TRUE, n_boot = 200,
               show_reference = FALSE, base_size = BASE_PT,
               legend.position = c(0.40, 0.18))

save_panel <- function(p, id) {
  ggsave(file.path(OUT_DIR, paste0(id, ".pdf")), p,
         width = PANEL_W, height = PANEL_H, units = "mm",
         device = cairo_pdf, bg = "white")
  cat(sprintf("  -> %s.pdf\n", id))
}

save_panel(pA, "Figure3A_ROC_Qingdao")
save_panel(pB, "Figure3B_PRC_Qingdao")
save_panel(pC, "Figure3C_ROC_Xiamen")
save_panel(pD, "Figure3D_PRC_Xiamen")

# 2×2 合并版，便于直接查看四个面板
combined <- ggarrange(pA, pB, pC, pD,
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow = 2,
                      font.label = list(size = 10, face = "bold"))
ggsave(file.path(OUT_DIR, "Figure3_ABCD_combined.pdf"), combined,
       width = 175, height = 150, units = "mm",
       device = cairo_pdf, bg = "white")
cat("  -> Figure3_ABCD_combined.pdf\n")
cat(sprintf("输出目录: %s\n", OUT_DIR))
