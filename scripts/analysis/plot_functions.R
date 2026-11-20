#!/usr/bin/env Rscript
#===========================================================
# 共享绘图函数库 (Shared Plot Function Library)
#
# 每个图类型一个函数，所有函数返回 ggplot/grob 对象，
# 内部不调用 ggsave、不写文件 —— 由调用方决定保存方式。
#
# 用途:
#   * module*_plots.R  —— source 本文件后调用函数，各自 ggsave 输出单张图
#   * figure_publication.R —— source 本文件，调用函数生成面板后用 ggarrange 排版
#
# 约定:
#   * 分组顺序 ALS → AD → DC → HC，配色 ALL_GROUP_COLORS
#   * 用 resolve_group_theme() 动态截取数据中实际存在的分组
#   * Seurat / enrichplot 等重包在依赖它们的函数内部懒加载
#===========================================================

# ========== 1. 包加载 ==========
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(ggrepel)
  library(ggsignif)
  library(scales)
  library(RColorBrewer)
  library(pheatmap)
  library(pROC)
  library(grid)
  library(ggplotify)
  library(showtext)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ========== 2. 路径/读取助手 ==========
# 自动推导项目根目录（复用各脚本的 --file= 模式）
get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg[1]))), "..", "..")))
  }
  getwd()
}

# 解析 --tag 参数
parse_cli_tag <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  i <- 1; while (i <= length(args)) {
    if (args[i] == "--tag" && i < length(args)) return(args[i + 1])
    i <- i + 1
  }
  NULL
}

# 带 memo 缓存的读取（figure 脚本加载 250MB Seurat 对象时只读一次）
.read_cache <- new.env(parent = emptyenv())
read_csv_cached <- function(path) {
  key <- paste0("csv:", path)
  if (exists(key, envir = .read_cache)) return(get(key, envir = .read_cache))
  v <- read.csv(path, stringsAsFactors = FALSE)
  assign(key, v, envir = .read_cache); v
}
read_rds_cached <- function(path) {
  key <- paste0("rds:", path)
  if (exists(key, envir = .read_cache)) return(get(key, envir = .read_cache))
  v <- readRDS(path)
  assign(key, v, envir = .read_cache); v
}

# ========== 3. 常量（集中各脚本复制的配色）==========
ALL_GROUP_COLORS  <- c("ALS" = "#E64B35", "AD" = "#00A087", "DC" = "#3C5488", "HC" = "#4DBBD5")
ALL_GROUP_LEVELS  <- c("ALS", "AD", "DC", "HC")
CENTER_COLORS     <- c("Jinan" = "#F39B7F", "Qingdao" = "#8491B4", "Xiamen" = "#91D1C2")

CT_COLORS <- c(
  "Neuron" = "#D62728", "Astrocyte" = "#FF7F0E", "Microglia" = "#2CA02C",
  "Oligodendrocyte" = "#17BECF", "BBB_Endothelial" = "#9467BD",
  "B_cell" = "#E377C2", "T_cell" = "#8C564B", "NK_cell" = "#BCBD22",
  "Neutrophil" = "#1F77B4", "Macrophage" = "#7F7F7F",
  "Platelet" = "#393B79", "Hepatocyte" = "#F7B6D2",
  "Myeloid_progenitor" = "#9EDAE5"
)
ORIGIN_GROUP_COLORS <- c("CNS" = "#E64B35", "Immune" = "#3C5488", "Other" = "#4DBBD5",
                         "Unannotated" = "grey75")   # 0 命中 cluster 不参与来源分组
# 蛋白驱动的 EV 功能来源配色（Phase 0 可行性评估 + 免疫 marker 扩充）
SOURCE_COLORS <- c(
  "Complement" = "#E64B35", "Platelet" = "#4DBBD5", "Endothelial" = "#00A087",
  "Neutrophil" = "#3C5488", "Monocyte" = "#FF7F0E", "T_cell" = "#8C564B",
  "NK_cell" = "#BCBD22", "B_cell" = "#E377C2", "Hepatic" = "#F39B7F",
  "Unannotated" = "grey75"   # 灰色在 10-80% 范围内
)
TISSUE_COLORS <- c(
  "CNS" = "#E64B35", "Immune" = "#3C5488", "Endothelial" = "#00A087",
  "Epithelial" = "#FF7F0E", "Stromal" = "#9467BD", "Blood_platelet" = "#E377C2",
  "Other" = "grey50", "Unannotated" = "grey80"   # 灰在 10-80% 之间且相隔 ≥20%
)

MODEL_COLORS <- c(LR = "#E64B35", RF = "#4DBBD5", SVM = "#00A087", XGB = "#3C5488")
FS_COLORS <- c("Markers" = "#E64B35", "Metadata" = "#4DBBD5", "MarkersPlusMetadata" = "#00A087")
NFL_BAR_COLORS <- c("NFL" = "#A6B1C2", "12-Marker RF" = "#E8846F")   # 柔和灰蓝/珊瑚色

# Nature 期刊风格亚群配色（自动扩展）
BASE_PALETTE <- c(
  "#D62728", "#1F77B4", "#2CA02C", "#9467BD", "#FF7F0E",
  "#17BECF", "#E377C2", "#8C564B", "#BCBD22", "#7F7F7F",
  "#393B79", "#637939", "#9EDAE5", "#C49C94", "#F7B6D2"
)
extend_palette <- function(cluster_ids) {
  n <- length(cluster_ids)
  colors <- if (n <= length(BASE_PALETTE)) BASE_PALETTE[1:n] else colorRampPalette(BASE_PALETTE)(n)
  stats::setNames(colors, as.character(cluster_ids))
}

# 动态截取数据中实际存在的分组/配色
resolve_group_theme <- function(x, levels = ALL_GROUP_LEVELS, colors = ALL_GROUP_COLORS) {
  present <- intersect(levels, as.character(x))
  present <- present[!is.na(present)]
  list(levels = present, colors = colors[present])
}

# ========== 4.5 出版级设计系统（pt / mm / DPI 对应关系） ==========
# 核心换算：1 in = 25.4 mm = 72 pt  →  1 pt = 0.3528 mm；1 mm = 2.8346 pt
#   * theme 的 element_text(size=) 单位是「磅(pt)」——即最终打印字号，直接写 6-8 即可
#   * geom 的 size / linewidth 单位是「毫米(mm)」——转磅需 ×0.3528：
#       geom_text(size = pt2mm(7))  ≈ 7pt 文字；geom_line(linewidth = pt2mm(1)) ≈ 1pt 线宽
#   * DPI 只影响位图密度，不影响打印字号：像素 = 物理尺寸(mm) × dpi / 25.4
#   * 设计原则：直接在最终打印尺寸(175mm)下设计；字号用 pt，线宽用 pt2mm()，位图按 dpi 导出
PUB_WIDTH_MM  <- 175      # 双栏期刊标准宽度
PUB_MAX_H_MM  <- 200      # 高度上限
PUB_DPI       <- 1000     # 位图分辨率（发布级）
PUB_FONT      <- "Arial"  # 需已安装到 ~/.fonts 并注册到 sysfonts
PUB_TEXT_PT   <- 5        # 正文基准字号（用户调整为 5pt）
PUB_MIN_PT    <- 4
PUB_MAX_PT    <- 6
PUB_LINE_PT   <- 1        # 线宽 0.5-1.5pt 取中值

# pt → mm（geom_text/geom_point size、geom_line/geom_path linewidth 用）
pt2mm <- function(pt) pt * 0.3528
# 把尺寸钳制到出版字号范围（避免 base_size 小时 `base+4` 等偏移越界）
clamp_pt <- function(x) pmin(PUB_MAX_PT, pmax(PUB_MIN_PT, x))

# 发布用：给单个 ggplot 面板统一应用 Arial + 把文字字号钳到 [6,8]pt
# 在 compose 阶段对每个面板调用即可，无需逐个函数传 family
apply_pub <- function(p, family = PUB_FONT, min_pt = PUB_MIN_PT, max_pt = PUB_MAX_PT) {
  if (!inherits(p, "ggplot")) return(p)
  th <- p$theme
  els <- c("text", "plot.title", "plot.subtitle", "axis.text", "axis.text.x", "axis.text.y",
           "axis.title", "axis.title.x", "axis.title.y",
           "legend.text", "legend.title", "strip.text")
  for (el in els) {
    cur <- th[[el]]
    if (inherits(cur, "element_text")) {
      sz <- cur$size
      if (is.numeric(sz)) sz <- min(max_pt, max(min_pt, sz))
      th[[el]] <- element_text(family = family, size = sz, face = cur$face,
                               colour = cur$colour, hjust = cur$hjust, vjust = cur$vjust,
                               angle = cur$angle, lineheight = cur$lineheight,
                               margin = cur$margin, debug = cur$debug)
    } else if (is.null(cur) && el == "text") {
      # 只有基础 text 为 NULL 时才创建；其余 NULL 元素保持继承
      # （继承自 element_blank 的如 venn/heatmap 的 axis 不会因此被重新点亮）
      th[["text"]] <- element_text(family = family)
    }
  }
  p$theme <- th
  p
}

# ========== 4. 主题 ==========
# UMAP/t-SNE 无坐标轴主题（参数化合并 module3 三个变体）
theme_clean <- function(title_size = 18, legend = c("none", "right"),
                        facet = FALSE, subtitle = FALSE,
                        key_size = 0.9, legend_spacing = 0.3, ...) {
  legend <- match.arg(legend)
  th <- theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = title_size),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = if (legend == "right") "right" else "none",
      ...)
  if (facet) {
    th <- th + theme(strip.text = element_text(size = title_size * 0.75, face = "bold"),
                     strip.background = element_blank())
  }
  if (legend == "right") {
    th <- th + theme(legend.title = element_blank(),
                     legend.text = element_text(size = max(10, title_size * 0.75)),
                     legend.key.size = unit(key_size, "cm"),
                     legend.key.width = unit(key_size, "cm"),
                     legend.spacing.y = unit(legend_spacing, "cm"))
  }
  if (subtitle) {
    th <- th + theme(plot.subtitle = element_text(hjust = 0.5, size = title_size * 0.7, color = "grey40"))
  }
  th
}

# 发表级面板主题（期刊字号；字体经 apply_pub 统一设 Arial）
theme_publication <- function(base_size = PUB_TEXT_PT) {
  theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = clamp_pt(base_size + 1)),
          plot.subtitle = element_text(hjust = 0.5, size = clamp_pt(base_size - 1), colour = "grey40"),
          axis.text = element_text(size = clamp_pt(base_size - 1)),
          axis.title = element_text(size = base_size),
          legend.text = element_text(size = clamp_pt(base_size - 1)),
          legend.title = element_text(size = base_size),
          strip.text = element_text(size = clamp_pt(base_size - 1), face = "bold"),
          legend.key.size = unit(0.35, "cm"),
          plot.margin = ggplot2::margin(3, 3, 3, 3))
}

# ========== 5. 非 ggplot 归一化 ==========
# pheatmap 对象 → ggplot（供 ggarrange 排版）
as_ggplot_pheatmap <- function(p) {
  ggplotify::as.ggplot(p)
}

# 动态富集气泡图尺寸（沿用 module1_pathway_enrichment.R 的宽度计算）
enrichment_panel_size <- function(res, n_show) {
  n_char <- max(nchar(res@result$Description[seq_len(n_show)]))
  list(w = max(10, 5 + n_char * 0.08),
       h = max(7, n_show * 0.42))
}

# ========== 6. Module 1 函数 ==========

# 火山图（limma 结果；logFC 翻转：正数 = ALS 上调）
# 图例默认置顶，key 显示圆点；虚线阈值线宽符合 0.5-1.5pt 规范
plot_volcano <- function(res, title = NULL, top_n = 5, fc_thr = 1.5, p_thr = 0.05,
                         base_size = 18, colors = c(Up = "#E64B35", Down = "#4DBBD5", NS = "grey60"),
                         legend.position = "top", linewidth = pt2mm(1), point_size = 0.8,
                         label_size = 6, ...) {
  df <- res %>%
    mutate(plot_logFC = -limma_logFC) %>%
    mutate(Direction = case_when(
      limma_sig & plot_logFC > 0 ~ "Up",
      limma_sig & plot_logFC < 0 ~ "Down",
      TRUE ~ "NS")) %>%
    mutate(Direction = factor(Direction, levels = c("Up", "Down", "NS")))
  p <- ggplot(df, aes(x = plot_logFC, y = -log10(limma_P.Value), color = Direction)) +
    geom_point(alpha = 0.5, size = point_size) +
    scale_color_manual(values = colors, name = NULL) +
    geom_hline(yintercept = -log10(p_thr), linetype = "dashed", alpha = 0.4, linewidth = linewidth) +
    geom_vline(xintercept = c(-log2(fc_thr), log2(fc_thr)), linetype = "dashed", alpha = 0.3, linewidth = linewidth) +
    labs(title = title %||% "Volcano Plot",
         x = "log2 Fold Change", y = "-log10(P-value)") +
    theme_minimal(base_size = base_size) +
    theme(legend.position = legend.position,
          legend.title = element_blank(),
          legend.text = element_text(size = base_size - 4),
          plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 2)) +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 2, alpha = 1)))
  up_sig <- df %>% filter(Direction == "Up") %>% arrange(limma_P.Value) %>% head(top_n)
  down_sig <- df %>% filter(Direction == "Down") %>% arrange(limma_P.Value) %>% head(top_n)
  label_data <- bind_rows(up_sig, down_sig)
  if (nrow(label_data) > 0) {
    # 不画连接线（segment.color=NA），小 padding 让标签贴近对应点
    p <- p + geom_text_repel(data = label_data, aes(label = protein),
                             max.overlaps = 30, size = pt2mm(label_size),
                             segment.color = NA,
                             box.padding = 0.3, point.padding = 0.2)
  }
  p
}

# 三集合 Venn 图（DE 共识 ∩ LASSO 并集 ∩ RF 并集，返回 ggplot）
plot_venn <- function(venn_data, set_labels = NULL, set_colors = NULL, base_size = 18,
                      label_size = 6, detail_size = 5, ...) {
  s1 <- venn_data$de_consensus
  s2 <- venn_data$lasso_union
  s3 <- venn_data$rf_union
  core <- venn_data$core_set
  set_sizes <- c(length(s1), length(s2), length(s3))
  if (is.null(set_labels)) {
    set_labels <- c("DE Consensus", "LASSO", "Random Forest")
  }
  set_details <- c(
    sprintf("%d candidate proteins\n(≥1 comparison)", set_sizes[1]),
    sprintf("%d candidate proteins\n(union of 2 comparisons)", set_sizes[2]),
    sprintf("%d candidate proteins\n(Top30 by MDA, union)", set_sizes[3])
  )
  if (is.null(set_colors)) set_colors <- c("#E64B35", "#4DBBD5", "#00A087")

  cx <- 5; cy <- 4.8; r <- 2.6
  angles_deg <- c(90, 210, 330)
  centers <- t(sapply(angles_deg, function(a) c(cx + r * 0.65 * cos(a * pi / 180),
                                                 cy + r * 0.65 * sin(a * pi / 180))))

  circ_pts <- function(cx, cy, r, n = 200) {
    th <- seq(0, 2 * pi, length.out = n)
    data.frame(x = cx + r * cos(th), y = cy + r * sin(th))
  }
  circles <- do.call(rbind, lapply(seq_len(3), function(i) {
    circ_pts(centers[i, 1], centers[i, 2], r) %>%
      mutate(group = i, fill_col = adjustcolor(set_colors[i], alpha.f = 0.12))
  }))

  text_positions <- list(
    c(centers[1, 1], centers[1, 2] + r * 0.45),
    c(centers[2, 1] - r * 0.42, centers[2, 2] - r * 0.20),
    c(centers[3, 1] + r * 0.42, centers[3, 2] - r * 0.20)
  )

  p <- ggplot() +
    geom_polygon(data = circles, aes(x = x, y = y, group = group, fill = fill_col),
                 color = NA) +
    scale_fill_identity() +
    coord_equal(xlim = c(0, 10), ylim = c(0, 10)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          # 显式删除坐标轴文字/标题/刻度
          axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
  for (i in seq_len(3)) {
    p <- p +
      annotate("text", x = text_positions[[i]][1], y = text_positions[[i]][2] + 0.30,
               label = set_labels[i], size = pt2mm(label_size), fontface = "bold", hjust = 0.5) +
      annotate("text", x = text_positions[[i]][1], y = text_positions[[i]][2] - 0.55,
               label = set_details[i], size = pt2mm(detail_size), hjust = 0.5)
  }
  p + annotate("text", x = cx, y = cy + 0.15,
               label = sprintf("%d proteins", length(core)),
               size = base_size * 0.6, fontface = "bold", color = "#E41A1C", hjust = 0.5)
}

# ============================================================
# ggplot 版热图（参照 plots.R 的 plotHeatmap 逻辑重写）
#  geom_raster 画色块 + geom_segment 画聚类树（条带厚度可控，不占过多宽高）
#  + geom_raster 画行列分组注释条 + 连续色条图例
#  替代 pheatmap，聚类树条带大小可调；返回 ggplot 可被 ggarrange/apply_pub 组合
# ============================================================
plot_heatmap_ggplot <- function(mat,
                                scale = "row",
                                cluster_rows = TRUE, cluster_cols = FALSE,
                                clustering_method = "ward.D2",
                                dendro_row_w = 2, dendro_col_w = 2,   # 聚类树条带厚度(相对 tile 数)
                                row_anno = NULL, col_anno = NULL,     # 分组注释 data.frame（行为行/列名，列为分组）
                                row_anno_colors = NULL, col_anno_colors = NULL,  # 命名颜色向量
                                title = NULL,
                                colors = colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100),
                                show_rownames = TRUE, show_colnames = FALSE,
                                fontsize = 7, value_name = "Scaled expression",
                                color_limits = NULL, ...) {
  if (is.null(dimnames(mat))) dimnames(mat) <- list(paste0("R", seq_len(nrow(mat))),
                                                    paste0("C", seq_len(ncol(mat))))
  if (scale == "row" && nrow(mat) > 1) mat <- t(scale(t(mat)))

  row_clust <- if (cluster_rows && nrow(mat) > 1) hclust(stats::dist(mat), method = clustering_method) else NULL
  col_clust <- if (cluster_cols && ncol(mat) > 1) hclust(stats::dist(t(mat)), method = clustering_method) else NULL
  row_order <- if (!is.null(row_clust)) row_clust$order else seq_len(nrow(mat))
  col_order <- if (!is.null(col_clust)) col_clust$order else seq_len(ncol(mat))

  # 主色块长表（按聚类顺序）
  # 注意：expand.grid 返回 factor 列，cbind(factor,factor) 会变成整数层级码
  #       导致 mat[矩阵] 按位置而非按名称取值 → 必须用 match() 转数值下标
  df <- expand.grid(row = rownames(mat)[row_order], col = colnames(mat)[col_order],
                    KEEP.OUT.ATTRS = FALSE) %>%
    mutate(row_pos = match(row, rownames(mat)[row_order]),
           col_pos = match(col, colnames(mat)[col_order])) %>%
    mutate(val = mat[cbind(match(row, rownames(mat)), match(col, colnames(mat)))])

  rw <- if (!is.null(row_clust)) dendro_row_w else 0
  cw <- if (!is.null(col_clust)) dendro_col_w else 0
  raw <- if (!is.null(row_anno)) ncol(row_anno) else 0
  caw <- if (!is.null(col_anno)) ncol(col_anno) else 0
  total_anno <- raw + caw   # 注释列总数；legend 顺序：Center→...→Group→(row_anno)→Scaled expression

  p <- ggplot(df, aes(x = col_pos, y = row_pos, fill = val)) +
    geom_raster() +
    scale_fill_gradientn(colors = colors, name = value_name,
                         limits = color_limits, oob = scales::squish,
                         guide = guide_colorbar(barheight = unit(1, "cm"), barwidth = unit(0.25, "cm"),
                                                order = total_anno + 1)) +
    scale_y_reverse(breaks = seq_len(nrow(mat)),
                    labels = if (show_rownames) rownames(mat)[row_order] else NULL,
                    expand = c(0, 0), position = "right") +
    scale_x_continuous(breaks = seq_len(ncol(mat)),
                       labels = if (show_colnames) colnames(mat)[col_order] else NULL,
                       expand = c(0, 0)) +
    theme_void() +
    theme(axis.text.y = element_text(size = fontsize - 2, hjust = 0, margin = ggplot2::margin(0, 2, 0, 2, "pt")),
          axis.text.x = element_text(size = fontsize, angle = 90, hjust = 1, vjust = 0.5,
                                     margin = ggplot2::margin(2, 0, 2, 0, "pt")),
          plot.margin = ggplot2::margin(4, 4, 4, 4, "pt"))

  # 行聚类树（左侧条带；dendro 高度映射到 x 负方向）
  if (!is.null(row_clust)) {
    seg <- ggdendro::dendro_data(as.dendrogram(row_clust))$segments
    hmax <- max(c(seg$y, seg$yend)); if (hmax <= 0) hmax <- 1
    seg <- seg %>%
      mutate(y = y / hmax * rw, yend = yend / hmax * rw) %>%
      mutate(y = -y + 0.5 - raw, yend = -yend + 0.5 - raw)
    p <- p + geom_segment(data = seg, aes(x = y, y = x, xend = yend, yend = xend),
                          color = "grey30", linewidth = pt2mm(0.5), inherit.aes = FALSE)
  }
  # 列聚类树（上方条带；整体上移，避免与最后一个注释条如 Center 重叠遮挡）
  if (!is.null(col_clust)) {
    seg <- ggdendro::dendro_data(as.dendrogram(col_clust))$segments
    hmax <- max(c(seg$y, seg$yend)); if (hmax <= 0) hmax <- 1
    seg <- seg %>%
      mutate(y = y / hmax * cw, yend = yend / hmax * cw) %>%
      mutate(y = -y - caw, yend = -yend - caw)
    p <- p + geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend),
                          color = "grey30", linewidth = pt2mm(0.5), inherit.aes = FALSE)
  }
  # 注释条图例（用 ggnewscale 建立独立 fill 尺度，每个注释列一个图例）
  p <- p + ggnewscale::new_scale_fill()
  # 行注释条（左侧，颜色映射到 group 生成图例；legend 顺序排在列注释之后）
  if (!is.null(row_anno)) {
    for (k in seq_len(ncol(row_anno))) {
      grp <- row_anno[[k]]
      cols <- if (!is.null(row_anno_colors)) row_anno_colors[[names(row_anno)[k]]] else
        extend_palette(sort(unique(as.character(grp))))
      ctb <- tibble(row = rownames(row_anno), group = factor(as.character(grp), levels = names(cols))) %>%
        mutate(id = match(row, rownames(mat)[row_order])) %>%
        mutate(x = -k + 0.5)
      p <- p + geom_raster(data = ctb, aes(x = x, y = id, fill = group), inherit.aes = FALSE) +
        scale_fill_manual(name = names(row_anno)[k], values = cols,
                          guide = guide_legend(order = caw + k,
                                               keyheight = unit(0.22, "cm"), keywidth = unit(0.22, "cm")))
      if (k < ncol(row_anno)) p <- p + ggnewscale::new_scale_fill()
    }
  }
  # 列注释条（上方；legend 顺序倒置，使最后一个注释列如 Center 在最上）
  if (!is.null(col_anno)) {
    for (k in seq_len(ncol(col_anno))) {
      grp <- col_anno[[k]]
      cols <- if (!is.null(col_anno_colors)) col_anno_colors[[names(col_anno)[k]]] else
        extend_palette(sort(unique(as.character(grp))))
      ctb <- tibble(col = rownames(col_anno), group = factor(as.character(grp), levels = names(cols))) %>%
        mutate(id = match(col, colnames(mat)[col_order])) %>%
        mutate(y = -k + 0.5)
      p <- p + geom_raster(data = ctb, aes(x = id, y = y, fill = group), inherit.aes = FALSE) +
        scale_fill_manual(name = names(col_anno)[k], values = cols,
                          guide = guide_legend(order = caw - k + 1,
                                               keyheight = unit(0.22, "cm"), keywidth = unit(0.22, "cm")))
      if (k < ncol(col_anno)) p <- p + ggnewscale::new_scale_fill()
    }
  }
  if (!is.null(title)) p <- p + ggtitle(title)
  p
}

# 核心标志物热图（pheatmap → ggplot；支持 Group 或 Group+Center 注释）
plot_core_heatmap <- function(mat, anno_col = NULL, anno_colors = NULL, title = NULL,
                              scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
                              clustering_method = "ward.D2", fontsize = 7,
                              colors = colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100),
                              show_colnames = FALSE, show_rownames = TRUE, blank_axes = FALSE,
                              dendro_col_w = 2, dendro_row_w = 2, color_limits = c(-3, 3), ...) {
  # 用 ggplot 版热图替代 pheatmap（聚类树条带可控，不占过多宽高）
  g <- plot_heatmap_ggplot(mat, scale = scale,
                           cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                           clustering_method = clustering_method,
                           row_anno = NULL, col_anno = anno_col,
                           row_anno_colors = NULL, col_anno_colors = anno_colors,
                           title = title, colors = colors,
                           show_rownames = show_rownames, show_colnames = show_colnames,
                           fontsize = fontsize, dendro_col_w = dendro_col_w, dendro_row_w = dendro_row_w,
                           color_limits = color_limits)
  if (blank_axes) {
    g <- g + theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
  }
  g
}

# GO/KEGG/Reactome 富集气泡图（enrichplot::dotplot，返回 ggplot）
plot_go_dotplot <- function(res, ont = NULL, comp_label = "", n_show = 20, title = NULL, ...) {
  if (!requireNamespace("enrichplot", quietly = TRUE)) stop("enrichplot not installed")
  library(enrichplot)
  if (is.null(title)) {
    title <- if (is.null(ont)) paste0("KEGG - ", comp_label) else paste0("GO ", ont, " - ", comp_label)
  }
  p <- dotplot(res, showCategory = n_show, title = title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.margin = ggplot2::margin(10, 10, 10, 60),
          axis.text.y = element_text(size = 11, lineheight = 0.85),
          panel.spacing = unit(1.5, "lines"))
  p
}

# 长文本自动换行（在词边界按 width 字符折行，不加省略号）
# 用足够大的 width 使大多数通路名落在 2 行内；避免纵轴过宽把图压成细缝
wrap_labels <- function(x, width = 40) {
  vapply(as.character(x), function(s) {
    words <- strsplit(s, " ")[[1]]
    out <- ""; cur <- ""
    for (w in words) {
      if (nchar(cur) + nchar(w) + 1 <= width) {
        cur <- if (cur == "") w else paste(cur, w)
      } else {
        out <- if (out == "") cur else paste(out, cur, sep = "\n")
        cur <- w
      }
    }
    if (cur != "") out <- if (out == "") cur else paste(out, cur, sep = "\n")
    out
  }, character(1))
}

# 富集气泡图（直接从富集结果 CSV 绘制，figure 脚本避免重跑 enrichGO）
# legend 默认置底、key 小尺寸，避免挤压绘图区；y 轴标题显式删除；网格线辅助阅读
plot_enrichment_dotplot_csv <- function(res_df, title = NULL, n_show = 20, base_size = 14,
                                        max_desc_chars = 85, wrap_width = 40,
                                        legend.position = "bottom", legend_shift_pt = 100,
                                        ...) {
  df <- as.data.frame(res_df) %>% arrange(.data[["p.adjust"]]) %>% head(n_show)
  df$Description <- wrap_labels(df$Description, width = wrap_width)
  gr_split <- strsplit(as.character(df$GeneRatio), "/")
  df$GeneRatio_num <- round(vapply(gr_split, function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1)), 3)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  p <- ggplot(df, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
    geom_point(alpha = 0.7) +
    # 用 guide 的 order 参数强制 legend 顺序：Count(order=1) 在前，p.adjust(order=2) 在后
    scale_size_continuous(range = c(0.6, 2.0), name = "Count",
                          breaks = scales::pretty_breaks(n = 4),
                          guide = guide_legend(nrow = 1,
                                               keywidth = unit(0.3, "cm"), keyheight = unit(0.3, "cm"),
                                               label.position = "bottom", order = 1)) +
    scale_color_gradient(low = "#B2182B", high = "#2166AC", name = "p.adjust",
                         # 限 3 个断点：取 bar 25%/50%/75% 相对位置，不贴两端
                         breaks = function(limits) {
                           r <- range(limits, na.rm = TRUE)
                           if (r[1] == r[2]) r[1] else r[1] + (r[2] - r[1]) * c(0.25, 0.5, 0.75)
                         },
                         labels = function(x) sprintf("%.0e", x),   # p 值科学记数法（短格式）
                         guide = guide_colorbar(barwidth = unit(2, "cm"), barheight = unit(0.25, "cm"),
                                                order = 2)) +
    labs(title = title %||% "Pathway Enrichment",
         x = "GeneRatio", y = NULL, size = "Count") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          plot.margin = ggplot2::margin(4, 4, 4, 4),
          legend.position = legend.position,
          legend.title = element_text(size = base_size - 1),
          legend.text = element_text(size = base_size - 1),
          legend.spacing.x = unit(0.05, "cm"),   # Count key 更紧凑
          # 底部 legend 居中在"含 y 轴标签的面板"下方导致偏右：加右侧盒边距把整条 legend 左移
          legend.box.margin = ggplot2::margin(0, legend_shift_pt, 0, 0),
          axis.title.y = element_blank(),        # 显式删除 y 轴标题
          axis.text.y = element_text(size = base_size, lineheight = 0.9),
          panel.grid.major = element_line(colour = "grey90", linewidth = pt2mm(0.25)),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1.5, "lines"))
  p
}

# 敏感性分析热图（阈值网格 → pheatmap → ggplot）
plot_sensitivity_heatmap <- function(grid_df, fc_thresholds, p_thresholds, metric = "n_GO_BP",
                                     title = NULL, base_size = 14, ...) {
  heat_mat <- matrix(grid_df[[metric]],
                     nrow = length(fc_thresholds), ncol = length(p_thresholds),
                     byrow = TRUE,   # grid 为 fc 外层循环（fc-major），byrow 保证行=FC、列=p 标注正确
                     dimnames = list(paste0("FC≥", fc_thresholds), paste0("p≤", p_thresholds)))
  # 常数矩阵（如 KEGG 全为 0）会令 pheatmap 的 cut() 报错，扩一个小范围
  rng <- range(as.numeric(heat_mat), na.rm = TRUE)
  if (rng[1] == rng[2]) rng <- c(rng[1], rng[1] + 1)
  breaks <- seq(rng[1], rng[2], length.out = 101)
  p <- pheatmap(heat_mat, cluster_rows = FALSE, cluster_cols = FALSE,
                display_numbers = TRUE, main = title %||% metric,
                color = colorRampPalette(c("white", "#FFD4D4", "#E64B35"))(100),
                breaks = breaks, fontsize = base_size + 2, fontsize_number = base_size, silent = TRUE)
  as_ggplot_pheatmap(p)
}

# 稳定通路柱状图（facet 按比较）
plot_stability_barplot <- function(stable_pathways, n_top = 30, base_size = 16, ...) {
  stable_pathways$PathwayID <- paste(stable_pathways$Type, stable_pathways$Pathway, sep = ": ")
  top_stable <- stable_pathways %>% group_by(Comparison) %>%
    arrange(desc(stability_ratio)) %>% slice_head(n = n_top) %>% ungroup()
  top_stable$PathwayID <- factor(top_stable$PathwayID,
                                 levels = unique(top_stable$PathwayID[order(top_stable$stability_ratio)]))
  ggplot(top_stable, aes(x = stability_ratio, y = PathwayID, fill = Type)) +
    geom_bar(stat = "identity", width = 0.7) +
    facet_wrap(~ Comparison, scales = "free_y") +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(title = "Stable Pathways (≥30% threshold combinations)",
         x = "Stability Ratio", y = NULL, fill = "Type") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2))
}

# ========== 7. Module 2 函数 ==========

# PRC AUC 计算 + Bootstrap CI（沿用 module2_plots.R 逻辑）
calc_prc_auc <- function(pred_prob, true_label, pos_class = "ALS") {
  ord <- order(pred_prob, decreasing = TRUE)
  tp <- cumsum(true_label[ord] == pos_class)
  fp <- cumsum(true_label[ord] != pos_class)
  n_pos <- sum(true_label == pos_class)
  rec <- tp / n_pos
  prec <- tp / (tp + fp)
  prec[is.nan(prec)] <- 0
  rec <- c(0, rec)
  prec <- c(prec[1], prec)
  auc <- sum(diff(rec) * (prec[-1] + prec[-length(prec)]) / 2)
  list(Recall = rec, Precision = prec, AUC = auc)
}

bootstrap_prc_ci <- function(pred_prob, true_label, pos_class = "ALS", n_boot = 500, seed = 42) {
  set.seed(seed)
  n <- length(pred_prob)
  boot_aucs <- replicate(n_boot, {
    idx <- sample(n, replace = TRUE)
    calc_prc_auc(pred_prob[idx], true_label[idx], pos_class)$AUC
  })
  c(quantile(boot_aucs, 0.025, names = FALSE), quantile(boot_aucs, 0.975, names = FALSE))
}

# 从多分类预测构建 ALS-vs-Other 二进制预测表
.as_als_binary <- function(multi_pred, pos_class = "ALS") {
  multi_pred %>%
    mutate(TrueLabel = ifelse(TrueLabel == pos_class, pos_class, "Other"),
           Prob_pos = .data[[pos_class]]) %>%
    select(Model, TestSet, TrueLabel, Prob_pos)
}

# ROC 曲线（每模型一条，AUC+CI 图例）
plot_roc <- function(multi_pred, test_set = NULL, pos_class = "ALS", ci = TRUE,
                     n_boot = 500, seed = 42, colors = MODEL_COLORS, base_size = 18,
                     title = NULL, legend.position = c(0.6, 0.18),
                     show_reference = TRUE, ...) {
  pred_df <- .as_als_binary(multi_pred, pos_class)
  if (!is.null(test_set)) pred_df <- pred_df %>% filter(TestSet == test_set)
  models <- unique(pred_df$Model)
  roc_list <- list(); roc_ann <- data.frame()
  for (m in models) {
    sub_m <- pred_df %>% filter(Model == m)
    y_bin <- ifelse(sub_m$TrueLabel == pos_class, 1, 0)
    roc_obj <- suppressMessages(roc(y_bin, sub_m$Prob_pos, quiet = TRUE, direction = "<"))
    if (ci) {
      roc_ci <- pROC::ci(roc_obj, of = "auc", method = "bootstrap", boot.n = n_boot, quiet = TRUE)
    } else {
      roc_ci <- rep(auc(roc_obj), 3)
    }
    roc_coords <- coords(roc_obj, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
    roc_list[[m]] <- data.frame(Model = m, FPR = 1 - roc_coords$specificity, TPR = roc_coords$sensitivity)
    roc_ann <- rbind(roc_ann, data.frame(Model = m,
                                         AUC = round(auc(roc_obj), 3),
                                         CI_low = round(roc_ci[1], 3),
                                         CI_high = round(roc_ci[3], 3)))
  }
  roc_df <- bind_rows(roc_list)
  roc_labels <- if (ci) {
    setNames(sprintf("%s  AUC=%.3f (%.3f-%.3f)", roc_ann$Model, roc_ann$AUC,
                     roc_ann$CI_low, roc_ann$CI_high), roc_ann$Model)
  } else {
    setNames(sprintf("%s  AUC=%.3f", roc_ann$Model, roc_ann$AUC), roc_ann$Model)
  }
  p <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
    geom_path(linewidth = pt2mm(PUB_LINE_PT)) +
    scale_color_manual(values = colors, labels = roc_labels)
  if (isTRUE(show_reference)) {
    p <- p + geom_abline(slope = 1, intercept = 0, color = "grey50",
                         linewidth = pt2mm(0.5), linetype = "dashed")
  }
  p <- p +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title %||% paste("ROC (from multiclass ALS prob) -", test_set %||% "all"),
         x = "False Positive Rate",
         y = "True Positive Rate", color = NULL) +
    theme_few(base_size = base_size) +
    theme(legend.position = legend.position,
          legend.justification = c(0.5, 0.5),   # 指定坐标为图例中心，同类图例位置一致
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.title = element_blank(),
          legend.text = element_text(size = base_size * 0.8, lineheight = 0.8),
          legend.key.size = unit(0.25, "cm"),
          legend.spacing.y = unit(0.03, "cm"),
          axis.title = element_text(size = base_size + 2),
          axis.text = element_text(size = base_size - 1))
  p
}

# PRC 曲线（每模型一条，AUC+CI 图例）
plot_prc <- function(multi_pred, test_set = NULL, pos_class = "ALS", ci = TRUE,
                     n_boot = 500, seed = 42, colors = MODEL_COLORS, base_size = 18,
                     title = NULL, legend.position = c(0.4, 0.18),
                     show_reference = TRUE, ...) {
  pred_df <- .as_als_binary(multi_pred, pos_class)
  if (!is.null(test_set)) pred_df <- pred_df %>% filter(TestSet == test_set)
  models <- unique(pred_df$Model)
  prc_list <- list(); prc_ann <- data.frame()
  for (m in models) {
    sub_m <- pred_df %>% filter(Model == m)
    prc_res <- calc_prc_auc(sub_m$Prob_pos, sub_m$TrueLabel, pos_class)
    prc_ci <- if (ci) bootstrap_prc_ci(sub_m$Prob_pos, sub_m$TrueLabel, pos_class, n_boot, seed) else c(prc_res$AUC, prc_res$AUC)
    prc_list[[m]] <- data.frame(Model = m, Recall = prc_res$Recall, Precision = prc_res$Precision)
    prc_ann <- rbind(prc_ann, data.frame(Model = m, AUC = round(prc_res$AUC, 3),
                                         CI_low = round(prc_ci[1], 3), CI_high = round(prc_ci[2], 3)))
  }
  prc_df <- bind_rows(prc_list)
  prc_df <- prc_df[is.finite(prc_df$Precision) & is.finite(prc_df$Recall), ]
  prc_labels <- if (ci) {
    setNames(sprintf("%s  AUC=%.3f (%.3f-%.3f)", prc_ann$Model, prc_ann$AUC,
                     prc_ann$CI_low, prc_ann$CI_high), prc_ann$Model)
  } else {
    setNames(sprintf("%s  AUC=%.3f", prc_ann$Model, prc_ann$AUC), prc_ann$Model)
  }
  p <- ggplot(prc_df, aes(x = Recall, y = Precision, color = Model)) +
    geom_path(linewidth = pt2mm(PUB_LINE_PT)) +
    scale_color_manual(values = colors, labels = prc_labels)
  if (isTRUE(show_reference)) {
    p <- p + geom_abline(slope = -1, intercept = 1, color = "grey50", linetype = "dashed",
                         linewidth = pt2mm(0.5))
  }
  p <- p +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title %||% paste("PRC Curve -", test_set %||% "all"),
         x = "Recall", y = "Precision", color = NULL) +
    theme_few(base_size = base_size) +
    theme(legend.position = legend.position,
          legend.justification = c(0.5, 0.5),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.title = element_blank(),
          legend.text = element_text(size = base_size * 0.8, lineheight = 0.8),
          legend.key.size = unit(0.25, "cm"),
          legend.spacing.y = unit(0.03, "cm"),
          axis.title = element_text(size = base_size + 2),
          axis.text = element_text(size = base_size - 1))
  p
}

# 混淆矩阵（自动选 Accuracy 最高的模型）
plot_confusion_matrix <- function(multi_pred, test_set = NULL, model = NULL, base_size = 18, title = NULL, ...) {
  cm_acc <- multi_pred %>% group_by(Model) %>%
    summarise(Accuracy = mean(TrueLabel == PredClass), .groups = "drop")
  best_cm <- if (is.null(model)) cm_acc$Model[which.max(cm_acc$Accuracy)] else model
  sub_best <- multi_pred %>% filter(Model == best_cm)
  if (!is.null(test_set)) sub_best <- sub_best %>% filter(TestSet == test_set)
  all_classes <- sort(unique(c(as.character(sub_best$PredClass), as.character(sub_best$TrueLabel))))
  cm_tbl <- table(Predicted = sub_best$PredClass, Reference = sub_best$TrueLabel)
  df_cm <- as.data.frame(as.table(cm_tbl))
  colnames(df_cm) <- c("Predicted", "Reference", "Freq")
  df_cm$Predicted <- factor(df_cm$Predicted, levels = all_classes)
  df_cm$Reference <- factor(df_cm$Reference, levels = all_classes)
  df_cm <- df_cm %>%
    group_by(Reference) %>%
    mutate(Pct = Freq / sum(Freq) * 100,
           Label = sprintf("%d\n(%.1f%%)", Freq, Pct))
  acc_val <- mean(sub_best$TrueLabel == sub_best$PredClass)
  p <- ggplot(df_cm, aes(x = Reference, y = Predicted, fill = Freq)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Label), size = pt2mm(base_size * 0.96), fontface = "bold") +
    scale_fill_gradient(low = "white", high = "#2c7bb6", name = "Count",
                        guide = guide_colorbar(barwidth = unit(0.25, "cm"), barheight = unit(3, "cm"))) +
    coord_fixed() +   # 保持方形
    labs(title = title %||% paste("Confusion Matrix -", best_cm, "-", test_set %||% "all"),
         subtitle = sprintf("Accuracy = %.3f", acc_val)) +
    theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          plot.subtitle = element_text(hjust = 0.5, size = base_size),
          axis.title = element_text(size = base_size + 2),
          axis.text = element_text(size = base_size - 1, face = "bold"),
          legend.title = element_text(size = base_size - 2),
          legend.text = element_text(size = base_size - 4))
  p
}

# SHAP 特征重要性柱状图
plot_shap_importance <- function(shap_imp, title = NULL, base_size = 20, ...) {
  shap_imp <- shap_imp %>% arrange(MeanAbsSHAP)
  shap_imp$Feature <- factor(shap_imp$Feature, levels = shap_imp$Feature)
  p <- ggplot(shap_imp, aes(x = MeanAbsSHAP, y = Feature, fill = MeanAbsSHAP)) +
    geom_bar(stat = "identity", color = NA, width = 0.7) +
    scale_fill_gradient(low = "#F7F7F7", high = "#B2182B", guide = "none") +
    labs(title = title %||% "SHAP Feature Importance", x = "Mean |SHAP Value|", y = NULL) +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 6),
          axis.title = element_text(size = base_size + 2),
          axis.text = element_text(size = base_size - 1))
  p
}

# SHAP Beeswarm（Top N 特征；feature 值归一化配色）
plot_shap_beeswarm <- function(shap_imp, shap_vals, feat_vals, top_n = 15, title = NULL,
                               subtitle = NULL, base_size = 20, ...) {
  top_n_beeswarm <- min(top_n, nrow(shap_imp))
  shap_imp_ord <- shap_imp %>% arrange(MeanAbsSHAP)
  shap_imp_ord$Feature <- factor(shap_imp_ord$Feature, levels = shap_imp_ord$Feature)
  top_feat <- tail(levels(shap_imp_ord$Feature), top_n_beeswarm)

  shap_bee <- shap_vals %>%
    select(all_of(top_feat)) %>%
    mutate(SampleID = row_number()) %>%
    pivot_longer(-SampleID, names_to = "Feature", values_to = "SHAP")

  feat_bee <- feat_vals %>%
    select(all_of(top_feat)) %>%
    mutate(SampleID = row_number()) %>%
    pivot_longer(-SampleID, names_to = "Feature", values_to = "RawValue") %>%
    group_by(Feature) %>%
    mutate(MinVal = min(RawValue, na.rm = TRUE),
           MaxVal = max(RawValue, na.rm = TRUE),
           Feature_Value = ifelse(MaxVal > MinVal, (RawValue - MinVal) / (MaxVal - MinVal), 0.5)) %>%
    ungroup() %>%
    select(SampleID, Feature, Feature_Value)

  bee_data <- shap_bee %>%
    inner_join(feat_bee, by = c("SampleID", "Feature")) %>%
    mutate(Feature = factor(Feature, levels = top_feat))

  p <- ggplot(bee_data, aes(x = SHAP, y = Feature, color = Feature_Value)) +
    geom_point(alpha = 0.6, size = 0.8, shape = 16, position = position_jitter(height = 0.2)) +
    scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                          midpoint = 0.5, name = "Feature Value\n(normalized)",
                          guide = guide_colorbar(barwidth = unit(0.25, "cm"), barheight = unit(3, "cm"))) +
    labs(title = title %||% "SHAP Beeswarm Plot",
         subtitle = subtitle, x = "SHAP Value (Impact on Model Output)", y = NULL) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 6),
          plot.subtitle = element_text(hjust = 0.5, size = base_size - 3, color = "grey40"),
          axis.title = element_text(size = base_size + 2),
          axis.text = element_text(size = base_size - 1),
          legend.position = "right",
          legend.title = element_text(size = base_size - 4),
          plot.margin = ggplot2::margin(5, 5, 5, 15))   # 左侧留白，避免蛋白名被裁剪
  p
}

# 多分类 SHAP 特征重要性（按类别 dodged bar）
plot_shap_importance_multi <- function(shap_m_imp, title = NULL, base_size = 20,
                                       class_colors = c("ALS" = "#E64B35", "DC" = "#3C5488", "HC" = "#4DBBD5"), ...) {
  desired_order <- intersect(names(class_colors), unique(shap_m_imp$Class))
  shap_m_imp$Class <- factor(shap_m_imp$Class, levels = desired_order)
  cc <- class_colors[intersect(names(class_colors), desired_order)]
  p <- ggplot(shap_m_imp, aes(x = MeanAbsSHAP, y = reorder(Feature, MeanAbsSHAP), fill = Class)) +
    geom_bar(stat = "identity", position = "dodge", color = NA, width = 0.7) +
    scale_fill_manual(values = cc) +
    labs(title = title %||% "Multiclass SHAP Feature Importance",
         x = "Mean |SHAP Value|", y = NULL, fill = "Class") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 6),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 2),
          legend.position = "right")
  p
}

# 单类 SHAP 特征重要性
plot_shap_importance_class <- function(shap_m_imp, class, class_colors = c("ALS" = "#E64B35", "DC" = "#3C5488", "HC" = "#4DBBD5"),
                                       base_size = 20, title = NULL, ...) {
  class_data <- shap_m_imp %>% filter(Class == class) %>% arrange(desc(MeanAbsSHAP))
  class_data$Feature <- factor(class_data$Feature, levels = rev(class_data$Feature))
  p <- ggplot(class_data, aes(x = MeanAbsSHAP, y = Feature, fill = MeanAbsSHAP)) +
    geom_bar(stat = "identity", color = NA, width = 0.7) +
    scale_fill_gradient(low = "#F7F7F7", high = class_colors[class], guide = "none") +
    labs(title = title %||% sprintf("SHAP Feature Importance - %s", class),
         x = "Mean |SHAP Value|", y = NULL) +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 6),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 2))
  p
}

# Benchmark AUC 柱状图（3 特征集 × 4 模型）
plot_benchmark_bar <- function(bm, title = NULL, subtitle = NULL, base_size = 20, ...) {
  bm$Avg_AUC <- (bm$Test1_AUC + bm$Test2_AUC) / 2
  bm$FeatureSet <- factor(bm$FeatureSet, levels = c("Markers", "Metadata", "MarkersPlusMetadata"))
  p <- ggplot(bm, aes(x = Model, y = Avg_AUC, fill = FeatureSet)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.65, color = "grey30", linewidth = 0.3) +
    geom_text(aes(label = round(Avg_AUC, 3)), position = position_dodge(0.65),
              vjust = -0.3, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = FS_COLORS) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    labs(title = title %||% "Benchmark: 3 Feature Sets × 4 Models",
         subtitle = subtitle %||% "Average AUC across validation cohorts",
         x = "Model", y = "Average AUC", fill = "Feature Set") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          plot.subtitle = element_text(hjust = 0.5, size = base_size - 4, color = "grey40"),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 2),
          legend.position = "right",
          legend.title = element_text(size = base_size - 4))
  p
}

# Benchmark ROC（从 benchmark_data.rds$predictions 绘制，facet 按 TestSet）
plot_benchmark_roc <- function(bmd, best_model = NULL, title = NULL, base_size = 20, ...) {
  preds <- bmd$predictions
  ts1 <- bmd$test1_name; ts2 <- bmd$test2_name
  if (is.null(best_model) || length(best_model) == 0) {
    best_model <- unique(sub(".* ", "", unique(preds$ModelSet)))[1]
  }
  roc_list <- list()
  for (fs_name in c("Markers", "Metadata", "MarkersPlusMetadata")) {
    key_pattern <- paste0(fs_name, " ", best_model)
    sub <- preds[grepl(key_pattern, preds$ModelSet), ]
    if (nrow(sub) == 0) next
    for (tn in c(ts1, ts2)) {
      s <- sub[sub$Test == tn, ]
      if (nrow(s) < 5) next
      r <- suppressMessages(roc(s$True, s$Pred, levels = c("Other", "ALS"), direction = "<", quiet = TRUE))
      roc_list[[paste(fs_name, tn)]] <- data.frame(
        FeatureSet = fs_name, TestSet = tn,
        FPR = 1 - r$specificities, TPR = r$sensitivities,
        AUC = round(auc(r), 4), stringsAsFactors = FALSE)
    }
  }
  roc_df <- bind_rows(roc_list)
  if (nrow(roc_df) == 0) {
    stop("plot_benchmark_roc: 没有可绘制的 ROC 数据（ModelSet 未匹配到 best_model）")
  }
  roc_df$FeatureSet <- factor(roc_df$FeatureSet, levels = c("Markers", "Metadata", "MarkersPlusMetadata"))
  ggplot(roc_df, aes(x = FPR, y = TPR, color = FeatureSet)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50",
                linewidth = pt2mm(0.5)) +
    geom_path(linewidth = pt2mm(PUB_LINE_PT)) +
    scale_color_manual(values = FS_COLORS) +
    facet_wrap(~ TestSet) +
    coord_fixed() +
    labs(title = title %||% paste("ROC Comparison -", best_model),
         x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", color = "Feature Set") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          axis.title = element_text(size = base_size - 2),
          axis.text = element_text(size = base_size - 6),
          legend.position = "right")
}

# ALS 亚组分析柱状图
plot_subgroup <- function(sub_df, title = NULL, subtitle = NULL, base_size = 20, ...) {
  sub_df$VarSub <- paste(sub_df$Variable, sub_df$Subgroup, sep = ": ")
  p <- ggplot(sub_df, aes(x = Subgroup, y = MeanProb, fill = TestSet)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6, color = "grey30") +
    geom_errorbar(aes(ymin = MeanProb - SD_Prob, ymax = MeanProb + SD_Prob),
                  position = position_dodge(0.6), width = 0.2) +
    facet_wrap(~ Variable, scales = "free_x", ncol = 3) +
    labs(title = title %||% "ALS Subgroup Analysis",
         subtitle = subtitle %||% "Predicted ALS probability by subgroup",
         x = NULL, y = "Mean Predicted ALS Probability", fill = "Test Set") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          axis.title = element_text(size = base_size - 2),
          axis.text.x = element_text(size = max(9, base_size - 8), angle = 30, hjust = 1),
          axis.text.y = element_text(size = base_size - 6),
          legend.position = "bottom",
          strip.text = element_text(size = base_size - 4, face = "bold"))
  p
}

# NFL vs 12-Marker 面板 AUC 柱状图（数据 = nfl_vs_marker_auc.csv）
# p_values：命名向量/列表（按 TestSet）给出 NFL vs 12-Marker 的组内 AUC 差异 p 值 → 标注显著性
plot_nfl_barplot <- function(nfl_auc, title = NULL, subtitle = NULL, base_size = 16,
                             colors = NFL_BAR_COLORS, p_values = NULL, ...) {
  nfl_auc$Model <- factor(nfl_auc$Model, levels = names(colors))
  nfl_auc$TestSet <- factor(nfl_auc$TestSet)
  p <- ggplot(nfl_auc, aes(x = TestSet, y = AUC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.5, color = "grey30", linewidth = 0.3) +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                  position = position_dodge(0.75), width = 0.15) +
    # AUC 数值标签放在 error bar 顶端之上，避免重叠（字号 0.5×）
    geom_text(aes(label = round(AUC, 3), y = pmax(AUC, CI_high) + 0.03),
              position = position_dodge(0.75), vjust = 0, size = pt2mm(3.5), fontface = "bold") +
    scale_fill_manual(values = colors,
                      guide = guide_legend(keywidth = unit(0.25, "cm"), keyheight = unit(0.25, "cm"))) +
    scale_y_continuous(limits = c(0, 1.35), labels = scales::percent_format()) +
    labs(title = title %||% "12-Marker Panel vs NFL",
         subtitle = subtitle %||% "AUC with 95% CI",
         x = NULL, y = "AUC", fill = "Model") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size),
          plot.subtitle = element_text(hjust = 0.5, size = base_size, color = "grey40"),
          axis.title = element_text(size = base_size + 2),
          axis.text = element_text(size = base_size - 1),
          legend.position = "right",
          legend.title = element_text(size = base_size),
          legend.text = element_text(size = base_size))
  # 组内（每队列）NFL vs 12-Marker 显著性：画一条横跨两柱的括号线 + 星号
  if (!is.null(p_values)) {
    max_top <- nfl_auc %>% group_by(TestSet) %>%
      summarise(m = max(AUC, CI_high, na.rm = TRUE), .groups = "drop")
    ts_levels <- levels(nfl_auc$TestSet) %||% unique(as.character(nfl_auc$TestSet))
    for (i in seq_along(ts_levels)) {
      ts <- ts_levels[i]
      pv <- p_values[[ts]]
      if (is.null(pv) || length(pv) == 0 || is.na(pv)) next
      star <- if (pv < 0.001) "***" else if (pv < 0.01) "**" else if (pv < 0.05) "*" else "n.s."
      top <- max_top$m[max_top$TestSet == ts] + 0.10   # 括号线高度（max_top 已是最高的柱顶/CI 上界）
      x1 <- i - 0.375; x2 <- i + 0.375                  # dodge(0.75) 下两柱的中心位置
      p <- p +
        annotate("segment", x = x1, xend = x2, y = top, yend = top,
                 color = "grey30", linewidth = pt2mm(0.3)) +
        annotate("segment", x = x1, xend = x1, y = top - 0.02, yend = top,
                 color = "grey30", linewidth = pt2mm(0.3)) +
        annotate("segment", x = x2, xend = x2, y = top - 0.02, yend = top,
                 color = "grey30", linewidth = pt2mm(0.3)) +
        annotate("text", x = i, y = top + 0.04, label = star, size = pt2mm(7))
    }
  }
  p
}

# NFL ROC 与 12-Marker RF ROC 叠加（新写；nfl_vals = 每样本 NFL 长表，可 NULL）
#   nfl_vals: data.frame(Sample, TestSet, TrueLabel, NFL)
#   若 nfl_vals 为 NULL，则在 RF ROC 上标注 NFL AUC 虚线（降级路径）
plot_nfl_roc <- function(bin_pred, nfl_vals = NULL, test_set = NULL, nfl_auc_row = NULL,
                         base_size = 16, colors = c("12-Marker RF" = "#E64B35", "NFL" = "#7F7F7F"), ...) {
  # 12-Marker RF ROC（最优模型 RF）
  rf <- bin_pred %>% filter(Model == "RF")
  if (!is.null(test_set)) rf <- rf %>% filter(TestSet == test_set)
  y_bin <- ifelse(rf$TrueLabel == "ALS", 1, 0)
  roc_rf <- suppressMessages(roc(y_bin, rf$Prob_ALS, quiet = TRUE, direction = "<"))
  rf_df <- data.frame(FPR = 1 - roc_rf$specificities, TPR = roc_rf$sensitivities, Curve = "12-Marker RF")

  p <- ggplot(rf_df, aes(x = FPR, y = TPR, color = Curve)) +
    geom_path(linewidth = pt2mm(PUB_LINE_PT)) +
    scale_color_manual(values = colors) +
    geom_abline(slope = 1, intercept = 0, color = "grey50",
                linewidth = pt2mm(0.5), linetype = "dashed") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "False Positive Rate", y = "True Positive Rate", color = NULL) +
    theme_few(base_size = base_size) +
    theme(legend.position = c(0.72, 0.28),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA))

  if (!is.null(nfl_vals)) {
    nfl <- nfl_vals
    if (!is.null(test_set)) nfl <- nfl %>% filter(TestSet == test_set)
    if (nrow(nfl) >= 5) {
      nfl_roc <- suppressMessages(roc(ifelse(nfl$TrueLabel == "ALS", 1, 0), nfl$NFL, quiet = TRUE, direction = "<"))
      nfl_df <- data.frame(FPR = 1 - nfl_roc$specificities, TPR = nfl_roc$sensitivities, Curve = "NFL")
      p <- p + geom_path(data = nfl_df, aes(x = FPR, y = TPR, color = Curve), linewidth = pt2mm(PUB_LINE_PT))
    }
  }
  # AUC 标注
  auc_text <- ""
  if (!is.null(nfl_auc_row)) {
    auc_text <- sprintf("12-Marker RF AUC=%.3f | NFL AUC=%.3f",
                        round(auc(roc_rf), 3), nfl_auc_row$AUC)
  } else {
    auc_text <- sprintf("12-Marker RF AUC=%.3f", round(auc(roc_rf), 3))
  }
  p + annotate("text", x = 0.04, y = 0.04,
               label = sub(" \\| ", "\n", auc_text),
               hjust = 0, vjust = 0, size = 3.0,
               color = "grey20", fontface = "bold")
}

# ========== 8. Module 3 函数 ==========

# Seurat 懒加载
ensure_seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat not installed")
  if (!"package:Seurat" %in% search()) library(Seurat)
}

# PCA Elbow
plot_pca_elbow <- function(stdev, n_pcs = 30, base_size = 18, ...) {
  n <- min(n_pcs, length(stdev))
  df <- data.frame(PC = seq_len(n), Stdev = stdev[seq_len(n)])
  p <- ggplot(df, aes(x = PC, y = Stdev)) +
    geom_point(size = 2.5, color = "#E64B35") +
    geom_line(linewidth = 0.8, color = "#E64B35") +
    scale_x_continuous(breaks = seq(0, n, by = 5)) +
    labs(title = "PCA Elbow Plot", x = "PC", y = "Stdev") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 3),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  p
}

# PCA 方差解释（双轴 bar + line）
plot_pca_variance <- function(stdev, n_pcs = 30, base_size = 18, ...) {
  var_exp <- stdev^2 / sum(stdev^2)
  cum_var <- cumsum(var_exp)
  n <- min(n_pcs, length(var_exp))
  df <- data.frame(PC = seq_len(n), Variance = var_exp[seq_len(n)] * 100, Cumulative = cum_var[seq_len(n)] * 100)
  p <- ggplot(df, aes(x = PC)) +
    geom_bar(aes(y = Variance, fill = "单PC方差"), stat = "identity", alpha = 0.5, width = 0.6) +
    geom_point(aes(y = Cumulative, color = "累积方差"), size = 2.5) +
    geom_line(aes(y = Cumulative, color = "累积方差", group = 1), linewidth = 0.9) +
    scale_x_continuous(breaks = seq(0, n, by = 5)) +
    scale_y_continuous(name = "方差解释率 (%)", limits = c(0, max(df$Cumulative) * 1.1)) +
    scale_fill_manual(name = NULL, values = c("单PC方差" = "#4DBBD5")) +
    scale_color_manual(name = NULL, values = c("累积方差" = "#E64B35")) +
    labs(title = "PCA 方差解释率与累积方差", x = "PC") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 4),
          legend.position = "bottom",
          legend.text = element_text(size = base_size - 3),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  p
}

# 通用降维散点图（合并 UMAP/t-SNE 的 cluster/group/facet/origin/origin_group/HPA 变体）
plot_embedding <- function(coords, x = "UMAP_X", y = "UMAP_Y", color_by = NULL,
                           palette = NULL, facet_by = NULL, label_by = NULL,
                           label_type = c("text_repel", "label_repel", "none"),
                           point_size = 0.3, point_alpha = 0.05, coord_fixed = TRUE,
                           title = NULL, subtitle = NULL, theme = NULL,
                           label_size = 8, label_repel = TRUE, legend_pt_size = 6, ...) {
  label_type <- match.arg(label_type)
  p <- ggplot(coords, aes(x = .data[[x]], y = .data[[y]]))
  if (!is.null(color_by)) p <- p + aes(color = .data[[color_by]])
  p <- p + geom_point(size = point_size, alpha = point_alpha)
  if (!is.null(palette)) p <- p + scale_color_manual(values = palette)

  if (!is.null(label_by)) {
    cen <- coords %>%
      group_by(.data[[label_by]]) %>%
      summarise(!!sym(x) := median(.data[[x]]), !!sym(y) := median(.data[[y]]), .groups = "drop")
    if (!label_repel) {
      # 数字直接位于 cluster 质心，无连接线
      p <- p + geom_text(data = cen, aes(label = .data[[label_by]]),
                         size = pt2mm(label_size), fontface = "bold", color = "black")
    } else if (label_type == "text_repel") {
      p <- p + geom_text_repel(data = cen, aes(label = .data[[label_by]]),
                               size = pt2mm(label_size), fontface = "bold", color = "black",
                               box.padding = 0.3, point.padding = NA,
                               min.segment.length = 0, segment.color = "grey50",
                               segment.size = 0.3, max.overlaps = 50)
    } else if (label_type == "label_repel") {
      p <- p + geom_label_repel(data = cen, aes(label = .data[[label_by]]),
                                size = pt2mm(label_size), fontface = "bold", fill = "white",
                                alpha = 0.85, color = "black", box.padding = 0.5,
                                min.segment.length = 0, segment.color = "grey50",
                                segment.size = 0.3, max.overlaps = 40)
    }
  }
  if (!is.null(color_by)) {
    p <- p + guides(color = guide_legend(override.aes = list(size = legend_pt_size, alpha = 1)))
  }
  if (!is.null(facet_by)) p <- p + facet_wrap(stats::reformulate(facet_by), ncol = 2)
  if (coord_fixed) p <- p + ggplot2::coord_fixed()   # 显式命名空间，避免形参 coord_fixed 遮蔽
  p <- p + labs(title = title, subtitle = subtitle)
  p <- p + (theme %||% theme_clean())
  p
}

# Seurat DotPlot（coord_flip）
plot_marker_dotplot <- function(seu, features, assay = "EV", base_size = 18, title = NULL, ...) {
  ensure_seurat()
  p <- DotPlot(seu, features = features, assay = assay) +
    coord_flip() +
    labs(title = title %||% "各亚群 Top5 标记蛋白表达", x = "Gene", y = "Cluster") +
    theme_few(base_size = base_size) +
    theme(axis.text.x = element_text(size = base_size - 3, face = "bold"),
          axis.text.y = element_text(size = base_size - 4, face = "bold"),
          axis.title.x = element_text(size = base_size, face = "bold"),
          axis.title.y = element_text(size = base_size, face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          legend.title = element_text(size = base_size - 4),
          legend.text = element_text(size = base_size - 6))
  p
}

# Seurat DoHeatmap（每亚群抽 cells_per_cluster 个 EV）
plot_marker_heatmap <- function(seu, features, group_colors = NULL, cells_per_cluster = 200,
                                seed = 42, assay = "EV", base_size = 18, title = NULL, ...) {
  ensure_seurat()
  set.seed(seed)
  cells_use <- unlist(sapply(levels(Idents(seu)), function(cl) {
    cells <- WhichCells(seu, idents = cl)
    if (length(cells) > cells_per_cluster) sample(cells, cells_per_cluster) else cells
  }))
  p <- DoHeatmap(subset(seu, cells = cells_use), features = features, size = 4,
                 assay = assay, group.colors = group_colors) +
    labs(title = title %||% "各亚群 Top5 标记蛋白热图") +
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          axis.text.y = element_text(size = base_size - 4, face = "bold"))
  p
}

# 堆叠占比柱状图（cluster 或 EV origin）
plot_proportion_stack <- function(prop_df, x = "Group", fill = "Cluster", colors = NULL,
                                  title = NULL, base_size = 18, ...) {
  p <- ggplot(prop_df, aes(x = .data[[x]], y = Proportion, fill = .data[[fill]])) +
    geom_col(position = "fill", width = 0.6, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = title %||% "Proportion by Group", x = NULL, y = "Proportion", fill = fill) +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 2),
          legend.title = element_text(size = base_size - 2),
          legend.text = element_text(size = base_size - 4))
  p
}

# 分组箱线图（合并 EV counts 箱线图 + CNS Score 箱线图）
plot_group_boxplot <- function(df, x = "Group", y = NULL, fill = "Group", colors = NULL,
                               title = NULL, subtitle = NULL, add_jitter = TRUE,
                               base_size = 18, ...) {
  stopifnot(!is.null(y))
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA)
  if (add_jitter) p <- p + geom_jitter(width = 0.15, size = 1.5, alpha = 0.4, color = "grey30")
  p <- p + scale_fill_manual(values = colors) +
    labs(title = title %||% "", subtitle = subtitle, x = NULL, y = y) +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
          plot.subtitle = element_text(hjust = 0.5, size = base_size - 4, color = "grey40"),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 2),
          legend.position = "none")
  p
}

# ALS-high 亚群 EV 数量提琴图（返回 ggplot 列表，调用方 wrap_plots / ggarrange 组合）
plot_cluster_violins <- function(prop_df, als_high_clusters, diff_df,
                                 group_levels = ALL_GROUP_LEVELS, group_colors = ALL_GROUP_COLORS,
                                 base_size = 20, ...) {
  plot_list <- list()
  for (cl in als_high_clusters) {
    sub <- prop_df %>% filter(Cluster == cl)
    sample_ev <- sub %>%
      group_by(Sample, Group) %>%
      summarise(Count = sum(Count), .groups = "drop") %>%
      mutate(Group = factor(Group, levels = group_levels))
    test_info <- diff_df %>% filter(Cluster == cl)
    label_text <- sprintf("KW adj.p = %.2e", test_info$KW_padj)
    other_groups <- setdiff(group_levels, "ALS")
    comparisons <- lapply(other_groups, function(g) c("ALS", g))
    p_violin <- ggplot(sample_ev, aes(x = Group, y = Count, fill = Group)) +
      geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.4, color = "grey30") +
      geom_boxplot(width = 0.12, fill = "white", alpha = 0.7, outlier.shape = NA, color = "grey40") +
      scale_fill_manual(values = group_colors) +
      labs(title = paste("Cluster", cl, "- EV Counts"), subtitle = label_text,
           x = NULL, y = "EV Count") +
      theme_few(base_size = base_size) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
            plot.subtitle = element_text(hjust = 0.5, size = base_size - 4, color = "grey40"),
            axis.title = element_text(size = base_size),
            axis.text = element_text(size = base_size - 2),
            legend.position = "none") +
      ggsignif::geom_signif(comparisons = comparisons, map_signif_level = TRUE,
                            test = "wilcox.test", step_increase = 0.1, tip_length = 0.02,
                            size = 0.7, textsize = 5.5, color = "grey20")
    plot_list[[length(plot_list) + 1]] <- p_violin
  }
  plot_list
}

# 来源 marker 阳性率热图（pheatmap → ggplot）
plot_origin_heatmap <- function(pos_mat, cluster_origin_vec, ct_colors = CT_COLORS,
                                title = NULL, ...) {
  cluster_order <- as.character(sort(as.numeric(rownames(pos_mat))))
  origin_vec <- cluster_origin_vec[cluster_order]
  annotation_row <- data.frame(Origin = factor(origin_vec, levels = names(ct_colors)),
                               row.names = cluster_order)
  ann_colors <- list(Origin = ct_colors)
  p <- pheatmap(pos_mat, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
                annotation_row = annotation_row, annotation_colors = ann_colors,
                color = colorRampPalette(c("white", "#FFD4D4", "#E64B35"))(100),
                main = title %||% "亚群 × 来源标志物 阳性率",
                fontsize = 10, fontsize_row = 8, fontsize_col = 8, silent = TRUE)
  as_ggplot_pheatmap(p)
}

# cluster top-marker 气泡图（DotPlot 风格：color=平均表达 z-score, size=阳性率非线性, PosRate 从 0 映射）
# marker_data 需含 Cluster, Protein, PosRate, avg_expr；每 (cluster, marker) 组合都要有一行（PosRate=0 也显示小圈）
# 颜色：avg_expr -> log1p（压缩右偏）-> 全局 z-score -> 颜色范围按数据分位数 [z_lo, z_hi]（数据驱动，
#   避免固定 [-2,2] 把非对称分布压到渐变中段显紫红），oob=squish 让极端值落到端点色
# size：PosRate^size_exp 非线性（接近 1 的圈更突出），整体偏小、0 最小
plot_cluster_marker_bubble <- function(marker_data, base_size = 14, z_lo = 0.01, z_hi = 0.99,
                                       size_exp = 1.8, size_range = c(0.4, 3.0), ...) {
  d <- marker_data
  color_lim <- NULL
  # 颜色标准化：log1p 后全局 z-score（若无预计算列 avg_expr_norm）
  if ("avg_expr" %in% names(d) && !"avg_expr_norm" %in% names(d)) {
    ex <- log1p(pmax(d$avg_expr, 0))
    z <- (ex - mean(ex, na.rm = TRUE)) / stats::sd(ex, na.rm = TRUE)
    lo <- unname(quantile(z, z_lo, na.rm = TRUE)); hi <- unname(quantile(z, z_hi, na.rm = TRUE))
    d$avg_expr_norm <- if (hi > lo) pmin(pmax(z, lo), hi) else z
    color_lim <- c(lo, hi)
  }
  color_col <- if ("avg_expr_norm" %in% names(d)) "avg_expr_norm" else "avg_expr"
  if (is.null(color_lim)) color_lim <- range(d[[color_col]], na.rm = TRUE)
  # 色标刻度（整数 z 值，落在范围内才显示）
  cb_breaks <- seq(ceiling(color_lim[1]), floor(color_lim[2]), by = 1)
  # size 非线性映射：PosRate^size_exp（p>1 使接近 1 的圈更大更突出、0 更小）
  d$size_map <- d$PosRate ^ size_exp
  size_pos <- c(0, 0.25, 0.5, 0.75, 1)
  p <- ggplot(d, aes(x = Cluster, y = Protein, size = size_map, color = .data[[color_col]])) +
    geom_point(alpha = 0.85) +
    # 蓝红连续渐变（低=蓝 #2166AC，高=红 #B2182B），风格同 Figure1 富集气泡图；红=高表达、蓝=低表达
    scale_color_gradient(low = "#2166AC", high = "#B2182B", name = "Expression\n(z-score)",
                         limits = color_lim, oob = scales::squish,
                         breaks = if (length(cb_breaks) > 0) cb_breaks else waiver(),
                         guide = guide_colorbar(barwidth = unit(0.25, "cm"),   # 与 Figure3 beeswarm 色条一致
                                                barheight = unit(1.5, "cm"))) +
    scale_size_continuous(range = size_range, name = "PosRate",
                          breaks = size_pos ^ size_exp, labels = as.character(size_pos),
                          guide = guide_legend(keywidth = unit(0.25, "cm"), keyheight = unit(0.25, "cm"))) +
    labs(x = "EV Cluster", y = NULL) +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
          axis.text.x = element_text(size = base_size - 2, angle = 45, hjust = 1),
          axis.text.y = element_text(size = base_size - 1),
          plot.margin = ggplot2::margin(5, 5, 5, 12),   # 左侧留白避免蛋白名被裁剪
          legend.spacing.y = unit(0.1, "cm"),   # PosRate 图例 key 间距更小
          panel.grid.major = element_line(colour = "grey90", linewidth = pt2mm(0.25)),   # 淡网格线
          panel.grid.minor = element_blank(),
          legend.position = "right")
  p
}

# cluster 比例按 Group 箱线图（x=cluster, y=Proportion, fill=Group）
# x_labels：命名向量（cluster → 标签），替换 x 轴刻度；y_lower：线性 y 轴下边界（设负数让箱线图上升）
# sig_padj：命名向量（cluster → ALS vs 其它 的 BH-padj），在每 cluster 顶部加显著性星标
# y_log2=TRUE：log2 变换 y 轴（0 映射到最小正值/2 避免 -Inf，保留 0 观测为底部点），显示 outlier；
#   log2 压缩高离群值、放大低比例 cluster 间差异，轴刻度用原始百分比
plot_cluster_prop_boxplot <- function(prop_df, clusters, colors = ALL_GROUP_COLORS, base_size = 14,
                                      x_labels = NULL, y_lower = NULL, sig_padj = NULL,
                                      y_log2 = FALSE, ...) {
  d <- prop_df %>% filter(Cluster %in% clusters)
  d$Cluster <- factor(d$Cluster, levels = clusters)
  gth <- resolve_group_theme(d$Group, colors = colors)
  d$Group <- factor(d$Group, levels = gth$levels)
  # y 值：线性比例 或 log2（0 → 最小正值/2）
  if (y_log2) {
    pos <- d$Proportion[d$Proportion > 0]
    eps <- if (length(pos) > 0) min(pos) / 2 else 1e-4
    d$y_val <- log2(pmax(d$Proportion, eps))
  } else {
    d$y_val <- d$Proportion
  }
  p <- ggplot(d, aes(x = Cluster, y = y_val, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = pt2mm(1), outlier.alpha = 0.5) +
    scale_fill_manual(values = gth$colors, name = "Group") +
    labs(x = "EV Cluster", y = if (y_log2) "Proportion of EVs (log2)" else "Proportion of EVs") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
          axis.text.x = element_text(size = base_size - 2),
          plot.margin = ggplot2::margin(8, 5, 5, 18),   # 左侧留白避免面板标签遮挡纵坐标数字
          legend.position = "right")
  if (!is.null(x_labels)) {
    p <- p + scale_x_discrete(labels = x_labels)
  }
  # y 轴刻度与范围
  if (y_log2) {
    brks_prop <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
    brks_lab <- ifelse(brks_prop >= 0.01, sprintf("%g%%", 100 * brks_prop),
                       sprintf("%.1f%%", 100 * brks_prop))
    p <- p + scale_y_continuous(breaks = log2(brks_prop), labels = brks_lab,
                                expand = expansion(mult = c(0.03, 0.15)))
  } else if (!is.null(y_lower)) {
    p <- p + scale_y_continuous(limits = c(y_lower, NA),
                                expand = expansion(mult = c(0, if (!is.null(sig_padj)) 0.15 else 0.05)))
  }
  # 显著性星标（放每 cluster 最高观测之上，含 outlier）
  if (!is.null(sig_padj)) {
    star <- ifelse(sig_padj < 1e-3, "***", ifelse(sig_padj < 1e-2, "**",
                   ifelse(sig_padj < 0.05, "*", "ns")))
    ymax <- d %>% group_by(Cluster) %>% summarise(ym = max(y_val, na.rm = TRUE), .groups = "drop")
    off <- if (y_log2) 0.5 else (max(d$y_val, na.rm = TRUE) - min(d$y_val, na.rm = TRUE)) * 0.04
    sig_df <- data.frame(cl_id = names(sig_padj), star = unname(star), stringsAsFactors = FALSE)
    sig_df$y_pos <- ymax$ym[match(sig_df$cl_id, as.character(ymax$Cluster))] + off
    sig_df$Cluster <- factor(sig_df$cl_id, levels = clusters)   # 与 x 轴 cluster 层级一致，保证对齐
    p <- p + geom_text(data = sig_df, aes(x = Cluster, y = y_pos, label = star),
                       inherit.aes = FALSE, size = pt2mm(5), vjust = 0)
  }
  p
}

# ========== 9. Module 4 函数 ==========

# DE 蛋白 × EV cluster 阳性率热图
plot_de_cluster_heatmap <- function(de_mat_sel, annotation_row = NULL, ann_colors = NULL,
                                    title = NULL, ...) {
  p <- pheatmap(de_mat_sel, scale = "none", cluster_rows = TRUE, cluster_cols = FALSE,
                annotation_row = annotation_row, annotation_colors = ann_colors,
                color = colorRampPalette(c("white", "#FFD4D4", "#E64B35"))(100),
                main = title %||% "DE Protein Positive Rate across EV Clusters",
                fontsize = 8, fontsize_row = 7, fontsize_col = 9, silent = TRUE)
  as_ggplot_pheatmap(p)
}

# 核心标志物 × cluster 条形图（facet 按 Protein；x_col 可切换 Cluster 或 HPA label）
plot_core_cluster_bar <- function(core_pos, fill = "origin_group", colors = ORIGIN_GROUP_COLORS,
                                  x_col = "Cluster", reorder_x = FALSE,
                                  facet_ncol = 4, title = NULL, subtitle = NULL,
                                  base_size = 14, ...) {
  p <- ggplot(core_pos, aes(x = .data[[x_col]], y = PosRate, fill = .data[[fill]])) +
    geom_bar(stat = "identity", width = 0.7, color = "grey30", linewidth = 0.2) +
    facet_wrap(~Protein, ncol = facet_ncol, scales = "free_y") +
    scale_fill_manual(values = colors) +
    labs(title = title %||% "Core Markers: Positive Rate across EV Clusters",
         subtitle = subtitle, x = "EV Cluster", y = "Positive Rate", fill = "Origin") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          plot.subtitle = element_text(hjust = 0.5, size = base_size - 1, color = "grey40"),
          strip.text = element_text(face = "bold", size = base_size - 3),
          axis.text.x = element_text(size = base_size - 6, angle = 45, hjust = 1))
  if (reorder_x) {
    p <- p + aes(x = reorder(.data[[x_col]], PosRate))
  }
  p
}

# DE × cluster 组间差异气泡热图
plot_diff_bubble <- function(sig_diff, title = NULL, subtitle = NULL, base_size = 14, ...) {
  p <- ggplot(sig_diff, aes(x = Cluster, y = reorder(Protein, logOR),
                            size = -log10(p_adjusted), color = logOR)) +
    geom_point(alpha = 0.7) +
    scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", name = "log2(OR)") +
    labs(title = title %||% "DE Protein × EV Cluster: ALS vs HC (single EV)",
         subtitle = subtitle %||% "Fisher exact test, BH adjusted p < 0.05",
         x = "EV Cluster", y = "DE Protein", size = "-log10(p.adj)") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          axis.text.x = element_text(size = base_size - 4, angle = 45, hjust = 1))
  p
}

# 跨模块一致性散点（bulk log2FC vs 单 EV 阳性率）
plot_consistency_scatter <- function(consistency, cor_val = NULL, title = NULL,
                                     subtitle = NULL, base_size = 16, ...) {
  if (is.null(cor_val)) cor_val <- cor(consistency$log2FC, consistency$EV_pos_rate, method = "spearman")
  p <- ggplot(consistency, aes(x = log2FC, y = EV_pos_rate, color = is_core)) +
    geom_point(size = 2.5, alpha = 0.6) +
    geom_smooth(method = "lm", color = "grey50", linewidth = 0.5, se = TRUE) +
    geom_text(data = subset(consistency, is_core), aes(label = Protein),
              size = 3.5, fontface = "bold", vjust = -0.8) +
    scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey60")) +
    annotate("text", x = min(consistency$log2FC) * 0.9, y = max(consistency$EV_pos_rate) * 0.95,
             label = paste("Spearman r =", round(cor_val, 4)), size = 5, hjust = 0) +
    labs(title = title %||% "Cross-Module Consistency: Bulk DE vs Single EV",
         subtitle = subtitle %||% "log2FC (ALS vs HC, Jinan) vs PosRate in ALS-high clusters",
         x = "Bulk log2FC (Module 1)", y = "Single EV PosRate (Module 3)",
         color = "Core Marker") +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4))
  p
}

# DE 蛋白来源饼图
plot_origin_pie <- function(origin_summary, title = NULL, subtitle = NULL,
                            colors = ORIGIN_GROUP_COLORS, base_size = 16, ...) {
  origin_summary$fraction <- origin_summary$fraction / sum(origin_summary$fraction)
  origin_summary$label <- paste0(origin_summary$origin_group, "\n",
                                 round(origin_summary$fraction * 100, 1), "%")
  p <- ggplot(origin_summary, aes(x = "", y = fraction, fill = origin_group)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y") +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 6, fontface = "bold") +
    labs(title = title %||% "DE Protein Distribution by EV Origin",
         subtitle = subtitle) +
    theme_void(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          plot.subtitle = element_text(hjust = 0.5, size = base_size - 3, color = "grey40"),
          legend.position = "none")
  p
}

# ========== 10. 新增面板 ==========

# 12 核心标志物分组箱线图（Figure 2C；数据 = TMM 表达矩阵 + sample info + core_markers）
# 只显示箱线图 + 离群点（pch=16 实心圆，无描边）
plot_marker_boxplots <- function(expr, sample_info, proteins, group_cols = ALL_GROUP_COLORS,
                                 center_filter = NULL, facet_ncol = 4, base_size = 14,
                                 title = NULL, ylab = NULL, exclude = "AD", log2_transform = TRUE,
                                 add_signif = TRUE, ...) {
  keep <- intersect(proteins, rownames(expr))
  if (length(keep) == 0) stop("没有核心标志物在表达矩阵中")
  long <- expr[keep, , drop = FALSE] %>%
    as.data.frame() %>%
    mutate(Protein = rownames(.)) %>%
    pivot_longer(-Protein, names_to = "sample_name", values_to = "value")
  meta <- sample_info %>% select(any_of(c("sample_name", "group", "center")))
  long <- long %>% inner_join(meta, by = "sample_name")
  if (!is.null(center_filter)) long <- long %>% filter(center %in% center_filter)
  long <- long %>% filter(!(group %in% exclude)) %>% filter(group %in% names(group_cols))
  if (log2_transform) long <- long %>% mutate(value = log2(value + 1))
  ylab <- ylab %||% if (log2_transform) "log2(TMM + 1)" else "TMM expression"
  gth <- resolve_group_theme(long$group, colors = group_cols)
  long$group <- factor(long$group, levels = gth$levels)
  p <- ggplot(long, aes(x = group, y = value, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = pt2mm(1.5), outlier.alpha = 0.6) +
    scale_fill_manual(values = gth$colors) +
    facet_wrap(~ Protein, scales = "free_y", ncol = facet_ncol) +
    labs(title = title %||% "Core Markers: Expression by Group",
         x = NULL, y = ylab, fill = NULL) +
    theme_few(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(size = base_size - 2),
          legend.position = "none")
  if (add_signif && length(gth$levels) >= 2) {
    comps <- utils::combn(gth$levels, 2, simplify = FALSE)
    # "ns" → "n.s."；textsize 比默认小 2pt
    sig_label <- function(pv) {
      ifelse(pv < 0.001, "***", ifelse(pv < 0.01, "**", ifelse(pv < 0.05, "*", "n.s.")))
    }
    p <- p + ggsignif::geom_signif(comparisons = comps, map_signif_level = sig_label,
                                   test = "wilcox.test", step_increase = 0.1,
                                   tip_length = 0.02, size = pt2mm(0.5), textsize = pt2mm(4))
  }
  p
}

# ========== 11. 组合助手（figure_publication.R 使用） ==========

# 面板规格构造器（pub_max：该面板的 apply_pub 字号上限，可让某些面板文字更大）
panel <- function(fun, label = NULL, ..., pub_max = NULL) {
  structure(list(fun = fun, label = label, args = list(...), pub_max = pub_max), class = "panel_spec")
}

# ggarrange 排版（扩展自 plots.R 的 arrangeFixPlot）
# 支持两种布局：
#   1) 扁平网格：spec$panels + spec$ncol/nrow/widths/heights
#   2) 多行布局：spec$rows（每行一个 list(ncol/widths/labels/panels)，行间用 spec$heights 分配高度）
compose_figure <- function(spec, panels = NULL) {
  if (!requireNamespace("ggpubr", quietly = TRUE)) stop("ggpubr not installed")

  build_panel <- function(p) {
    f <- get(p$fun, envir = parent.frame())
    apply_pub(do.call(f, p$args), family = spec$family %||% PUB_FONT,
              max_pt = p$pub_max %||% PUB_MAX_PT)   # 统一 Arial + 字号钳制（pub_max 可放宽）
  }
  arrange_one <- function(row) {
    rpl <- lapply(row$panels, build_panel)
    ra <- list(plotlist = rpl, ncol = row$ncol %||% length(rpl),
               labels = row$labels %||% NULL,
               font.label = list(size = 10, face = "bold"),
               common.legend = FALSE, legend = NULL)
    if (!is.null(row$widths)) ra$widths <- row$widths
    if (!is.null(row$heights)) ra$heights <- row$heights
    do.call(ggpubr::ggarrange, ra)
  }

  # ---- 多行布局 ----
  if (!is.null(spec$rows)) {
    row_plots <- lapply(spec$rows, arrange_one)
    ga <- list(plotlist = row_plots, nrow = length(row_plots), ncol = 1,
               labels = NULL, common.legend = FALSE, legend = NULL)
    if (!is.null(spec$heights)) ga$heights <- spec$heights
    return(do.call(ggpubr::ggarrange, ga))
  }

  # ---- 扁平网格 ----
  pl <- lapply(panels, build_panel)
  ncol <- spec$ncol %||% 1
  nrow <- spec$nrow %||% ceiling(length(pl) / ncol)
  labels <- if (isTRUE(spec$labels %||% TRUE)) LETTERS[seq_along(pl)] else NULL
  ga_args <- list(
    plotlist = pl, ncol = ncol, nrow = nrow,
    labels = labels, font.label = list(size = 10, face = "bold"),
    common.legend = isTRUE(spec$common_legend %||% FALSE),
    legend = spec$legend %||% NULL,   # 默认不强制图例位置，各面板保留自身 legend（如火山图置顶）
    align = spec$align %||% NULL)     # align="hv" 在 dotplot 面板上会产生零尺寸视口，默认关闭
  # ggarrange 显式传 NULL widths/heights 会报错，仅在提供时传入
  if (!is.null(spec$widths)) ga_args$widths <- spec$widths
  if (!is.null(spec$heights)) ga_args$heights <- spec$heights
  do.call(ggpubr::ggarrange, ga_args)
}

# 保存 figure（仅输出矢量 PDF，出版尺寸规范）
# width_mm / height_mm 为物理打印尺寸(mm)；不再生成 PNG（按用户要求）
save_figure <- function(g, id, width_mm = PUB_WIDTH_MM, height_mm = PUB_MAX_H_MM, out_dir) {
  showtext_auto(enable = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (height_mm > PUB_MAX_H_MM) {
    warning(sprintf("%s 高度 %.0fmm 超过上限 200mm", id, height_mm))
  }
  ggsave(file.path(out_dir, paste0(id, ".pdf")), g,
         width = width_mm, height = height_mm, units = "mm", device = cairo_pdf,
         limitsize = FALSE, bg = "white")
  cat(sprintf("  -> %s.pdf (%dx%d mm)\n", id, width_mm, height_mm))
}

# ========== 12. 数据整形助手（figure_publication.R 使用） ==========

# Seurat → UMAP 长表
umap_df <- function(seu) {
  ensure_seurat()
  coords <- as.data.frame(Embeddings(seu, "umap"))
  colnames(coords) <- c("UMAP_X", "UMAP_Y")
  coords$cluster <- as.character(Idents(seu))
  coords$group <- factor(seu$group, levels = resolve_group_theme(seu$group)$levels)
  coords
}

# 实际 cluster 标签
clusters <- function(seu) {
  ensure_seurat()
  sort(as.character(unique(Idents(seu))))
}

# 从 top5 markers CSV 派生 dotplot/heatmap 特征（沿用 module3 逻辑）
top_features <- function(markers_df, n = 5) {
  if (is.null(markers_df) || nrow(markers_df) == 0) return(character(0))
  if (!"avg_log2FC" %in% names(markers_df)) return(unique(markers_df$gene))
  topn <- markers_df %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = n, with_ties = FALSE) %>%
    mutate(cluster_num = as.numeric(as.character(cluster))) %>%
    arrange(cluster_num)
  unique(topn$gene, fromLast = FALSE)
}

# 每样本 CNS Score 表（cluster_proportions + ev_origin_annotation）
cns_score_df <- function(prop_df, origin_df) {
  cns_clusters <- origin_df$Cluster[origin_df$origin_group == "CNS"]
  gth <- resolve_group_theme(prop_df$Group)
  cns_score <- prop_df %>%
    group_by(Sample, Group) %>%
    summarise(CNS_Score = sum(Proportion[Cluster %in% cns_clusters], na.rm = TRUE), .groups = "drop")
  cns_score$Group <- factor(cns_score$Group, levels = gth$levels)
  cns_score
}

# CNS Score 的 ALS vs Others Wilcoxon 检验文字
cns_pval_text <- function(cns_score) {
  als <- cns_score$CNS_Score[cns_score$Group == "ALS"]
  others <- cns_score$CNS_Score[cns_score$Group != "ALS"]
  wt <- wilcox.test(als, others, exact = FALSE)
  sprintf("ALS vs Others\nWilcoxon p = %.2e", wt$p.value)
}

# 每样本 NFL 值长表（把 nfl.xlsx 的 B-code 样本与 multiclass_predictions 对齐）
assemble_nfl_values <- function(xlsx_path, multi_pred, base_dir = NULL, ...) {
  if (!file.exists(xlsx_path)) return(NULL)
  if (!requireNamespace("readxl", quietly = TRUE)) return(NULL)
  nfl <- readxl::read_excel(xlsx_path, sheet = 1)
  # 找到样本列（编号/B-code）与 NFL 列
  samp_col <- NULL; nfl_col <- NULL
  for (nm in names(nfl)) {
    v <- as.character(nfl[[nm]])
    if (grepl("^B", v[!is.na(v)][1])) samp_col <- nm
    if (grepl("NFL", nm)) nfl_col <- nm
  }
  if (is.null(samp_col) || is.null(nfl_col)) return(NULL)
  nfl_df <- data.frame(Sample = as.character(nfl[[samp_col]]),
                       NFL = as.numeric(nfl[[nfl_col]]),
                       stringsAsFactors = FALSE) %>%
    filter(!is.na(Sample) & !is.na(NFL)) %>%
    group_by(Sample) %>% summarise(NFL = mean(NFL), .groups = "drop")
  multi_pred %>%
    select(Sample, TestSet, TrueLabel) %>%
    distinct() %>%
    inner_join(nfl_df, by = "Sample", relationship = "many-to-many")
}
