#!/usr/bin/env Rscript
# 整合 PBA 样本队列基线数据 与 Project sample info
#
# 输入:
#   - data/raw/PBA样本队列基线数据.xlsx (4 个 sheet: DC / HC / NDC / ALS)
#   - data/csv/01-2.Project sample info.csv (样本分组信息)
# 输出:
#   - data/csv/Sample info merged.csv (清洗后的样本信息)
#   - data/csv/Sample info merged_column_map.csv (表头新旧对照表)
#
# 匹配方式: Excel 各 sheet 的 编号 <-> CSV$Original.Name (left join)
#
# 数据清洗:
#   - 列名: 中文 → 全小写英文下划线 (cohort → center)
#   - center 值: 中文 → 英文 (济南→Jinan, 青岛→Qingdao, 厦门→Xiamen)
#   - gender 值: 中文 → 英文 (男→Male, 女→Female)
#   - onset_site 值: 中文 → 英文 (球部→Bulbar, ...)
#   - age / kings_stage / functional_score / disease_duration → 数值型
#   - disease_duration 去除 "m" 后缀，单位: 月

# --- 自动获取项目根目录 ---
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}
project_root <- normalizePath(file.path(script_dir, "../.."))

# --- 文件路径 ---
baseline_file <- file.path(project_root, "data/raw/PBA样本队列基线数据.xlsx")
sample_info_file <- file.path(project_root, "data/csv/01-2.Project sample info.csv")
output_file <- file.path(project_root, "data/csv/Sample info merged.csv")
map_file <- file.path(project_root, "data/csv/Sample info merged_column_map.csv")

# --- 读取数据 ---
library(readxl)
sample_info <- read.csv(sample_info_file, stringsAsFactors = FALSE)

# 读取所有 sheet 并合并
sheets <- excel_sheets(baseline_file)
baseline_list <- list()

for (s in sheets) {
  df <- read_excel(baseline_file, sheet = s)
  # ALS sheet 用 "组别" 代替 "疾病分类"，统一列名
  if ("组别" %in% colnames(df)) {
    colnames(df)[colnames(df) == "组别"] <- "疾病分类"
  }
  # 统一所有列为字符型，避免跨 sheet 类型不一致导致 bind 失败
  df[] <- lapply(df, as.character)
  baseline_list[[s]] <- df
}

# 行合并（不同 sheet 列数不同时自动补 NA）
baseline <- do.call(dplyr::bind_rows, baseline_list)

cat(sprintf("[merge_sample_info] 基线数据: %d × %d (%d sheets)\n",
            nrow(baseline), ncol(baseline), length(sheets)))
cat(sprintf("[merge_sample_info] 样本信息: %d × %d\n",
            nrow(sample_info), ncol(sample_info)))

# ---- 1. 合并 ----
merged <- merge(
  sample_info,
  baseline,
  by.x = "Original.Name",
  by.y = "编号",
  all.x = TRUE
)

# ---- 2. 列名映射：中文 → 全小写英文下划线 ----
# 使用 merge 后实际列名，King's分期中的引号为 U+2018，需动态匹配
kings_col <- grep("分期", colnames(merged), value = TRUE)

col_map <- c(
  "Original.Name" = "original_name",
  "Sample.Name"   = "sample_name",
  "Group"         = "group",
  "队列"          = "center",
  "性别"          = "gender",
  "年龄"          = "age",
  "疾病分类"      = "disease",
  "发病部位"      = "onset_site",
  "功能评分"      = "functional_score",
  "就诊时病程"    = "disease_duration",
  "就诊时间"      = "visit_date"
)

if (length(kings_col) == 1) {
  col_map <- c(col_map, setNames("kings_stage", kings_col))
}

# 保存列名对照表
map_df <- data.frame(
  original_name = names(col_map),
  new_name      = unname(col_map),
  stringsAsFactors = FALSE
)
write.csv(map_df, map_file, row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("[merge_sample_info] 列名对照表: %s\n", map_file))

# 重命名列
colnames(merged) <- ifelse(colnames(merged) %in% names(col_map),
                           col_map[colnames(merged)],
                           colnames(merged))

# ---- 3. 数据清洗 ----

# 3a. center: 中文中心名 → 英文
center_map <- c(
  "济南"     = "Jinan",
  "青岛"     = "Qingdao",
  "厦门"     = "Xiamen",
  "厦门队列" = "Xiamen",
  "青岛队列" = "Qingdao"
)
merged$center <- ifelse(merged$center %in% names(center_map),
                        center_map[merged$center],
                        merged$center)

# 3b. gender: 中文 → 英文
gender_map <- c("男" = "Male", "女" = "Female")
merged$gender <- ifelse(merged$gender %in% names(gender_map),
                        gender_map[merged$gender],
                        merged$gender)

# 3c. onset_site: 中文发病部位 → 英文
onset_site_map <- c(
  "球部"       = "Bulbar",
  "四肢"       = "Quadrilateral",
  "双上肢"     = "Bilateral_upper_limbs",
  "双下肢"     = "Bilateral_lower_limbs",
  "右下肢"     = "Right_lower_limb",
  "左下肢"     = "Left_lower_limb",
  "右上肢"     = "Right_upper_limb",
  "左上肢"     = "Left_upper_limb",
  "左侧肢体"   = "Left_limb",
  "右上下肢"   = "Right_upper_and_lower_limbs",
  "呼吸肌"     = "Respiratory"
)
merged$onset_site <- ifelse(merged$onset_site %in% names(onset_site_map),
                            onset_site_map[merged$onset_site],
                            merged$onset_site)

# 3d. age → 数值型
merged$age <- as.numeric(merged$age)

# 3e. kings_stage: 去除"期"后缀 → 数值型 (如 "1期" → 1)
merged$kings_stage <- as.numeric(gsub("期", "", merged$kings_stage))

# 3f. functional_score → 数值型
merged$functional_score <- as.numeric(merged$functional_score)

# 3g. disease_duration: 去除 "m" 后缀 → 数值型，单位: 月
merged$disease_duration <- as.numeric(gsub("m", "", merged$disease_duration, ignore.case = TRUE))

# ---- 4. 输出 ----
write.csv(merged, output_file, row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("[merge_sample_info] 输出文件: %s\n", output_file))
cat(sprintf("[merge_sample_info] 整合完成: %d × %d\n", nrow(merged), ncol(merged)))
cat(sprintf("[merge_sample_info] 列名: %s\n", paste(colnames(merged), collapse = ", ")))
