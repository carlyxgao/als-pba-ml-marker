#!/bin/bash
# 从 raw data zip 中解压 CSV 文件到 data/csv，不保留目录结构
# 自动获取脚本所在路径，无需硬编码

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

ZIP_FILE="$PROJECT_ROOT/data/raw/RepReqPBA_COL031_Qilu.pla-02.reBE.report1.zip"
OUTPUT_DIR="$PROJECT_ROOT/data/csv"

mkdir -p "$OUTPUT_DIR"

FILES=(
  "report/src/Result/Folder01.Data file/01-2.Project sample info.csv"
  "report/src/Result/Folder03.Protein expression/03-1.Protein expression raw data.csv"
  "report/src/Result/Folder03.Protein expression/03-2.Protein expression TMM nomalized.csv"
)

for file in "${FILES[@]}"; do
  echo "解压: $file"
  unzip -j -o "$ZIP_FILE" "$file" -d "$OUTPUT_DIR"
done

echo ""
echo "文件已解压到: $OUTPUT_DIR"
ls -lh "$OUTPUT_DIR/"*.csv 2>/dev/null