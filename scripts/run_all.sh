#!/usr/bin/env bash
# 一键跑完整个分析流程
#
#   cp config_template.py my_config.py   # 先改好配置
#   bash run_all.sh my_config.py
#
# 依赖：numpy pandas scipy pillow openpyxl matplotlib
#       （LaTeX 报告另需 TeX Live + ctex + xeCJK）

set -euo pipefail
CONFIG="${1:-my_config.py}"
HERE="$(cd "$(dirname "$0")" && pwd)"

if [[ ! -f "$CONFIG" ]]; then
  echo "找不到配置文件 $CONFIG —— 请先 cp config_template.py my_config.py 并修改"
  exit 1
fi

echo "==> [1/5] EVOS 定量"
python3 "$HERE/01_quant_evos.py" --config "$CONFIG"

echo "==> [2/5] Spark 定量"
python3 "$HERE/02_quant_spark.py" --config "$CONFIG"

echo "==> [3/5] 统计与排序"
python3 "$HERE/03_analyze.py" --config "$CONFIG"

echo "==> [4/5] 发表级图"
python3 "$HERE/04_figures_pub.py" --config "$CONFIG"

echo "==> [5/5] LaTeX 表格"
python3 "$HERE/05_latex_tables.py" --config "$CONFIG"

echo
echo "完成。结果在 results/，报告模板在 templates/report.tex"
echo "下一步：cp templates/report.tex report/ && cd report && latexmk -xelatex report.tex"
