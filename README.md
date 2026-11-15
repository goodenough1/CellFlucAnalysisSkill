# CellFlucAnalysisSkill

细胞荧光报告系统（EGFP 报告 / mCherry 内参）突变体筛选的分析 skill：
从原始图像目录 → 逐孔定量 → 相对活性排序 → 统计 → 发表级图 → LaTeX 报告。

来自 2026-11-15「EL222 / PhyB 光遗传元件突变体」筛查实战，脚本已参数化。

## 安装

把本目录放到 Claude Code 的 skills 目录即可：

```bash
cp -r CellFlucAnalysisSkill ~/.claude/skills/
```

之后在对话中提到「荧光板筛选 / GFP-mCherry 比值 / 突变体活性 / EVOS + Spark 数据」时，
Claude 会自动读取 `SKILL.md` 并按流程执行。

## 单独使用（不经过 Claude）

```bash
pip install numpy pandas scipy pillow openpyxl matplotlib

cp scripts/config_template.py my_config.py   # 改路径、孔位、通道名
bash scripts/run_all.sh my_config.py
```

需要 LaTeX 报告时：

```bash
cp templates/report.tex report/
cd report && latexmk -xelatex report.tex
```

## 目录

```
.
├── SKILL.md                    Skill 主文件（流程、判据、检查清单）
├── references/
│   ├── pitfalls.md             25 个踩过的坑（附当时的数据证据）
│   └── methods.md              定量公式、统计方法、仪器选择判据
├── scripts/
│   ├── config_template.py      配置模板（从这里开始）
│   ├── common.py               读图/背景扣除/分割/孔位遍历
│   ├── 01_quant_evos.py        EVOS 10× 逐孔定量
│   ├── 02_quant_spark.py       Spark 2× 定量（Processed 图）
│   ├── 03_analyze.py           相对活性、ANOVA、Welch t、排序
│   ├── 04_figures_pub.py       发表级图（PDF + 600 dpi PNG）
│   ├── 05_latex_tables.py      LaTeX 表格
│   └── run_all.sh              一键跑完
└── templates/
    └── report.tex              中文报告模板（ctexart + XeLaTeX）
```

## 输出

```
results/
├── evos_well_level.csv         EVOS 逐孔（n_gfp / int_gfp / reporter_index …）
├── spark_well_level.csv        Spark 逐孔（Green/Red 积分、对象数）
├── activity_summary.csv        各指标 × 各构建：mean / SD / 相对活性 / p
├── final_ranking.csv           推荐排序
├── qc_signal.csv               信号质量
└── figures/pub/*.pdf|png       发表级图
```

## 三条核心原则

1. **比值优先**：报告基因 / 内参，消掉转染效率与细胞数差异。
2. **先看信号再看排序**：对象数 < 5 的仪器不能用；至少三个指标同向才下结论。
3. **原始数据只读**：所有输出写到 `results/`，不碰原始目录。

## 环境

Python ≥ 3.9，依赖 numpy / pandas / scipy / pillow / openpyxl / matplotlib。
LaTeX 部分另需 TeX Live（ctex + xeCJK），用 XeLaTeX 编译。

## License

内部实验工具，未设 license。使用前请确认原始数据不出现在仓库中。
