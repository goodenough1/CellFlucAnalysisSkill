---
name: fluo-plate-mutant-screen
description: 分析「双荧光报告系统」的板式突变体筛选实验（EVOS 荧光显微图 + Spark/Tecan 板式全孔图），以 EGFP/mCherry（报告基因/内参）比值计算各突变体的相对活性、排序并出发表级图表与 LaTeX 报告。当用户提到荧光板筛选、转染后激发结果、GFP/mCherry 比值、突变体活性比较、EVOS+Spark 数据、或者要把这类原始目录整理成统计结果与论文图时使用。触发词：荧光筛选、报告基因、内参、mCherry、EGFP、EVOS、Spark、突变体活性、relative activity、plate screen、reporter assay。
---

# 双荧光报告系统的突变体筛选分析

把一个「原始数据目录」变成「可发表的结论」：定量 → 排序 → 统计 → 图 → 报告。
本 skill 来自 2026-11-15 OT2 磁珠法提取质粒 / 转染 20 h 的 EL222-PhyB 突变体筛查，
脚本已参数化，换一批数据只需改 `config`。

## 什么时候用

- 有 EVOS（或同类显微镜）拍的多孔位 GFP + mCherry 图像，和/或
- 有 Spark/Tecan 的荧光板读数或全孔成像，
- 想回答「哪个突变体活性强、强多少、统计上是否显著」。

**不适用**：流式数据、qPCR、纯生化酶活；这些没有图像通道的配对关系。

## 分析前的四件事（必做）

1. **读孔位图，不要猜。** 让用户提供板图或实验记录表格（xlsx/tsv），确认：
   哪些列是实验组、哪些是重复、WT/阴性对照在哪、每孔转染了哪几个质粒。
   本类实验常见布局：`列 2-4 / 5-7 / 8-10` 各为一组的三重复，
   最后三列放单质粒对照（EGFP-only、mCherry-only）。
2. **确认哪台仪器可信。** 不要默认两台都能用。先看信号：
   - 对象数 < 5/孔 → 该仪器对这组不可用；
   - 2× 物镜对弱信号不友好（单个细胞只有几个像素）。
   用 `03_analyze` 里输出的「EVOS vs Spark Spearman」作为客观判据。
3. **确认通道定义。** 报告基因 = 需要被诱导的那个（通常是 EGFP 在 Green 通道），
   内参 = 组成型表达（通常是 mCherry 在 Red 通道）。
4. **先看提取后质粒 QC。** 浓度 + A260/280 + A260/230，并算每个质粒取 100 ng
   需要多少 µL。低浓度（<20 ng/µL）会让取样体积变成 10 µL 级，稀释 DMEM、
   带入胍盐/蛋白残留，该孔的低信号是**转染失败造成的假阴性**，不能判为失活。
   判据：A260/280 1.7–2.0；A260/230 ≥1.5 可接受、<1.5 边缘、
   <0.5 或浓度 <20 ng/µL 为不合格。

## 流程

```
原始目录（只读）
   │
   ├─ 01_quant_evos.py    → results/evos_well_level.csv      逐孔 EVOS 定量
   ├─ 02_quant_spark.py   → results/spark_well_level.csv     逐孔 Spark 定量
   ├─ 03_analyze.py       → results/activity_summary.csv     相对活性 + 统计
   │                        results/final_ranking.csv        最终排序
   │                        results/qc_signal.csv            信号质量
   ├─ 04_figures_pub.py   → results/figures/pub/Fig*.pdf     发表级图
   └─ 05_latex_tables.py  → report/well_table.tex            附录表
                            templates/report.tex → 编译 PDF
```

### 0. 建 config

复制 `scripts/config_template.py`，按真实实验改：
孔位映射、路径、通道名、像素尺寸、空孔行（用于 Spark 阈值）。

关键字段：
- `groups.<名>.columns` / `.labels`（行号→构建名）/ `.control` / `.primary_instrument`
- `evos_glob`：本实验是 `{well}_*_{channel}.tiff`（中间那串是图像流水号）
- `spark_image_glob`：指向 **Processed** 图，见下

### 1. EVOS 定量

```bash
python scripts/01_quant_evos.py --config my_config.py
```

算法：大窗口中值背景 → 高斯平滑 → 4σ 分割 → 积分。
**先检查曝光/增益是否全孔一致**（脚本会打印），不一致时脚本会按
`image/(exposure*gain)` 归一化——但这只是补救，最好的做法是拍摄时固定参数。

指标：
- `reporter_index = ΣGFP/ΣmCherry`（主指标）
- `gfp_in_rfp_frac`：mCherry⁺ 细胞里被点亮的面积比（稳健，抗单点伪影）
- `n_gfp` / `n_rfp`：细胞数（判断统计力够不够）

### 2. Spark 定量

```bash
python scripts/02_quant_spark.py --config my_config.py
```

**必须用 `ProcessingPlugin/Images/*_Processed.tiff`**，那是做过串色矫正＋背景扣除的；
原始 `Images/` 里 Red 会漏进 Green（实测可达 20%），直接算比值会把
「只表达 mCherry 的细胞」判成 GFP⁺。

Processed 图仍有正噪声本底，所以阈值取**空孔分布 p99.99**，只积分超阈值像素。
同时导出仪器自带 Counting 的对象数/荧光和做交叉验证。

### 3. 统计与排序

```bash
python scripts/03_analyze.py --config my_config.py
```

- 相对活性：`(mean - ctrl)/(best - ctrl) × 100`
- 单因素 ANOVA + Welch t vs 同组对照
- **多指标一致性检查**：若三个指标排序一致，结论才可信
- 仪器间 Spearman：ρ > 0.9 才算互证；ρ < 0.5 说明其中一台不可用

### 4. 发表级图

```bash
python scripts/04_figures_pub.py --config my_config.py
```

输出矢量 PDF + 600 dpi PNG，规格见 `references/pitfalls.md` 的「出图」节。
每张图含：柱（mean ± s.d.）、复孔散点、显著性星号、代表性视野。

**出图三条铁律：**
1. **代表视野必须来自本组的主判据仪器。** 用 Spark 的积分下结论，就显示 Spark
   的 2× 全孔图；用 EVOS 的 10× 下结论，就显示 EVOS 视野。`montage()` 默认读
   `primary_instrument`，不要手动混搭。
2. **2× 全孔图必须去噪后再展示。** 原始 Processed 图散粒噪声很大，直接拉伸
   看不清。用 top-hat（核 ~13 px）+ 4σ 掩膜叠加在压暗的原图上；图注写明
   "仅用于展示，定量用原始积分"。2× 没有可信像素标定，不要画 100 µm 标尺。
3. **对照按生物学身份命名。** `pBApo-EL222 → EL222-WT`，`pPK-351 → PhyB-WT`；
   不要写 "Vector" 或 "Apo-EL222"。表格里保留质粒编号时，在首次出现处注明映射。

### 5. LaTeX 报告

```bash
python scripts/05_latex_tables.py --config my_config.py
cp templates/report.tex report/
cd report && latexmk -xelatex report.tex
```

`templates/report.tex` 已经处理好了中文排版、浮动体、字体缺字、表格宽度等坑，
直接改文字即可。

## 结果怎么读（经验判据）

| 现象 | 含义 | 处理 |
|---|---|---|
| 对照孔比值≈0，突变体高 | 正常，有动态范围 | 用倍数 + 相对活性报告 |
| 所有孔都≈0 | 转染失败 / 没诱导 / 报告质粒不对 | 查单质粒 EGFP 对照是否有信号 |
| 对照组本身有明显信号 | 报告质粒基础渗漏 | 以它为基线归一化，并在讨论中说明 |
| 突变体低于对照 | 表达/功能受损 | 复核序列，别直接判「无活性」 |
| 复孔 CV > 50% | 动态范围太小或细胞太少 | 结论只作方向性，补充实验 |
| 两台仪器排序矛盾 | 先看对象数，弱的那台通常不可信 | 用强信号仪器做判据 |

## 常见陷阱（详见 references/pitfalls.md）

1. Spark 不用 Processed 图 → 串色污染。
2. 用全图积分而不是阈值积分 → 被噪声本底淹没。
3. 2× 物镜对象数 0–2 就宣布「无活性」→ 其实是分辨率不够。
4. 单个亮斑伪影把对照均值抬高 → 用中位数指标或面积比交叉验证。
5. 只报一个指标 → 至少三个独立指标同向才下结论。
6. 图片展示时的对比度拉伸被当成定量 → 定量永远用原始值，展示时才拉伸并注明。

## 交付物清单

```
results/
├── evos_well_level.csv          逐孔 EVOS 原始量
├── spark_well_level.csv         逐孔 Spark 原始量
├── activity_summary.csv         各指标 × 各构建的均值/SD/相对活性/p
├── final_ranking.csv            推荐排序
├── qc_signal.csv                信号质量（判断可信度）
├── 活性分析汇总.xlsx            汇总工作簿
└── figures/pub/Fig*.pdf|png     发表级图
report/report.tex + report.pdf   报告
```

## 参考文件

- `references/pitfalls.md` — 所有踩过的坑、参数取值、出图与 LaTeX 规格
- `references/methods.md` — 定量公式、统计方法、为什么这样选
- `scripts/config_template.py` — 配置模板（从它开始）
- `templates/report.tex` — 中文报告模板
