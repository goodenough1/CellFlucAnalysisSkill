#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
由结果 CSV 生成 LaTeX 附录表格（逐孔原始数据，longtable 跨页）。

用法
    python 05_latex_tables.py --config my_config.py

输出
    <config.report_dir>/well_table.tex
    <config.report_dir>/summary_tables.tex   （各组汇总表）
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common


def sci(x, digits=2):
    """7.98e-05 → $7.98\\times10^{-5}$"""
    if pd.isna(x):
        return "—"
    if x == 0:
        return "0"
    mant, exp = f"{x:.{digits}e}".split("e")
    return rf"${mant}\times10^{{{int(exp)}}}$"


def esc(s):
    return str(s).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


def well_table(df, cfg):
    """逐孔原始数据。两组各插一行组名。"""
    body = []
    first = True
    for group in list(cfg["groups"]) + ["Control"]:
        sub = df[df["group"] == group]
        if sub.empty:
            continue
        if not first:
            body.append(r"\midrule")
        first = False
        body.append(rf"\multicolumn{{6}}{{l}}{{\textit{{{esc(group)} 组}}}}\\")
        for _, r in sub.iterrows():
            body.append(
                f"{r['well']} & {esc(r['label'])} & {sci(r.get('reporter_index'))} & "
                f"{sci(r.get('spark_reporter_index'))} & "
                f"{int(r['n_gfp']) if pd.notna(r.get('n_gfp')) else '—'} & "
                f"{int(r['n_rfp']) if pd.notna(r.get('n_rfp')) else '—'} \\\\")

    tex = rf"""% 由 05_latex_tables.py 自动生成，请勿手改
{{\small\setlength{{\tabcolsep}}{{5pt}}%
\begin{{longtable}}{{@{{}}l l r r r r@{{}}}}
\caption{{逐孔原始定量数据}}\label{{tab:wells}}\\
\toprule
孔位 & 构建 & EVOS 报告/内参 & Spark 报告/内参 & 报告$^+$ 细胞 & 内参$^+$ 细胞\\
\midrule
\endfirsthead
\toprule
孔位 & 构建 & EVOS 报告/内参 & Spark 报告/内参 & 报告$^+$ 细胞 & 内参$^+$ 细胞\\
\midrule
\endhead
\bottomrule
\endfoot
{chr(10).join(body)}
\end{{longtable}}}}
"""
    return tex


def summary_tables(summ, cfg):
    """每组一张「构建 / 均值±SD / 相对活性 / 倍数 / n」表。"""
    out = []
    for group in cfg["groups"]:
        primary = {"EVOS": "EVOS reporter_index",
                   "Spark": "Spark reporter_index"}[
                       cfg["groups"][group].get("primary_instrument", "EVOS")]
        t = summ[(summ["group"] == group) & (summ["metric"] == primary)]
        if t.empty:
            continue
        rows = []
        for _, r in t.iterrows():
            rows.append(
                f"{esc(r['label'])} & {sci(r['mean'])} $\\pm$ {sci(r['sd'])} & "
                f"{r['relative_activity_pct']:.1f}\\% & "
                f"{r['fold_over_control']:.2f} & {int(r['n'])} \\\\")
        out.append(rf"""
\begin{{center}}
\captionof{{table}}{{{esc(group)} 变体活性（$n=\text{{rep}}$）。}}
\label{{tab:{group}}}
\small
\begin{{tabular}}{{@{{}}lcccc@{{}}}}
\toprule
构建 & 均值 $\pm$ 标准差 & 相对活性 & 倍数 vs 对照 & n\\
\midrule
{chr(10).join(rows)}
\bottomrule
\end{{tabular}}
\end{{center}}
""")
    return "\n".join(out)


def main():
    args = common.parse_args(__doc__)
    cfg = common.load_config(args.config)
    report_dir = cfg.get("report_dir", "./report")
    os.makedirs(report_dir, exist_ok=True)

    e = pd.read_csv(os.path.join(cfg["output_dir"], "evos_well_level.csv"))
    s = pd.read_csv(os.path.join(cfg["output_dir"], "spark_well_level.csv"))
    summ = pd.read_csv(os.path.join(cfg["output_dir"], "activity_summary.csv"))
    df = e.merge(s.drop(columns=["group", "label", "column"]),
                 on="well", how="outer")

    with open(os.path.join(report_dir, "well_table.tex"), "w",
              encoding="utf-8") as f:
        f.write(well_table(df, cfg))
    with open(os.path.join(report_dir, "summary_tables.tex"), "w",
              encoding="utf-8") as f:
        f.write(summary_tables(summ, cfg))
    print(f"已写入 {report_dir}/well_table.tex 与 summary_tables.tex")


if __name__ == "__main__":
    main()
