#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
合并两台仪器的定量结果，计算相对活性、做统计、给出推荐排序。

用法
    python 03_analyze.py --config my_config.py

输出
    results/activity_by_replicate.csv   每个复孔的活性值（长表）
    results/activity_summary.csv        各指标 × 各构建：mean/SD/相对活性/p
    results/final_ranking.csv           推荐排序（按各组 primary_instrument）
    results/qc_signal.csv               信号质量

定义
----
relative activity (%) = (mean_x - mean_ctrl) / (mean_best - mean_ctrl) × 100
  · 对照 = 同组阴性对照
  · best = 组内均值最高的变体
统计：单因素 ANOVA；各变体 vs 对照用 Welch t 检验（不假设等方差）。
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

# (显示名, 来源, 列名)
METRICS = [
    ("EVOS reporter_index", "EVOS", "reporter_index"),
    ("EVOS gfp_in_rfp_frac", "EVOS", "gfp_in_rfp_frac"),
    ("EVOS gfp_pos_frac", "EVOS", "gfp_pos_frac"),
    ("Spark reporter_index", "Spark", "spark_reporter_index"),
    ("Spark gfp_in_red_area", "Spark", "spark_gfp_in_red_frac"),
]


def load(cfg):
    e = pd.read_csv(os.path.join(cfg["output_dir"], "evos_well_level.csv"))
    s = pd.read_csv(os.path.join(cfg["output_dir"], "spark_well_level.csv"))
    s["spark_gfp_in_red_frac"] = (
        np.minimum(s["spark_green_area_frac"], s["spark_red_area_frac"])
        / s["spark_red_area_frac"].replace(0, np.nan))
    df = e.merge(s, on=["well", "group", "label", "column"], how="outer")
    df["gfp_pos_frac"] = df["n_gfp"] / df["n_rfp"].replace(0, np.nan)
    return df


def summarize(df, cfg, group, metric):
    g = cfg["groups"][group]
    ctrl = g["control"]
    labels = common.ordered_labels(cfg, group)
    sub = df[df["group"] == group]

    by = {}
    for lab in labels:
        v = sub.loc[sub["label"] == lab, metric].astype(float).dropna().values
        if v.size:
            by[lab] = v
    if ctrl not in by:
        return pd.DataFrame()

    present = [l for l in labels if l in by]
    anova_p = np.nan
    if len(present) >= 2:
        _, anova_p = stats.f_oneway(*[by[l] for l in present])

    ctrl_mean = by[ctrl].mean()
    mutants = [l for l in present if l != ctrl]
    best = max(mutants, key=lambda l: by[l].mean())
    denom = by[best].mean() - ctrl_mean

    rows = []
    for lab in present:
        v = by[lab]
        p = np.nan
        if lab != ctrl and len(v) > 1 and len(by[ctrl]) > 1:
            _, p = stats.ttest_ind(v, by[ctrl], equal_var=False)
        rows.append({
            "group": group, "label": lab, "n": len(v),
            "mean": v.mean(),
            "sd": v.std(ddof=1) if len(v) > 1 else np.nan,
            "cv_pct": 100 * v.std(ddof=1) / v.mean() if len(v) > 1 and v.mean() else np.nan,
            "fold_over_control": v.mean() / ctrl_mean if ctrl_mean else np.nan,
            "relative_activity_pct":
                100 * (v.mean() - ctrl_mean) / denom if denom > 0 else np.nan,
            "t_p_vs_control": p,
            "anova_p": anova_p,
            "control_label": ctrl, "control_mean": ctrl_mean,
            "best_label": best,
        })
    out = pd.DataFrame(rows)
    rankable = out["label"] != ctrl
    out.loc[rankable, "rank"] = out.loc[rankable, "mean"] \
        .rank(ascending=False).astype("Int64")
    return out


def main():
    args = common.parse_args(__doc__)
    cfg = common.load_config(args.config)
    df = load(cfg)

    # 逐复孔长表
    rec = []
    for _, r in df.iterrows():
        for name, src, col in METRICS:
            rec.append({"group": r["group"], "label": r["label"], "well": r["well"],
                        "source": src, "metric": name, "value": r.get(col)})
    pd.DataFrame(rec).to_csv(
        os.path.join(cfg["output_dir"], "activity_by_replicate.csv"),
        index=False, encoding="utf-8-sig")

    # 汇总
    frames = []
    for group in cfg["groups"]:
        for name, src, col in METRICS:
            if col not in df.columns:
                continue
            t = summarize(df, cfg, group, col)
            if not t.empty:
                t.insert(0, "metric", name)
                t.insert(1, "source", src)
                frames.append(t)
    summ = pd.concat(frames, ignore_index=True)
    summ.to_csv(os.path.join(cfg["output_dir"], "activity_summary.csv"),
                index=False, encoding="utf-8-sig")

    # 推荐排序：每组用 primary_instrument 的主力指标
    primary_metric = {"EVOS": "EVOS reporter_index",
                      "Spark": "Spark reporter_index"}
    recs = []
    for group, g in cfg["groups"].items():
        m = primary_metric[g.get("primary_instrument", "EVOS")]
        recs.append(summ[(summ["group"] == group) & (summ["metric"] == m)])
    pd.concat(recs, ignore_index=True).to_csv(
        os.path.join(cfg["output_dir"], "final_ranking.csv"),
        index=False, encoding="utf-8-sig")

    # QC
    qc_cols = ["well", "group", "label", "n_gfp", "n_rfp",
               "mean_rfp_per_cell", "spark_green_count", "spark_red_count",
               "spark_green_area_frac", "spark_red_area_frac"]
    df[[c for c in qc_cols if c in df.columns]].to_csv(
        os.path.join(cfg["output_dir"], "qc_signal.csv"),
        index=False, encoding="utf-8-sig")

    # 打印
    for group in cfg["groups"]:
        print(f"\n========== {group} ==========")
        for name, src, col in METRICS:
            t = summ[(summ["group"] == group) & (summ["metric"] == name)]
            if t.empty:
                continue
            print(f"\n[{name}]  ANOVA p={t['anova_p'].iloc[0]:.3g}")
            for _, r in t.iterrows():
                p = r["t_p_vs_control"]
                p_str = f"p={p:.3f}" if pd.notna(p) else "对照"
                rel = r["relative_activity_pct"]
                rel_str = f"{rel:6.1f}%" if pd.notna(rel) else "   n/a "
                print(f"  {r['label']:16s} n={r['n']}  "
                      f"{r['mean']:.4g} ± {r['sd']:.3g}  "
                      f"fold={r['fold_over_control']:.2f}  "
                      f"rel={rel_str}  {p_str}")

    # 仪器一致性
    print("\n========== 仪器一致性 ==========")
    for group in cfg["groups"]:
        s = df[df["group"] == group].dropna(
            subset=["reporter_index", "spark_reporter_index"])
        if len(s) > 3:
            rho, p = stats.spearmanr(s["reporter_index"],
                                     s["spark_reporter_index"])
            verdict = "互证可靠" if rho > 0.9 else (
                "一致性一般，需人工判断" if rho > 0.5 else "其中一台不可用")
            print(f"{group}: Spearman ρ={rho:.3f}, P={p:.2g} (n={len(s)}) → {verdict}")

    print(f"\n已写入 {cfg['output_dir']}")


if __name__ == "__main__":
    main()
