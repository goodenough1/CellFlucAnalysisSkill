#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
发表级图表：每组一张「三指标 + 代表性视野」主图，加一张一致性与一张总览。

用法
    python 04_figures_pub.py --config my_config.py

输出
    results/figures/pub/{组}_screen.pdf|png   每组主图
    results/figures/pub/validation.pdf|png    仪器一致性
    results/figures/pub/summary.pdf|png       相对活性总览

规格（Nature 参考）：双栏 183 mm，Arial 7 pt，线宽 0.5 pt，
矢量 PDF + 600 dpi PNG，Okabe–Ito 配色，误差线 = s.d.，显微图带标尺。
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import ndimage as ndi
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

MM = 1 / 25.4
W_DOUBLE = 183 * MM

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7,
    "axes.linewidth": 0.5, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.direction": "out", "ytick.direction": "out",
    "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.unicode_minus": False, "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

# 数据列名 → activity_summary.csv 里的 metric 名（用于取 p 值）
METRIC_ALIAS = {
    "reporter_index": "EVOS reporter_index",
    "gfp_in_rfp_frac": "EVOS gfp_in_rfp_frac",
    "spark_reporter_index": "Spark reporter_index",
    "spark_green_area_frac": "Spark reporter_index",
    "spark_green_area_x1e4": "Spark reporter_index",
}


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(pad=2)


def panel_letter(ax, letter, x=-0.22, y=1.10):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="left")


def p_stars(p):
    if pd.isna(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


# ---------------------------------------------------------------- 柱状图
def bars(ax, df, cfg, group, metric, ylabel, letter, summ=None):
    g = cfg["groups"][group]
    labels = common.ordered_labels(cfg, group)
    ctrl = g["control"]
    accent = cfg["group_colors"].get(group, "#0072B2")

    means, sds, pts = [], [], []
    for lab in labels:
        v = df.loc[(df["group"] == group) & (df["label"] == lab),
                   metric].astype(float).dropna().values
        means.append(v.mean() if v.size else np.nan)
        sds.append(v.std(ddof=1) if v.size > 1 else 0)
        pts.append(v)

    xs = np.arange(len(labels))
    colors = [cfg["control_color"] if l == ctrl else accent for l in labels]
    ax.bar(xs, means, width=0.62, color=colors, edgecolor="none",
           yerr=sds, capsize=2,
           error_kw={"elinewidth": 0.6, "ecolor": "#333333", "capthick": 0.6})
    for i, v in enumerate(pts):
        if len(v):
            j = np.linspace(-0.10, 0.10, len(v))
            ax.plot(np.full(len(v), xs[i]) + j, v, "o", ms=2.6,
                    mfc="white", mec="#222222", mew=0.5, zorder=3)

    ax.set_xticks(xs)
    ax.set_xticklabels([common.short(cfg, l) for l in labels], rotation=32,
                       ha="right", rotation_mode="anchor")
    ax.set_ylabel(ylabel)
    ax.yaxis.set_major_locator(MaxNLocator(4, prune="lower"))
    style(ax)

    top = max(m + (s if np.isfinite(s) else 0) for m, s in zip(means, sds))
    ax.set_ylim(0, top * 1.28)

    alias = METRIC_ALIAS.get(metric)
    if summ is not None and alias:
        sub = summ[(summ["group"] == group) & (summ["metric"] == alias)]
        pmap = dict(zip(sub["label"], sub["t_p_vs_control"]))
        for i, lab in enumerate(labels):
            s = p_stars(pmap.get(lab, np.nan))
            if s:
                ax.text(xs[i], means[i] + sds[i] + top * 0.045, s,
                        ha="center", va="bottom", fontsize=6.5)
    panel_letter(ax, letter)


# ---------------------------------------------------------------- 显微图
def display_img(cfg, well, channel):
    """去照明梯度 + 百分位拉伸；仅用于展示，定量用原始值。"""
    a, _ = common.load_evos(cfg, well, channel)
    small = a[::4, ::4]
    bg = ndi.gaussian_filter(small, 10, mode="nearest")
    a = a - ndi.zoom(bg, 4, order=1)[:a.shape[0], :a.shape[1]]
    v0, v1 = np.percentile(a, [1.0, 99.9])
    return np.clip((a - v0) / (v1 - v0 + 1e-6), 0, 1)


def add_scale_bar(ax, cfg):
    um = cfg.get("scale_bar_um", 100)
    px = um / cfg["pixel_size_um"]
    y = ax.get_ylim()[0] * 0.93
    x0 = ax.get_xlim()[1] * 0.08
    ax.plot([x0, x0 + px], [y, y], color="white", lw=1.4,
            solid_capstyle="butt")
    ax.text(x0 + px / 2, y - ax.get_ylim()[0] * 0.03, f"{um} µm",
            color="white", fontsize=5.5, ha="center", va="bottom")


def montage(fig, sub_gs, cfg, group):
    g = cfg["groups"][group]
    col0 = cfg.get("montage_columns", {}).get(group, g["columns"][0])
    rows = [r for r in cfg["rows"] if r in g["labels"]]
    wells = [f"{r}{col0}" for r in rows]
    y0, y1, x0, x1 = cfg.get("montage_crop", [300, 1200, 380, 1580])
    ch_rep = cfg["evos_channels"]["reporter"]
    ch_int = cfg["evos_channels"]["internal"]

    first = None
    for c, well in enumerate(wells):
        for r, ch in enumerate([ch_rep, ch_int]):
            ax = fig.add_subplot(sub_gs[r, c])
            ax.imshow(display_img(cfg, well, ch)[y0:y1, x0:x1],
                      cmap="gray", vmin=0, vmax=1, aspect="equal")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            if r == 0:
                ax.set_title(common.short(cfg, g["labels"][rows[c]]),
                             fontsize=6.5, pad=1.5)
            if c == 0:
                ax.set_ylabel("Reporter" if r == 0 else "Internal",
                              fontsize=6.5, labelpad=1)
            if first is None:
                first = ax
    add_scale_bar(first, cfg)


# ---------------------------------------------------------------- 各组主图
def group_figure(df, summ, cfg, group, panels, out):
    fig = plt.figure(figsize=(W_DOUBLE, 4.9))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.05], hspace=0.35)
    top = outer[0].subgridspec(1, len(panels), wspace=0.55)
    for i, (metric, ylabel, title) in enumerate(panels):
        ax = fig.add_subplot(top[0, i])
        bars(ax, df, cfg, group, metric, ylabel, "abc"[i], summ=summ)
        ax.set_title(title, fontsize=7, pad=4)
    bottom = outer[1].subgridspec(2, len(cfg["rows"]), wspace=0.04, hspace=0.04)
    montage(fig, bottom, cfg, group)
    fig.text(0.5, 0.012,
             "Bars, mean; error bars, s.d.; n = replicate wells; "
             "*P < 0.05, **P < 0.01, ***P < 0.001 (Welch t vs control).",
             ha="center", fontsize=6, color="#333333")
    save(fig, out)


# ---------------------------------------------------------------- 一致性
def validation_figure(df, cfg, out):
    fig = plt.figure(figsize=(W_DOUBLE * 0.82, 2.5))
    gs = fig.add_gridspec(1, 2, wspace=0.38)

    ax = fig.add_subplot(gs[0, 0])
    best = None
    for group in cfg["groups"]:
        d = df[df["group"] == group].dropna(
            subset=["reporter_index", "spark_reporter_index"]).copy()
        d = d[(d["reporter_index"] > 0) & (d["spark_reporter_index"] > 0)]
        if len(d) < 3:
            continue
        ax.scatter(d["reporter_index"], d["spark_reporter_index"], s=22,
                   color=cfg["group_colors"].get(group, "#0072B2"),
                   edgecolors="white", linewidths=0.4, zorder=3)
        rho, _ = stats.spearmanr(d["reporter_index"], d["spark_reporter_index"])
        if best is None or rho > best[0]:
            best = (rho, group, d)
    if best:
        rho, group, d = best
        lx = np.log10(d["reporter_index"])
        ly = np.log10(d["spark_reporter_index"])
        k, b = np.polyfit(lx, ly, 1)
        xx = np.linspace(lx.min(), lx.max(), 50)
        ax.plot(10 ** xx, 10 ** (k * xx + b), color="#333333", lw=0.7, ls="--")
        ax.text(0.04, 0.95, f"{group}: ρ = {rho:.2f}", transform=ax.transAxes,
                va="top", ha="left", fontsize=6.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("EVOS  reporter / internal")
    ax.set_ylabel("Spark  reporter / internal")
    ax.set_title("Quantification agreement", fontsize=7, pad=4)
    style(ax); panel_letter(ax, "a", x=-0.28)

    ax = fig.add_subplot(gs[0, 1])
    groups = list(cfg["groups"])
    data = [df.loc[df["group"] == g, "spark_green_count"].dropna().values
            for g in groups]
    bp = ax.boxplot(data, widths=0.42, showfliers=False, patch_artist=True,
                    medianprops={"color": "#222222", "lw": 0.7},
                    boxprops={"lw": 0.5, "edgecolor": "#333333"},
                    whiskerprops={"lw": 0.5}, capprops={"lw": 0.5})
    for patch, g in zip(bp["boxes"], groups):
        patch.set_facecolor(cfg["group_colors"].get(g, "#0072B2"))
        patch.set_alpha(0.35)
    for i, (g, v) in enumerate(zip(groups, data), start=1):
        if len(v):
            j = np.linspace(-0.12, 0.12, len(v))
            ax.plot(np.full(len(v), i) + j, v, "o", ms=2.4,
                    mfc=cfg["group_colors"].get(g, "#0072B2"),
                    mec="white", mew=0.3, zorder=3)
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels([f"{g}\n(n={len(v)})" for g, v in zip(groups, data)])
    ax.set_ylabel("Spark reporter$^+$ objects per well")
    ax.set_title("Spark sensitivity by group", fontsize=7, pad=4)
    style(ax); panel_letter(ax, "b", x=-0.30)
    save(fig, out)


# ---------------------------------------------------------------- 总览
def summary_figure(df, cfg, out):
    groups = list(cfg["groups"])
    fig, axes = plt.subplots(1, len(groups), figsize=(W_DOUBLE * 0.86, 2.2),
                             gridspec_kw={"wspace": 0.62}, squeeze=False)
    for idx, (ax, (group, g)) in enumerate(zip(axes[0], cfg["groups"].items())):
        metric = ("reporter_index" if g.get("primary_instrument") == "EVOS"
                  else "spark_reporter_index")
        ctrl = g["control"]
        ctrl_mean = df.loc[(df["group"] == group) & (df["label"] == ctrl),
                           metric].astype(float).mean()
        means = df[df["group"] == group].groupby("label")[metric].mean()
        best_mean = means.drop(ctrl).max()

        labs = [l for l in common.ordered_labels(cfg, group)
                if l != ctrl][::-1]
        mu, sd, pts = [], [], []
        for lab in labs:
            v = df.loc[(df["group"] == group) & (df["label"] == lab),
                       metric].astype(float).values
            rel = 100 * (v - ctrl_mean) / (best_mean - ctrl_mean)
            mu.append(rel.mean()); sd.append(rel.std(ddof=1)); pts.append(rel)

        y = np.arange(len(labs))
        lo = np.clip(mu, 0, None)
        ax.barh(y, lo, height=0.55, color=cfg["group_colors"].get(group),
                xerr=[np.minimum(sd, lo), sd],
                error_kw={"elinewidth": 0.6, "ecolor": "#333333",
                          "capthick": 0.6, "capsize": 2})
        for i, v in enumerate(pts):
            j = np.linspace(-0.13, 0.13, len(v))
            ax.plot(np.clip(v, 0, None), np.full(len(v), y[i]) + j, "o", ms=2.6,
                    mfc="white", mec="#222222", mew=0.5, zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels([common.short(cfg, l) for l in labs])
        ax.set_xlim(0, 112)
        ax.set_xlabel("Relative activity (%)")
        ax.set_title(f"{group} ({g.get('primary_instrument', 'EVOS')})",
                     fontsize=7, pad=4)
        style(ax)
        panel_letter(ax, "abcd"[idx], x=-0.34)
    save(fig, out)


# ---------------------------------------------------------------- 保存
def save(fig, out_base):
    os.makedirs(os.path.dirname(out_base), exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_base}.{ext}", dpi=600)
    plt.close(fig)
    print("saved", out_base)


def main():
    args = common.parse_args(__doc__)
    cfg = common.load_config(args.config)
    out_dir = os.path.join(cfg["output_dir"], "figures", "pub")

    e = pd.read_csv(os.path.join(cfg["output_dir"], "evos_well_level.csv"))
    s = pd.read_csv(os.path.join(cfg["output_dir"], "spark_well_level.csv"))
    summ = pd.read_csv(os.path.join(cfg["output_dir"], "activity_summary.csv"))
    s["spark_green_area_x1e4"] = s["spark_green_area_frac"] * 1e4
    df = e.merge(s.drop(columns=["group", "label", "column"]),
                 on="well", how="outer")

    panels = {
        "EVOS": [("reporter_index", "Reporter / internal",
                  "Reporter / internal control"),
                 ("gfp_in_rfp_frac", "Reporter$^+$ area / internal$^+$ area",
                  "Fraction of transfected cells activated"),
                 ("n_gfp", "Reporter$^+$ cells per field",
                  "Reporter$^+$ cell count")],
        "Spark": [("spark_reporter_index", "Spark reporter / internal",
                   "Spark, reporter / internal control"),
                  ("reporter_index", "EVOS reporter / internal",
                   "EVOS, reporter / internal control"),
                  ("spark_green_area_x1e4",
                   "Reporter$^+$ area (×10$^{-4}$)",
                   "Spark positive area")],
    }
    for group, g in cfg["groups"].items():
        inst = g.get("primary_instrument", "EVOS")
        group_figure(df, summ, cfg, group, panels[inst],
                     os.path.join(out_dir, f"{group}_screen"))

    validation_figure(df, cfg, os.path.join(out_dir, "validation"))
    summary_figure(df, cfg, os.path.join(out_dir, "summary"))


if __name__ == "__main__":
    main()
