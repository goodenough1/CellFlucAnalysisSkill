#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EVOS 荧光图像逐孔定量（报告基因 / 内参）。

用法
    python 01_quant_evos.py --config my_config.py

输出
    results/evos_well_level.csv  每孔一行

要点
----
1. 先核对全部孔位的曝光时间与增益。一致 → 孔间直接可比；
   不一致 → 用 image/(exposure*gain) 归一化到参考值（脚本自动做并提示）。
2. 照明不均：大窗口中值滤波估计背景后相减，等效于地道的 rolling-ball。
3. 分割：高斯平滑 + 4σ 阈值 + 开运算 + 最小面积过滤。
4. 指标（比值用于消除转染效率差异）：
     reporter_index   = ΣGFP / ΣmCherry
     gfp_in_rfp_frac  = GFP⁺ 面积 / mCherry⁺（膨胀后）面积
     n_gfp / n_rfp    = 阳性细胞数
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common


def quantify(cfg, well):
    ch_rep = cfg["evos_channels"]["reporter"]
    ch_int = cfg["evos_channels"]["internal"]

    rep, m_rep = common.load_evos(cfg, well, ch_rep)
    inte, m_int = common.load_evos(cfg, well, ch_int)

    g = common.bg_subtract(rep, cfg["evos_bg_downsample"], cfg["evos_bg_kernel"])
    r = common.bg_subtract(inte, cfg["evos_bg_downsample"], cfg["evos_bg_kernel"])

    g_mask, n_g, g_lab, _ = common.segment(
        g, cfg["evos_smooth_sigma"], cfg["evos_k_sigma"], cfg["evos_min_area"])
    r_mask, n_r, r_lab, _ = common.segment(
        r, cfg["evos_smooth_sigma"], cfg["evos_k_sigma"], cfg["evos_min_area"])

    int_g = float(np.clip(g, 0, None)[g_mask].sum())
    int_r = float(np.clip(r, 0, None)[r_mask].sum())

    r_dil = ndi.binary_dilation(r_mask, iterations=5)   # 容忍通道位移
    r_area = float(r_dil.sum())
    gfp_in_rfp = float((g_mask & r_dil).sum()) / r_area if r_area else np.nan

    per_cell_g = common.per_object_means(np.clip(g, 0, None), g_lab,
                                         np.arange(1, n_g + 1))
    per_cell_r = common.per_object_means(np.clip(r, 0, None), r_lab,
                                         np.arange(1, n_r + 1))

    group, label, col = common.group_of(cfg, well)
    return {
        "well": well, "group": group, "label": label, "column": col,
        "n_gfp": n_g, "n_rfp": n_r,
        "int_gfp": int_g, "int_rfp": int_r,
        "reporter_index": int_g / int_r if int_r > 0 else np.nan,
        "gfp_in_rfp_frac": gfp_in_rfp,
        "mean_gfp_per_cell": float(per_cell_g.mean()) if per_cell_g.size else 0.0,
        "median_gfp_per_cell": float(np.median(per_cell_g)) if per_cell_g.size else 0.0,
        "mean_rfp_per_cell": float(per_cell_r.mean()) if per_cell_r.size else 0.0,
        "exposure_reporter": m_rep["exposure_ms"], "gain_reporter": m_rep["gain"],
        "exposure_internal": m_int["exposure_ms"], "gain_internal": m_int["gain"],
    }


def check_settings(df, cfg):
    """核对曝光/增益一致性；不一致时提示并给出归一化因子。"""
    for ch in ["reporter", "internal"]:
        e = df[f"exposure_{ch}"].astype(float)
        g = df[f"gain_{ch}"].astype(float)
        if e.isna().all():
            print(f"[QC] {ch}: 元数据里没有曝光/增益，无法核对")
            continue
        if e.nunique() == 1 and g.nunique() == 1:
            print(f"[QC] {ch}: 曝光/增益全孔一致 "
                  f"({e.iloc[0]:.1f} ms, gain {g.iloc[0]:.1f}) —— 可直接比较")
        else:
            print(f"[QC] {ch}: ⚠ 曝光/增益不一致 "
                  f"(exposure {e.min():.1f}–{e.max():.1f}, "
                  f"gain {g.min():.1f}–{g.max():.1f})。"
                  "如未归一化，孔间不可直接比较！")


def main():
    args = common.parse_args(__doc__)
    cfg = common.load_config(args.config)

    rows = []
    for well, _, _ in common.all_wells(cfg):
        try:
            rows.append(quantify(cfg, well))
        except FileNotFoundError:
            print(f"[skip] {well}: 缺少图像")
    df = pd.DataFrame(rows)

    # 曝光/增益不一致时按 exposure*gain 归一化到中位数乘积
    for ch, col in [("reporter", "reporter_index"), ("internal", "reporter_index")]:
        if df[f"exposure_{ch}"].notna().any():
            k = (df[f"exposure_{ch}"] * df[f"gain_{ch}"])
            ref = k.median()
            if k.nunique() > 1:
                df[f"scale_{ch}"] = ref / k
            else:
                df[f"scale_{ch}"] = 1.0
    if (df.get("scale_reporter", pd.Series([1])).nunique() > 1
            or df.get("scale_internal", pd.Series([1])).nunique() > 1):
        df["reporter_index"] = (df["int_gfp"] * df["scale_reporter"]) / \
                               (df["int_rfp"] * df["scale_internal"])
        df["gfp_in_rfp_frac"] = df["gfp_in_rfp_frac"]     # 比值本身无需缩放
        print("[QC] 已按曝光×增益归一化 reporter_index")

    check_settings(df, cfg)

    out = os.path.join(cfg["output_dir"], "evos_well_level.csv")
    df.to_csv(out, index=False, encoding="utf-8-sig")
    cols = ["well", "group", "label", "n_gfp", "n_rfp",
            "int_gfp", "int_rfp", "reporter_index"]
    print(df[cols].to_string(index=False))
    print(f"\n已写入 {out}")
    print("提示：n_gfp < 10 的孔统计力很弱，解读时需谨慎。")


if __name__ == "__main__":
    main()
