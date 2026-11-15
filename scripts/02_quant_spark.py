#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spark 全孔荧光定量（报告基因 / 内参）。

用法
    python 02_quant_spark.py --config my_config.py

输出
    results/spark_well_level.csv  每孔一行

关键
----
1. 必须使用 ProcessingPlugin/Images/*_Processed.tiff：
   SparkControl 已在这些图上做了串色矫正（Green 扣 20% Red）和背景扣除。
   原始 Images/ 未矫正，Red 会漏进 Green，直接算比值会严重高估 GFP。
2. Processed 图仍有正噪声本底（约一半像素非零），不能整图积分；
   以空孔分布的 99.99 百分位为阈值，只积分超阈值像素。
3. 同时解析仪器自带 Counting 结果（对象数 / 荧光和）作交叉验证。
"""

import os
import sys

import numpy as np
import pandas as pd
import openpyxl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common


def read_instrument_counts(cfg):
    """解析 SparkControl 导出的 xlsx 中 Counting 各区块。"""
    path = os.path.join(cfg["spark_dir"], cfg["spark_xlsx"])
    if not os.path.exists(path):
        print(f"[warn] 找不到 {path}，跳过仪器自带计数")
        return {}

    ws = openpyxl.load_workbook(path, data_only=True).active
    titles = {}
    for r in range(1, ws.max_row + 1):
        v = ws.cell(r, 1).value
        if isinstance(v, str):
            titles[v.strip()] = r

    wanted = {
        "spark_green_count": "Counting 1 - Green object count",
        "spark_red_count": "Counting 1 - Red object count",
        "spark_green_sum": "Counting 1 - Green fluorescence - Sum [RFU]",
        "spark_red_sum": "Counting 1 - Red fluorescence - Sum [RFU]",
        "spark_green_avg": "Counting 1 - Green fluorescence - Average [RFU]",
        "spark_red_avg": "Counting 1 - Red fluorescence - Average [RFU]",
    }
    out = {}
    for key, title in wanted.items():
        r0 = titles.get(title)
        if r0 is None:
            continue
        header = r0 + 1                     # 表头行：<>,1,2,...
        vals = {}
        for rr in range(header + 1, header + 9):
            row = ws.cell(rr, 1).value
            if row not in list("ABCDEFGH"):
                continue
            for c in range(2, 14):
                col = ws.cell(header, c).value
                try:
                    col = int(col)
                except (TypeError, ValueError):
                    continue
                v = ws.cell(rr, c).value
                if isinstance(v, (int, float)):
                    vals[f"{row}{col}"] = float(v)
        out[key] = vals
    return out


def main():
    args = common.parse_args(__doc__)
    cfg = common.load_config(args.config)
    thr = common.spark_blank_threshold(cfg)
    print("阳性阈值（空孔分位）：" +
          "  ".join(f"{k}={v:.0f}" for k, v in thr.items()))

    inst = read_instrument_counts(cfg)
    ch_rep = cfg["spark_channels"]["reporter"]
    ch_int = cfg["spark_channels"]["internal"]

    rows = []
    for well, _, _ in common.all_wells(cfg):
        try:
            g = common.load_spark_processed(cfg, well, ch_rep)
            r = common.load_spark_processed(cfg, well, ch_int)
        except FileNotFoundError:
            print(f"[skip] {well}: 缺少 Processed 图")
            continue

        g_mask = g > thr["reporter"]
        r_mask = r > thr["internal"]
        sum_g = float(g[g_mask].sum())
        sum_r = float(r[r_mask].sum())
        group, label, col = common.group_of(cfg, well)

        row = {
            "well": well, "group": group, "label": label, "column": col,
            "spark_green_area_frac": float(g_mask.mean()),
            "spark_red_area_frac": float(r_mask.mean()),
            "spark_int_green": sum_g,
            "spark_int_red": sum_r,
            "spark_reporter_index": sum_g / sum_r if sum_r > 0 else np.nan,
        }
        for k, d in inst.items():
            row[k] = d.get(well, np.nan)
        rows.append(row)

    df = pd.DataFrame(rows)
    out = os.path.join(cfg["output_dir"], "spark_well_level.csv")
    df.to_csv(out, index=False, encoding="utf-8-sig")

    cols = ["well", "group", "label", "spark_green_area_frac",
            "spark_int_green", "spark_int_red", "spark_reporter_index",
            "spark_green_count"]
    cols = [c for c in cols if c in df.columns]
    print(df[cols].to_string(index=False))
    print(f"\n已写入 {out}")
    print("提示：若 spark_green_count 普遍 < 5，说明 2× 物镜分辨不了这组信号，"
          "该仪器对本组不可用，改用 EVOS 10×。")


if __name__ == "__main__":
    main()
