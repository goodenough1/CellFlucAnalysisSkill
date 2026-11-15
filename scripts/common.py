#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""各脚本共用的工具：读配置、读图、背景扣除、分割、孔位遍历。"""

import argparse
import glob
import importlib.util
import json
import os
import re

import numpy as np
from PIL import Image
from scipy import ndimage as ndi


# ---------------------------------------------------------------- 配置
def load_config(path):
    spec = importlib.util.spec_from_file_location("expcfg", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    cfg = mod.CONFIG
    os.makedirs(cfg["output_dir"], exist_ok=True)
    return cfg


def parse_args(desc):
    ap = argparse.ArgumentParser(description=desc)
    ap.add_argument("--config", default="my_config.py",
                    help="实验配置文件（见 config_template.py）")
    return ap.parse_args()


# ---------------------------------------------------------------- 孔位
def all_wells(cfg):
    """实验组 + 对照组的全部孔位。"""
    wells = []
    for rows in cfg["groups"].values():
        for r in cfg["rows"]:
            for c in rows["columns"]:
                wells.append((f"{r}{c}", r, c))
    for r in cfg["rows"]:
        for c in cfg["controls_columns"]:
            wells.append((f"{r}{c}", r, c))
    seen, out = set(), []
    for w in wells:                       # 去重保序
        if w[0] not in seen:
            seen.add(w[0]); out.append(w)
    return out


def group_of(cfg, well):
    """返回 (组名, 构建名, 列号)；对照孔返回 ("Control", 名字, 列)。"""
    r, c = well[0], int(well[1:])
    for name, g in cfg["groups"].items():
        if c in g["columns"]:
            return name, g["labels"].get(r, f"{r}{c}"), c
    if c in cfg["controls_columns"]:
        return "Control", cfg["controls_labels"].get(r, f"{r}{c}"), c
    return None, None, c


def ordered_labels(cfg, group):
    """组内构建顺序：对照在前的行顺序。"""
    g = cfg["groups"][group]
    labels = [g["labels"][r] for r in cfg["rows"] if r in g["labels"]]
    ctrl = g["control"]
    return [ctrl] + [l for l in labels if l != ctrl]


def short(cfg, label):
    return cfg.get("short_names", {}).get(label, label)


# ---------------------------------------------------------------- EVOS 读图
def evos_path(cfg, well, channel):
    pat = os.path.join(cfg["evos_dir"], cfg["evos_glob"]).format(
        well=well, channel=channel)
    hits = sorted(glob.glob(pat))
    if not hits:
        raise FileNotFoundError(pat)
    return hits[0]


def load_evos(cfg, well, channel):
    """返回 (float32 图像, 元数据)。元数据从 OME ImageDescription 里抠。"""
    path = evos_path(cfg, well, channel)
    im = Image.open(path)
    arr = np.asarray(im.convert("L"), dtype=np.float32)
    desc = im.tag_v2.get(270, "") or ""
    m = re.search(r"\{.*\}", desc, re.S)
    meta = {}
    if m:
        try:
            meta = json.loads(m.group(0))
        except json.JSONDecodeError:
            meta = {}
    return arr, {
        "exposure_ms": meta.get("Exposure Time"),
        "gain": meta.get("Gain"),
        "path": path,
    }


# ---------------------------------------------------------------- 图像处理
def bg_subtract(a, downsample=4, kernel=25):
    """大窗口中值滤波估计照明背景并相减（比逐像素滚动球快很多）。"""
    small = a[::downsample, ::downsample]
    bg = ndi.median_filter(small, size=kernel, mode="nearest")
    bg = ndi.zoom(bg, downsample, order=1)[: a.shape[0], : a.shape[1]]
    return a - bg


def robust_sigma(x):
    med = np.median(x)
    return 1.4826 * np.median(np.abs(x - med)), med


def segment(sig, sigma=2.0, k=4.0, min_area=40, min_threshold=2.5):
    """高斯平滑 → kσ 阈值 → 开运算 → 去小连通域。

    返回 (mask, n_objects, labels, keep_ids)。
    """
    sm = ndi.gaussian_filter(sig, sigma)
    s, _ = robust_sigma(sm)
    thr = max(k * s, min_threshold)
    mask = ndi.binary_opening(sm > thr, iterations=1)
    lab, n = ndi.label(mask)
    if n == 0:
        return mask, 0, lab, np.array([], dtype=int)
    sizes = np.bincount(lab.ravel())[1:]
    keep = np.where(sizes >= min_area)[0] + 1
    if len(keep) == 0:
        return mask, 0, lab, np.array([], dtype=int)
    return np.isin(lab, keep), len(keep), lab, sizes[keep - 1]


def per_object_means(sig, lab, keep_ids):
    if len(keep_ids) == 0:
        return np.array([])
    sums = ndi.sum(sig, lab, keep_ids)
    area = ndi.sum(np.ones_like(sig), lab, keep_ids)
    return sums / np.maximum(area, 1)


# ---------------------------------------------------------------- Spark 读图
def spark_processed_path(cfg, well, channel):
    pat = os.path.join(cfg["spark_dir"], cfg["spark_processed_subdir"],
                       cfg["spark_image_glob"]).format(well=well, channel=channel)
    hits = sorted(glob.glob(pat))
    if not hits:
        raise FileNotFoundError(
            f"{pat}\n"
            "  → 确认用的是 ProcessingPlugin/Images 下的 *_Processed.tiff；"
            "原始 Images/ 未做串色矫正，不能用。")
    return hits[0]


def load_spark_processed(cfg, well, channel):
    return np.asarray(Image.open(spark_processed_path(cfg, well, channel)),
                      dtype=np.float32)


def spark_blank_threshold(cfg):
    """用空孔 Processed 图分布的分位数作阳性阈值。"""
    out = {}
    for role, ch in cfg["spark_channels"].items():
        vals = []
        for c in cfg["spark_blank_cols"]:
            try:
                a = load_spark_processed(cfg, f"{cfg['spark_blank_row']}{c}", ch)
                vals.append(a.ravel())
            except FileNotFoundError:
                pass
        if not vals:
            raise RuntimeError("找不到空孔图像，无法确定 Spark 阈值")
        out[role] = float(np.percentile(np.concatenate(vals),
                                        cfg.get("spark_threshold_pct", 99.99)))
    return out
