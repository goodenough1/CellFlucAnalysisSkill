#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
实验配置模板 —— 复制成 my_config.py 后按真实实验修改。

所有脚本都通过 `--config my_config.py` 读取本文件的 CONFIG 字典。
唯一必须改的是 CONFIG；下面的注释是最好的填写说明。

本示例即 2026-11-15 EL222/PhyB 突变体筛查（可直接跑通）。
"""

CONFIG = {
    # ------------------------------------------------------------ 基本信息
    "name": "20261115 EL222/PhyB 突变体筛查",
    "date": "2026-11-15",
    "output_dir": "./results",          # 相对脚本运行目录

    # ------------------------------------------------------------ EVOS 图像
    "evos_dir": ("/Volumes/原始数据/EVOS/"
                 "20261115-OT2磁珠法质粒提取转染20h激发结果"),
    # {well} 会被替换成 "B2"，{channel} 替换成 GFP / RFP。
    # 很多 EVOS 导出的文件名形如 B2_9083_GFP.tiff，中间是流水号，用 * 吃掉。
    "evos_glob": "{well}_*_{channel}.tiff",
    "evos_channels": {"reporter": "GFP", "internal": "RFP"},
    # 若同一孔有多张图（多个视野），用列表 glob；只会取第一个匹配。
    # 多视野时建议改造成对 filename 列表循环再平均。

    # 图像处理参数（一般不用改）
    "evos_bg_downsample": 4,       # 背景估计降采样倍数
    "evos_bg_kernel": 25,          # 降采样后的中值窗口 → 100 px ≈ 62 µm
    "evos_smooth_sigma": 2.0,      # 分割前高斯平滑 (px)
    "evos_k_sigma": 4.0,           # 阈值 = k × 噪声 σ
    "evos_min_area": 40,           # 最小阳性面积 (px²)
    "pixel_size_um": 0.619225,     # EVOS 10× 的 µm/px，取自 OME 元数据

    # ------------------------------------------------------------ Spark
    "spark_dir": ("/Volumes/原始数据/Spark/"
                  "20261115-OT2提取质粒转染后激发20h/"
                  "dy-311 (Modified)_20261115_160024"),
    # 相对 spark_dir
    "spark_processed_subdir": "ProcessingPlugin/Images",
    "spark_xlsx": "Export/xlsx/dy-311 (Modified)_20261115_160820.xlsx",
    "spark_image_glob": "Measurement 1_{well}_R_{channel}_*_Processed.tiff",
    "spark_channels": {"reporter": "G", "internal": "R"},
    # 未铺细胞的空孔所在行 —— 用它们估计 Processed 图的噪声本底
    "spark_blank_row": "G",
    "spark_blank_cols": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
    "spark_threshold_pct": 99.99,   # 阳性阈值取的百分位

    # ------------------------------------------------------------ 孔位布局
    # rows：有数据的行；controls_columns：单质粒对照列
    "rows": ["B", "C", "D", "E", "F"],
    "controls_columns": [8, 9, 10],
    "controls_labels": {
        "B": "pBApo",
        "C": "EGFP-Kan-Bsd",
        "D": "mCherry-Kan-Bsd",
        "E": "pN1-016",
        "F": "pN1-017",
    },

    # 实验组：每组 = 若干列（重复）+ 行→构建名 + 阴性对照 + 主判据仪器
    "groups": {
        "EL222": {
            "columns": [2, 3, 4],
            "control": "pBApo-EL222",
            "primary_instrument": "EVOS",     # 依据：Spark 绿色对象数 0–2
            "labels": {
                "B": "pBApo-EL222",           # 阴性对照（无功能 EL222）
                "C": "EL222-AQTrip",
                "D": "EL222-E87D",
                "E": "EL222-L87I",
                "F": "EL222-G85R",
            },
        },
        "PhyB": {
            "columns": [5, 6, 7],
            "control": "pPK-351",             # 空载体
            "primary_instrument": "Spark",
            "labels": {
                "B": "pPK-351",
                "C": "PhyB-S32A",
                "D": "PhyB-E22P",
                "E": "PhyB-S73T",
                "F": "PhyB-E101L",
            },
        },
    },

    # ------------------------------------------------------------ 出图/报告
    # 图上用的短名（不填则用原名）
    "short_names": {
        "pBApo-EL222": "Apo-EL222", "EL222-AQTrip": "AQTrip",
        "EL222-E87D": "E87D", "EL222-L87I": "L87I", "EL222-G85R": "G85R",
        "pPK-351": "Vector", "PhyB-S32A": "S32A", "PhyB-E22P": "E22P",
        "PhyB-S73T": "S73T", "PhyB-E101L": "E101L",
    },
    # 每组配色（Okabe–Ito 色盲友好）
    "group_colors": {"EL222": "#D55E00", "PhyB": "#0072B2"},
    "control_color": "#999999",
    # 每组代表性视野用的孔（默认取该组第一列 × 每一行）
    "montage_columns": {"EL222": 2, "PhyB": 5},
    "montage_crop": [300, 1200, 380, 1580],   # y0, y1, x0, x1（像素）
    "scale_bar_um": 100,

    # 报告标题
    "report_title": "EL222 与 PhyB 光遗传元件突变体活性筛查",
    "report_subtitle": "OT2 磁珠法提取质粒、转染 20 h 激发实验数据分析",
}
