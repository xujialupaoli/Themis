#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Themis helper: mix two abundance predictions (DBG/ggcat vs Ganon).

输入格式（两者相同）：
- 无表头，两列 TSV
- 第 1 列：物种 ID（字符串或数字都行）
- 第 2 列：数值（相对丰度，建议当做浮点读取；若为整数也能自动转浮点）
"""

import argparse
from pathlib import Path
import pandas as pd


# ======================
# I/O：固定两列无表头
# ======================
def load_abundance(filepath: str) -> pd.Series:
    """
    读取两列、无表头的 TSV：
      col0 = 物种ID（保留为字符串）
      col1 = abundance 数值（转 float）
    返回：Series(index=id[str], values=float)，已去除 NaN 并裁掉负值为 0
    """
    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=["id", "val"],
        dtype={"id": str},         # id 一律按字符串读取，避免数值化丢前导零等
        engine="python"
    )
    s = pd.to_numeric(df["val"], errors="coerce").fillna(0.0).astype(float)
    s = s.clip(lower=0.0)         # 防止极小负数（数值噪声）
    # 统一 id 去空白
    idx = df["id"].astype(str).str.strip()
    ser = pd.Series(s.values, index=idx)

    # （可选）若总和为 0，保持原样；归一化在后续统一处理
    return ser


# ======================
# 归一化与对齐
# ======================
def normalize_series(s: pd.Series) -> pd.Series:
    tot = float(s.sum())
    if tot > 0:
        return (s / tot).astype(float)
    return s.astype(float)


def adjust_ganon(ganon: pd.Series, dbg_series: pd.Series) -> pd.Series:
    """
    仅保留与 DBG 的交集，并在交集上归一化
    """
    common = ganon.index.intersection(dbg_series.index)
    if len(common) == 0:
        raise ValueError("错误：Ganon 和 DBG 的物种集合没有交集，无法混合。")
    g_sub = ganon[common].astype(float)
    return normalize_series(g_sub)


# ======================
# 混合：线性加权 + 归一化
# ======================
def hybrid_predict(dbg_series: pd.Series,
                   adjusted_ganon: pd.Series,
                   dbg_weight: float) -> pd.Series:
    """
    final = dbg_weight * dbg + (1 - dbg_weight) * adjusted_ganon
    然后在交集上归一化
    """
    if not 0.0 <= dbg_weight <= 1.0:
        raise ValueError("权重必须在 0 到 1 之间。")

    common = dbg_series.index.intersection(adjusted_ganon.index)
    if len(common) == 0:
        raise ValueError("错误：混合时找不到交集物种。")

    d_sub = dbg_series[common].astype(float)
    g_sub = adjusted_ganon[common].astype(float)

    hybrid = dbg_weight * d_sub + (1.0 - dbg_weight) * g_sub
    hybrid = hybrid.clip(lower=0.0)
    hybrid = normalize_series(hybrid)

    if hybrid.sum() <= 0:
        raise ValueError("错误：混合后丰度总和为 0。")

    return hybrid


# ======================
# Themis 内部调用接口
# ======================
def run(dbg_file: str, ganon_file: str, weight: float, output: str):
    """
    参数：
      dbg_file   两列无表头 TSV（id \t abundance）
      ganon_file 两列无表头 TSV（id \t abundance）
      weight     DBG 权重（0-1）
      output     输出文件路径（两列 TSV，无表头）
    """
    if not 0.0 <= weight <= 1.0:
        raise ValueError("权重必须在 0 到 1 之间。")

    dbg = load_abundance(dbg_file)
    gan = load_abundance(ganon_file)

    # 归一化（全局）后再求交集加权，避免极端值影响
    dbg = normalize_series(dbg)
    gan = normalize_series(gan)

    gan_adj = adjust_ganon(gan, dbg)
    final = hybrid_predict(dbg, gan_adj, weight)

    # 按丰度降序输出，更利于核查
    final = final.sort_values(ascending=False)

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(output, sep="\t", header=False)


# ======================
# 命令行
# ======================
def parse_args():
    p = argparse.ArgumentParser(description="混合两份丰度预测 (DBG/ggcat 与 Ganon)")
    p.add_argument("-d", "--dbgtax", required=True, help="DBG/ggcat 丰度文件（两列无表头 TSV）")
    p.add_argument("-g", "--ganon", required=True, help="Ganon 丰度文件（两列无表头 TSV）")
    p.add_argument("-w", "--weight", type=float, required=True, help="DBG 权重 (0-1)")
    p.add_argument("-o", "--output", default="hybrid_prediction.tsv", help="输出路径")
    return p.parse_args()


def main():
    args = parse_args()
    run(args.dbgtax, args.ganon, args.weight, args.output)
    print(f"[mix] 完成：{args.output}")


if __name__ == "__main__":
    main()
