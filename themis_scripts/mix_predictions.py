#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import argparse
from pathlib import Path
import pandas as pd


# ======================

# ======================
def load_abundance(filepath: str) -> pd.Series:
    
    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=["id", "val"],
        dtype={"id": str},         
    )
    s = pd.to_numeric(df["val"], errors="coerce").fillna(0.0).astype(float)
    s = s.clip(lower=0.0)         
    # 
    idx = df["id"].astype(str).str.strip()
    ser = pd.Series(s.values, index=idx)

    # 
    return ser


# ======================

# ======================
def normalize_series(s: pd.Series) -> pd.Series:
    tot = float(s.sum())
    if tot > 0:
        return (s / tot).astype(float)
    return s.astype(float)


def adjust_ganon(ganon: pd.Series, dbg_series: pd.Series) -> pd.Series:
    common = ganon.index.intersection(dbg_series.index)
    if len(common) == 0:
        raise ValueError("Error: The species sets have no overlap and cannot be mixed.")
    g_sub = ganon[common].astype(float)
    return normalize_series(g_sub)


# ======================

# ======================
def hybrid_predict(dbg_series: pd.Series,
                   adjusted_ganon: pd.Series,
                   dbg_weight: float) -> pd.Series:
    """
    final = dbg_weight * dbg + (1 - dbg_weight) * adjusted_ganon
    """
    if not 0.0 <= dbg_weight <= 1.0:
        raise ValueError("The weights must be between 0 and 1.")

    common = dbg_series.index.intersection(adjusted_ganon.index)
    if len(common) == 0:
        raise ValueError("Error: No intersecting species found during mixing.")

    d_sub = dbg_series[common].astype(float)
    g_sub = adjusted_ganon[common].astype(float)

    hybrid = dbg_weight * d_sub + (1.0 - dbg_weight) * g_sub
    hybrid = hybrid.clip(lower=0.0)
    hybrid = normalize_series(hybrid)

    if hybrid.sum() <= 0:
        raise ValueError("Error: The sum of abundances after mixing is 0.")

    return hybrid


# ======================

# ======================
def run(dbg_file: str, ganon_file: str, weight: float, output: str):
    if not 0.0 <= weight <= 1.0:
        raise ValueError("The weights must be between 0 and 1.")

    dbg = load_abundance(dbg_file)
    gan = load_abundance(ganon_file)

    # 
    dbg = normalize_series(dbg)
    gan = normalize_series(gan)

    gan_adj = adjust_ganon(gan, dbg)
    final = hybrid_predict(dbg, gan_adj, weight)

    # 
    final = final.sort_values(ascending=False)

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(output, sep="\t", header=False)


# ======================

# ======================
def parse_args():
    p = argparse.ArgumentParser(description="Mixed abundance prediction")
    p.add_argument("-d", "--dbgtax", required=True, help="Abundance file (two columns without header TSV)")
    p.add_argument("-g", "--ganon", required=True, help="Abundance file (two columns without header TSV)")
    p.add_argument("-w", "--weight", type=float, required=True, help="DBG weights (0-1)")
    p.add_argument("-o", "--output", default="hybrid_prediction.tsv", help="Output path")
    return p.parse_args()


def main():
    args = parse_args()
    run(args.dbgtax, args.ganon, args.weight, args.output)
    print(f"[mix] Completed:{args.output}")


if __name__ == "__main__":
    main()
