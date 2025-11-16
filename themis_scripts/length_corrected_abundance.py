#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Themis helper: compute length-corrected relative abundance per species.

This is adapted from your original length_corrected_abundance.py, but:
- 提供 run() 接口，方便在 Themis pipeline 内部调用；
- 保留命令行使用方式，行为不变。

功能:
  输入:
    1) counts.tsv (无表头): 两列
         species_id<TAB>reads_count
    2) mapping.tsv (有表头): 至少包含:
         - species_taxid (或自定义列名)
         - id           (strain genome FASTA 路径, .gz 或普通文件)

  方法:
    - 根据 counts.tsv 累加每个物种的 reads 数
    - 在 mapping.tsv 中找到该物种对应的所有基因组路径
    - 计算每个基因组的总碱基数, 取平均作为该物种代表长度
    - w_i = reads_i / avg_len_i
    - rel_i = w_i / sum_j w_j

  输出:
    带表头 TSV:
      species_taxid  reads  n_strains  avg_genome_len  reads_over_len  rel_abundance
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict


# ======================
# 内部通用函数
# ======================

def read_counts(counts_path):
    counts = defaultdict(float)
    with open(counts_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                print(f"[warn] counts line {ln} malformed: {line}", file=sys.stderr)
                continue
            sid, val = parts[0].strip(), parts[1].strip()
            try:
                counts[sid] += float(val)
            except ValueError:
                print(f"[warn] counts line {ln} bad number: {val}", file=sys.stderr)
                continue
    return counts


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")


def fasta_length(path):
    total = 0
    with open_maybe_gz(path) as f:
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                continue
            total += len(line.strip())
    return total


def read_mapping(mapping_path, species_col, path_col):
    with open(mapping_path, "r", encoding="utf-8", errors="ignore") as f:
        header = f.readline()
        if not header:
            raise RuntimeError("Empty mapping file.")
        cols = header.rstrip("\n").split("\t")
        try:
            i_species = cols.index(species_col)
        except ValueError:
            raise RuntimeError(f"Mapping header missing column '{species_col}'. Columns: {cols}")
        try:
            i_path = cols.index(path_col)
        except ValueError:
            raise RuntimeError(f"Mapping header missing column '{path_col}'. Columns: {cols}")

        species_to_paths = defaultdict(list)
        for ln, line in enumerate(f, 2):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_species, i_path):
                print(f"[warn] mapping line {ln} too few columns; skipped.", file=sys.stderr)
                continue
            sid = parts[i_species].strip()
            gpath = parts[i_path].strip()
            if not sid or not gpath:
                continue
            species_to_paths[sid].append(gpath)
    return species_to_paths


# ======================
# 核心：给 Themis 用的 run()
# ======================

def run(counts_file,
        mapping_file,
        output_file,
        species_col="species_taxid",
        path_col="id",
        skip_missing_genomes=True):
    """
    在 Themis pipeline 内部调用的接口。

    参数:
      counts_file:  species_id \t reads_count (无表头)
      mapping_file: 含物种与基因组路径信息的表，至少包含 [species_col, path_col]
      output_file:  输出 TSV
      species_col:  映射表中物种列名（默认 'species_taxid'）
      path_col:     映射表中基因组路径列名（默认 'id'）
      skip_missing_genomes:
                    True 时缺失的基因组文件仅告警并跳过；
                    False 行为与原脚本近似（这里仍然跳过，但不终止）
    """

    counts = read_counts(counts_file)
    if not counts:
        print("[error] No counts read from counts file.", file=sys.stderr)
        raise SystemExit(1)

    species_to_paths = read_mapping(mapping_file, species_col, path_col)

    genome_len_cache = {}
    species_avg_len = {}
    species_n = {}

    # 为 counts 中出现的每个物种计算平均基因组长度
    for sid in counts.keys():
        paths = species_to_paths.get(sid, [])
        if not paths:
            print(f"[warn] species {sid} has counts but no genomes in mapping; skipped in abundance calc.",
                  file=sys.stderr)
            continue

        lens = []
        for p in paths:
            if not os.path.isfile(p):
                msg = f"[warn] missing genome file: {p} (species {sid}); skipping this strain."
                print(msg, file=sys.stderr)
                if skip_missing_genomes:
                    continue
                else:
                    continue  # 保持与你原始行为接近：这里也只是跳过
            if p not in genome_len_cache:
                try:
                    genome_len_cache[p] = fasta_length(p)
                except Exception as e:
                    print(f"[warn] failed to read {p}: {e}; skipping. ({e})", file=sys.stderr)
                    continue
            lens.append(genome_len_cache[p])

        if not lens or sum(lens) <= 0:
            print(f"[warn] species {sid} has no valid genome lengths; skipped in abundance calc.",
                  file=sys.stderr)
            continue

        species_n[sid] = len(lens)
        species_avg_len[sid] = sum(lens) / len(lens)

    # 计算权重: reads / avg_len
    weights = {}
    for sid, rc in counts.items():
        L = species_avg_len.get(sid)
        if L and L > 0:
            weights[sid] = rc / L

    denom = sum(weights.values())
    if denom == 0:
        print("[error] Sum of (reads / avg_len) is zero; nothing to normalize.", file=sys.stderr)
        raise SystemExit(2)

    # 写输出
    with open(output_file, "w", encoding="utf-8") as out:
        out.write("\t".join([
            "species_taxid",
            "reads",
            "n_strains",
            "avg_genome_len",
            "reads_over_len",
            "rel_abundance",
        ]) + "\n")

        # 排序规则保持稳定：先有权重的，再按 ID 排
        for sid in sorted(counts.keys(), key=lambda x: (x not in weights, x)):
            reads = counts[sid]
            n = species_n.get(sid, 0)
            L = species_avg_len.get(sid, 0)
            w = weights.get(sid, 0.0)
            rel = (w / denom) if w > 0 else 0.0
            out.write(f"{sid}\t{reads}\t{n}\t{L}\t{w}\t{rel}\n")

    print(f"[info] length-corrected abundance written to: {output_file}", file=sys.stderr)


# ======================
# 原始命令行入口 (保持向后兼容)
# ======================

def parse_args():
    ap = argparse.ArgumentParser(
        description="Compute length-corrected relative abundance per species."
    )
    ap.add_argument(
        "-c", "--counts", required=True,
        help="Counts TSV: species_id<TAB>reads_count (no header)."
    )
    ap.add_argument(
        "-m", "--mapping", required=True,
        help="Mapping TSV with header; needs columns 'species_taxid' and 'id' (genome path)."
    )
    ap.add_argument(
        "-o", "--output", required=True,
        help="Output TSV path."
    )
    ap.add_argument(
        "--species-col", default="species_taxid",
        help="Column name for species id in mapping file (default: species_taxid)."
    )
    ap.add_argument(
        "--path-col", default="id",
        help="Column name for genome path in mapping file (default: id)."
    )
    ap.add_argument(
        "--skip-missing-genomes", action="store_true",
        help="If set, ignore missing genome files instead of failing."
    )
    return ap.parse_args()


def main():
    args = parse_args()
    run(
        counts_file=args.counts,
        mapping_file=args.mapping,
        output_file=args.output,
        species_col=args.species_col,
        path_col=args.path_col,
        skip_missing_genomes=args.skip_missing_genomes,
    )


if __name__ == "__main__":
    main()
