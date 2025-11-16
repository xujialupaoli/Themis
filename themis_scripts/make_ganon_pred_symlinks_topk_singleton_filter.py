#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Themis helper: select top-K strains per species based on ganon predictions,
with singleton abundance filtering.

保持原逻辑不变；仅确保所有 ID（species_taxid / strain_taxid / genome_ID）一律按字符串处理，
只有 abundance 转为 float。
"""

import csv
import sys
import argparse
from pathlib import Path
from collections import defaultdict


def read_species_set(ganon_species_file: Path):
    """读取 ganon 预测的物种ID集合（第一列为 species_taxid）。"""
    species_set = set()
    with open(ganon_species_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        if header and (
            "species_taxid" in header[0].lower()
            or "species" in header[0].lower()
        ):
            # 有表头，从下一行开始
            pass
        else:
            # 无表头，则把首行当作数据行
            if header and header[0].strip():
                # NOTE: 强制字符串
                species_set.add(str(header[0]).strip())

        for row in reader:
            if not row:
                continue
            # NOTE: 强制字符串
            sid = str(row[0]).strip()
            if sid:
                species_set.add(sid)

    return species_set


def read_ref_mapping(ref_info_file: Path):
    """
    从 ref_info 建立:
        genome_ID -> (species_taxid, id_path)
    需要列：genome_ID, species_taxid, id
    """
    with open(ref_info_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            sys.exit("[错误] Ref 信息文件为空或缺少表头。")

        cols = {name: i for i, name in enumerate(header)}
        for need in ("genome_ID", "species_taxid", "id"):
            if need not in cols:
                sys.exit(f"[错误] Ref 信息文件缺少必要列：{need}")

        i_genome = cols["genome_ID"]
        i_species = cols["species_taxid"]
        i_idpath = cols["id"]

        mapping = {}
        for row in reader:
            if not row or len(row) <= i_idpath:
                continue
            # NOTE: 三个关键列均强制字符串
            genome_id = str(row[i_genome]).strip()
            species_taxid = str(row[i_species]).strip()
            id_path = str(row[i_idpath]).strip()
            if not genome_id or not species_taxid or not id_path:
                continue
            mapping[genome_id] = (species_taxid, id_path)

    return mapping


def read_ganon_strains_with_abund(ganon_strain_file: Path):
    """
    读取 ganon 预测到的菌株及其相对丰度：
      第一列：strain_taxid（与 Ref 的 genome_ID 对应，字符串）
      第二列：abundance（浮点）
    返回列表 [(strain_taxid, abundance)]，按出现顺序。
    """
    items = []
    with open(ganon_strain_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        # 判断首行是不是表头
        if header and len(header) >= 2 and (
            "strain" in header[0].lower()
            or "abund" in header[1].lower()
        ):
            # 有表头，从下一行读数据
            pass
        else:
            # 无表头：尝试把首行当作数据
            if header and len(header) >= 2:
                # NOTE: ID 强制字符串
                st = str(header[0]).strip()
                try:
                    abund = float(header[1])  # 仅 abundance 转 float
                except Exception:
                    abund = None
                if st and abund is not None:
                    items.append((st, abund))

        for row in reader:
            if not row or len(row) < 2:
                continue
            # NOTE: ID 强制字符串
            st = str(row[0]).strip()
            try:
                abund = float(row[1])  # 仅 abundance 转 float
            except Exception:
                continue
            if st:
                items.append((st, abund))

    return items


# ======================
# Themis 用的核心接口
# ======================

def run(ref_info,
        ganon_species,
        ganon_strain,
        top_k,
        singleton_min_abund,
        out_tsv,
        tmp_dir=None):
    """
    根据 ganon 物种 & 菌株丰度 + ref_info，筛选 Top-K 菌株并输出映射表。

    参数:
      ref_info:            Ref 信息表路径 (含 genome_ID, species_taxid, id)
      ganon_species:       ganon 物种丰度
      ganon_strain:        ganon 菌株丰度 (需含 abundance)
      top_k:               每个物种保留菌株数
      singleton_min_abund: 若某物种仅1个菌株预测，则其 abundance 必须 > 此阈值才保留
      out_tsv:             输出 TSV，列:
                             species_taxid  strain_taxid  new_strain_taxid  id
      tmp_dir:             仅为兼容参数，Themis 中不再创建软链接，可为 None

    返回:
      out_tsv 的 Path
    """
    if top_k <= 0:
        raise SystemExit("[错误] top_k 必须为正整数。")

    ref_info_file = Path(ref_info)
    ganon_species_file = Path(ganon_species)
    ganon_strain_file = Path(ganon_strain)
    out_tsv_path = Path(out_tsv)

    # 1) 物种集合（ID 皆为字符串）
    species_set = read_species_set(ganon_species_file)
    if not species_set:
        raise SystemExit("[错误] 未从 ganon 物种文件读到任何 species_taxid。")

    # 2) Ref 映射（键与值中的 ID 皆为字符串）
    ref_map = read_ref_mapping(ref_info_file)
    if not ref_map:
        raise SystemExit("[错误] Ref 映射为空，请检查 ref_info 文件。")

    # 3) 菌株 + abundance（ID=字符串；abundance=float）
    predicted_items = read_ganon_strains_with_abund(ganon_strain_file)
    if not predicted_items:
        raise SystemExit("[错误] 未从 ganon 菌株文件读到任何 (strain_taxid, abundance)。")

    # 4) 映射到物种，仅保留物种在 species_set 中的条目
    per_species_strains = defaultdict(list)  # sp -> [(strain_taxid, abund, id_path)]
    not_in_ref = 0

    for strain_taxid, abund in predicted_items:
        # NOTE: 这里的 strain_taxid = genome_ID（字符串）
        if strain_taxid not in ref_map:
            not_in_ref += 1
            continue
        species_taxid, id_path = ref_map[strain_taxid]  # 两者都是字符串
        if species_taxid not in species_set:
            continue
        per_species_strains[species_taxid].append((strain_taxid, abund, id_path))

    # 5) 单菌株物种阈值过滤
    filtered_species = {}
    thr = float(singleton_min_abund)
    dropped_singletons = 0

    for sp, arr in per_species_strains.items():
        if len(arr) == 1:
            (_, abund, _) = arr[0]
            if not (abund > thr):
                dropped_singletons += 1
                continue
        filtered_species[sp] = arr

    if not filtered_species:
        raise SystemExit("[错误] 过滤后没有可用的物种。单菌株阈值可能过高或输入不匹配。")

    # 6) 每个物种按 abundance 降序取 top_k
    selected_records = []  # (species_taxid, strain_taxid, new_strain_taxid, id_path)
    for sp, arr in filtered_species.items():
        # NOTE: x[0] 是字符串 ID；x[1] 是 float abundance
        arr_sorted = sorted(arr, key=lambda x: (-x[1], x[0]))  # abundance desc, then id
        top_arr = arr_sorted[:top_k]
        for idx, (strain_taxid, abund, id_path) in enumerate(top_arr, start=1):
            # NOTE: sp 是字符串（species_taxid），字符串拼接无歧义
            new_name = f"{sp}_{idx}"
            selected_records.append((sp, strain_taxid, new_name, id_path))

    if not selected_records:
        raise SystemExit("[错误] Top-K 选择后没有记录，请检查 top_k 与输入数据。")

    # 7) 写出 TSV（全为字符串 + abundance未输出）
    out_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["species_taxid", "strain_taxid", "new_strain_taxid", "id"])
        w.writerows(selected_records)

    # 不再创建软链接，只提示一下兼容行为
    if tmp_dir:
        print(f"[提示] Themis: 已忽略 tmp_dir={tmp_dir}，当前版本不再创建软链接。", file=sys.stderr)

    print(f"[完成] 已写出 Top-K 映射表: {out_tsv_path}")
    print(f"[统计] 通过筛选的物种数: {len(filtered_species)}；最终菌株条目: {len(selected_records)}")
    if dropped_singletons:
        print(f"[提示] 因单菌株 abundance ≤ 阈值而被丢弃的物种数：{dropped_singletons}", file=sys.stderr)
    if not_in_ref:
        print(f"[提示] 有 {not_in_ref} 个预测菌株在 Ref 中未找到，对应条目已跳过。", file=sys.stderr)

    return out_tsv_path


# ======================
# 命令行入口（兼容原用法）
# ======================

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "基于 ganon 物种/菌株丰度和 Ref 信息，为每个物种选择 Top-K 菌株，"
            "并对单菌株物种应用 abundance 阈值过滤（不再创建软链接）。"
        )
    )
    p.add_argument("--ref_info", required=True,
                   help="Ref 信息表路径（需含 genome_ID, species_taxid, id 列）。")
    p.add_argument("--ganon_species", required=True,
                   help="ganon_species_abundance.txt 路径。")
    p.add_argument("--ganon_strain", required=True,
                   help="ganon_strain_abundance.txt 路径（需含 abundance 列）。")
    p.add_argument("--top_k", type=int, required=True,
                   help="每个物种保留的 Top-K 菌株数。")
    p.add_argument("--singleton_min_abund", type=float,
                   default=1.0000000000000001e-07,
                   help="单菌株物种的最小 abundance 阈值（严格大于此值才保留）。")
    p.add_argument("--out_tsv", default="ganon_species_strain_topk.tsv",
                   help="输出 TSV 文件名（默认：ganon_species_strain_topk.tsv）。")
    p.add_argument("--tmp_dir", default=None,
                   help="兼容参数，Themis 版本中不再创建软链接，可忽略。")
    return p.parse_args()


def main():
    args = parse_args()
    run(
        ref_info=args.ref_info,
        ganon_species=args.ganon_species,
        ganon_strain=args.ganon_strain,
        top_k=args.top_k,
        singleton_min_abund=args.singleton_min_abund,
        out_tsv=args.out_tsv,
        tmp_dir=args.tmp_dir,
    )


if __name__ == "__main__":
    main()
