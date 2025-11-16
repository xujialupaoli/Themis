#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Themis helper: select species-strain mapping based on ganon predictions.
(最小修改版：ID 一律按字符串；表头判断更稳)
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

        # 统一规范化首列内容
        if header:
            h0 = str(header[0]).strip().lower()
        else:
            h0 = ""

        # 判断是否有表头（更稳）
        if h0 and ("species_taxid" in h0 or "species" in h0):
            # 有表头，从下一行开始读数据
            pass
        else:
            # 无表头：把第一行当数据
            if header and str(header[0]).strip():
                species_set.add(str(header[0]).strip())

        for row in reader:
            if not row:
                continue
            sid = str(row[0]).strip()  # ID 一律字符串
            if sid:
                species_set.add(sid)

    return species_set


def read_ref_mapping(ref_info_file: Path):
    """
    从 ref_info 建立:
        genome_ID -> (species_taxid, id_path)
    要求 ref_info 至少包含列:
        genome_ID, species_taxid, id
    """
    mapping = {}
    with open(ref_info_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            sys.exit("[错误] Ref 信息文件为空或缺少表头。")

        # 列名做 strip，避免尾随空白导致 KeyError
        cols = {str(name).strip(): i for i, name in enumerate(header)}
        try:
            i_genome = cols["genome_ID"]
            i_species = cols["species_taxid"]
            i_idpath = cols["id"]
        except KeyError as e:
            sys.exit(f"[错误] Ref 信息文件缺少必要列：{e} (需要: genome_ID, species_taxid, id)")

        for row in reader:
            if not row or len(row) <= i_idpath:
                continue
            genome_id = str(row[i_genome]).strip()   # 强制字符串
            species_taxid = str(row[i_species]).strip()
            id_path = str(row[i_idpath]).strip()
            if not genome_id or not species_taxid or not id_path:
                continue
            mapping[genome_id] = (species_taxid, id_path)

    return mapping


def read_ganon_strains(ganon_strain_file: Path):
    """
    读取 ganon 预测到的菌株列表（第一列为 strain_taxid / genome_ID）。
    按出现顺序保留。
    """
    strains = []
    with open(ganon_strain_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        # 判断是否有表头（更稳）
        if header:
            h0 = str(header[0]).strip().lower()
        else:
            h0 = ""

        if not h0 or ("strain_taxid" not in h0 and "strain" not in h0):
            # 无表头：把第一行当数据
            if header and str(header[0]).strip():
                strains.append(str(header[0]).strip())

        for row in reader:
            if not row:
                continue
            sid = str(row[0]).strip()  # ID 一律字符串
            if sid:
                strains.append(sid)

    return strains


# ======================
# Themis pipeline 用的 run()
# ======================

def run(ref_info,
        ganon_species,
        ganon_strain,
        out_tsv,
        tmp_dir=None):
    """
    根据 ganon 物种 & 菌株预测 + ref_info，生成筛选后的 species-strain 映射 TSV。

    参数:
      ref_info:      Ref 信息表路径 (含 genome_ID, species_taxid, id 列)
      ganon_species: ganon 物种丰度文件
      ganon_strain:  ganon 菌株丰度文件
      out_tsv:       输出 TSV 路径，列为:
                       species_taxid  strain_taxid  new_strain_taxid  id
      tmp_dir:       为兼容保留，不再创建软链接；可传 None 或忽略

    返回:
      输出文件的 Path 对象
    """
    ref_info_file = Path(ref_info)
    ganon_species_file = Path(ganon_species)
    ganon_strain_file = Path(ganon_strain)
    out_tsv_path = Path(out_tsv)

    # 1) 读取 ganon 预测物种集合（ID 均为字符串）
    species_set = read_species_set(ganon_species_file)
    if not species_set:
        raise SystemExit("[错误] 未从 ganon 物种文件读到任何 species_taxid。")

    # 2) 读取 Ref 映射（键值 ID 全是字符串）
    ref_map = read_ref_mapping(ref_info_file)
    if not ref_map:
        raise SystemExit("[错误] Ref 映射为空，请检查 ref_info 文件。")

    # 3) 读取 ganon 预测菌株（ID 字符串）
    predicted_strains = read_ganon_strains(ganon_strain_file)
    if not predicted_strains:
        raise SystemExit("[错误] 未从 ganon 菌株文件读到任何 strain_taxid。")

    # 4) 过滤并生成记录
    per_species_counter = defaultdict(int)
    records = []  # (species_taxid, strain_taxid, new_strain_taxid, id_path)

    for strain_taxid in predicted_strains:
        if strain_taxid not in ref_map:
            # 不在 ref 表里的菌株，提示但不终止
            print(f"[提示] 预测菌株在 Ref 中未找到：{strain_taxid}", file=sys.stderr)
            continue

        species_taxid, id_path = ref_map[strain_taxid]
        if species_taxid not in species_set:
            # 该菌株所属物种不在 ganon 物种预测集合中 -> 丢弃
            continue

        per_species_counter[species_taxid] += 1
        new_strain_taxid = f"{species_taxid}_{per_species_counter[species_taxid]}"
        records.append((species_taxid, strain_taxid, new_strain_taxid, id_path))

    if not records:
        raise SystemExit("[错误] 过滤后没有可用记录。请检查 ganon 预测与 ref_info 是否匹配。")

    # 5) 写出结果 TSV
    out_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["species_taxid", "strain_taxid", "new_strain_taxid", "id"])
        writer.writerows(records)

    # 不再创建软链接；仅提示一次（如果调用方传了 tmp_dir）
    if tmp_dir:
        print(f"[提示] Themis: 已忽略 tmp_dir={tmp_dir}，当前版本不再创建软链接。", file=sys.stderr)

    print(f"[完成] 已写出筛选映射表: {out_tsv_path}")
    print(f"[统计] 物种数: {len(per_species_counter)}；菌株条目: {len(records)}")

    return out_tsv_path


# ======================
# 命令行入口 (兼容)
# ======================

def parse_args():
    p = argparse.ArgumentParser(
        description="基于 Ref 数据库与 ganon 预测结果，输出 species-strain 映射表（不创建软链接）。"
    )
    p.add_argument("--ref_info", required=True,
                   help="Ref 信息表路径（含 genome_ID, species_taxid, id 列）。")
    p.add_argument("--ganon_species", required=True,
                   help="ganon_species_abundance.txt 路径。")
    p.add_argument("--ganon_strain", required=True,
                   help="ganon_strain_abundance.txt 路径。")
    p.add_argument("--out_tsv", default="ganon_species_strain_selected.tsv",
                   help="输出 TSV 文件名（默认：ganon_species_strain_selected.tsv）。")
    p.add_argument("--tmp_dir", default=None,
                   help="(可选) 旧版用于软链接的目录，Themis 中已不再使用。")
    return p.parse_args()


def main():
    args = parse_args()
    run(
        ref_info=args.ref_info,
        ganon_species=args.ganon_species,
        ganon_strain=args.ganon_strain,
        out_tsv=args.out_tsv,
        tmp_dir=args.tmp_dir,
    )


if __name__ == "__main__":
    main()
