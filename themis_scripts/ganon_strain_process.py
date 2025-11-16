#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Themis helper:
从 ganon tax_profile.tre 中抽取 strain/genome 层面的丰度，生成 strain_abundance.txt。

行为基于你提供的原始脚本:
  - 解析以 "strain" 开头的行
  - tokens[3] 中查找 "GCA"/"GCF"，截取到行尾作为 genome_ID
  - abundance = tokens[8] / 100
  - 若提供 genomes_info_file，则仅保留 genome_ID 在其中者，并在此集合内重新归一化
  - 输出:
      strain_taxid \t abundance
"""

import sys
from pathlib import Path

import pandas as pd


def run(tre_file, output_path, genomes_info_file=None):
    # 1) 可选 genomes 列表
    genomesID = None
    if genomes_info_file and genomes_info_file != "-":
        df = pd.read_csv(genomes_info_file, sep="\t")
        if "genome_ID" in df.columns:
            genomesID = set(df["genome_ID"].astype(str).tolist())

    tax_profile_dict = {}

    # 2) 从 tre_file 提取 strain abundance
    with open(tre_file, "r", encoding="utf-8") as f_in:
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith("strain"):
                tokens = line.split("\t")
                if len(tokens) < 9:
                    continue

                name = tokens[3]  # 原脚本用 tokens[3][3:], 这里行为保持尽量接近
                # 在 name 中寻找 GCA/GCF
                start = name.find("GCA")
                if start < 0:
                    start = name.find("GCF")
                if start < 0:
                    # 找不到 GCA/GCF，就跳过
                    continue
                genome_ID = name[start:].strip()
                if not genome_ID:
                    continue

                try:
                    abundance = float(tokens[8]) / 100.0
                except ValueError:
                    continue

                # 可选过滤: 只保留在 genomesID 中的
                if genomesID is not None:
                    if genome_ID in genomesID:
                        tax_profile_dict[genome_ID] = abundance
                else:
                    tax_profile_dict[genome_ID] = abundance

    # 3) 若提供 genomesID，则重新归一化
    if genomesID is not None and tax_profile_dict:
        s = sum(tax_profile_dict.values())
        if s > 0:
            for k in list(tax_profile_dict.keys()):
                tax_profile_dict[k] = tax_profile_dict[k] / s

    # 4) 排序 & 输出
    sorted_items = sorted(
        tax_profile_dict.items(),
        key=lambda item: item[1],
        reverse=True,
    )

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f_out:
        f_out.write("strain_taxid\tabundance\n")
        for k, v in sorted_items:
            f_out.write(f"{k}\t{v}\n")


def main():
    # 命令行兼容模式:
    # python ganon_strain_process.py tax_profile.tre [genomes_info.txt]
    if len(sys.argv) not in (3, 4):
        print("Usage: ganon_strain_process.py tax_profile.tre strain_abundance.txt [genomes_info.txt]", file=sys.stderr)
        sys.exit(1)

    tre_file = sys.argv[1]
    out_file = sys.argv[2]
    genomes_info = sys.argv[3] if len(sys.argv) == 4 else None

    run(tre_file, out_file, genomes_info_file=genomes_info)


if __name__ == "__main__":
    main()
