#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Themis helper:
从 ganon 的 tax_profile.tre 中抽取 species 层级丰度，生成 species_abundance.txt 风格表格。

对应你原始脚本的行为：
  - 只读取以 "species" 开头的行
  - tokens[2] 的最后一个 '|' 之后作为 species_taxid
  - tokens[8] (百分比) / 100 作为 predicted_abundance
  - 按丰度从大到小排序
  - 输出列名: species_taxid, predicted_abundance
"""

import sys
from pathlib import Path


def run(tre_file, output_path):
    tre_file = Path(tre_file)
    output_path = Path(output_path)

    tax_profile_dict = {}

    with tre_file.open("r", encoding="utf-8") as f_in:
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith("species"):
                tokens = line.split("\t")
                if len(tokens) < 9:
                    continue
                # 第3列形如: something|...|species_taxid
                species_taxid = tokens[2].split("|")[-1]
                try:
                    abundance = float(tokens[8]) / 100.0
                except ValueError:
                    continue
                tax_profile_dict[species_taxid] = abundance

    # 排序（从大到小）
    sorted_tax_profile = sorted(
        tax_profile_dict.items(),
        key=lambda item: item[1],
        reverse=True,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f_out:
        f_out.write("species_taxid\tpredicted_abundance\n")
        for sid, abund in sorted_tax_profile:
            f_out.write(f"{sid}\t{abund}\n")


def main():
    # 命令行兼容用法：
    # python ganon_species_process.py tax_profile.tre
    if len(sys.argv) != 2:
        print("Usage: python ganon_species_process.py tax_profile.tre", file=sys.stderr)
        sys.exit(1)

    tre_file = sys.argv[1]
    run(tre_file, "species_abundance.txt")


if __name__ == "__main__":
    main()
