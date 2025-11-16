#!/usr/bin/env python3
# -*- coding: utf-8 -*-



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
                # 
                species_taxid = tokens[2].split("|")[-1]
                try:
                    abundance = float(tokens[8]) / 100.0
                except ValueError:
                    continue
                tax_profile_dict[species_taxid] = abundance

    # 
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
    # 
    # 
    if len(sys.argv) != 2:
        print("Usage: python ganon_species_process.py tax_profile.tre", file=sys.stderr)
        sys.exit(1)

    tre_file = sys.argv[1]
    run(tre_file, "species_abundance.txt")


if __name__ == "__main__":
    main()
