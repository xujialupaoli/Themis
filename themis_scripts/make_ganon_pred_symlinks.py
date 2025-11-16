#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import csv
import sys
import argparse
from pathlib import Path
from collections import defaultdict


def read_species_set(ganon_species_file: Path):
    
    species_set = set()
    with open(ganon_species_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        # 
        if header:
            h0 = str(header[0]).strip().lower()
        else:
            h0 = ""

        # 
        if h0 and ("species_taxid" in h0 or "species" in h0):
            # 
            pass
        else:
            # 
            if header and str(header[0]).strip():
                species_set.add(str(header[0]).strip())

        for row in reader:
            if not row:
                continue
            sid = str(row[0]).strip()  
            if sid:
                species_set.add(sid)

    return species_set


def read_ref_mapping(ref_info_file: Path):
    mapping = {}
    with open(ref_info_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            sys.exit("[Error] The Ref information file is empty or missing a header.")

        # 
        cols = {str(name).strip(): i for i, name in enumerate(header)}
        try:
            i_genome = cols["genome_ID"]
            i_species = cols["species_taxid"]
            i_idpath = cols["id"]
        except KeyError as e:
            sys.exit(f"[Error] Ref information file is missing a required column: {e} (Required: genome_ID, species_taxid, id)")

        for row in reader:
            if not row or len(row) <= i_idpath:
                continue
            genome_id = str(row[i_genome]).strip()   
            species_taxid = str(row[i_species]).strip()
            id_path = str(row[i_idpath]).strip()
            if not genome_id or not species_taxid or not id_path:
                continue
            mapping[genome_id] = (species_taxid, id_path)

    return mapping


def read_ganon_strains(ganon_strain_file: Path):
    strains = []
    with open(ganon_strain_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        # 
        if header:
            h0 = str(header[0]).strip().lower()
        else:
            h0 = ""

        if not h0 or ("strain_taxid" not in h0 and "strain" not in h0):
            # 
            if header and str(header[0]).strip():
                strains.append(str(header[0]).strip())

        for row in reader:
            if not row:
                continue
            sid = str(row[0]).strip()  
            if sid:
                strains.append(sid)

    return strains


# ======================

# ======================

def run(ref_info,
        ganon_species,
        ganon_strain,
        out_tsv,
        tmp_dir=None):
    
    ref_info_file = Path(ref_info)
    ganon_species_file = Path(ganon_species)
    ganon_strain_file = Path(ganon_strain)
    out_tsv_path = Path(out_tsv)

    # 1) 
    species_set = read_species_set(ganon_species_file)
    if not species_set:
        raise SystemExit("[Error] No species_taxid read from ganon species file.")

    # 2) 
    ref_map = read_ref_mapping(ref_info_file)
    if not ref_map:
        raise SystemExit("[Error] Ref mapping is empty. Please check the ref_info file.")

    # 3) 
    predicted_strains = read_ganon_strains(ganon_strain_file)
    if not predicted_strains:
        raise SystemExit("[Error] No strain_taxid read from the ganon strain file.")

    # 4) 
    per_species_counter = defaultdict(int)
    records = []  # (species_taxid, strain_taxid, new_strain_taxid, id_path)

    for strain_taxid in predicted_strains:
        if strain_taxid not in ref_map:
            # 
            print(f"[Hint] The predicted strain was not found in the Ref:{strain_taxid}", file=sys.stderr)
            continue

        species_taxid, id_path = ref_map[strain_taxid]
        if species_taxid not in species_set:
            # 
            continue

        per_species_counter[species_taxid] += 1
        new_strain_taxid = f"{species_taxid}_{per_species_counter[species_taxid]}"
        records.append((species_taxid, strain_taxid, new_strain_taxid, id_path))

    if not records:
        raise SystemExit("[Error] No records were found after filtering. ")

    # 5) 写出结果 TSV
    out_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["species_taxid", "strain_taxid", "new_strain_taxid", "id"])
        writer.writerows(records)

    # No more symbolic links will be created; only one prompt will be given (if the caller provides tmp_dir).
    if tmp_dir:
        print(f"[Note] Themis: tmp_dir={tmp_dir} has been ignored; symbolic links are no longer created in the current version.", file=sys.stderr)

    print(f"[Completed] The filter mapping table has been written:{out_tsv_path}")
    print(f"[Statistics] Number of species: {len(per_species_counter)}；Strains entry: {len(records)}")

    return out_tsv_path


# ======================

# ======================

def parse_args():
    p = argparse.ArgumentParser(
        description="Based on the Ref database and  prediction results, output a species-strain mapping table (without creating soft links)."
    )
    p.add_argument("--ref_info", required=True,
                   help="Ref information table path (including genome_ID, species_taxid, id columns).")
    p.add_argument("--ganon_species", required=True,
                   help="The path to species_abundance.txt.")
    p.add_argument("--ganon_strain", required=True,
                   help="The path to strain_abundance.txt.")
    p.add_argument("--out_tsv", default="ganon_species_strain_selected.tsv",
                   help="Output TSV filename (default: ganon_species_strain_selected.tsv).")
    p.add_argument("--tmp_dir", default=None,
                   help="(Optional) The directory used for symbolic links in the old version is no longer used in Themis.")
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
