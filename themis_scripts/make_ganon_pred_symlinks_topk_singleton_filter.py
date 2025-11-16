#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Themis helper: select top-K strains per species based on ganon predictions,
with singleton abundance filtering.
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

        if header and (
            "species_taxid" in header[0].lower()
            or "species" in header[0].lower()
        ):
            # 
            pass
        else:
            # 
            if header and header[0].strip():
                # NOTE: 
                species_set.add(str(header[0]).strip())

        for row in reader:
            if not row:
                continue
            # NOTE: 
            sid = str(row[0]).strip()
            if sid:
                species_set.add(sid)

    return species_set


def read_ref_mapping(ref_info_file: Path):
    with open(ref_info_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            sys.exit("[Error] The Ref information file is empty or missing a header.")

        cols = {name: i for i, name in enumerate(header)}
        for need in ("genome_ID", "species_taxid", "id"):
            if need not in cols:
                sys.exit(f"[Error] Ref information file is missing a necessary column:：{need}")

        i_genome = cols["genome_ID"]
        i_species = cols["species_taxid"]
        i_idpath = cols["id"]

        mapping = {}
        for row in reader:
            if not row or len(row) <= i_idpath:
                continue
            # NOTE: All three key columns must be forced to be strings.
            genome_id = str(row[i_genome]).strip()
            species_taxid = str(row[i_species]).strip()
            id_path = str(row[i_idpath]).strip()
            if not genome_id or not species_taxid or not id_path:
                continue
            mapping[genome_id] = (species_taxid, id_path)

    return mapping


def read_ganon_strains_with_abund(ganon_strain_file: Path):
    items = []
    with open(ganon_strain_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        # 
        if header and len(header) >= 2 and (
            "strain" in header[0].lower()
            or "abund" in header[1].lower()
        ):
            # 
            pass
        else:
            # 
            if header and len(header) >= 2:
                # NOTE: ID Forced string
                st = str(header[0]).strip()
                try:
                    abund = float(header[1])  
                except Exception:
                    abund = None
                if st and abund is not None:
                    items.append((st, abund))

        for row in reader:
            if not row or len(row) < 2:
                continue
            # NOTE: ID Forced string
            st = str(row[0]).strip()
            try:
                abund = float(row[1])  
            except Exception:
                continue
            if st:
                items.append((st, abund))

    return items


# ======================

# ======================

def run(ref_info,
        ganon_species,
        ganon_strain,
        top_k,
        singleton_min_abund,
        out_tsv,
        tmp_dir=None):
    
    if top_k <= 0:
        raise SystemExit("[Error] top_k must be a positive integer.")

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
    predicted_items = read_ganon_strains_with_abund(ganon_strain_file)
    if not predicted_items:
        raise SystemExit("[Error] No (strain_taxid, abundance) was read from the ganon strain file.")

    # 4) 
    per_species_strains = defaultdict(list)  
    not_in_ref = 0

    for strain_taxid, abund in predicted_items:
        # NOTE: Here, strain_taxid = genome_ID (string)
        if strain_taxid not in ref_map:
            not_in_ref += 1
            continue
        species_taxid, id_path = ref_map[strain_taxid]  
        if species_taxid not in species_set:
            continue
        per_species_strains[species_taxid].append((strain_taxid, abund, id_path))

    # 5) 
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
        raise SystemExit("[Error] No species were found after filtering. The single strain threshold may be too high or the input may not match.")

    # 6) For each species, select the top_k in descending order of abundance.
    selected_records = []  # (species_taxid, strain_taxid, new_strain_taxid, id_path)
    for sp, arr in filtered_species.items():
        # NOTE: x[0] is a string ID; x[1] is a float abundance
        arr_sorted = sorted(arr, key=lambda x: (-x[1], x[0]))  # abundance desc, then id
        top_arr = arr_sorted[:top_k]
        for idx, (strain_taxid, abund, id_path) in enumerate(top_arr, start=1):
            # NOTE: sp is a string (species_taxid), and string concatenation is unambiguous.
            new_name = f"{sp}_{idx}"
            selected_records.append((sp, strain_taxid, new_name, id_path))

    if not selected_records:
        raise SystemExit("[Error] No records were found after selecting Top-K. Please check top_k against the input data.")

    # 7) 
    out_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["species_taxid", "strain_taxid", "new_strain_taxid", "id"])
        w.writerows(selected_records)

    # 
    if tmp_dir:
        print(f"[Note] Themis: tmp_dir={tmp_dir} has been ignored; symbolic links are no longer created in the current version.", file=sys.stderr)

    print(f"[Completed] The Top-K mapping table has been written: {out_tsv_path}")
    print(f"[Statistics] Number of species that passed the screening: {len(filtered_species)}; Final number of strain entries: {len(selected_records)}")
    if dropped_singletons:
        print(f"[Hint] Number of species discarded because single-strain abundance ≤ threshold:{dropped_singletons}", file=sys.stderr)
    if not_in_ref:
        print(f"[Note] {not_in_ref} predicted strains were not found in the Ref; the corresponding entries have been skipped.", file=sys.stderr)

    return out_tsv_path


# ======================

# ======================

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Top-K strains were selected for each species."
            "And apply abundance threshold filtering to single strain species (no more soft links are created)."
        )
    )
    p.add_argument("--ref_info", required=True,
                   help="Ref information table path (must include genome_ID, species_taxid, id columns).")
    p.add_argument("--ganon_species", required=True,
                   help="The path to ganon_species_abundance.txt.")
    p.add_argument("--ganon_strain", required=True,
                   help="The path to ganon_strain_abundance.txt (must include the abundance column).")
    p.add_argument("--top_k", type=int, required=True,
                   help="The number of Top-K strains retained for each species.")
    p.add_argument("--singleton_min_abund", type=float,
                   default=1.0000000000000001e-07,
                   help="Minimum abundance threshold for a single strain species (only those strictly greater than this value are retained).")
    p.add_argument("--out_tsv", default="ganon_species_strain_topk.tsv",
                   help="Output TSV filename (default: ganon_species_strain_topk.tsv).")
    p.add_argument("--tmp_dir", default=None,
                   help="For compatibility reasons, symbolic links are no longer created in the Themis version and can be ignored.")
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
