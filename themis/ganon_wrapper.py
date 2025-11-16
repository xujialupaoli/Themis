#!/usr/bin/env python3
# -*- coding: utf-8 -*-



from __future__ import annotations
import os
import sys
import subprocess
from pathlib import Path
from typing import List, Optional


try:
    from .utils import run_cmd, ensure_dir
except Exception:
    def run_cmd(cmd: List[str], log_prefix: str = ""):
        if log_prefix:
            print(f"{log_prefix} $", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
    def ensure_dir(p: Path):
        Path(p).mkdir(parents=True, exist_ok=True)


_gsp_run = None
_gst_run = None
try:
    from themis_scripts import ganon_species_process as _gsp
    from themis_scripts import ganon_strain_process as _gst
    _gsp_run = getattr(_gsp, "run", None)
    _gst_run = getattr(_gst, "run", None)
except Exception:
    pass

GANON = os.environ.get("THEMIS_GANON_BIN", "ganon")  

def _check_ganon_available():
    from shutil import which
    if which(GANON) is None:
        raise SystemExit(
            f"[Themis][error] cannot find '{GANON}' in PATH. "
            "Please `pip install ganon` in the same environment."
        )

# =========================
# classify
# =========================

def run_classify_paired(
    db_prefix: str,
    read1: str,
    read2: str,
    out_prefix: str,
    threads: int,
    report_type: str,
):
    """
    双端： ganon classify -d DB -p R1 R2 --output-prefix OUT --report-type TYPE -t N
    """
    _check_ganon_available()
    ensure_dir(Path(out_prefix).parent)
    cmd = [
        GANON, "classify",
        "-d", db_prefix,                 
        "-p", read1, read2,
        "--output-prefix", out_prefix,
        "--report-type", report_type,
        "-t", str(threads),
        "--quiet",
    ]
    run_cmd(cmd, echo=False, silence=True)

def run_classify_single(
    db_prefix: str,
    read: str,
    out_prefix: str,
    threads: int,
    report_type: str,
):
    
    _check_ganon_available()
    ensure_dir(Path(out_prefix).parent)
    cmd = [
        GANON, "classify",
        "-d", db_prefix,                 
        "-s", read,
        "--output-prefix", out_prefix,
        "--report-type", report_type,
        "-t", str(threads),
        "--quiet", 
    ]
    run_cmd(cmd, echo=False, silence=True)

# =========================

# =========================

def run_report_and_postprocess(
    db_prefix: str,
    classify_prefix: str,
    tre_out: str,
    species_out: str,
    strain_out: str,
    predict_spy_out: Optional[str] = None,
):
    
    _check_ganon_available()

    classify_prefix = Path(classify_prefix)
    tre_out = Path(tre_out)
    species_out = Path(species_out)
    strain_out = Path(strain_out)

    ensure_dir(tre_out.parent)
    ensure_dir(species_out.parent)
    ensure_dir(strain_out.parent)

    rep_in = f"{classify_prefix}.rep"
    if not Path(rep_in).exists():
        raise SystemExit(f"[Themis][error] cannot find REP: {rep_in}")

    # 1) ganon report
    report_prefix = tre_out.with_suffix("")  # /.../tax_profile
    cmd = [
        GANON, "report",
        "-i", rep_in,
        "--db-prefix", db_prefix,
        "--output-prefix", str(report_prefix),
        "--report-type", "abundance",
        "-r", "all",
        "--quiet",
    ]
    run_cmd(cmd, echo=False, silence=True)

    tre_file = f"{report_prefix}.tre"
    if not Path(tre_file).exists():
        raise SystemExit(f"[Themis][error] ganon report did not create: {tre_file}")

    # 2)  species / strain
    _parse_to_outputs_using_tre(
        tre_file=tre_file,
        species_out=str(species_out),
        strain_out=str(strain_out),
        prefer_scripts=True
    )

    # 3)  predict_spy
    _maybe_write_predict_spy(species_out, predict_spy_out)

# ---------  .tre ->  ---------

def _parse_to_outputs_using_tre(tre_file: str, species_out: str, strain_out: str, prefer_scripts: bool = True):
    if prefer_scripts and _gsp_run is not None:
        _gsp_run(tre_file=tre_file, output_path=species_out)
    else:
        _parse_species_from_tre(tre_file, species_out)

    if prefer_scripts and _gst_run is not None:
        _gst_run(tre_file=tre_file, output_path=strain_out)
    else:
        _parse_strain_from_tre(tre_file, strain_out)

def _maybe_write_predict_spy(species_out: Path, predict_spy_out: Optional[str]):
    if not predict_spy_out:
        return
    with open(species_out, "r", encoding="utf-8", errors="ignore") as fin, \
         open(predict_spy_out, "w", encoding="utf-8") as fout:
        fin.readline()  # skip header
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                fout.write(f"{parts[0]}\t{parts[1]}\n")

# ---------  .tre ---------

def _parse_species_from_tre(tre_file: str, out_path: str):
    items = []
    with open(tre_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("species"):
                continue
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 9:
                continue
            sid = toks[2].split("|")[-1].strip()
            try:
                abund = float(toks[8]) / 100.0
            except Exception:
                continue
            if sid:
                items.append((sid, abund))
    items.sort(key=lambda x: x[1], reverse=True)
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("speciesID\tabundance\n")
        for sid, abund in items:
            out.write(f"{sid}\t{abund}\n")

def _parse_strain_from_tre(tre_file: str, out_path: str):
    items = []
    with open(tre_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("strain"):
                continue
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 9:
                continue
            sid = toks[2].split("|")[-1].strip() if toks[2] else ""
            if not sid:
                sid = toks[3].strip()
            try:
                abund = float(toks[8]) / 100.0
            except Exception:
                continue
            if sid:
                items.append((sid, abund))
    items.sort(key=lambda x: x[1], reverse=True)
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("strain_taxid\tabundance\n")
        for sid, abund in items:
            out.write(f"{sid}\t{abund}\n")
