#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from pathlib import Path

from . import profile as profile_mod
from .utils import run_cmd



def _ganon():
    return os.environ.get("THEMIS_GANON_BIN", "ganon")

def subcmd_build_custom(argv):
    cmd = [_ganon(), "build-custom"] + list(argv)
    run_cmd(cmd, echo=False)


def subcmd_profile(args):
    if not args.reads:
        raise SystemExit("[Themis] --reads/-r .Paired: -r R1 -r R2；Single: -r R.")

    profile_mod.run(
        reads=args.reads,
        single=args.single,
        db_prefix=args.db_prefix,
        out_prefix=args.out,
        report_type="abundance",
        file_info=args.ref_info,
        threads=args.threads,
        k=args.kmer,
    )


def build_parser():
    p = argparse.ArgumentParser(
        prog="themis",
        description="Themis: a robust and accurate species-level metagenomic profiler."
    )
    sub = p.add_subparsers(dest="subcmd", metavar="<command>")

   
    p_build = sub.add_parser(
        "build-custom",
        help="Build custom themis databases."
    )
    p_build.add_argument("ganon_args", nargs=argparse.REMAINDER,
                         help="Arguments passed directly to 'ganon build-cutstom'.")


    p_prof = sub.add_parser("profile", help="Profile reads against custom databases.")
    p_prof.add_argument("-r", "--reads", action="append", required=True,
                        help=("For paired-end data, specify mates consecutively: -r R1.fq -r R2.fq. "
                        "For single-end data, use --single and give one -r per file. "))
    p_prof.add_argument("--single", action="store_true", help="Treat input as single-end reads. ")
    p_prof.add_argument("--db-prefix", required=True,
                        help="Database input prefix.")
    p_prof.add_argument("--ref-info", required=True,
                        help=("Tab-separated reference metadata file. Fields: "
                        "strain_name <tab> strain_taxid <tab> species_taxid "
                        "<tab> species_name <tab> genome_path. "
                        "strain_name and strain_taxid must be unique."))
    p_prof.add_argument("--out", required=True, help="Output directory for profiling results.")
    p_prof.add_argument("--threads", type=int, default=8, help="Number of threads.")
    p_prof.add_argument("-k", "--kmer", type=int, default=31, help="k-mer size used in the ccDBG-based profiling step.")

    return p


def main():
    p = build_parser()
    argv = sys.argv[1:]

    
    if not argv:
        p.print_help()
        sys.exit(0)

    
    if argv[0] in ("-h", "--help"):
        p.print_help()
        sys.exit(0)

    
    subcmd = argv[0]

    if subcmd == "build-custom":
        
        subcmd_build_custom(argv[1:])
    elif subcmd == "profile":
        
        args = p.parse_args(argv)
        subcmd_profile(args)
    else:
        
        p.print_help()
        sys.exit(1)


def build_custom_main():
    subcmd_build_custom(sys.argv[1:])

def profile_main():
    p = build_parser()
    fake_argv = ["themis", "profile"] + sys.argv[1:]
    args = p.parse_args(fake_argv[1:])
    subcmd_profile(args)
