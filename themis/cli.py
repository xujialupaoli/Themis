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
        "themis",
        description="Themis:  metagenomic profiling ."
    )
    sub = p.add_subparsers(dest="subcmd", metavar="<command>")

   
    p_build = sub.add_parser(
        "build-custom",
        help="Wrap upstream 'ganon build-custom' as direct 'ganon-build' (pass-through)."
    )
    p_build.add_argument("ganon_args", nargs=argparse.REMAINDER,
                         help="Arguments passed directly to 'ganon-build'.")


    p_prof = sub.add_parser("profile", help="Run Themis profiling pipeline")
    p_prof.add_argument("-r", "--reads", action="append", required=True,
                        help="Input reads. Paired: -r R1 -r R2. Single: -r R.")
    p_prof.add_argument("--single", action="store_true", help="Single-end mode.")
    p_prof.add_argument("--db-prefix", required=True,
                        help="Ganon DB prefix built by 'themis build-custom'.")
    p_prof.add_argument("--ref-info", required=True,
                        help="Reference info TSV (RefDB_xxx_genomes_info.txt).")
    p_prof.add_argument("--out", required=True, help="Output directory.")
    p_prof.add_argument("--threads", type=int, default=8, help="Number of threads.")
    p_prof.add_argument("-k", "--kmer", type=int, default=31, help="k-mer size for ggcat.")

    return p


def main():
    p = build_parser()
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(0)

    subcmd = sys.argv[1]
    if subcmd == "build-custom":
        argv = sys.argv[2:]
        subcmd_build_custom(argv)
    elif subcmd == "profile":
        args = p.parse_args()
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
