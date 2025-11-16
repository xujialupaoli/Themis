# themis/build.py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import os
import sys
from pathlib import Path
from .utils import run_cmd

GANON_BIN = os.environ.get("THEMIS_GANON_BIN", "ganon")


def build_custom(args):
    cmd = [GANON_BIN, "build-custom"] + list(args)
    run_cmd(cmd, "[Themis][build-custom]")


def main():
    if len(sys.argv) < 2:
        print("Usage: themis-build-custom [ganon build-custom options...]", file=sys.stderr)
        sys.exit(1)
    cmd = [GANON_BIN, "build-custom"] + sys.argv[1:]
    run_cmd(cmd, "[Themis][build-custom]")


if __name__ == "__main__":
    main()
