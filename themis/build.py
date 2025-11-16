# themis/build.py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Themis build:
封装 ganon build-custom，给用户一个统一命令：
  themis build-custom [ganon 原来的参数]

你修改过的 ganon 仍然是核心，只是用户不必自己记一长串命令。
"""

import os
import sys
from pathlib import Path
from .utils import run_cmd

GANON_BIN = os.environ.get("THEMIS_GANON_BIN", "ganon")


def build_custom(args):
    """
    args: 直接是命令行参数列表，例如:
      ["-i", "genomes/", "-d", "themis_db", "-x", "ncbi", ...]
    内部就是：
      ganon build-custom + args
    """
    cmd = [GANON_BIN, "build-custom"] + list(args)
    run_cmd(cmd, "[Themis][build-custom]")


def main():
    # 直接把参数透传给 ganon build-custom
    # 用法:
    #   themis-build-custom [ganon build-custom 原来的所有参数]
    if len(sys.argv) < 2:
        print("Usage: themis-build-custom [ganon build-custom options...]", file=sys.stderr)
        sys.exit(1)
    cmd = [GANON_BIN, "build-custom"] + sys.argv[1:]
    run_cmd(cmd, "[Themis][build-custom]")


if __name__ == "__main__":
    main()
