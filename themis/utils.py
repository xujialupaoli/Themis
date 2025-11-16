# themis/utils.py
import subprocess
import sys
from pathlib import Path

def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)

def run_cmd(cmd, log_prefix="[Themis] ", *, echo=True, silence=False):
    """
    echo=False  -> 不打印我们自己的提示行（不暴露 '[Themis][...] ...'）
    silence=True -> 把外部程序的 stdout/stderr 全部丢到 /dev/null（彻底静音）
    两个参数都不传时，行为与旧版 100% 一致。
    """
    if echo and log_prefix:
        print(log_prefix + " " + " ".join(cmd), file=sys.stderr)
    stdout = subprocess.DEVNULL if silence else None
    stderr = subprocess.DEVNULL if silence else None
    subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)
