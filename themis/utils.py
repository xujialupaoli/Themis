# themis/utils.py
import subprocess
import sys
from pathlib import Path

def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)

def run_cmd(cmd, log_prefix="[Themis] ", *, echo=True, silence=False):
    if echo and log_prefix:
        print(log_prefix + " " + " ".join(cmd), file=sys.stderr)
    stdout = subprocess.DEVNULL if silence else None
    stderr = subprocess.DEVNULL if silence else None
    subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)
