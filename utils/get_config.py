#!/usr/bin/env python3
"""Print a value from config/config.yaml for use in bash scripts.

    cfg() { python3 "$REPO_DIR/utils/get_config.py" "$@"; }
    OUTDIR=$(cfg path gatk_dir)          # absolute path under `paths:`
    CPUS=$(cfg value gatk.cpus)          # any value, dotted key
    SAMPLES=($(cfg samples gatk))        # sample_id column of a sample sheet
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.config import cfg, cfg_path, sample_ids  # noqa: E402

kind, key = sys.argv[1], sys.argv[2]
if kind == "path":
    print(cfg_path(key))
elif kind == "samples":
    print(" ".join(sample_ids(key)))
elif kind == "value":
    v = cfg
    for k in key.split("."):
        v = v[k]
    print(" ".join(map(str, v)) if isinstance(v, list) else v)
else:
    sys.exit(f"unknown kind: {kind} (use path | value | samples)")
