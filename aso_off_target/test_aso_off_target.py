#!/usr/bin/env python3
import json
import sys
from pathlib import Path
import subprocess
from get_aso_off_target import get_aso_off_target, SCRIPT_DIR

def compare_json(dir1,dir2):
    for file1 in dir1.glob("*.json"):
        file2 = dir2 / file1.name
        if file2.exists():
            d1 = json.load(open(file1))
            d2 = json.load(open(file2))

            # Normalize: missing keys treated as empty list
            all_keys = set(d1.keys()) | set(d2.keys())
            nd1 = {k: d1.get(k, []) for k in all_keys}
            nd2 = {k: d2.get(k, []) for k in all_keys}

            diffs = []
            for k in all_keys:
                if nd1[k] != nd2[k]:
                    diffs.append(f"  Key '{k}': {nd1[k]} != {nd2[k]}")

            return len(diffs)==0

def test_script_json_outputs():
    test_output_dir = SCRIPT_DIR / "test"
    reference_dir = SCRIPT_DIR / "test/hamming_distance_2"
    get_aso_off_target("","",0,True)
    # Compare generated JSONs with reference
    assert compare_json(test_output_dir, reference_dir), "JSON outputs differ!"

if __name__ == "__main__":
    test_script_json_outputs()
