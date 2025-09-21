#!/usr/bin/env python3

import subprocess
from pathlib import Path
import argparse
import json

SCRIPT_DIR = Path(__file__).parent


def get_aso_off_target(ref_fasta_file,aso_seqs_fasta_file,max_hamming_dis,test=False):
    script_path = SCRIPT_DIR / "get_aso_off_target.sh"
    json_path = SCRIPT_DIR / "all_transcripts.json"
    def get_aso_off_target_from_file(json_file=json_path):
        with open(json_file) as f:
            return json.load(f)

    def run_script(script_path: str, *args):
        """
        Run a script with given arguments, ignoring output.

        :param script_path: Path to the script (.sh, .py, etc.)
        :param args: Arguments to pass to the script
        """
        script_path = Path(script_path)
        if not script_path.exists():
            raise FileNotFoundError(f"Script not found: {script_path}")

        try:
            subprocess.run(
                [str(script_path), *args],
                check=True,
                stdout=subprocess.DEVNULL,  # suppress stdout
                stderr=subprocess.DEVNULL   # optionally suppress stderr
            )
        except subprocess.CalledProcessError as e:
            print(f"Script failed with code {e.returncode}")

    if test:
        print("Testing:")
        run_script(script_path,"-t")
        return()
    run_script(
        script_path,
        "-r", ref_fasta_file,
        "-q", aso_seqs_fasta_file,
        "-k", str(max_hamming_dis)
    )
    return get_aso_off_target_from_file()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ASO off-target prediction.")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA file path")
    parser.add_argument("-q", "--query", required=True, help="ASO sequences FASTA file path")
    parser.add_argument("-k", "--max_hamming", type=int, required=True, help="Maximum Hamming distance")

    args = parser.parse_args()

    res = get_aso_off_target(args.ref, args.query, args.max_hamming)
    print(res)
