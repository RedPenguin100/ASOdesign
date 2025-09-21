#!/usr/bin/env python3
import pysam
import argparse
import json
from collections import defaultdict
import os
from Bio.Seq import Seq


def main():
    parser = argparse.ArgumentParser(description="Extract ASO reads info per transcript from BAM")
    parser.add_argument("-b", "--bam", required=True, help="Input sorted BAM file")
    parser.add_argument("-s", "--strand", choices=["+", "-"], required=True, help="Which strand to process (+ or -)")
    parser.add_argument("-mm", "--max-mismatch", type=int, default=None, help="Maximum mismatch to catch; higher are ignored")
    parser.add_argument("-o", "--outdir", default="transcript_jsons", help="Output directory for transcript JSONs")
    parser.add_argument("--single-json", action="store_true", help="Save all transcripts into one JSON file instead of separate files")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    bamfile = pysam.AlignmentFile(args.bam, "rb")

    # transcript -> ASO -> list of [s, e, mm]
    transcript_aso_reads = defaultdict(lambda: defaultdict(list))

    for read in bamfile.fetch():
        # Determine strand
        strand = "-" if read.is_reverse else "+"
        if strand != args.strand:
            continue

        aso_seq = read.query_sequence
        if args.strand == "-" and args.single_json:
            aso_seq = str(Seq(aso_seq).reverse_complement())
        transcript = read.reference_name
        start = read.reference_start
        end = read.reference_end

        mismatches = read.get_tag("NM") if read.has_tag("NM") else 0
        if args.max_mismatch is not None and mismatches > args.max_mismatch:
            continue

        transcript_aso_reads[transcript][aso_seq].append([start, end, mismatches])

    if args.single_json:
        # Invert structure to {aso_seq: {transcript: [[s, e, mm], ...]}}
        aso_dict = defaultdict(lambda: defaultdict(list))
        for transcript, aso_reads in transcript_aso_reads.items():
            for aso_seq, records in aso_reads.items():
                aso_dict[aso_seq][transcript].extend(records)

        out_file = os.path.join(args.outdir, "all_transcripts.json")
        with open(out_file, "w") as f:
            json.dump(aso_dict, f, indent=2)
        print(f"Saved all transcripts into single JSON: {out_file}")
    else:
        for transcript, aso_dict in transcript_aso_reads.items():
            out_file = os.path.join(args.outdir, f"{transcript}.json")
            with open(out_file, "w") as f:
                json.dump(aso_dict, f, indent=2)
        print(f"Saved {len(transcript_aso_reads)} transcript JSONs in '{args.outdir}'")

if __name__ == "__main__":
    main()
