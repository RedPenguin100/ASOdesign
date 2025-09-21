#!/bin/bash
set -e

# ---------------------------
# Default parameters
# ---------------------------
TEST_MODE=0  # 0 = off, 1 = on

# ---------------------------
# Usage / Help
# ---------------------------
usage() {
    echo "Usage: $0 -r <reference.fa> -q <query.fa> -k <max_hamming_distance> [-t]"
    echo
    echo "  -r    Reference file (FASTA)"
    echo "  -q    Query file (FASTA)"
    echo "  -k    Maximum Hamming distance allowed"
    echo "  -t    Test mode (does not execute main commands)"
    echo "  -h    Show this help message and exit"
    exit 1
}

# ---------------------------
# Parse arguments
# ---------------------------
while getopts "r:q:k:th" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        q) QUERY="$OPTARG" ;;
        k) K="$OPTARG" ;;
        t) TEST_MODE=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done

# ---------------------------
# Check required args
# ---------------------------
if [[ $TEST_MODE -eq 0 ]]; then
    # Only check required arguments in normal mode
    if [[ -z "$REF" || -z "$QUERY" || -z "$K" ]]; then
        echo "Error: Missing required arguments."
        usage
    fi
fi

# ---------------------------
# Run main commands
# ---------------------------
if [[ $TEST_MODE -eq 1 ]]; then
    echo "pass_0"
    ./create_bam.sh -r "fasta_files/ref.fasta" -q "fasta_files/aso_reads.fasta" -o _res_  -m 2
    echo "pass_1"
    ./find_matches_in_bam.py -b _res_.sorted.bam -s - -mm 2 -o ./test
    rm _res_*
else
    ./create_bam.sh -r $REF -q $QUERY -o _res_ -m $K
    ./find_matches_in_bam.py -b _res_.sorted.bam -s - -mm $K -o . --single-json
    rm _res_*
fi
