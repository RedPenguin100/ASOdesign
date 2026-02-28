#!/bin/bash
set -e

METHOD="h"
MAX_MISMATCHES=3
OUTPUT_PREFIX="res"

# ---------------------------
# Usage / Help
# ---------------------------
usage() {
    cat <<EOF
Usage: $0 [options]

This script performs alignment using Bowtie1 with Hamming distance.
(Bowtie2 / Levenshtein distance is currently disabled.)

Options:
  -r <ref.fa>        Reference FASTA file (required)
  -q <reads.fa>      Query FASTA file (required)
  -o <prefix>        Output prefix (default: res)
  -m <int>           Max mismatches (default: 3)
  -h                 Show this help message and exit

Example:
  $0 -r ref.fa -q reads.fa -o results -m 2
EOF
    exit 0
}

# ---------------------------
# Parse command-line arguments
# ---------------------------
while getopts "r:q:o:m:h" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        q) QUERY="$OPTARG" ;;
        o) OUTPUT_PREFIX="$OPTARG" ;;
        m) MAX_MISMATCHES="$OPTARG" ;;
        h) usage ;;
        *) echo "Invalid option"; usage ;;
    esac
done

# ---------------------------
# Validate required arguments
# ---------------------------
for var in REF QUERY OUTPUT_PREFIX; do
    if [ -z "${!var}" ]; then
        echo "Missing argument: $var"
        usage
    fi
done

# ---------------------------
# Bowtie1 / Bowtie2 index functions
# ---------------------------
build_index_bowtie1() {
    echo "Building Bowtie1 index..."
    bowtie-build "$REF" "${OUTPUT_PREFIX}_index" \
        || { echo "bowtie-build failed"; exit 1; }
}

# (disabled: works poorly)
# build_index_bowtie2() {
#     echo "Building Bowtie2 index..."
#     bowtie2-build "$REF" "${OUTPUT_PREFIX}_index" \
#         || { echo "bowtie2-build failed"; exit 1; }
# }

# ---------------------------
# Alignment functions
# ---------------------------
align_bowtie1() {
    SAM_FILE="${OUTPUT_PREFIX}.sam"
    echo "Running Bowtie1..."
    bowtie -v "$MAX_MISMATCHES" -a --best -f -x "${OUTPUT_PREFIX}_index" -S "$QUERY" > "$SAM_FILE" \
        || { echo "Bowtie1 alignment failed"; exit 1; }
    echo "SAM saved for debugging: $SAM_FILE"
}

# (disabled: works poorly)
# align_bowtie2() {
#     SAM_FILE="${OUTPUT_PREFIX}.sam"
#     echo "Running Bowtie2..."
#     bowtie2 -a --very-sensitive-local -x "${OUTPUT_PREFIX}_index" -f -U "$QUERY" \
#         -S "$SAM_FILE" --mp 1,1 --rdg 1,1 --rfg 1,1 --score-min G,20,8
#     echo "SAM saved for debugging: $SAM_FILE"
# }

# ---------------------------
# Convert SAM → sorted BAM
# ---------------------------
sam_to_bam() {
    SAM_FILE="${OUTPUT_PREFIX}.sam"
    BAM_FILE="${OUTPUT_PREFIX}.sorted.bam"
    echo "Converting SAM → BAM..."
    samtools view -bS "$SAM_FILE" | samtools sort -o "$BAM_FILE" - \
        || { echo "BAM conversion failed"; exit 1; }
    echo "BAM ready: $BAM_FILE"
}

# ---------------------------
# Index BAM
# ---------------------------
index_bam() {
    BAM_FILE="${OUTPUT_PREFIX}.sorted.bam"
    if command -v samtools &>/dev/null; then
        echo "Indexing BAM file..."
        samtools index "$BAM_FILE" || { echo "BAM indexing failed"; exit 1; }
        echo "BAM index ready: ${BAM_FILE}.bai"
    else
        echo "samtools not found. Skipping indexing."
    fi
}

# ---------------------------
# Main workflow
# ---------------------------
if [[ "$METHOD" == "h" || "$METHOD" == "hamming" ]]; then
    build_index_bowtie1
    align_bowtie1
elif [[ "$METHOD" == "l" || "$METHOD" == "levenstein" ]]; then
    build_index_bowtie2
    align_bowtie2
else
    echo "Unknown method: $METHOD"
    exit 1
fi

sam_to_bam
index_bam
