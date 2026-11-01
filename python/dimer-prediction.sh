#!/bin/bash
# ── settings — edit here ─────────────────────────────────────────────────
INPUT="../primers.xlsx"    # input file (.xlsx or .csv)
MODE="multiplex"                    # analysis mode: multiplex or pairs
FORMAT="xlsx"                   # output table format: xlsx  or  tsv
THRESHOLD="-1.00"               # dimer threshold in kcal/mol (lower = dimer)
SORT="asc"                      # sort order: asc | desc | input
TYPE="heterodimers"               # interactions to report: heterodimers | homodimers | both
# ─────────────────────────────────────────────────────────────────────────

# allow passing the input file as a command-line argument instead:
#   ./dimer-prediction.sh ../my-primers.xlsx
INPUT="${1:-$INPUT}"

# derive output paths from the input file name
STEM="${INPUT%.*}"                      # strip extension  →  ../primers
OUTPUT="${STEM}_results.${FORMAT}"      # e.g. ../primers_results.xlsx
ALIGNMENTS="${STEM}_alignments.txt"     # e.g. ../primers_alignments.txt

python make_primer_alignments.py \
    --mode       "$MODE"       \
    --input      "$INPUT"      \
    --output     "$OUTPUT"     \
    --alignments "$ALIGNMENTS" \
    --threshold  "$THRESHOLD"  \
    --sort       "$SORT"       \
    --type       "$TYPE"