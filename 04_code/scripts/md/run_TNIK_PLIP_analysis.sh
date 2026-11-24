#!/bin/bash
set -euo pipefail

ROOT="$HOME/Desktop/TNIK_IPF_Project"
INPUT="$ROOT/06_docking/interaction_analysis/PLIP_complexes"
OUTPUT="$ROOT/06_docking/interaction_analysis/PLIP_results"

mkdir -p "$OUTPUT"

run_plip() {
    local NAME="$1"
    local COMPLEX="$2"
    local OUTDIR="$OUTPUT/$NAME"

    mkdir -p "$OUTDIR"

    echo ""
    echo "=================================================="
    echo "Running PLIP: $NAME"
    echo "=================================================="

    plip \
        --file "$COMPLEX" \
        --out "$OUTDIR" \
        --name "$NAME" \
        --xml \
        --txt \
        --maxthreads 1 \
        --verbose

    echo "Completed: $NAME"
}

run_plip \
    "ponatinib_seed11" \
    "$INPUT/ponatinib_seed11_complex.pdb"

run_plip \
    "tepotinib_seed22" \
    "$INPUT/tepotinib_seed22_complex.pdb"

echo ""
echo "=================================================="
echo "PLIP ANALYSIS COMPLETE"
echo "=================================================="

find "$OUTPUT" -maxdepth 3 -type f -print | sort

echo ""
echo "Report sizes:"
du -sh "$OUTPUT"/*
