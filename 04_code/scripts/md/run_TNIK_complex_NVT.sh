#!/bin/bash

set -u

PROJECT="$HOME/Desktop/TNIK_IPF_Project"
BASE="$PROJECT/06_docking/molecular_dynamics/energy_minimization"
LOGDIR="$PROJECT/06_docking/molecular_dynamics/logs"
REPORTDIR="$PROJECT/06_docking/molecular_dynamics/analysis"
SYSTEMS=(
    "ponatinib_N32_H"
    "ponatinib_N37_H"
    "tepotinib_pH74_plus1"
)

mkdir -p "$LOGDIR" "$REPORTDIR"

REPORT="$REPORTDIR/complex_NVT_run_summary.txt"

echo "TNIK COMPLEX NVT SIMULATION SUMMARY" > "$REPORT"
echo "Started: $(date)" >> "$REPORT"
echo "========================================" >> "$REPORT"

if ! command -v gmx >/dev/null 2>&1; then
    echo "ERROR: GROMACS is not available. Activate tnik_md."
    exit 1
fi

for SYSTEM in "${SYSTEMS[@]}"; do

    echo
    echo "=============================================="
    echo "SYSTEM: $SYSTEM"
    echo "=============================================="

    DIR="$BASE/$SYSTEM"
    OUT="$DIR/nvt_500ps"
    LOG="$LOGDIR/${SYSTEM}_nvt_500ps.log"

    if [ ! -d "$DIR" ]; then
        echo "ERROR: System directory missing: $DIR"
        echo "$SYSTEM: FAILED - directory missing" >> "$REPORT"
        exit 1
    fi

    # Find the NVT TPR produced by the successful preprocessing step.
    TPR=""
    for CANDIDATE in \
        "$DIR/nvt_validation.tpr" \
        "$DIR/nvt.tpr" \
        "$DIR/nvt_posres.tpr"
    do
        if [ -s "$CANDIDATE" ]; then
            TPR="$CANDIDATE"
            break
        fi
    done

    if [ -z "$TPR" ]; then
        echo "ERROR: No expected NVT TPR found in $DIR"
        echo "$SYSTEM: FAILED - NVT TPR missing" >> "$REPORT"
        exit 1
    fi

    # Never overwrite an earlier trajectory.
    if [ -e "${OUT}.xtc" ] || [ -e "${OUT}.gro" ] || \
       [ -e "${OUT}.edr" ] || [ -e "${OUT}.log" ] || \
       [ -e "${OUT}.cpt" ]; then
        echo "ERROR: Existing NVT output detected for $SYSTEM."
        echo "Move or archive it before rerunning; nothing was overwritten."
        echo "$SYSTEM: NOT RUN - existing output detected" >> "$REPORT"
        exit 1
    fi

    echo "Input TPR: $TPR"
    echo "Starting NVT simulation..."

    gmx mdrun \
        -s "$TPR" \
        -deffnm "$OUT" \
        -ntmpi 1 \
        -v \
        > "$LOG" 2>&1

    RC=$?

    if [ "$RC" -ne 0 ]; then
        echo "ERROR: NVT failed for $SYSTEM. Inspect:"
        echo "$LOG"
        echo "$SYSTEM: FAILED - mdrun exit code $RC" >> "$REPORT"
        exit 1
    fi

    # Confirm essential outputs exist and are non-empty.
    MISSING=0
    for EXT in gro xtc edr log cpt; do
        if [ ! -s "${OUT}.${EXT}" ]; then
            echo "Missing output: ${OUT}.${EXT}"
            MISSING=1
        fi
    done

    if [ "$MISSING" -ne 0 ]; then
        echo "$SYSTEM: FAILED - incomplete output files" >> "$REPORT"
        exit 1
    fi

    echo "Checking trajectory..."
    gmx check -f "${OUT}.xtc" \
        > "$LOGDIR/${SYSTEM}_nvt_trajectory_check.log" 2>&1

    if [ "$?" -ne 0 ]; then
        echo "WARNING: trajectory check requires inspection."
        echo "$SYSTEM: SIMULATION FINISHED - trajectory check needs review" >> "$REPORT"
    else
        echo "$SYSTEM: NVT COMPLETED; trajectory file check passed" >> "$REPORT"
    fi

    echo "NVT completed for $SYSTEM"
    echo "Trajectory: ${OUT}.xtc"
    echo "Run log: $LOG"

done

echo "========================================" >> "$REPORT"
echo "Finished: $(date)" >> "$REPORT"

echo
echo "ALL THREE NVT RUNS HAVE FINISHED."
echo "Summary: $REPORT"
echo "Logs: $LOGDIR"
