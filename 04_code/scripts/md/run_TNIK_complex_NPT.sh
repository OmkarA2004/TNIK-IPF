#!/bin/bash

cd ~/Desktop/TNIK_IPF_Project || exit 1

for SYS in ponatinib_N32_H ponatinib_N37_H tepotinib_pH74_plus1
do
    BASE="06_docking/molecular_dynamics/energy_minimization/$SYS"
    OUT="$BASE/npt_1ns"
    TPR="$OUT/npt_1ns.tpr"

    echo "========================================"
    echo "STARTING NPT: $SYS"
    echo "========================================"

    if [ ! -f "$TPR" ]; then
        echo "ERROR: TPR missing for $SYS"
        exit 1
    fi

    gmx mdrun -s "$TPR" -deffnm "$OUT/npt_1ns" -ntmpi 1 -v

    STATUS=$?

    if [ "$STATUS" -ne 0 ]; then
        echo "NPT FAILED: $SYS"
        exit 1
    fi

    echo "NPT COMPLETED: $SYS"
done

echo "ALL THREE NPT SIMULATIONS COMPLETED"
