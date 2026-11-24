#!/bin/bash

cd ~/Desktop/TNIK_IPF_Project || exit 1

OUTDIR="06_docking/molecular_dynamics/analysis/npt_stability"
mkdir -p "$OUTDIR"

for SYS in ponatinib_N32_H ponatinib_N37_H tepotinib_pH74_plus1
do
    EDR="06_docking/molecular_dynamics/energy_minimization/$SYS/npt_1ns/npt_1ns.edr"

    echo "========================================"
    echo "EXTRACTING NPT ENERGY DATA: $SYS"
    echo "========================================"

    if [ ! -f "$EDR" ]; then
        echo "ERROR: Energy file missing: $EDR"
        exit 1
    fi

    gmx energy \
        -f "$EDR" \
        -o "$OUTDIR/${SYS}_npt_energy.xvg" \
        <<'SELECTION'
Temperature
Pressure
Density
Volume
Potential
Total-Energy
0
SELECTION

    if [ $? -ne 0 ]; then
        echo "ERROR extracting energy data for $SYS"
        exit 1
    fi

    echo "EXTRACTION COMPLETED: $SYS"
done

echo "ALL THREE NPT ENERGY EXTRACTIONS COMPLETED"
