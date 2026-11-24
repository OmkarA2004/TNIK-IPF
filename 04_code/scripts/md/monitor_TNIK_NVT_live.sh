#!/bin/bash

# ============================================================
# TNIK-IPF PROJECT
# LIVE GROMACS NVT SIMULATION MONITOR
# Monitoring only: does not modify or control simulations
# ============================================================

PROJECT="$HOME/Desktop/TNIK_IPF_Project"

MD="$PROJECT/06_docking/molecular_dynamics"

SYSTEMS=(
    "ponatinib_N32_H"
    "ponatinib_N37_H"
    "tepotinib_pH74_plus1"
)

TARGET_PS=500
REFRESH=5

get_value() {
    local file="$1"
    local pattern="$2"

    if [ -f "$file" ]; then
        grep -E "$pattern" "$file" | tail -n 1
    fi
}

while true; do

    clear

    echo "============================================================"
    echo "           TNIK-IPF MOLECULAR DYNAMICS"
    echo "                  LIVE NVT MONITOR"
    echo "============================================================"

    echo
    echo "Date and time: $(date)"
    echo "Target simulation duration: ${TARGET_PS} ps"
    echo "Refresh interval: ${REFRESH} seconds"
    echo
    echo "------------------------------------------------------------"

    ACTIVE_PROCESS=$(pgrep -f '[g]mx mdrun' | tr '\n' ' ')

    if [ -n "$ACTIVE_PROCESS" ]; then
        echo "GROMACS PROCESS: ACTIVE"
        echo "Process ID(s): $ACTIVE_PROCESS"
    else
        echo "GROMACS PROCESS: NO ACTIVE MATCHING PROCESS"
    fi

    echo
    echo "============================================================"

    for SYSTEM in "${SYSTEMS[@]}"; do

        DIR="$MD/energy_minimization/$SYSTEM"

        LOG="$DIR/nvt_500ps.log"
        GRO="$DIR/nvt_500ps.gro"
        CPT="$DIR/nvt_500ps.cpt"
        XTC="$DIR/nvt_500ps.xtc"
        EDR="$DIR/nvt_500ps.edr"

        echo
        echo "############################################################"
        echo "SYSTEM: $SYSTEM"
        echo "############################################################"

        if [ ! -f "$LOG" ]; then
            echo "STATUS: WAITING / LOG NOT CREATED"
            continue
        fi

        # ----------------------------------------------------
        # Latest recorded simulation time
        # ----------------------------------------------------

        TIME_LINE=$(grep -E '^[[:space:]]*Time[[:space:]]' "$LOG" | tail -n 1)

        STEP_LINE=$(grep -E '^[[:space:]]*Step[[:space:]]' "$LOG" | tail -n 1)

        SIM_TIME=$(echo "$TIME_LINE" | awk '{print $2}')
        SIM_STEP=$(echo "$STEP_LINE" | awk '{print $2}')

        echo
        echo "LATEST RECORDED PROGRESS"

        if [ -n "$SIM_TIME" ]; then

            PERCENT=$(awk -v t="$SIM_TIME" -v target="$TARGET_PS" \
                'BEGIN {
                    p=(t/target)*100
                    if (p>100) p=100
                    if (p<0) p=0
                    printf "%.2f",p
                }')

            REMAINING=$(awk -v t="$SIM_TIME" -v target="$TARGET_PS" \
                'BEGIN {
                    r=target-t
                    if (r<0) r=0
                    printf "%.2f",r
                }')

            echo "Simulation time : $SIM_TIME ps"
            echo "Simulation step : ${SIM_STEP:-Not recorded}"
            echo "Target time     : $TARGET_PS ps"
            echo "Completion      : $PERCENT %"
            echo "Remaining time  : $REMAINING ps"

            BAR=$(awk -v p="$PERCENT" 'BEGIN {
                n=int(p/2)
                for(i=0;i<n;i++) printf "#"
                for(i=n;i<50;i++) printf "-"
            }')

            echo
            echo "Progress: [$BAR] $PERCENT%"

        else
            echo "No simulation time recorded yet."
        fi

        # ----------------------------------------------------
        # Determine whether this system has finished
        # ----------------------------------------------------

        if grep -q "Finished mdrun" "$LOG"; then

            STATUS="COMPLETED"

        elif pgrep -f "gmx mdrun.*$SYSTEM" >/dev/null 2>&1; then

            STATUS="RUNNING"

        elif [ -f "$CPT" ] || [ -f "$XTC" ]; then

            STATUS="NOT CURRENTLY RUNNING / OUTPUT EXISTS"

        else

            STATUS="WAITING OR NOT STARTED"

        fi

        echo
        echo "SIMULATION STATUS: $STATUS"

        # ----------------------------------------------------
        # Output files
        # ----------------------------------------------------

        echo
        echo "OUTPUT FILES"

        for FILE in "$LOG" "$GRO" "$CPT" "$XTC" "$EDR"; do

            if [ -f "$FILE" ]; then

                SIZE=$(du -h "$FILE" | awk '{print $1}')
                MODIFIED=$(stat -f "%Sm" -t "%H:%M:%S" "$FILE")

                echo "$(basename "$FILE") : $SIZE | Updated: $MODIFIED"

            else

                echo "$(basename "$FILE") : Not created yet"

            fi

        done

        # ----------------------------------------------------
        # Latest temperature and energy records
        # ----------------------------------------------------

        echo
        echo "LATEST LOG ACTIVITY"

        tail -n 6 "$LOG"

        echo
        echo "------------------------------------------------------------"

    done

    echo
    echo "============================================================"
    echo "Automatic refresh in $REFRESH seconds."
    echo "Press CONTROL + C to exit this monitor only."
    echo "Your GROMACS simulations will continue running."
    echo "============================================================"

    sleep "$REFRESH"

done
