#!/bin/bash

set -e

while :
do
    cryspy -n
    LOG_LASTLINE=$(tail -n 1 log_cryspy)

    # 途中経過を表示
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') Last log: $LOG_LASTLINE"

    if [ "$LOG_LASTLINE" = "Done all structures!" ]; then
        echo "[INFO] All structures processed. Exiting."
        exit 0
    # ---------- for EA
    elif [[ "${LOG_LASTLINE:0:17}" = "Reached maxgen_ea" ]]; then
        echo "[INFO] EA max generation reached. Exiting."
        exit 0
    elif [ "$LOG_LASTLINE" = "EA is ready" ]; then
        echo "[INFO] EA ready. Running EA steps..."
        cryspy -n
        LOG_LASTLINE=$(tail -n 1 log_cryspy)
        echo "[INFO] After EA step: $LOG_LASTLINE"
        if [[ "${LOG_LASTLINE:0:17}" = "Reached maxgen_ea" ]]; then
            exit 0
        fi
        cryspy -n
    # ---------- for BO
    elif [[ "${LOG_LASTLINE:0:21}" = "Reached max_select_bo" ]]; then
        echo "[INFO] BO max selection reached. Exiting."
        exit 0
    elif [ "$LOG_LASTLINE" = "BO is ready" ]; then
        echo "[INFO] BO ready. Running BO steps..."
        cryspy -n
        LOG_LASTLINE=$(tail -n 1 log_cryspy)
        echo "[INFO] After BO step: $LOG_LASTLINE"
        if [[ "${LOG_LASTLINE:0:21}" = "Reached max_select_bo" ]]; then
            exit 0
        fi
        cryspy -n
    # ---------- for LAQA
    elif [ "$LOG_LASTLINE" = "LAQA is ready" ]; then
        echo "[INFO] LAQA ready. Running LAQA steps..."
        cryspy -n
        cryspy -n
    fi

    sleep 300  # seconds
done
