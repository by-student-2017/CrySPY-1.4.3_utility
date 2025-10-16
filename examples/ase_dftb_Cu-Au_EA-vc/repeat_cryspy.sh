#!/bin/bash

set -e

while :
do
    cryspy -n

    echo "======================================"
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') CrySPY log:"
    cat log_cryspy
    echo "======================================"

    LOG_LASTLINE=$(tail -n 1 log_cryspy)

    if [ "$LOG_LASTLINE" = "Done all structures!" ]; then
        echo "[INFO] All structures processed. Exiting."
        exit 0
    elif [[ "${LOG_LASTLINE:0:17}" = "Reached maxgen_ea" ]]; then
        echo "[INFO] EA max generation reached. Exiting."
        exit 0
    elif [ "$LOG_LASTLINE" = "EA is ready" ]; then
        echo "[INFO] EA ready. Running EA steps..."
        cryspy -n
        cat log_cryspy
        LOG_LASTLINE=$(tail -n 1 log_cryspy)
        if [[ "${LOG_LASTLINE:0:17}" = "Reached maxgen_ea" ]]; then
            exit 0
        fi
        cryspy -n
    elif [[ "${LOG_LASTLINE:0:21}" = "Reached max_select_bo" ]]; then
        echo "[INFO] BO max selection reached. Exiting."
        exit 0
    elif [ "$LOG_LASTLINE" = "BO is ready" ]; then
        echo "[INFO] BO ready. Running BO steps..."
        cryspy -n
        cat log_cryspy
        LOG_LASTLINE=$(tail -n 1 log_cryspy)
        if [[ "${LOG_LASTLINE:0:21}" = "Reached max_select_bo" ]]; then
            exit 0
        fi
        cryspy -n
    elif [ "$LOG_LASTLINE" = "LAQA is ready" ]; then
        echo "[INFO] LAQA ready. Running LAQA steps..."
        cryspy -n
        cryspy -n
    fi

    sleep 300
done
