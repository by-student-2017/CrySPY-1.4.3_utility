#!/bin/bash

MAX_REPEAT=10
declare -A id_count
SKIP_LOG="skipped_ids.log"

echo "[INFO] Monitoring CrySPY jobs. MAX_REPEAT=$MAX_REPEAT"
echo "[INFO] Skipped IDs will be logged in $SKIP_LOG"

while :
do
    cryspy -n || true

    echo "======================================"
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') CrySPY log:"
    cat log_cryspy
    echo "======================================"

    LOG_LASTLINE=$(tail -n 1 log_cryspy)

    # ID抽出: 「ID」の次のフィールドを取得、末尾のコロン削除、重複除外
    IDS=$(grep "still queueing or running" log_cryspy | \
          awk '{for(i=1;i<=NF;i++){if($i=="ID"){print $(i+1)}}}' | \
          sed 's/://' | sort -u)

    for ID in $IDS; do
        ((id_count[$ID]++))
        echo "[DEBUG] ID=$ID count=${id_count[$ID]}"
        if [ "${id_count[$ID]}" -ge "$MAX_REPEAT" ]; then
            STAT_FILE="work/${ID}/stat_job"
            if [ -f "$STAT_FILE" ]; then
                echo "[INFO] ID $ID exceeded $MAX_REPEAT times. Marking as skip."
                sed -i '3s/.*/skip/' "$STAT_FILE"
                echo "$(date '+%Y-%m-%d %H:%M:%S') ID $ID skipped" >> "$SKIP_LOG"
            else
                echo "[WARN] stat_job file not found for ID $ID"
            fi
            unset id_count[$ID]
        fi
    done

    # 終了条件
    if [ "$LOG_LASTLINE" = "Done all structures!" ]; then
        echo "[INFO] All structures processed. Exiting."
        exit 0
    elif [[ "$LOG_LASTLINE" =~ ^Reached\ maxgen_ea ]]; then
        echo "[INFO] EA max generation reached. Exiting."
        exit 0
    elif [ "$LOG_LASTLINE" = "EA is ready" ]; then
        echo "[INFO] EA ready. Running EA steps..."
        cryspy -n || true
        cat log_cryspy
        LOG_LASTLINE=$(tail -n 1 log_cryspy)
        if [[ "$LOG_LASTLINE" =~ ^Reached\ maxgen_ea ]]; then
            exit 0
        fi
        cryspy -n || true
    elif [[ "$LOG_LASTLINE" =~ ^Reached\ max_select_bo ]]; then
        echo "[INFO] BO max selection reached. Exiting."
        exit 0
    elif [ "$LOG_LASTLINE" = "BO is ready" ]; then
        echo "[INFO] BO ready. Running BO steps..."
        cryspy -n || true
        cat log_cryspy
        LOG_LASTLINE=$(tail -n 1 log_cryspy)
        if [[ "$LOG_LASTLINE" =~ ^Reached\ max_select_bo ]]; then
            exit 0
        fi
        cryspy -n || true
    elif [ "$LOG_LASTLINE" = "LAQA is ready" ]; then
        echo "[INFO] LAQA ready. Running LAQA steps..."
        cryspy -n || true
        cryspy -n || true
    fi

    sleep 150
done
