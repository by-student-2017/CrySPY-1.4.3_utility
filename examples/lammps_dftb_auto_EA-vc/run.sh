#!/bin/bash
set -e

run_stage() {
  echo "Monitoring CrySPY progress..."
  tail -f log_cryspy &
  TAIL_PID=$!
  while :
  do
    cryspy -n || { echo "CrySPY failed"; kill $TAIL_PID; exit 1; }
    LOG_LASTLINE=$(tail -n 1 log_cryspy)
    if [ "$LOG_LASTLINE" = "Done all structures!" ]; then
      echo "Stage completed."
      break
    elif [ "${LOG_LASTLINE:0:17}" = "Reached maxgen_ea" ]; then
      echo "[INFO] EA max generation reached. Exiting."
      break
    elif [ "$LOG_LASTLINE" = "EA is ready" ]; then
      echo "EA stage detected, continuing..."
      cryspy -n
    fi
    sleep 60
  done
  kill $TAIL_PID
}

restart_if_needed() {
  if grep -q "Unhandled exception" log_cryspy; then
    echo "CrySPY error detected. Restarting..."
    cryspy_restart || { echo "CrySPY restart failed"; exit 1; }
  fi
}

# --- Stage 1: Lammps 計算 ---
echo "===== Stage 1: Lammps calculation START ====="
cp -r calc_in_lammps calc_in
cp -r Xx_tmp_lammps Xx_tmp
python3 make_input_lammps.py "$@"
cryspy -n
run_stage
restart_if_needed
rm -rf cryspy.in calc_in
echo "===== Stage 1: Lammps calculation END ====="

# --- Stage 2: DFTB+ 計算 ---
echo "===== Stage 2: DFTB+ calculation START ====="
rm -fr calc_in
rm -fr Xx_tmp
cp -r calc_in_dftb calc_in
cp -r Xx_tmp_dftb Xx_tmp
python3 make_input_dftb.py "$@"
cryspy -n
run_stage
restart_if_needed
rm -rf calc_in
echo "===== Stage 2: DFTB+ calculation END ====="

