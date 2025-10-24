#!/bin/bash
set -e

run_stage() {
  echo "Monitoring CrySPY progress..."
  tail -f log_cryspy &
  TAIL_PID=$!
  trap "kill $TAIL_PID" EXIT
  while :
  do
    cryspy -n || { echo "CrySPY failed"; exit 1; }
    if grep -q "Done all structures!" log_cryspy; then
      echo "Stage completed."
      break
    elif grep -q "Reached maxgen_ea" log_cryspy; then
      echo "[INFO] EA max generation reached. Exiting."
      break
    fi
    sleep 60
  done
}

restart_if_needed() {
  if grep -q "Unhandled exception" log_cryspy; then
    echo "CrySPY error detected. Restarting..."
    cryspy_restart || { echo "CrySPY restart failed"; exit 1; }
  fi
}

# --- Stage 1: CHGNet計算 ---
echo "===== Stage 1: CHGNet calculation START ====="
cp -r calc_in_chgnet calc_in
cp -r Xx_tmp_chgnet Xx_tmp
python3 make_input_chgnet.py "$@"
cryspy -n
run_stage
restart_if_needed
rm -rf cryspy.in calc_in
echo "===== Stage 1: CHGNet calculation END ====="

echo "===== Clean up previous files and reset environment ====="
find ./data -type f ! -name 'init_struc_data.pkl' ! -name 'opt_struc_data.pkl' -delete
find ./data -type d -empty -delete
rm cryspy.stat err_cryspy log_cryspy
rm -fr work

# --- Stage 2: QE計算 ---
echo "===== Stage 2: QE calculation START ====="
rm -fr calc_in
cp -r calc_in_qe calc_in
rm -fr Xx_tmp
cp -r Xx_tmp_qe Xx_tmp
python3 make_input_qe.py "$@"
cryspy -n
run_stage
restart_if_needed
rm -rf calc_in
echo "===== Stage 2: QE calculation END ====="

