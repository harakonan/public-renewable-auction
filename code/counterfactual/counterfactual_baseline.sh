#!/bin/bash

set -euo pipefail

SCRIPT_DIRECTORY=$(cd "$(dirname "$0")" && pwd)
PROJECT_DIRECTORY=$(cd "$SCRIPT_DIRECTORY/../.." && pwd)
DATA_DIRECTORY=$(<"$PROJECT_DIRECTORY/data.txt")
DATA_DIRECTORY=${DATA_DIRECTORY%$'\r'}
TASK_DIRECTORY="$DATA_DIRECTORY/intermediate/counterfactual_task"
LOG_DIRECTORY="$DATA_DIRECTORY/log/counterfactual"
BASELINE_FILE="$TASK_DIRECTORY/counterfactual_baseline.jls"
ANCHOR_FILE="$TASK_DIRECTORY/counterfactual_marginal_participant_anchor.csv"

mkdir -p "$TASK_DIRECTORY" "$LOG_DIRECTORY"
if [[ -f "$BASELINE_FILE" && -f "$ANCHOR_FILE" ]]; then
	printf 'Reusing baseline common-support calibration: %s\n' "$BASELINE_FILE"
	exit 0
fi
if [[ -f "$BASELINE_FILE" ]]; then
	printf 'Replacing a baseline cache without marginal-participant anchor diagnostics.\n'
	rm -f "$BASELINE_FILE"
fi

cd "$SCRIPT_DIRECTORY"
env \
	JULIA_NUM_THREADS=1 \
	OPENBLAS_NUM_THREADS=1 \
	CF_MODE=fine \
	CF_N_AUX=${CF_N_AUX:-5000} \
	CF_N_WINPROB=${CF_N_WINPROB:-1000} \
	CF_N_ALLOC=${CF_N_ALLOC:-2000} \
	CF_N_OBJECTIVE_DRAWS=${CF_N_OBJECTIVE_DRAWS:-2000} \
	CF_MARGINAL_BID_DRAWS=${CF_MARGINAL_BID_DRAWS:-5000} \
	CF_REGISTERED_MAXEVAL=${CF_REGISTERED_MAXEVAL:-2000} \
	CF_BASELINE_FILE="$BASELINE_FILE" \
	julia --threads 1 counterfactual_baseline.jl \
	2>&1 | tee "$LOG_DIRECTORY/baseline.log"
