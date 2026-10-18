#!/bin/bash

# Local production entry point for the counterfactual analysis.

set -euo pipefail

SCRIPT_DIRECTORY=$(cd "$(dirname "$0")" && pwd)
PROJECT_DIRECTORY=$(cd "$SCRIPT_DIRECTORY/../.." && pwd)
DATA_DIRECTORY=$(<"$PROJECT_DIRECTORY/data.txt")
DATA_DIRECTORY=${DATA_DIRECTORY%$'\r'}
TASK_DIRECTORY="$DATA_DIRECTORY/intermediate/counterfactual_task"
OUTPUT_DIRECTORY="$DATA_DIRECTORY/intermediate/counterfactual_output"
LOG_DIRECTORY="$DATA_DIRECTORY/log/counterfactual"
FRESH=${CF_FRESH:-0}
RHO_GRID=${CF_RHO_GRID:--0.658,0,0.13284679368741764,0.454}
LOGGAMMA_SD_GRID=${CF_LOGGAMMA_SD_GRID:-0.0383370135,0.0884369676}
COST_SHIFT_GRID=${CF_COST_SHIFT_GRID:-0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10}

[[ "$FRESH" == "0" || "$FRESH" == "1" ]] || {
	printf 'CF_FRESH must be 0 or 1.\n' >&2
	exit 2
}
[[ -n "$DATA_DIRECTORY" && "$DATA_DIRECTORY" != "/" ]] || {
	printf 'Refusing an empty or root data directory.\n' >&2
	exit 2
}

mkdir -p "$TASK_DIRECTORY" "$OUTPUT_DIRECTORY" "$LOG_DIRECTORY"
REQUIRED_INPUTS=(
	"$TASK_DIRECTORY/participant_primitives.csv"
	"$DATA_DIRECTORY/processed/generation_auction_wind20.csv"
	"$PROJECT_DIRECTORY/output/result/gamma_c_dist_result.csv"
	"$PROJECT_DIRECTORY/output/result/params1_hat.csv"
	"$PROJECT_DIRECTORY/output/result/params2_dist_result.csv"
	"$PROJECT_DIRECTORY/output/result/params_specified.csv"
	"$PROJECT_DIRECTORY/output/result/spot_result.csv"
)
for input_file in "${REQUIRED_INPUTS[@]}"; do
	[[ -f "$input_file" ]] || {
		printf 'Required counterfactual input is missing: %s\n' "$input_file" >&2
		exit 1
	}
done

if [[ "$FRESH" == "1" ]]; then
	printf 'Removing the canonical common-support cache and task outputs.\n'
	rm -f "$TASK_DIRECTORY/counterfactual_baseline.jls"
	rm -f "$TASK_DIRECTORY/counterfactual_marginal_participant_anchor.csv"
	find "$TASK_DIRECTORY" -maxdepth 1 -type d \
		-name 'counterfactual_task_[0-9][0-9][0-9]' -exec rm -rf {} +
	find "$OUTPUT_DIRECTORY" -maxdepth 1 -type f \
		-name 'counterfactual_*.csv' -delete
	find "$LOG_DIRECTORY" -maxdepth 1 -type f \
		-name '*.log' -delete
fi

cat >"$LOG_DIRECTORY/run.env" <<EOF
CF_MODEL=two_type_common_upper_support
CF_N_AUX=${CF_N_AUX:-5000}
CF_N_WINPROB=${CF_N_WINPROB:-1000}
CF_N_ALLOC=${CF_N_ALLOC:-2000}
CF_N_OBJECTIVE_DRAWS=${CF_N_OBJECTIVE_DRAWS:-2000}
CF_MARGINAL_BID_DRAWS=${CF_MARGINAL_BID_DRAWS:-5000}
CF_COST_SHIFT_GRID=$COST_SHIFT_GRID
CF_LAMBDA_GRID=${CF_LAMBDA_GRID:-0,0.05,...,1.00}
CF_RHO_GRID=$RHO_GRID
CF_LOGGAMMA_SD_GRID=$LOGGAMMA_SD_GRID
CF_ODE_MAXIMUM_BID_STEP=${CF_ODE_MAXIMUM_BID_STEP:-0.10}
CF_ODE_MAXIMUM_PROBABILITY_STEP=${CF_ODE_MAXIMUM_PROBABILITY_STEP:-0.002}
EOF

cd "$SCRIPT_DIRECTORY"
IFS=',' read -r -a RHO_VALUES <<< "$RHO_GRID"
IFS=',' read -r -a LOGGAMMA_SD_VALUES <<< "$LOGGAMMA_SD_GRID"
TASK_COUNT=$((${#RHO_VALUES[@]} * ${#LOGGAMMA_SD_VALUES[@]}))
WORKERS=${CF_WORKERS:-4}
./counterfactual_baseline.sh
printf 'Launching %s primitive scenarios with at most %s workers.\n' \
	"$TASK_COUNT" "$WORKERS"
seq 1 "$TASK_COUNT" | xargs -P "$WORKERS" -I {} \
	env CF_TASK_ID={} ./counterfactual_task.sh
./counterfactual_result.sh
