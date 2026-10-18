#!/usr/bin/env bash

# Generate the paper-facing tables, figures, and sentence files.

set -euo pipefail

SCRIPT_DIRECTORY=$(cd "$(dirname "$0")" && pwd)
PROJECT_DIRECTORY=$(cd "$SCRIPT_DIRECTORY/../.." && pwd)
DATA_FILE="$PROJECT_DIRECTORY/data.txt"

[[ -f "$DATA_FILE" ]] || {
	printf 'Missing data-directory configuration: %s\n' "$DATA_FILE" >&2
	exit 1
}
DATA_DIRECTORY=$(<"$DATA_FILE")
DATA_DIRECTORY=${DATA_DIRECTORY%$'\r'}
[[ -n "$DATA_DIRECTORY" && "$DATA_DIRECTORY" != "/" ]] || {
	printf 'data.txt must contain a nonempty data-directory path.\n' >&2
	exit 2
}
LOG_DIRECTORY="$DATA_DIRECTORY/log/replication"
LOG_FILE="$LOG_DIRECTORY/result.log"
RUNTIME_FILE="$LOG_DIRECTORY/runtime.tsv"
mkdir -p "$LOG_DIRECTORY"

REQUIRED_INPUTS=(
	"$PROJECT_DIRECTORY/output/result/params1_result.csv"
	"$PROJECT_DIRECTORY/output/result/gamma_c_dist_result.csv"
	"$PROJECT_DIRECTORY/output/result/counterfactual_primary.csv"
	"$DATA_DIRECTORY/processed/generation_auction_wind20.csv"
)
for input_file in "${REQUIRED_INPUTS[@]}"; do
	[[ -f "$input_file" ]] || {
		printf 'Missing result-stage input: %s\n' "$input_file" >&2
		exit 1
	}
done

START_EPOCH=$(date +%s)
START_ISO=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
STATUS=failed
record_runtime() {
	end_epoch=$(date +%s)
	end_iso=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
	elapsed=$((end_epoch - START_EPOCH))
	if [[ ! -f "$RUNTIME_FILE" ]]; then
		printf 'stage\tenvironment\tstart_utc\tend_utc\telapsed_seconds\tstatus\n' \
			>"$RUNTIME_FILE"
	fi
	printf 'result\tlocal\t%s\t%s\t%s\t%s\n' \
		"$START_ISO" "$end_iso" "$elapsed" "$STATUS" >>"$RUNTIME_FILE"
}
trap record_runtime EXIT

: >"$LOG_FILE"
run() {
	printf '\n+ %s\n' "$*" | tee -a "$LOG_FILE"
	"$@" 2>&1 | tee -a "$LOG_FILE"
}

cd "$PROJECT_DIRECTORY"
run julia --threads 1 code/result/bne_check.jl
run Rscript code/result/descriptive.R
run Rscript code/result/estimation.R
run Rscript code/result/counterfactual.R

STATUS=completed
printf '\nPaper outputs completed. Log: %s\n' "$LOG_FILE"
