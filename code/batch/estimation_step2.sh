#!/usr/bin/env bash

# Submit the complete Step 2 dependency chain on a Slurm cluster.

set -euo pipefail

SCRIPT_DIRECTORY=$(cd "$(dirname "$0")" && pwd)
PROJECT_DIRECTORY=$(cd "$SCRIPT_DIRECTORY/../.." && pwd)
DATA_FILE="$PROJECT_DIRECTORY/data.txt"

command -v sbatch >/dev/null 2>&1 || {
	printf 'estimation_step2 must be launched on a Slurm login node.\n' >&2
	exit 1
}
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
REQUIRED_INPUTS=(
	"$PROJECT_DIRECTORY/output/result/params1_hat.csv"
	"$PROJECT_DIRECTORY/output/result/params1_result.csv"
	"$PROJECT_DIRECTORY/output/result/params_specified.csv"
	"$PROJECT_DIRECTORY/output/result/spot_result.csv"
)
for input_file in "${REQUIRED_INPUTS[@]}"; do
	[[ -f "$input_file" ]] || {
		printf 'Missing completed Step 1 output: %s\n' "$input_file" >&2
		exit 1
	}
done

LOG_DIRECTORY="$DATA_DIRECTORY/log/replication/estimation_step2"
SUBMISSION_LOG="$LOG_DIRECTORY/submission_$(date -u '+%Y%m%dT%H%M%SZ').log"
mkdir -p "$LOG_DIRECTORY"

cd "$PROJECT_DIRECTORY/code/estimation"
export REPLICATION_SLURM_LOG_DIRECTORY="$LOG_DIRECTORY"
./step2.sh | tee "$SUBMISSION_LOG"

FINAL_JOB_ID=$(awk '/Step 2 result job:/ {print $NF}' "$SUBMISSION_LOG")
[[ -n "$FINAL_JOB_ID" ]] || {
	printf 'Could not identify the final Step 2 job in %s\n' "$SUBMISSION_LOG" >&2
	exit 1
}
printf '\nStep 2 was submitted. Final dependent job: %s\n' "$FINAL_JOB_ID"
printf 'Wait for this job to complete successfully before synchronizing outputs locally.\n'
printf 'Submission record: %s\n' "$SUBMISSION_LOG"
