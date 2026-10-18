#!/usr/bin/env bash

# Submit the complete Step 1 dependency chain on a Slurm cluster.

set -euo pipefail

SCRIPT_DIRECTORY=$(cd "$(dirname "$0")" && pwd)
PROJECT_DIRECTORY=$(cd "$SCRIPT_DIRECTORY/../.." && pwd)
DATA_FILE="$PROJECT_DIRECTORY/data.txt"

command -v sbatch >/dev/null 2>&1 || {
	printf 'estimation_step1 must be launched on a Slurm login node.\n' >&2
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
INPUT_FILE="$DATA_DIRECTORY/processed/generation_auction_wind20.csv"
[[ -f "$INPUT_FILE" ]] || {
	printf 'Missing Step 1 input: %s\n' "$INPUT_FILE" >&2
	exit 1
}

# A clean HPCC replication starts without generated output or intermediate
# directories. Create the two roots before submitting jobs that write to them.
mkdir -p "$DATA_DIRECTORY/intermediate" "$PROJECT_DIRECTORY/output/result"

LOG_DIRECTORY="$DATA_DIRECTORY/log/replication/estimation_step1"
SUBMISSION_LOG="$LOG_DIRECTORY/submission_$(date -u '+%Y%m%dT%H%M%SZ').log"
mkdir -p "$LOG_DIRECTORY"

cd "$PROJECT_DIRECTORY/code/estimation"
export REPLICATION_SLURM_LOG_DIRECTORY="$LOG_DIRECTORY"
./step1.sh | tee "$SUBMISSION_LOG"

FINAL_JOB_ID=$(awk '/Step 1 result job:/ {print $NF}' "$SUBMISSION_LOG")
[[ -n "$FINAL_JOB_ID" ]] || {
	printf 'Could not identify the final Step 1 job in %s\n' "$SUBMISSION_LOG" >&2
	exit 1
}
printf '\nStep 1 was submitted. Final dependent job: %s\n' "$FINAL_JOB_ID"
printf 'Wait for this job to complete successfully before submitting Step 2.\n'
printf 'Submission record: %s\n' "$SUBMISSION_LOG"
