#!/bin/bash

# Submit the complete Step 1 pipeline on Slurm. After the original estimation,
# the bootstrap and sensitivity arrays run in parallel. The result job waits for
# both arrays to finish successfully.
# Usage: Run from `code/estimation` with `bash step1.sh` (or `./step1.sh`).

set -euo pipefail

LOG_DIRECTORY=${REPLICATION_SLURM_LOG_DIRECTORY:-$PWD}
mkdir -p "$LOG_DIRECTORY"

J1=$(sbatch --parsable \
	--output="$LOG_DIRECTORY/step1_%j.out" \
	--error="$LOG_DIRECTORY/step1_%j.err" \
	step1.sbatch)
J2=$(sbatch --parsable \
	--dependency=afterok:$J1 \
	--output="$LOG_DIRECTORY/step1_bs_%A_%a.out" \
	--error="$LOG_DIRECTORY/step1_bs_%A_%a.err" \
	step1_bs.sbatch)
J3=$(sbatch --parsable \
	--dependency=afterok:$J1 \
	--output="$LOG_DIRECTORY/step1_sa_%A_%a.out" \
	--error="$LOG_DIRECTORY/step1_sa_%A_%a.err" \
	step1_sa.sbatch)
J4=$(sbatch --parsable \
	--dependency=afterok:$J2:$J3 \
	--output="$LOG_DIRECTORY/step1_result_%j.out" \
	--error="$LOG_DIRECTORY/step1_result_%j.err" \
	step1_result.sbatch)

printf 'Original estimation job: %s\n' "$J1"
printf 'Bootstrap array job: %s\n' "$J2"
printf 'Sensitivity array job: %s\n' "$J3"
printf 'Step 1 result job: %s\n' "$J4"
