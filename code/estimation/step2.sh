#!/bin/bash

set -euo pipefail

# Submit the Step 2 estimation pipeline on Slurm:
# `step2_dist.sbatch`, then `step2.sbatch`, `step2_bs.sbatch`, and
# `step2_sa.sbatch` in parallel, then `step2_result.sbatch`.
# Usage: Run from `code/estimation` with `bash step2.sh` (or `./step2.sh`).

LOG_DIRECTORY=${REPLICATION_SLURM_LOG_DIRECTORY:-$PWD}
mkdir -p "$LOG_DIRECTORY"

J1=$(sbatch --parsable \
	--output="$LOG_DIRECTORY/step2_dist_%j.out" \
	--error="$LOG_DIRECTORY/step2_dist_%j.err" \
	step2_dist.sbatch)
J2=$(sbatch --parsable \
	--dependency=afterok:$J1 \
	--output="$LOG_DIRECTORY/step2_%j.out" \
	--error="$LOG_DIRECTORY/step2_%j.err" \
	step2.sbatch)
J3=$(sbatch --parsable \
	--dependency=afterok:$J1 \
	--output="$LOG_DIRECTORY/step2_bs_%A_%a.out" \
	--error="$LOG_DIRECTORY/step2_bs_%A_%a.err" \
	step2_bs.sbatch)
J4=$(sbatch --parsable \
	--dependency=afterok:$J1 \
	--output="$LOG_DIRECTORY/step2_sa_%A_%a.out" \
	--error="$LOG_DIRECTORY/step2_sa_%A_%a.err" \
	step2_sa.sbatch)
J5=$(sbatch --parsable \
	--dependency=afterok:$J2:$J3:$J4 \
	--output="$LOG_DIRECTORY/step2_result_%j.out" \
	--error="$LOG_DIRECTORY/step2_result_%j.err" \
	step2_result.sbatch)

printf 'Step 2 distribution job: %s\n' "$J1"
printf 'Step 2 main job: %s\n' "$J2"
printf 'Step 2 bootstrap array job: %s\n' "$J3"
printf 'Step 2 sensitivity array job: %s\n' "$J4"
printf 'Step 2 result job: %s\n' "$J5"
