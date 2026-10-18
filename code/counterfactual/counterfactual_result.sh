#!/bin/bash

set -euo pipefail

SCRIPT_DIRECTORY=$(cd "$(dirname "$0")" && pwd)
PROJECT_DIRECTORY=$(cd "$SCRIPT_DIRECTORY/../.." && pwd)
DATA_DIRECTORY=$(<"$PROJECT_DIRECTORY/data.txt")
DATA_DIRECTORY=${DATA_DIRECTORY%$'\r'}
LOG_DIRECTORY="$DATA_DIRECTORY/log/counterfactual"
mkdir -p "$LOG_DIRECTORY"
cd "$SCRIPT_DIRECTORY"

env \
	JULIA_NUM_THREADS=1 \
	OPENBLAS_NUM_THREADS=1 \
	CF_RHO_GRID=${CF_RHO_GRID:--0.658,0,0.13284679368741764,0.454} \
	CF_LOGGAMMA_SD_GRID=${CF_LOGGAMMA_SD_GRID:-0.0383370135,0.0884369676} \
	julia --threads 1 counterfactual_result.jl \
	2>&1 | tee "$LOG_DIRECTORY/aggregation.log"

printf 'Generating compact repository outputs.\n'
julia --threads 1 "$PROJECT_DIRECTORY/code/result/counterfactual.jl" \
	2>&1 | tee "$LOG_DIRECTORY/repository_output.log"
