#!/usr/bin/env bash

# Construct the processed datasets used by estimation and paper-output code.

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
[[ -d "$DATA_DIRECTORY/original" ]] || {
	printf 'Missing original-data directory: %s\n' "$DATA_DIRECTORY/original" >&2
	exit 1
}

PROCESSED_DIRECTORY="$DATA_DIRECTORY/processed"
LOG_DIRECTORY="$DATA_DIRECTORY/log/replication"
LOG_FILE="$LOG_DIRECTORY/data_processing.log"
RUNTIME_FILE="$LOG_DIRECTORY/runtime.tsv"
mkdir -p "$PROCESSED_DIRECTORY/step1" "$LOG_DIRECTORY"

# These hand-coded crosswalks are supplied replication inputs. The firm script
# validates them and constructs generation_auction_represent.csv.
MANUAL_INPUTS=(
	"$PROCESSED_DIRECTORY/generation_auction_firm.csv"
	"$PROCESSED_DIRECTORY/generation_auction_parent.csv"
	"$PROCESSED_DIRECTORY/generation_auction_represent_other.csv"
)
for input_file in "${MANUAL_INPUTS[@]}"; do
	[[ -f "$input_file" ]] || {
		printf 'Missing supplied firm crosswalk: %s\n' "$input_file" >&2
		exit 1
	}
done

# Check command-line dependencies before beginning a stage that overwrites
# processed outputs. In particular, spot_us_step1.sh requires Gnumeric to
# convert the legacy ICE .xls workbooks.
for required_command in Rscript python3 ssconvert; do
	command -v "$required_command" >/dev/null 2>&1 || {
		printf 'Missing required command: %s\n' "$required_command" >&2
		if [[ "$required_command" == "ssconvert" ]]; then
			printf 'On macOS, install it with: brew install gnumeric\n' >&2
		fi
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
	printf 'data_processing\tlocal\t%s\t%s\t%s\t%s\n' \
		"$START_ISO" "$end_iso" "$elapsed" "$STATUS" >>"$RUNTIME_FILE"
}
trap record_runtime EXIT

: >"$LOG_FILE"
run() {
	printf '\n+ %s\n' "$*" | tee -a "$LOG_FILE"
	"$@" 2>&1 | tee -a "$LOG_FILE"
}

cd "$PROJECT_DIRECTORY"
run Rscript code/data_processing/generation_auction.R
run Rscript code/data_processing/generation_auction_participation.R
run Rscript code/data_processing/corporate_structure.R
run Rscript code/data_processing/generation_data.R
run Rscript code/data_processing/generation_auction_firm.R
run Rscript code/data_processing/generation_auction_lee_technology.R
run Rscript code/data_processing/spot_average.R
run bash code/data_processing/spot_us_step1.sh "$DATA_DIRECTORY" "$PROJECT_DIRECTORY/code"
run Rscript code/data_processing/spot_us_step2.R
run Rscript code/data_processing/spot_us_average.R
run Rscript code/data_processing/epe_cost_step2.R
run Rscript code/data_processing/base_interest_rate.R
run Rscript code/data_processing/estimation_prep.R

STATUS=completed
printf '\nData processing completed. Log: %s\n' "$LOG_FILE"
