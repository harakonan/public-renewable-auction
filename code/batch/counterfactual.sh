#!/usr/bin/env bash

# Run the local production counterfactual and create compact repository CSVs.

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
RUNTIME_FILE="$LOG_DIRECTORY/runtime.tsv"
mkdir -p "$LOG_DIRECTORY"

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
	printf 'counterfactual\tlocal\t%s\t%s\t%s\t%s\n' \
		"$START_ISO" "$end_iso" "$elapsed" "$STATUS" >>"$RUNTIME_FILE"
}
trap record_runtime EXIT

cd "$PROJECT_DIRECTORY"
CF_FRESH=${CF_FRESH:-1} \
CF_WORKERS=${CF_WORKERS:-4} \
	code/counterfactual/counterfactual.sh

STATUS=completed
printf 'Counterfactual stage completed. Runtime record: %s\n' "$RUNTIME_FILE"
