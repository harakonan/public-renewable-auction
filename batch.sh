#!/usr/bin/env bash

# Stage dispatcher for the replication package. The complete replication spans
# a local machine and a Slurm cluster, so stages must be launched separately.

set -euo pipefail

SCRIPT_DIRECTORY=$(cd "$(dirname "$0")" && pwd)
BATCH_DIRECTORY="$SCRIPT_DIRECTORY/code/batch"

usage() {
	cat <<'EOF'
Usage: ./batch.sh STAGE

Stages, in execution order:
  data_processing   Local: construct the processed analysis datasets
  estimation_step1  HPCC: submit Step 1, bootstrap, and sensitivity jobs
  estimation_step2  HPCC: submit Step 2, bootstrap, and sensitivity jobs
  counterfactual    Local: solve and aggregate the counterfactual equilibria
  result            Local: generate paper tables, figures, and sentence files

The two estimation stages submit asynchronous Slurm dependency chains. Wait for
the final result job from estimation_step1 to finish before submitting
estimation_step2. After estimation_step2 finishes, synchronize its outputs back
to the local machine before running counterfactual and result.
EOF
}

[[ $# -eq 1 ]] || {
	usage
	exit 2
}

case "$1" in
	data_processing|estimation_step1|estimation_step2|counterfactual|result)
		exec "$BATCH_DIRECTORY/$1.sh"
		;;
	-h|--help|help)
		usage
		;;
	*)
		printf 'Unknown replication stage: %s\n\n' "$1" >&2
		usage >&2
		exit 2
		;;
esac
