#!/bin/bash --login

#SBATCH --job-name=step1
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --constraint=amd24
#SBATCH --time=0-2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

# completed in ~ 10 min

set -euo pipefail

module load Julia/1.11.3-linux-x86_64

: "${SLURM_SUBMIT_DIR:?Submit this file from code/estimation.}"
cd "$SLURM_SUBMIT_DIR"

data_root=$(<../../data.txt)
data_file="$data_root/processed/generation_auction_wind20.csv"
header=$(head -n 1 "$data_file")
for column in capacity bid_quantity prior_started_wind_capacity; do
    if [[ ",$header," != *",$column,"* ]]; then
        echo "Required column '$column' is missing from $data_file" >&2
        echo "Synchronize the current processed dataset before submitting Step 1." >&2
        exit 1
    fi
done

srun julia --threads 1 step1_spot.jl
srun julia --threads 1 step1_prespecified.jl
srun julia --threads 1 step1_cdfdc.jl
srun julia --threads 1 step1.jl
