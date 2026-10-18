#!/bin/bash --login

#SBATCH --job-name=s1result
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --constraint=amd24
#SBATCH --time=0-0:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

# completed in < 1 min

module load Julia/1.11.3-linux-x86_64

: "${SLURM_SUBMIT_DIR:?Submit this file from code/estimation.}"
cd "$SLURM_SUBMIT_DIR"
srun julia --threads 1 step1_result.jl
