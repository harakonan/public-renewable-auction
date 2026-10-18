#!/bin/bash --login

#SBATCH --job-name=s1_bs
#SBATCH --output=step1_bs/%x_%a.out
#SBATCH --error=step1_bs/%x_%a.err
#SBATCH --array=1-200
#SBATCH --constraint=amd24
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

module load Julia/1.11.3-linux-x86_64

: "${SLURM_SUBMIT_DIR:?Submit this file from code/estimation.}"
cd "$SLURM_SUBMIT_DIR"
srun julia --threads 1 step1.jl
