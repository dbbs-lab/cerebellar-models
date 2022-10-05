#!/bin/bash -l
#SBATCH --job-name="cereb recon"
#SBATCH --account="ich027"
#SBATCH --time=24:00:00
#SBATCH --nodes=20
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --hint=nomultithread

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
source $HOME/bsbenv/bin/activate

srun bsb compile $HOME/cerebellum/cerebellum.json
