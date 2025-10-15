#!/bin/bash
#SBATCH --mail-user=mi542876@ucf.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=10:00:00
#SBATCH --output=myjobresults-%J.out
#SBATCH --job-name=scRNA-Seq

echo "Slurm nodes assigned: $SLURM_JOB_NODELIST"

# Load modules
module load anaconda
module load libpng
conda activate immune_r

Rscript samples/GSE298578_SLE.R


