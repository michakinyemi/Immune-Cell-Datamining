#!/bin/bash
# #SBATCH --account=hnguyen
# #SBATCH --mail-user=mi542876@ucf.edu
# #SBATCH --mail-type=FAIL
# #SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=10:00:00
# #SBATCH --error=myjobresults-%J.err
#SBATCH --output=myjobresults-%J.out
#SBATCH --job-name=scRNA-Seq

echo "Slurm nodes assigned: $SLURM_JOB_NODELIST"

# Load modules
module load anaconda
module load samtools
module load bowtie2

conda activate compbio

bamCoverage --version
macs3 --version

bash src/atac/atac_pipeline.sh Sample4
