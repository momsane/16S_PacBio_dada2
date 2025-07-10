#!/bin/bash

#SBATCH --account pengel_general_data
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 1000
#SBATCH --partition cpu
#SBATCH --time 00:20:00
#SBATCH --array=2-232
#SBATCH --error /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees/logs/03_fastqc_postproc/%a.log
#SBATCH --output /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees/logs/03_fastqc_postproc/%a.log

echo -e "$(date) job $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"

module purge # Make sure nothing is already loaded

# modify the path to your conda installation, or use the instructions from the curnagl wiki if using the cluster conda
CONDA_HOME=/work/FAC/FBM/DMF/pengel/general_data/mgarci14/miniforge3 # Path to Conda installation
source $CONDA_HOME/etc/profile.d/conda.sh # Source Conda initialization script
conda activate qc # Activate Conda env

# Variable to modify
root=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees

# do not modify below this line
reads="$root"/results/preprocessing/trimmed_filtered_reads
metadata="$root"/workflow/config/metadata.tsv

sample=$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} 'NR==ArrayTaskID {print $1}' "$metadata")
out="$root"/results/fastqc_postproc/"$sample"

# Execute fastqc

mkdir -p "$out"
fastqc -o "$out" "$reads"/"$sample".fastq.gz