#!/bin/bash

#SBATCH --account pengel_general_data
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 16000
#SBATCH --partition cpu
#SBATCH --time 05:00:00
#SBATCH --error /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/logs/04_denoising.log
#SBATCH --output /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/logs/04_denoising.log

echo -e "$(date) job $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"

module purge # Make sure nothing is already loaded

# modify the path to your conda installation, or use the instructions from the curnagl wiki if using the cluster conda
CONDA_HOME=/work/FAC/FBM/DMF/pengel/general_data/mgarci14/miniforge3 # Path to Conda installation
source $CONDA_HOME/etc/profile.d/conda.sh # Source Conda initialization script
conda activate R # Activate Conda env

# Variables: modify these paths to your own
root=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis
script="$root"/workflow/scripts/04_denoising.R
processed_fastq_dir="$root"/results/preprocessing/trimmed_filtered_reads
readcounts="$root"/results/preprocessing/read_count_before_after.tsv
maxReads=1000000
errModel=binnedQualErrfun # use `binnedQualErrfun` for Revio data, or `PacBioErrfun` for Sequel data
maxBases=10000000000 # 1E10 strongly recommended
removeSingletons=F # option for the dada function
maxraref=5000 # use the multiqc output to set this value close to the max number of reads in a sample
out_denois="$root"/results/denoising
out_plots="$root"/plots

# Execute the R script

echo "Parameters:"
echo input.reads: "$processed_fastq_dir"
echo input.readcounts: "$readcounts"
echo maxReads: "$maxReads"
echo errModel: "$errModel"
echo maxBases: "$maxBases"
echo removeSingletons: "$removeSingletons"
echo maxraref: "$maxraref"
echo out.denois: "$out_denois"
echo out.plots: "$out_plots"

echo "Refer to the Rscript for information on the parameters"

Rscript --vanilla "$script" \
    "$processed_fastq_dir" \
    "$readcounts" \
    "$maxReads" \
    "$errModel" \
    "$maxBases" \
    "$removeSingletons" \
    "$maxraref" \
    "$out_denois" \
    "$out_plots"

echo -e "$(date)"