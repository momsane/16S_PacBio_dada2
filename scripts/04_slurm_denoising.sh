#!/bin/bash

#SBATCH --account pengel_general_data
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 32000
#SBATCH --partition cpu
#SBATCH --time 05:00:00
#SBATCH --error /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees/logs/04_denoising.log
#SBATCH --output /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees/logs/04_denoising.log

echo -e "$(date) job $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"

module purge # Make sure nothing is already loaded

# modify the path to your conda installation, or use the instructions from the curnagl wiki if using the cluster conda
CONDA_HOME=/work/FAC/FBM/DMF/pengel/general_data/mgarci14/miniforge3 # Path to Conda installation
source $CONDA_HOME/etc/profile.d/conda.sh # Source Conda initialization script
conda activate R # Activate Conda env

# Variables to modify
root=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees
errModel=binnedQualErrfun # use `binnedQualErrfun` if you have binned quality score, or else `PacBioErrfun`
maxReads=1000000 # change only if memory issues arise
maxBases=10000000000 # 1E10 strongly recommended
db2="$root"/data/databases/syncom_custom_db_addSpecies.fa # give dada a set of expected ASVs, or set to ""
pool=F # if T, increases sensitivity in detecting rare taxa -> risk of keeping errors
maxraref=8000 # use the multiqc output to set this value close to the max number of reads in a sample

# do not modify below this line
script="$root"/workflow/scripts/04_denoising.R
out_denois="$root"/results/denoising
out_plots="$root"/plots
processed_fastq_dir="$root"/results/preprocessing/trimmed_filtered_reads
readcounts="$root"/results/preprocessing/read_count_before_after.tsv

# Execute the R script

echo "Parameters:"
echo input.reads: "$processed_fastq_dir"
echo input.readcounts: "$readcounts"
echo maxReads: "$maxReads"
echo errModel: "$errModel"
echo maxBases: "$maxBases"
echo db2: "$db2"
echo pool: "$pool"
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
    "$db2" \
    "$pool" \
    "$maxraref" \
    "$out_denois" \
    "$out_plots"

echo -e "$(date)"