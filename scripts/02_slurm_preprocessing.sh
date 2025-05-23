#!/bin/bash

#SBATCH --account pengel_general_data
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 52000
#SBATCH --partition cpu
#SBATCH --time 02:00:00
#SBATCH --error /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/logs/02_preprocessing.log
#SBATCH --output /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/logs/02_preprocessing.log

echo -e "$(date) job $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"

module purge # Make sure nothing is already loaded

# modify the path to your conda installation, or use the instructions from the curnagl wiki if using the cluster conda
CONDA_HOME=/work/FAC/FBM/DMF/pengel/general_data/mgarci14/miniforge3 # Path to Conda installation
source $CONDA_HOME/etc/profile.d/conda.sh # Source Conda initialization script
conda activate R # Activate Conda env

# Variables: modify these paths to your own
root=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis
script="$root"/workflow/scripts/02_preprocessing.R

#reads="$root"/data/raw_reads # if you did not pre-rarefy
reads="$root"/results/prerarefied_reads # if you pre-rarefied

fwd_primer=AGRGTTYGATYMTGGCTCAG
rev_primer=RGYTACCTTGTTACGACTT
minLen=1000
maxLen=1600
out_preproc="$root"/results/preprocessing
out_plots="$root"/plots

# Execute the R script

echo "Parameters:"
echo input.raw: "$reads"
echo fwd.primer: "$fwd_primer"
echo rev.primer: "$rev_primer"
echo minLen: "$minLen"
echo maxLen: "$maxLen"
echo out.preproc: "$out_preproc"
echo out.plots: "$out_plots"

echo "Refer to the Rscript for information on the parameters"

Rscript --vanilla "$script" \
    "$reads" \
    "$fwd_primer" \
    "$rev_primer" \
    "$minLen" \
    "$maxLe" \
    "$out_preproc" \
    "$out_plots"

echo -e "$(date)"