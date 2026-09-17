#!/bin/bash

#SBATCH --account pengel_general_data
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 6
#SBATCH --mem 24000
#SBATCH --partition cpu
#SBATCH --time 02:00:00
#SBATCH --error /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees/logs/06_quantify_strains.log
#SBATCH --output /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees/logs/06_quantify_strains.log

echo -e "$(date) job $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"

module purge # Make sure nothing is already loaded

# modify the path to your conda installation, or use the instructions from the curnagl wiki if using the cluster conda
CONDA_HOME=/work/FAC/FBM/DMF/pengel/general_data/mgarci14/miniforge3 # Path to Conda installation
source $CONDA_HOME/etc/profile.d/conda.sh # Source Conda initialization script
conda activate R # Activate Conda env

#  modify the path to your project folder
root=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees

# do not modify below this line

script="$root"/workflow/scripts/06_quantify_strains.R
config="$root"/workflow/config/config.R

## execute the R script
Rscript --vanilla "$script" \
    "$root" \
    "$config"

echo -e "$(date)"