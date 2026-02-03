#!/bin/bash

#SBATCH --account pengel_general_data
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 4000
#SBATCH --partition cpu
#SBATCH --time 00:30:00
#SBATCH --error /work/FAC/FBM/DMF/pengel/general_data/20260129_syncom_invivo/logs/06_decontaminate.log
#SBATCH --output /work/FAC/FBM/DMF/pengel/general_data/20260129_syncom_invivo/logs/06_decontaminate.log

echo -e "$(date) job $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"

module purge # Make sure nothing is already loaded

# modify the path to your conda installation, or use the instructions from the curnagl wiki if using the cluster conda
CONDA_HOME=/work/FAC/FBM/DMF/pengel/general_data/mgarci14/miniforge3 # Path to Conda installation
source $CONDA_HOME/etc/profile.d/conda.sh # Source Conda initialization script
conda activate R # Activate Conda env

# Variables to modify
root=/work/FAC/FBM/DMF/pengel/general_data/20260129_syncom_invivo
group_var=SampleType
blank_name=blank

# do not modify below this line
script="$root"/workflow/scripts/06_decontaminate.R
ps="$root"/results/assign_taxonomy/phyloseq_object_filtered.RDS
out_decont="$root"/results/decontaminate
out_plots="$root"/plots

# Execute the R script

echo "Parameters:"
echo input.ps: "$asvs"
echo group_var: "$facet_var"
echo blank_name: "$blank_name"
echo out.decont: "$out_decont"
echo out.plots: "$out_plots"

echo "Refer to the Rscript for information on the parameters"

Rscript --vanilla "$script" \
    "$ps" \
    "$group_var" \
    "$blank_name" \
    "$out_decont" \
    "$out_plots"

echo -e "$(date)"