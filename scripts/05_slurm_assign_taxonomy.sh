#!/bin/bash

#SBATCH --account pengel_general_data
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32000
#SBATCH --partition cpu
#SBATCH --time 02:00:00
#SBATCH --error /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/logs/05_assign_taxonomy.log
#SBATCH --output /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/logs/05_assign_taxonomy.log

echo -e "$(date) job $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"

module purge # Make sure nothing is already loaded

# modify the path to your conda installation, or use the instructions from the curnagl wiki if using the cluster conda
CONDA_HOME=/work/FAC/FBM/DMF/pengel/general_data/mgarci14/miniforge3 # Path to Conda installation
source $CONDA_HOME/etc/profile.d/conda.sh # Source Conda initialization script
conda activate R # Activate Conda env

# Variables: modify these paths to your own
root=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis
script="$root"/workflow/scripts/05_assign_taxonomy.R
asvs="$root"/results/denoising/ASV_samples_table_noChim.rds
metadata="$root"/workflow/config/metadata.tsv
db1="$root"/data/databases/gg2_2024_09_toSpecies_trainset.fa
db2=""
rarefy_to=-1
facet_var=SampleType
out_tax="$root"/results/assign_taxonomy
out_plots="$root"/plots

# Execute the R script

echo "Parameters:"
echo input.asvs: "$asvs"
echo input.metadata: "$metadata"
echo db1: "$db1"
echo db2: "$db2"
echo rarefy_to: "$rarefy_to"
echo facet_var: "$facet_var"
echo out.tax: "$out_tax"
echo out.plots: "$out_plots"

echo "Refer to the Rscript for information on the parameters"

Rscript --vanilla "$script" \
    "$asvs" \
    "$metadata" \
    "$db1" \
    "$db2" \
    "$rarefy_to" \
    "$facet_var" \
    "$out_tax" \
    "$out_plots"

echo -e "$(date)"