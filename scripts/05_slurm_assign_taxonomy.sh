#!/bin/bash

#SBATCH --account pengel_general_data
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 16000
#SBATCH --partition cpu
#SBATCH --time 02:00:00
#SBATCH --error /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees/logs/05_assign_taxonomy.log
#SBATCH --output /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees/logs/05_assign_taxonomy.log

echo -e "$(date) job $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"

module purge # Make sure nothing is already loaded

# modify the path to your conda installation, or use the instructions from the curnagl wiki if using the cluster conda
CONDA_HOME=/work/FAC/FBM/DMF/pengel/general_data/mgarci14/miniforge3 # Path to Conda installation
source $CONDA_HOME/etc/profile.d/conda.sh # Source Conda initialization script
conda activate R # Activate Conda env

# Variables to modify
root=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees
db1="$root"/data/databases/syncom_custom_db_toSpecies_trainset.fa
db2="$root"/data/databases/syncom_custom_db_addSpecies.fa
clusters="$root"/workflow/config/all_16S_cd-hit_clusters_tax_full.tsv
rarefy_to=-1 # -1 means no rarefaction; use the rarefaction curves to set this value to a relevant number if needed
facet_var=SampleType

# do not modify below this line
script="$root"/workflow/scripts/05_assign_taxonomy.R
asvs="$root"/results/denoising/ASV_samples_table_noChim.rds
metadata="$root"/workflow/config/metadata.tsv
out_tax="$root"/results/assign_taxonomy
out_plots="$root"/plots

# Execute the R script

echo "Parameters:"
echo input.asvs: "$asvs"
echo input.metadata: "$metadata"
echo db1: "$db1"
echo db2: "$db2"
echo clusters_tax: "$clusters"
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
    "$clusters" \
    "$rarefy_to" \
    "$facet_var" \
    "$out_tax" \
    "$out_plots"

echo -e "$(date)"