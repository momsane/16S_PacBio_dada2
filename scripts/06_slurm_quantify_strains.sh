#!/bin/bash

#SBATCH --account pengel_general_data
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 4000
#SBATCH --partition cpu
#SBATCH --time 00:30:00
#SBATCH --error /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees/logs/06_quantify_strains.log
#SBATCH --output /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees/logs/06_quantify_strains.log

echo -e "$(date) job $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"

module purge # Make sure nothing is already loaded

# modify the path to your conda installation, or use the instructions from the curnagl wiki if using the cluster conda
CONDA_HOME=/work/FAC/FBM/DMF/pengel/general_data/mgarci14/miniforge3 # Path to Conda installation
source $CONDA_HOME/etc/profile.d/conda.sh # Source Conda initialization script
conda activate R # Activate Conda env

# Variables to modify
root=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees

## choose the correct option:
ps="$root"/results/assign_taxonomy/phyloseq_object_filtered_nonrarefied.RDS
# ps="$root"/results/assign_taxonomy/phyloseq_object_filtered_rarefied.RDS

clusters="$root"/workflow/config/all_16S_cd-hit_clusters_tax_full.tsv
qpcr="$root"/data/qPCR_results_analyzed.tsv # tab-delimited table containing qpcr data, or "" if none
abundance_col=normalized_16S_copies_gut # column in which the total abundance is reported, or "" if none
facet_var=SampleType
maxraref=-1 # set to ~4 times less the value you used in the denoising step, or -1 to skip rarefaction curves

# do not modify below this line
script="$root"/workflow/scripts/06_quantify_strains.R

out_quant="$root"/results/quantify_strains
out_plots="$root"/plots

# Execute the R script

echo "Parameters:"
echo input.ps: "$asvs"
echo input.clusters: "$clusters"
echo input.qpcr: "$qpcr"
echo abundance_col: "$abundance_col"
echo facet_var: "$facet_var"
echo maxraref: "$maxraref"
echo out.quant: "$out_quant"
echo out.plots: "$out_plots"

Rscript --vanilla "$script" \
    "$ps" \
    "$clusters" \
    "$qpcr" \
    "$abundance_col" \
    "$facet_var" \
    "$maxraref" \
    "$out_quant" \
    "$out_plots"

echo -e "$(date)"