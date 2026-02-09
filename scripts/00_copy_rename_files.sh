#!/bin/bash

# execute this script on the login node! (no access to nas from the compute nodes)

# modify the path to your conda installation, or use the instructions from the curnagl wiki if using the cluster conda
CONDA_HOME=/work/FAC/FBM/DMF/pengel/general_data/mgarci14/miniforge3 # Path to Conda installation
source $CONDA_HOME/etc/profile.d/conda.sh # Source Conda initialization script
conda activate bbmap # Activate Conda env containing dos2unix

# modify these variables
path_to_nas=/nas/FAC/FBM/DMF/pengel/general_data/D2c/datasets/NGS_data/20250318_Kinnex_Pesticom_JV_MG_AQ/raw_reads/PE_Kinnex16S_Apr24_A123_Circular_Consensus_Sequencing_Segmented_Reads
path_to_cluster=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/data/raw_reads
name_table=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/workflow/config/rename_files.tsv 
# tab-separated table with two columns only (no header): first column = original name, second column = new name
# make sure to add an empty line at the bottom, otherwise the last file won't be processed


# do not modify below this line
mkdir -p "$path_to_cluster"

dos2unix "$name_table"

# copy and rename files
while read -r oldname newname; do
    cp "$path_to_nas"/"$oldname" "$path_to_cluster"/"$newname"
done < "$name_table"

# check file names
echo "Please check that file names are correct and do not contain special characters:"
ls -la "$path_to_cluster"


