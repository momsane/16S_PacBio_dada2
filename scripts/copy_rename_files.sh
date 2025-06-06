#!/bin/bash

# execute this script on the login node! (no access to nas from the compute nodes)

path_to_nas=/nas/FAC/FBM/DMF/pengel/general_data/D2c/datasets/NGS_data/20250318_Kinnex_Pesticom_JV_MG_AQ/raw_reads/PE_Kinnex16S_Apr24_A123_Circular_Consensus_Sequencing_Segmented_Reads
path_to_cluster=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/data/raw_reads
name_table=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/workflow/config/rename_files.tsv 
# tab-separated table with two columns only (no header): first column = original name, second column = new name
# make sure to add an empty line at the bottom, otherwise the last file won't be renamed
# make sure the table is in Unix format! use dos2unix to convert it if necessary

# copy files
cp "$path_to_nas"/*.fastq.gz "$path_to_cluster"

# rename files
cd "$path_to_cluster"
while read -r oldname newname; do
    mv "$oldname" "$newname"
done < "$name_table"

# check file names
echo "Please check file names (first few files):"
ls "$path_to_cluster" | head


