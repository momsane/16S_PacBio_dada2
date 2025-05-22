#!/bin/bash

#SBATCH --account pengel_general_data
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 4000
#SBATCH --partition cpu
#SBATCH --time 00:10:00
#SBATCH --array=2-51
#SBATCH --error /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/logs/00_rarefy/%a.log
#SBATCH --output /work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/logs/00_rarefy/%a.log

echo -e "$(date) job $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"

module purge # Make sure nothing is already loaded

# modify the path to your conda installation, or use the instructions from the curnagl wiki if using the cluster conda
CONDA_HOME=/work/FAC/FBM/DMF/pengel/general_data/mgarci14/miniforge3 # Path to Conda installation
source $CONDA_HOME/etc/profile.d/conda.sh # Source Conda initialization script
conda activate bbmap-39.25 # Activate Conda env

# Variables: modify these paths to your own
root=/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis
raw_reads="$root"/data/raw_reads
rarefied_reads="$root"/results/prerarefied_reads
table="$root"/workflow/config/pre_rarefaction.tsv # tab-separated table of 2 columns: SampleId and rarefy_to, with header

# do not modify below this line
sample=$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} 'NR==ArrayTaskID {print $1}' "$table")
rarefy_to=$(awk -v ArrayTaskID=${SLURM_ARRAY_TASK_ID} 'NR==ArrayTaskID {print $2+0}' "$table")
file="$raw_reads"/"$sample".fastq.gz

mkdir -p "$rarefied_reads"
echo "$sample"

count=$(( $(zcat "$file" | wc -l) / 4 )) # count reads to see if there are more or less reads than rarefy_to
count=$(($count + 0)) # make sure it is numeric

if [[ "$count" -gt "$rarefy_to" ]]
    then
        echo "Rarefying to $rarefy_to reads"
        reformat.sh in="$file" out="$rarefied_reads"/"$sample".fastq.gz ow=t reads=$rarefy_to sampleseed=42
    else
        echo "Not rarefying since sample has only $count reads"
        cp "$file" "$rarefied_reads"/"$sample".fastq.gz
fi

echo -e "$(date)"