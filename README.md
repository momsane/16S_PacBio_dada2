# Analysis of 16S PacBio data with dada2



- Table of contents
    - [Running the dada2 pipeline](#running-the-dada2-pipeline)
        - [Overview](#adapting-the-scripts)
        - [Requirements](#adapting-the-scripts)
        - [Setting up the work environment](#setting-up-the-work-environment)
        - [Data preparation](#data-preparation)
        - [Adapting the scripts](#adapting-the-scripts)
        - [Running the pipeline](#running-the-pipeline)
    - [Custom database for defined communities](#creating-a-custom-database-for-defined-communities)
        - [Merging and dereplicating 16S sequences](#merging-and-dereplicating-16s-sequences)
        - [Formatting databases for dada2](#formatting-databases-for-dada2)
        - [Comparing and merging databases](#comparing-and-merging-databases)

---

## Running the dada2 pipeline

### Overview

The pipeline consists of a few main steps:

- *(optional step)* pre-rarefaction of the reads for large samples `00_rarefy.sh`
- quality check on the raw reads `01_fastqc_preproc.sh` and `01_multiqc_preproc.sh`
- pre-processing of the reads with dada2 (primer removal, length and quality trimming) `02_preprocessing.sh`
- quality check on the processed reads `03_fastqc_postproc.sh` and `03_multiqc_posteproc.sh`
- denoising into ASVs with dada2 `04_denoising.sh`
- taxonmy assignment with dada2 `05_assign_taxonomy.sh`
- *(only for defined communities)* compute strain abundance from ASV table `06_quantify_strains.sh`.
dada2 is implemented in Rscripts called by the bash scripts. **There should not be any need to modify the R scripts**.

### Requirements

- This pipeline is designed to be executed on a slurm cluster.
- Conda/mamba.
- R, fastQC, multiQC and bbmap (see envs).

### Setting up the work environment

#### Structure

First, you will create a working directory with the following structure:

```
.
├── data
    └── raw_reads
├── logs
├── plots
├── results
└── workflow <- this folder is downloaded with git clone https://github.com/momsane/16S_PacBio_dada2
    ├── config
    ├── envs
    └── scripts
```

#### Conda environments

You'll need access to a conda installation on your cluster. You can either have yours, installed for instance through miniforge3, or use the one provided by the cluster. This will change slightly how you activate environments at the beginning of the scripts. To use the cluster conda, follow the instructions [here](https://wiki.unil.ch/ci/books/high-performance-computing-hpc/page/using-conda-and-anaconda). To use your own conda installation, you just need to modify the `$CONDA_HOME` variable in the scripts.

Install the required conda environments using the .yaml files located in envs with the command `conda env create -f envs/<env>.yaml`.

### Data Preparation

Before running the pipeline, you need to prepare your data:

1.  **Metadata file:** modify `config/metadata.tsv` according to your samples. You do not need to keep the same columns except for the first one. This first column must contain the sample names (filenames without the `.fastq.gz` extension).
2.  **File naming Table:** the raw read files you got from the sequencing facility have long non-informative names. You will rename them with the SampleID. Create a table like `config/rename_files.tsv` where the first column is the current name of each file, and the second column is the new name. This table has no header. It is important that this table be in Unix format. If you modify it in excel for instance, it will not be in Unix format. You can use the command-line tool dos2unix to convert it. Make sure there is a line return after the last line in the table.
3.  **Read rarefaction table (optional):** if you have very uneven depth in your dataset with very large samples, you might want to consider rarefying the raw reads to limit unnecessary computation time and resources. Modify `config/pre_rarefaction.tsv` according to your samples. As with `config/rename_files.tsv`, make sure it is in Unix format. If you do not want to rarefy your reads, skip this step.
4.  **Raw reads:** you are now ready to copy them from the nas. Modify the script `copy_rename_files.sh` with the correct paths. Then execute it from the login node (i.e. use `bash` instead of `sbatch` to submit it), because the NAS can only be accessed from the login node.

### Adapting the scripts

Only the `.sh` scripts need to be modified. You will need to modify the beginning of these scripts:

- modify the commands to initialize conda according to the type of installation you are using
- modify the input variables, for instance the path to the root directory
- the fastQC scripts are array jobs (argument `--array` in the slurm header), so you need to modify the range of the arrays. `2-50` means you will process files described in lines 2 to 50 of `config/metadata.tsv`. We start at 2 to skip the header. So your array range should be `2-<number of samples>+1`
- you should not need to modify the resource requirements, unless your jobs get killed.

On top of this, each script requires some specific customization to match your specific data:

- `01_fastqc_preproc.sh` and `02_slurm_preprocessing.sh`: depending on whether you rarefied your raw reads or not, modify the `reads` variable
- `02_preprocessing.sh` for example:

| Variable           | Description                                       | Default Value                                          |
|--------------------|---------------------------------------------------|--------------------------------------------------------|
| `fwd_primer`        | Forward primer sequence                            | `AGRGTTYGATYMTGGCTCAG`                                  |
| `rev_primer`        | Reverse primer sequence                            | `RGYTACCTTGTTACGACTT`                                  |
| `minLen`            | Minimum read length for filtering                  | `1000`                                                  |
| `maxLen`            | Maximum read length for filtering                  | `1600`                                                  |

### Running the pipeline

To run each script:

1.  **Submit the job to the slurm scheduler:** use `sbatch <script_name>.sh`.
2.  **Monitor the job:** use `Squeue` to check the status of you jobs. After each step, check the logs for errors. 
    
| Script           | What to check                                       |
|--------------------|---------------------------------------------------|
| `00_rarefy.sh`        | Size of output fastQ files is > 0                            |
| `01_fastqc_preproc.sh` & `03_fastqc_postproc.sh`      | Combine all logs into a single file and check that "Analysis complete" appears as many times as samples: `cat ../../logs/01_fastqc_preproc/*.log > ../../logs/01_fastqc_preproc/combined_logs.log; grep -c "Analysis complete" ../../logs/01_fastqc_preproc/combined_logs.log`                            |
| `01_multiqc_preproc.sh` & `03_multiqc_postproc.sh`            | Look at HTML report                  |
| `02_slurm_preprocessing.sh`            | Check for errors in log; check number of fastQ files in results/preprocessing/primerfree_reads & results/preprocessing/trimmed_filtered_reads; check output plot                  |
| `04_slurm_denoising.sh`            | Check for errors in log; check output plots                  |
| `05_assign_taxonomy.sh`            | Check for errors in log; check output plots                  |
| `06_quantify_strains.sh`            | Check for errors in log; check output plots                  |

**To run `06_quantify_strains.sh`, you first need to create your custom database (see below) and run `05_assign_taxonomy.sh` with your database as `db2`.**

---


## Creating a custom database for defined communities

This can be done entirely on your local computer. Create a conda environment containing `seqkit`, `cd-hit` and `dos2unix`.
It is highly recommended to use 16S sequences from PacBio sequencing (WGS or amplicon) because 16S inference with Illumina or ONT genomes can be inaccurate.

### Merging and dereplicating 16S sequences

1. Put all 16S sequences in a single folder called `individual_16S`. Each sequence must have a unique ID that contains the strain name.

2. Concatenate the sequences using `cat individual_16S/*.fna >> all_16S.fna`.

3. Run **cd-hit-est** with 100% identity threshold: `cd-hit-est -i all_16S.fna -o all_16S_cd-hit -c 1 -n 10 -d 0`.

4. Parse the `.clstr` output into a table:
```
awk '
BEGIN { print "cluster\tsequence"}
/^>Cluster/ {
    cluster = $0
    gsub(">Cluster ", "", cluster)
    cluster = "Cluster" cluster
    next
}
/^[0-9]/ {
    match($0, />[^ ]+/)
    seq = substr($0, RSTART+1, RLENGTH-1)
    gsub(/\.{3}/, "", seq)
    print cluster "\t" seq
}
' all_16S_cd-hit.clstr > all_16S_cd-hit_clusters.tsv
```

5. Open `all_16S_cd-hit_clusters.tsv` in excel or a text editor and check that each cluster contains only sequences from the same strain. If it is not the case, this is fine. You will just need to keep in mind later that some ASVs cannot be used to quantify the abundance of your strains.

You are now in possession of a dereplicated database of all the 16S amplicons of your community. But we still need to **(1)** generate versions of this database with correctly formatted headers **(2)** merge these new versions with existing toSpecies and toGenus databases to use with `assignTaxonomy()` **(3)** generate a version of this database with only the species name to use with `addSpecies()`.

### Formatting databases for dada2

6. Prepare a file to change the sequence IDs in `all_16S_cd-hit` to reflect the cluster name. We just need to swap the columns in `all_16S_cd-hit_clusters.tsv` and remove the header:
```
awk 'BEGIN { FS = OFS = "\t" }; NR > 1 { print $2, $1 }' all_16S_cd-hit_clusters.tsv > all_16S_cd-hit_clusters_rename.tsv
```

7. Change the sequence IDs in `all_16S_cd-hit`: `seqkit replace -p '^(\S+)' -r '{kv}$2' -k all_16S_cd-hit_clusters_rename.tsv all_16S_cd-hit > all_16S_cd-hit_renamed.fna`.

8. Create a copy of `all_16S_cd-hit_clusters.tsv` named `all_16S_cd-hit_clusters_tax_full.tsv`. Open it in excel, and add a third column `taxonomy_full` with the full taxonomy of your strains. If you will merge it with GreenGenes2, it must be in GTDB-like taxonomy format, like *d__Bacteria;p__Bacillota_I;c__Bacilli_A;o__Lactobacillales;f__Lactobacillaceae;g__Bombilactobacillus;s__Bombilactobacillus mellifer*. If you are merging with SILVA, the format would be *Bacteria;Bacillota;Bacilli;Lactobacillales;Lactobacillaceae;Bombilactobacillus;mellifer;*. Then use `TEXTSPLIT` from excel to get a fourth column `taxonomy_genus` with only the taxonomy down to genus (keeping a semi-colon at the end), and a fifth column `genus_species` with the full species name. Append the strain name and cluster number to `genus_species` using "-" as a delimiter. Finally, create a sixth column `taxonomy_species` that is basically the same as `taxonomy_full` but you remove the genus in the species name. See the attached example for merging with GreenGenes2. Once you are done editing it, make sure it is **tab-delimited** and in **Unix format** `dos2unix all_16S_cd-hit_clusters_tax_full.tsv`.

9. Make a copy of `all_16S_cd-hit_clusters_tax_full.tsv` with only columns 1 and 4, remove the header, and remove redundant lines: `awk ' BEGIN { FS = OFS = "\t" }; NR > 1 {print $1, $4}' all_16S_cd-hit_clusters_tax_full.tsv | uniq > all_16S_cd-hit_clusters_tax_genus.txt`. Now you should end up with a table reporting the genus-level taxonomy of each cluster. You might have had several different strains in one cluster at **step 5**, but they should still all belong to the same genus, therefore each cluster should appear only once in the table.

10. Make a copy of `all_16S_cd-hit_clusters_tax_full.tsv` with only columns 1 and 6, remove the header, and remove redundant lines: `awk ' BEGIN { FS = OFS = "\t" }; NR > 1 {print $1, $6}' all_16S_cd-hit_clusters_tax_full.tsv | uniq > all_16S_cd-hit_clusters_tax_species.txt`. Now you should end up with a table reporting the species-level taxonomy of each cluster. You might have had several different strains in one cluster at **step 5**, but they should still all belong to the same species, therefore each cluster should appear only once in the table.

11. Make another copy with only columns 1 and 5, remove the header, and remove redundant lines: `awk ' BEGIN { FS = OFS = "\t" }; NR > 1 {print $1, $5}' all_16S_cd-hit_clusters_tax_full.tsv | uniq > all_16S_cd-hit_clusters_gs.txt`. Now you should end up with a table reporting the exact species name of each cluster. If a cluster had several different strains/species in it, it will appear more than once in this file. Open the file and remove lines so that each cluster appears only once. Then give a custom name to these mixed clusters. The new name does not matter since you will not use ASVs assigned to those clusters anyways.

12. Create a copy of the database with suitable headers for `assignTaxonomy()` at genus level: `seqkit replace -p '^(\S+)' -r '{kv}$2' -k all_16S_cd-hit_clusters_tax_genus.txt all_16S_cd-hit_renamed.fna > syncom_custom_db_toGenus.fa`.

13. Create a copy of the database with suitable headers for `assignTaxonomy()` at species level: `seqkit replace -p '^(\S+)' -r '{kv}$2' -k all_16S_cd-hit_clusters_tax_species.txt all_16S_cd-hit_renamed.fna > syncom_custom_db_toSpecies.fa`.

14. Create a copy of the database with suitable headers for `addSpecies()`: `seqkit replace -p '^(\S+)' -r '${1} {kv}' -k all_16S_cd-hit_clusters_gs.txt all_16S_cd-hit_renamed.fna > syncom_custom_db_addSpecies.fa`.

15. Open your databases in a text editor and check that the formatting corresponds to the requirements described [here](https://benjjneb.github.io/dada2/training.html#formatting-custom-databases).

### Comparing and merging databases

16. Copy the published databases you want to merge your custom databases to into the current folder. You can download them [here](https://benjjneb.github.io/dada2/training.html#dada2-formatted-reference-databases). Make sure to download the *_toGenus_trainset*  and the *_toSpecies_trainset* files.

17. Compare the two databases for `assignTaxonomy()`: we basically want to remove sequences in the published database that are identical to sequences in our custom database.
    - Use `cd-hit-est-2d -i syncom_custom_db_toSpecies.fa -i2 gg2_2024_09_toSpecies_trainset.fa -o compare_gg2_custom_toSpecies -c 1 -n 10 -d 0` for the species-level database.
    - Use `cd-hit-est-2d -i syncom_custom_db_toGenus.fa -i2 gg2_2024_09_toGenus_trainset.fa -o compare_gg2_custom_toGenus -c 1 -n 10 -d 0` for the genus-level database.
    - The outputs are: **(1)** a fasta file with all sequences from GreenGenes2 that are **not identical** to sequences in your custom db, and **(2)** a `.clstr` text file listing similar sequences between the two databases.

18. Merge the non-redundant databases:
```
# species-level
seqkit seq compare_gg2_custom_toSpecies > syncom_custom_db_toSpecies_trainset.fa
cat syncom_custom_db_toSpecies.fa >> syncom_custom_db_toSpecies_trainset.fa
# genus-level
seqkit seq compare_gg2_custom_toGenus > syncom_custom_db_toGenus_trainset.fa
cat syncom_custom_db_toGenus.fa >> syncom_custom_db_toGenus_trainset.fa
```

You are now ready to use the custom databases with dada2. You will also `all_16S_cd-hit_clusters_tax_full.tsv` to run the strain quantification script. I personally like to use the `toSpecies_trainset` with `assignTaxonomy()`, but you can instead use the `toGenus_trainset` with this function to limit memory usage if you are not interested in the species classification of 'contaminants'.

---

## Authors

This pipeline was written by Meline Garcia. Many thanks to [Malick N`Diaye](https://github.com/MalickNdiye) and [Aiswarya Prasad](https://github.com/Aiswarya-prasad) for their suggestions.