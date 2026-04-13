# Analysis of 16S PacBio data with dada2


Table of contents

- [Running the dada2 pipeline](#running-the-dada2-pipeline)
    - [Overview](#adapting-the-scripts)
    - [Requirements](#adapting-the-scripts)
    - [Setting up the work environment](#setting-up-the-work-environment)
    - [Data preparation](#data-preparation)
    - [Adapting the scripts](#adapting-the-scripts)
    - [Running the pipeline](#running-the-pipeline)
- [Custom database for defined communities](#creating-a-custom-database-for-defined-communities)
    - [Install required tools](#install-required-tools)
    - [Merging and dereplicating 16S sequences](#merging-and-dereplicating-16s-sequences)
    - [Formatting databases for dada2](#formatting-databases-for-dada2)
    - [Comparing and merging databases](#comparing-and-merging-databases)
- [Authors](#authors)
- [Citing](#citing)

---

## Running the dada2 pipeline

**Important note 1**: this pipeline should be executed independently for separate runs, as the error profile is unique to each run.

**Important note 2**: depending on the type of library you have and how the fastQ files were produced, you need to select between different options for the filtering and denoising steps. 
- Consensus sequences from a Kinnex library have less passes, therefore the quality is lower. I would recommend using maxEE=3-4 to avoid being too strict.
- The CSS software used to infer consensus sequences bins the quality scores to reduce file size, see [here](https://ccs.how/faq/qv-binning.html). If you have such binned quality scores, you need to choose the appropriate function for the **denoising** step. 

### Overview

Steps of the pipeline:

- *(optional step)* copy and rename raw reads files with sample ID if needed `00_copy_rename_files.sh`
- *(optional step)* pre-rarefaction of the reads `00_rarefy.sh`
- quality check on the raw reads `01_fastqc_preproc.sh` and `01_multiqc_preproc.sh`
- pre-processing of the reads with dada2 (primer removal, length and quality trimming) `02_slurm_preprocessing.sh`
- quality check on the processed reads `03_fastqc_postproc.sh` and `03_multiqc_posteproc.sh`
- denoising into ASVs with dada2 `04_slurm_denoising.sh`
- taxonomy assignment with dada2 `05_slurm_assign_taxonomy.sh`
- *(only for defined communities)* compute genome equivalents from ASV table `06_slurm_quantify_strains.sh`, with or without qPCR data.

dada2 is implemented in Rscripts called by the bash scripts. **There should not be any need to modify the R scripts**.

### Requirements

- This pipeline is designed to be executed on a slurm cluster.
- Conda/mamba.
- R, fastQC, multiQC and bbmap (see *envs*).

### Setting up the work environment

#### Structure

Set up your working directory:
1. Create a folder with the name of your project.
2. Enter this folder, clone this git using `git clone https://github.com/momsane/16S_PacBio_dada2` then rename the folder with `mv 16S_PacBio_dada2 workflow`.
3. Create additional folders to obtain the following tree:
```
.
├── data
    ├── raw_reads
    └── databases
├── logs
├── plots
├── results
└── workflow
    ├── config
    ├── envs
    └── scripts
```

#### Conda environments

You'll need access to a conda installation. You can either have yours, installed for instance through miniforge3, or use the one provided by the cluster. This will change slightly how you activate environments at the beginning of the scripts.
- To use the cluster conda, follow the instructions [here](https://wiki.unil.ch/ci/books/high-performance-computing-hpc/page/using-conda-and-anaconda).
- To use your own conda installation, you just need to modify the `$CONDA_HOME` variable in the scripts.

Install all the required conda environments using the .yaml files located in the *envs* folder with the command `conda env create -f workflow/envs/<env>.yaml`.

### Data Preparation

Before running the pipeline, you need to prepare some data. All files in `/config` should be tab-separated and in Unix format. The bash scripts include a **dos2unix** command to convert them. Make sure also there is a line return after the last row of the table otherwise it will not be read. Finally, the use of special characters (including spaces) other than _- in file names or tables must be avoided.

1.  **File naming table:** the raw read files you got from the sequencing facility have long non-informative names. If not done already, you will rename them with the SampleID. Create a table like `config/rename_files.tsv` where the first column is the current name of each file, and the second column is the new name. This table has no header. If you have samples from different pools, you will need to create one table per pool because some samples might have the same original name.
2.  **Metadata file:** modify `config/metadata.tsv` according to your samples. You do not need to keep the same columns except for the first one, `SampleID`. This first column must contain the sample names (final filenames without the `.fastq.gz` extension). Make sure there are no empty cells in this table - use NA values if necessary.
3.  **Read rarefaction table (optional):** if you have very uneven depth in your dataset, you might want to consider rarefying the raw reads to limit unnecessary computation time and resources for large samples. Modify `config/pre_rarefaction.tsv` according to your needs.
4.  **Raw reads:** you are now ready to copy them from the NAS. Modify the script `00_copy_rename_files.sh` with the correct paths. Then execute it from the login node (*i.e.* use `bash` instead of `sbatch` to submit it). If you have samples from different pools, you will need to execute this script independently for each pool.
5. **Databases:** you need to provide at least one database to assign taxonomy to your ASVs. Refer to [the dada2 website](https://benjjneb.github.io/dada2/training.html) for more information and links to download the databases.

### Adapting the scripts

Only the beginning of the `.sh` scripts needs to be modified:

- the commands to initialize conda according to the type of installation you are using
- the input variables, for instance the path to the root directory
- the pre-rarefaction and the fastQC scripts are array jobs (argument `--array` in the slurm header), so you need to modify the range of the arrays. `2-50` means you will process files described in lines 2 to 50 of `config/metadata.tsv`. We start at 2 to skip the header. So your array range should be `2-<number of samples + 1>`
- you should not need to modify the resource requirements, unless your jobs get killed. Before increasing memory and CPU requests, check the efficiency of your job using `seff <jobid>`.

### Running the pipeline

`00_copy_rename_files.sh` must be run from the login node using `bash` instead of `sbatch`.

To run the other scripts:

1.  **Submit the job to the slurm scheduler:** `sbatch <script_name>.sh`.
2.  **Monitor the job:** use `Squeue` to check the status of you jobs. **Check the log file after each step**: it will contain not only any error messages but also useful information. 
    
| Script           | What to check                                       |
|--------------------|---------------------------------------------------|
| `00_copy_rename_files.sh`        | All files have been copied; use `ls -la` to check that file names don't contain special characters                            |
| `00_rarefy.sh`        | Size of output fastQ files is > 0                            |
| `01_fastqc_preproc.sh` & `03_fastqc_postproc.sh`      | Output folders are not empty                            |
| `01_multiqc_preproc.sh` & `03_multiqc_postproc.sh`            | Log file indicates the expected of reports has been found; look at HTML report                  |
| `02_slurm_preprocessing.sh`            | Check log and plots; check number of fastQ files in results/preprocessing/primerfree_reads & results/preprocessing/trimmed_filtered_reads;                  |
| `04_slurm_denoising.sh`            | Check log and plots                  |
| `05_slurm_assign_taxonomy.sh`            | Check log and plots                  |
| `06_slurm_quantify_strains.sh`            | Check log and plots                  |


**Note 1:** the first time you run `02_slurm_preprocessing.sh` and `05_slurm_assign_taxonomy.sh`, some R packages will be installed. The execution might be halted with an error message after the last installation. This is because the R environment needs to be reloaded. Simply run the script again and it should work.

**Note 2:** to use `06_slurm_quantify_strains.sh`, you first need to create your custom database (see below) and run `05_slurm_assign_taxonomy.sh` with this custom database as `db2`.

---


## Creating a custom database for defined communities

This can be done entirely on a regular computer. It is highly recommended to use 16S sequences from PacBio sequencing (WGS or amplicon) because 16S inference with Illumina or ONT genomes can be inaccurate.

### Install required tools

Create and activate the following conda environment and define file names for the workflow:
```
conda create custom_db_dada2 bioconda::seqkit bioconda::cd-hit conda-forge::dos2unix
conda activate custom_db_dada2
prefix=all_16S_cd-hit
dbprefix=syncom_custom_db
```

### Merging and dereplicating 16S sequences

1. If needed:
  - put all 16S sequences in a single folder called `individual_16S`. **Each sequence must have a unique ID containing the strain name.**
  - concatenate the sequences: `cat individual_16S/*.fna >> all_16S.fna`.

2. Dereplicate sequences at 100% identity threshold: `cd-hit-est -i all_16S.fna -o "$prefix" -c 1 -n 10 -d 0`.

3. Parse the `.clstr` output into a table:
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
' "$prefix".clstr > "$prefix"_clusters.tsv
```

4. Open the `_clusters.tsv` table in Excel or a text editor and check that each cluster contains only sequences from the same strain. If it is not the case, this is fine. You will just need to keep in mind later that some ASVs cannot be used to quantify the abundance of your strains.

You are now in possession of a dereplicated database of all the 16S amplicons of your community. But we still need to **(1)** generate versions of this database with correctly formatted headers for dada2 **(2)** merge these new versions with existing databases to use with `assignTaxonomy()` **(3)** generate a version of this database with only the species name to use with `addSpecies()`.

### Formatting databases for dada2

5. Change the sequence IDs in the fasta file to reflect the cluster name:
```
awk 'BEGIN { FS = OFS = "\t" }; NR > 1 { print $2, $1 }' "$prefix"_clusters.tsv > "$prefix"_clusters_rename.tsv
seqkit replace -p '^(\S+)' -r '{kv}$2' -k "$prefix"_clusters_rename.tsv "$prefix" > "$prefix"_renamed.fna
```

6. 
    - Create a copy of the clusters table: `cp "$prefix"_clusters.tsv "$prefix"_clusters_tax_full.tsv`.
    - Open it in excel, and add a third column `taxonomy_full` with the full taxonomy of your strains. If you will merge it with GreenGenes2, it must be in GTDB-like taxonomy format, like *d__Bacteria;p__Bacillota_I;c__Bacilli_A;o__Lactobacillales;f__Lactobacillaceae;g__Bombilactobacillus;s__Bombilactobacillus mellifer*. If you are merging with SILVA, the format would be *Bacteria;Bacillota;Bacilli;Lactobacillales;Lactobacillaceae;Bombilactobacillus;mellifer;*.
    - Use `TEXTSPLIT` from excel to get a fourth column `taxonomy_genus` with only the taxonomy down to genus (keeping a semi-colon at the end), and a fifth column `genus_species` with the full species name.
    - Create a sixth column `strain` with the strain name.
    - Append the strain name and cluster number to the existing `genus_species` column using "-" as a delimiter.
    - Create a seventh column `taxonomy_species` that is basically the same as `taxonomy_full` but you remove the genus in the species name.
    - Create an eighth column `n_copies` indicating the number of copies of this specific ASV in the given strain. Put NA if you don't know.
    - See the example in `config` with GreenGenes2 format. Once you are done editing it, make sure it is **tab-delimited** and in **Unix format** `dos2unix "$prefix"_clusters_tax_full.tsv`.

7. Next we prepare some files to generate the toGenus and toSpecies databases:
```
awk ' BEGIN { FS = OFS = "\t" }; NR > 1 {print $1, $4}' "$prefix"_clusters_tax_full.tsv | sort -k1,1 | uniq > "$prefix"_clusters_tax_genus.txt 
awk ' BEGIN { FS = OFS = "\t" }; NR > 1 {print $1, $7}' "$prefix"_clusters_tax_full.tsv | sort -k1,1 | uniq > "$prefix"_clusters_tax_species.txt
```
In these tables, if a cluster appears twice but with different taxonomy (an ASV present in two different species for instance), you need to manually edit it so that this cluster appears only once. You can also edit the taxonomy column to replace the species epithet by a dummy value.

8. Next we prepare some files to generate the addSpecies database:
```
awk ' BEGIN { FS = OFS = "\t" }; NR > 1 {print $1, $5}' "$prefix"_clusters_tax_full.tsv | sort -k1,1 | uniq > "$prefix"_clusters_gs.txt
```
Similarly, if a cluster appears twice but with different taxonomy, you need to manually edit the table so that this cluster appears only once. You can also edit the taxonomy column to replace the species epithet and strain name by dummy values.

9. Now we create custom databases with suitable headers for each dada2 function:
```
# genus-level
seqkit replace -p '^(\S+)' -r '{kv}$2' -k "$prefix"_clusters_tax_genus.txt "$prefix"_renamed.fna > "$dbprefix"_toGenus.fa
# species-level
seqkit replace -p '^(\S+)' -r '{kv}$2' -k "$prefix"_clusters_tax_species.txt "$prefix"_renamed.fna > "$dbprefix"_toSpecies.fa
# addSpecies
seqkit replace -p '^(\S+)' -r '${1} {kv}' -k "$prefix"_clusters_gs.txt "$prefix"_renamed.fna > "$dbprefix"_addSpecies.fa
```

10. Open your databases in a text editor and check that the formatting corresponds to the requirements described [here](https://benjjneb.github.io/dada2/training.html#formatting-custom-databases).

### Comparing and merging databases

11. Copy the published databases you want to merge your custom databases to into the current folder. You can download them [here](https://benjjneb.github.io/dada2/training.html#dada2-formatted-reference-databases). Make sure to download the *_toGenus_trainset*  and the *_toSpecies_trainset* files.

12. Compare the two databases for `assignTaxonomy()`: we basically want to remove sequences in the published database that are identical to sequences in our custom database.
```
# species-level
cd-hit-est-2d -i "$dbprefix"_toSpecies.fa -i2 gg2_2024_09_toSpecies_trainset.fa -o compare_gg2_custom_toSpecies -c 1 -n 10 -d 0
# genus-level
cd-hit-est-2d -i "$dbprefix"_toGenus.fa -i2 gg2_2024_09_toGenus_trainset.fa -o compare_gg2_custom_toGenus -c 1 -n 10 -d 0
```
The outputs are: **(1)** a fasta file with all sequences from GreenGenes2 that are **not identical** to sequences in your custom db, and **(2)** a `.clstr` text file listing similar sequences between the two databases.

13. Merge the non-redundant databases:
```
# species-level
seqkit seq compare_gg2_custom_toSpecies > "$dbprefix"_toSpecies_trainset.fa
cat "$dbprefix"_toSpecies.fa >> "$dbprefix"_toSpecies_trainset.fa
# genus-level
seqkit seq compare_gg2_custom_toGenus > "$dbprefix"_toGenus_trainset.fa
cat "$dbprefix"_toGenus.fa >> "$dbprefix"_toGenus_trainset.fa
```

You are now ready to use the custom databases with dada2. You will also need the clusters table to run the strain quantification script. I personally like to use the `toSpecies_trainset` with `assignTaxonomy()`, but you can instead use the `toGenus_trainset` with this function to if you are not interested in the species-level classification of 'contaminants'.

---

## Authors

This pipeline was written by Meline Garcia. Many thanks to [Malick N`Diaye](https://github.com/MalickNdiye) and [Aiswarya Prasad](https://github.com/Aiswarya-prasad) for their suggestions.

---

## Citing

If you use this pipeline, please cite this repository as well as its main dependencies:

- [dada2](https://github.com/benjjneb/dada2): Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016). “DADA2: High-resolution sample inference from Illumina amplicon data.” Nature Methods, 13, 581-583. doi:10.1038/nmeth.3869. 
- [bbmap](https://github.com/bbushnell/BBTools): Bushnell, B. (2014) BBMap: A Fast, Accurate, Splice-Aware Aligner. 
9th Annual Genomics of Energy & Environment Meeting, Walnut Creek, CA.
- [fastQC](https://github.com/s-andrews/FastQC): https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- [multiQC](https://github.com/MultiQC/MultiQC): Philip Ewels, Måns Magnusson, Sverker Lundin, Max Käller, MultiQC: summarize analysis results for multiple tools and samples in a single report, Bioinformatics, Volume 32, Issue 19, October 2016, Pages 3047–3048, https://doi.org/10.1093/bioinformatics/btw354
- [iNEXT](https://github.com/AnneChao/iNEXT): Hsieh, T.C., Ma, K.H. & Chao, A. (2016) iNEXT: An R package for interpolation and extrapolation of species diversity (Hill numbers). Methods in Ecology and Evolution, 7, 1451-1456.
