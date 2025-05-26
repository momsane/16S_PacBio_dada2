# Analysis of 16S PacBio data with dada2

## Instructions

### Overview

The pipeline consists of a few main steps:

- *(optional step)* pre-rarefaction of the reads for large samples `00_rarefy.sh`
- quality check on the raw reads `01_fastqc_preproc.sh` and `01_multiqc_preproc.sh`
- pre-processing of the reads with dada2 (primer removal, length and quality trimming) `02_preprocessing.sh`
- quality check on the processed reads `03_fastqc_postproc.sh` and `03_multiqc_posteproc.sh`
- denoising into ASVs with dada2 `04_denoising.sh`
- taxonmy assignment with dada2 `05_assign_taxonomy.sh`.
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

You'll need access to a conda installation on your cluster. You can either have yours, installed for instance through miniforge3/miniconda3, or use the one provided by the cluster. This will change slightly how you activate environments at the beginning of the scripts. To use the cluster conda, follow the instructions [here](https://wiki.unil.ch/ci/books/high-performance-computing-hpc/page/using-conda-and-anaconda). To use your own conda installation, you just need to modify the `$CONDA_HOME` variable in the scripts.

Install the required conda environments using the .yaml files located in envs with the command `conda env create -f envs/<env>.yaml`.

### Data Preparation

Before running the pipeline, you need to prepare your data:

1.  **Metadata File:** modify `config/metadata.tsv` according to your samples. You do not need to keep the same columns except for the first one. This first column must contain the sample names (filenames without the `.fastq.gz` extension).
2.  **File Naming Table:** the raw read files you got from the sequencing facility have long non-informative names. You will rename them with the SampleID. Create a table like `config/rename_files.tsv` where the first column is the current name of each file, and the second column is the new name. This table has no header. It is important that this table be in Unix format. If you modify it in excel for instance, it will not be in Unix format. You can use the command-line tool dos2unix to convert it. Make sure there is a line return after the last line in the table.
3.  **Read Rarefaction Table (Optional):** if you have very uneven depth in your dataset with very large samples, you might want to consider rarefying the raw reads to limit unnecessary computation time and resources. Modify `config/pre_rarefaction.tsv` according to your samples. As with `config/rename_files.tsv`, make sure it is in Unix format. If you do not want to rarefy your reads, skip this step.
4.  **Raw Reads:** you are now ready to copy them from the nas. Modify the script `copy_rename_files.sh` with the correct paths. Then execute it from the login node (i.e. use `bash` instead of `sbatch` to submit it), because the NAS can only be accessed from the login node.

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

### Running the Pipeline

To run each script:

1.  **Submit the job to Slurm:** use `sbatch <script_name>.sh`.
2.  **Monitor the job:** use `Squeue` to check the status of you jobs. Check the logs for progress and errors. 


## Authors

This pipeline was written by Meline Garcia. Many thanks to Malick N`Diaye for his suggestions.
