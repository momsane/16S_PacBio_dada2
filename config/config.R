#file:config.R

facet_var = "SampleType" # one column in the metadata table to facet the taxonomy plot, or ""

# Pre-processing parameters
input.raw = "path/to/raw_reads" # folder with all raw read files or folder with pre-rarefied read files
fwd.primer = "AGRGTTYGATYMTGGCTCAG"
rev.primer = "RGYTACCTTGTTACGACTT"
minLen = 1200 # discard reads below this length
maxLen = 1700 # discard reads above this length
maxEE = 2 # discard reads with more than maxEE errors - use 3-4 for Kinnex
overwrite = "F" # whether to overwrite existing files ("T") or not ("F") - useful if your job did not finish in time

# Denoising parameters
errModel = "binnedQualErrfun" # dada2 function to estimate the error model
maxBases = 1E10 # max number of bases to use for error model inference
db2 = "path/to/_addSpecies.fa" # database of ASVs expected in the samples; or ""
pool = "F" # "T" or "pseudo" or "F", whether to pool samples for ASV inference
maxraref_denoising = 8000 # maximum number of reads to build rarefaction curves; -1 to disable rarefaction curves

# Taxonomic assignment parameters
db1 = "path/to/_trainset.fa" # taxonomy database (GreenGenes2, SILVA, or custom)
min_boot = 50 # numerical threshold to retain taxonomic assignment
rarefy_taxonomy = 37000 # number of reads to rarefy to; -1 to disable rarefaction

# Strain quantification parameters
input.ps = "path/to/phyloseq_object_filtered_rarefied.RDS" # change if you want to use the rarefied one or not
input.clusters = "path/to/all_16S_cd-hit_clusters_tax_full.tsv" # table of ASVs with their assigned cd-hit cluster and user-input taxonomy
input.qpcr = "path/to/qPCR_results_analyzed.tsv" # qpcr data, if available
abundance_col = "copies_16S_sample" # name of the column in qpcr containing the abundance in the sample
maxraref_strains = -1 # maximum number of 'cells' to build rarefaction curves; -1 to disable rarefaction curves