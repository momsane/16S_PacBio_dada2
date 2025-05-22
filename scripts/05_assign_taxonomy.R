# Load required libraries

if(!require(tidyverse)){
  install.packages(pkgs = 'tidyverse', repos = 'https://stat.ethz.ch/CRAN/')
  library(tidyverse)
}

if(!require(ggplot2)){
  install.packages(pkgs = 'ggplot2', repos = 'https://stat.ethz.ch/CRAN/')
  library(ggplot2)
}

if(!require(dada2)){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("dada2")
  library(dada2)
}

if(!require(GUniFrac)){
  install.packages(pkgs = 'GUniFrac', repos = 'https://stat.ethz.ch/CRAN/')
  library(GUniFrac)
}

if(!require(Biostrings)){
  install.packages(pkgs = 'Biostrings', repos = 'https://stat.ethz.ch/CRAN/')
  library(Biostrings)
}

if(!require(phyloseq)){
  install.packages(pkgs = 'phyloseq', repos = 'https://stat.ethz.ch/CRAN/')
  library(phyloseq)
}

if(!require(fantaxtic)){
  if(!require(devtools)){
    install.packages(pkgs = 'devtools', repos = 'https://stat.ethz.ch/CRAN/')
    library(devtools)
  }
  devtools::install_github("gmteunisse/fantaxtic")
  library(fantaxtic)
}

if(!require(ggnested)){
  if(!require(devtools)){
    install.packages(pkgs = 'devtools', repos = 'https://stat.ethz.ch/CRAN/')
    library(devtools)
  }
  devtools::install_github("gmteunisse/ggnested")
  library(ggnested)
}

# Collect arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8){
  stop(" Usage: 03_assign_taxonomy.R <ASV_table> <metadata_table.tsv> <db_tax> <db_species> <rarefy_to> <facet_var> <tax_results_dir> <plots_dir>", call.=FALSE)
} else {
  input.asvs <- args[1] # ASV table (no chimera)
  input.metadata <- args[2] # sample metadata table, tab-separated, first column is the the sample name
  db1 <- args[3] # taxonomy database for assignTaxonomy (GreenGenes2, SILVA, or custom)
  db2 <- args[4] # taxonomy database for addSpecies (SILVA or custom), put "" if not needed
  rarefy_to <- args[5] # number of reads to rarefy to; if equals -1, no rarefaction
  facet_var <- args[6] # a column in the metadata table to facet the taxonomy plot, put "" if not needed
  out.tax <- args[7] # folder to write denoising results
  out.plots <- args[8] # folder to write plots
}

# input.asvs <- "/Volumes/gr_Engel/mgarcia/SAGE_tuto_16S_FM_2024/1_Results/denoising/ASV_samples_table_noChim.rds"
# input.metadata <- "/Volumes/gr_Engel/mgarcia/SAGE_tuto_16S_FM_2024/workflow/config/metadata_SAGE23.tsv"
# db1 <- "/Volumes/gr_Engel/mgarcia/SAGE_tuto_16S_FM_2024/0_Data/databases/gg2_2024_09_toSpecies_trainset.fa"
# db2 <- ""
# rarefy_to <- 513
# facet_var <- "HostSpecies"
# out.tax <- "/Volumes/gr_Engel/mgarcia/SAGE_tuto_16S_FM_2024/1_Results/assign_taxonomy"
# out.plots <- "/Volumes/gr_Engel/mgarcia/SAGE_tuto_16S_FM_2024/2_Plots"

### Create outdirs ###

print("Creating directories")

dir.create(out.plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out.tax, recursive = TRUE, showWarnings = FALSE)

### Inputs ###

rarefy_to <- as.numeric(rarefy_to)

ASV_samples_table_noChim <- readRDS(input.asvs)
meta <- read.table(input.metadata, sep = "\t", header = T)
rownames(meta) = meta[ ,1]

### Rarefy reads if needed ###

if (rarefy_to == -1){
  print("No rarefaction")
  ASV_samples_table_noChim2 <- ASV_samples_table_noChim
} else {
  # rarefy reads
  print(paste0("Rarefaction to ", rarefy_to, " reads"))
  set.seed(42)
  ASV_samples_table_noChim2 <- Rarefy(ASV_samples_table_noChim, depth = rarefy_to)$otu.tab.rff
  
  # find ASVs with total abundance of 0
  ab <- apply(ASV_samples_table_noChim2, 2, sum)
  print(paste0("Rarefaction led to ", length(names(ab[ab == 0])), " ASVs being removed"))
  
  # save these ASVs somewhere
  seqs <- DNAStringSet(names(ab[ab == 0]))
  names(seqs) <- seq_along(seqs)
  writeXStringSet(seqs, file.path(out.tax, "ASV_removed_rarefaction.fna"), format="fasta")
  
  # remove the ASVs from table
  ASV_samples_table_noChim2 <- ASV_samples_table_noChim2[ ,setdiff(colnames(ASV_samples_table_noChim2), names(ab[ab == 0]))]
}

### Assign taxonomy ###

print(paste0("Using database ", db1, " to assign taxonomy"))

print("Note: no database for addSpecies with GreenGenes2, but the toSpecies trainset is good enough with full-length 16S for assignTaxonomy.")

ASV_taxonomy <- assignTaxonomy(seqs = ASV_samples_table_noChim2, refFasta = db1, multithread = TRUE, verbose = T)

if (db2 != ""){
  print(paste0("Using database ", db2, " to assign species"))
  ASV_taxonomy <- addSpecies(taxtab = ASV_taxonomy, refFasta = db2, verbose = T)
}

saveRDS(ASV_taxonomy, file.path(out.tax, "ASV_Taxonomy.RDS"))

print("Finished assigning taxonomy")


### Reformat species name to classic nomenclature (with Genus + Species or next assigned taxonomic level) ###

print("Reformatting species names")
print("For ASVs with no assigned genus, the lowest assigned taxonomic rank is used as species name")

tax_df2 <- as.data.frame(ASV_taxonomy)
tax_df2$ASV <- rownames(tax_df2)

tax_df2 <- tax_df2 %>% 
  mutate(across(where(is.character), str_trim)) %>% # removing leading and trailing whitespace
  mutate(lowest_assign_taxon = case_when(
    is.na(Class) ~ Phylum,
    is.na(Order) ~ Class,
    is.na(Family) ~ Order,
    is.na(Genus) ~ Family,
    is.na(Species) ~ Genus
  )) %>%
  separate(Species, into = c(NA, "clean_species"), sep = "s__", remove = F) %>% 
  mutate(Species_full = case_when(
    !is.na(clean_species) ~ paste(Genus, clean_species),
    !is.na(Genus) ~ paste(Genus, "sp."),
    is.na(Genus) ~ paste(lowest_assign_taxon, "sp.")
  )) %>% 
  mutate(Species_full = gsub("g__","s__", Species_full)) %>% 
  select(-c("clean_species", "lowest_assign_taxon", "ASV"))



### Create phyloseq object ###

print("Creating phyloseq object")

ps <- phyloseq(otu_table(ASV_samples_table_noChim2, taxa_are_rows=F), 
               sample_data(meta),
               tax_table(as.matrix(tax_df2)))


### Remove non-bacterial ASVs if present ###

print("Filtering out non-bacterial and chloroplast ASVs")

#check for non-bacterial ASVs
print(paste0("Represented domains: ", paste0(unique(tax_df2$Kingdom), collapse = ", ")))
print(paste0("Represented classes: ", paste0(sort(unique(tax_df2$Class)), collapse = ", ")))

# remove non-bacterial ASVs
ps <- subset_taxa(ps, Kingdom == "d__Bacteria" | is.na(Kingdom))
ps <- subset_taxa(ps, Class != "c__Chloroplast" | is.na(Class))
ps <- subset_taxa(ps, Family != "f__Mitochondria" | is.na(Family))

print(paste0("Found and removed ", length(which(ASV_taxonomy[,1] != "d__Bacteria")), " non-bacterial ASV(s)"))
print(paste0("Found and removed ", length(which(ASV_taxonomy[,3] == "c__Chloroplast")), " chloroplast ASV(s)"))
print(paste0("Found and removed ", length(which(ASV_taxonomy[,5] == "f__Mitochondria")), " mitochondrial ASV(s)"))

# save sequences as refseq and give new names to ASVs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

saveRDS(object = ps, file = file.path(out.tax, "phyloseq_object.RDS"))

print("Finished creating and filtering phyloseq object")


### Check for single- and doubletons ###

print("Checking for singletons and doubletons - note that singletons may be introduced by rarefaction")

tab <- otu_table(ps) # samples are rows and ASV are columns
class(tab) <- "matrix" # warning but it's fine

tot <- sort(apply(X=tab, MARGIN=2, sum)) # total abundance of each ASV
# show ASVs with total abundance of 1 or 2
sdtons <- names(tot[tot < 3])
print("The following ASVs have a total abundance of 1 or 2:")
print(tax_table(ps)[sdtons, "Species_full"])
extr <- tab[,sdtons]
extr_sum <- apply(extr, 1, sum) # show only samples in which they are present
print("Samples in which these ASVs were found:")
print(tab[names(extr_sum[extr_sum != 0]),sdtons])


### Plotting most abundant taxa ###

print("Plotting most abundant taxa")

top_nested <- nested_top_taxa(
  ps,
  top_tax_level = "Order", # most abundant order
  nested_tax_level = "Genus", # most abundant genera within those orders
  n_top_taxa = 5, # top 5 most abundant orders
  n_nested_taxa = 3 # top 5 most abundant genera
)

tax_plot = plot_nested_bar(ps_obj = top_nested$ps_obj, 
                           top_level = "Order", nested_level = "Genus",
                           legend_title = "Taxonomy") + 
  theme_nested(theme_classic) + ylab("Relative abundance") + 
  theme(strip.text.x = element_text(angle=90), 
        strip.background=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  guides(fill=guide_legend(ncol =1))

if (facet_var != ""){
  tax_plot <- tax_plot + facet_wrap(as.formula(paste0("~", facet_var)), drop=T)
}

ggsave(file.path(out.plots, "03_taxonomy_barplot.pdf"), tax_plot, device="pdf", width = 8, height = 6)

print("Taxonomy assignment done")
