# Load required libraries

if(!require(dplyr)){
  install.packages(pkgs = 'dplyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(dplyr)
}

if(!require(tidyr)){
  install.packages(pkgs = 'tidyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(tidyr)
}

if(!require(stringr)){
  install.packages(pkgs = 'stringr', repos = 'https://stat.ethz.ch/CRAN/')
  library(stringr)
}

if(!require(ggplot2)){
  install.packages(pkgs = 'ggplot2', repos = 'https://stat.ethz.ch/CRAN/')
  library(ggplot2)
}

if(!require(dada2)){
  if(!require(devtools)){
    install.packages(pkgs = 'devtools', repos = 'https://stat.ethz.ch/CRAN/')
  }
  devtools::install_github("benjjneb/dada2") # installing through GitHub to get latest updates
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
  if(!requireNamespace("BiocManager")){
    install.packages("BiocManager")
  }
  BiocManager::install("phyloseq")
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

if (length(args) != 9){
  stop(" Usage: 05_assign_taxonomy.R <ASV_table> <metadata_table.tsv> <db_tax> <db_species> <clusters_tax> <rarefy_to> <facet_var> <tax_results_dir> <plots_dir>", call.=FALSE)
} else {
  input.asvs <- args[1] # ASV table (no chimera)
  input.metadata <- args[2] # sample metadata table, tab-separated, first column is the the sample name
  db1 <- args[3] # taxonomy database for assignTaxonomy (GreenGenes2, SILVA, or custom)
  db2 <- args[4] # taxonomy database for addSpecies (SILVA or custom), put "" if not needed
  input.clusters <- args[5] # table of ASVs with their assigned cd-hit cluster and user-input taxonomy
  rarefy_to <- args[6] # number of reads to rarefy to; if equals -1, no rarefaction
  facet_var <- args[7] # one column in the metadata table to facet the taxonomy plot, put "" if not needed
  out.tax <- args[8] # folder to write denoising results
  out.plots <- args[9] # folder to write plots
}

# # root <- "/Volumes/D2c/mgarcia/20240708_mgarcia_syncom_invivo/exp01_inoculation_methods/pacbio_analysis"
# root <- "/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis/run1_bees"
# input.asvs <- file.path(root, "results", "denoising", "ASV_samples_table_noChim.rds")
# input.metadata <- file.path(root, "workflow", "config", "metadata.tsv")
# db1 <- file.path(root, "data", "databases", "syncom_custom_db_toSpecies_trainset.fa")
# db2 <- file.path(root, "data", "databases", "syncom_custom_db_addSpecies.fa")
# input.clusters <- file.path(root, "workflow", "config", "all_16S_cd-hit_clusters_tax_full.tsv")
# rarefy_to <- -1
# facet_var <- "SampleType"
# out.tax <- file.path(root, "results", "assign_taxonomy")
# out.plots <- file.path(root, "plots")

### Create outdirs ###

cat("\nCreating directories\n")

dir.create(out.plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out.tax, recursive = TRUE, showWarnings = FALSE)

### Inputs ###

rarefy_to <- as.numeric(rarefy_to)

ASV_samples_table_noChim <- readRDS(input.asvs)
meta <- read.table(input.metadata, sep = "\t", header = T)
rownames(meta) = meta[ ,1]

clusters <- read.table(input.clusters, sep = "\t", header = T)

### Rarefy reads if needed ###

if (rarefy_to == -1){
  cat("No rarefaction\n")
  ASV_samples_table_noChim2 <- ASV_samples_table_noChim
} else {
  # rarefy reads
  cat(paste0("Rarefaction to ", rarefy_to, " reads; sample with fewer reads will be kept as is; ASVs with a subsequent total abundance of 0 will be removed\n"))
  set.seed(42)
  raref <- Rarefy(ASV_samples_table_noChim, depth = rarefy_to)
  ASV_samples_table_noChim2 <- raref$otu.tab.rff 
  # add samples that did not have enough reads because Rarefy() removes them
  ASV_samples_table_noChim2 <- rbind(ASV_samples_table_noChim2, ASV_samples_table_noChim[raref$discard, ])
  
  # find ASVs with total abundance of 0
  ab <- apply(ASV_samples_table_noChim2, 2, sum)
  cat(paste0("Rarefaction led to ", length(names(ab[ab == 0])), " ASVs being removed\n"))
  
  # save these ASVs somewhere
  seqs <- DNAStringSet(names(ab[ab == 0]))
  names(seqs) <- seq_along(seqs)
  writeXStringSet(seqs, file.path(out.tax, "ASV_removed_rarefaction.fna"), format="fasta")
  
  # remove the ASVs from table
  ASV_samples_table_noChim2 <- ASV_samples_table_noChim2[ ,setdiff(colnames(ASV_samples_table_noChim2), names(ab[ab == 0]))]
}

### Assign taxonomy ###

cat(paste0("Using database ", db1, " to assign taxonomy\n"))

ASV_taxonomy <- assignTaxonomy(
  seqs = ASV_samples_table_noChim2,
  refFasta = db1,
  minBoot = 80,
  multithread = TRUE,
  verbose = T)

##### Custom function to assign species

# addSpecies uses exact matches, including length matching
# but in some cases either the reference sequence or the ASV could be a few nt shorter

addSpecies_custom <- function(seqs, refFasta) {
  # read inputs as DNAstrings
  asvs <- DNAStringSet(colnames(seqs))
  refs <- readDNAStringSet(refFasta)
  # convert to character for comparison
  asv_seqs <- as.character(asvs)
  ref_seqs <- as.character(refs)
  
  # extract only the description of the ref sequences
  full_headers <- names(refs)
  ref_descriptions <- sub("^[^ ]+\\s+", "", full_headers)  # keep only text after first space character
  
  results <- data.frame(ASV = character(0), ASV_index= character(0), ref_match = character(0), stringsAsFactors = FALSE)
  
  for (i in seq_along(asv_seqs)) {
    asv <- asv_seqs[i]
    
    # match if ASV is substring of reference or reference is substring of ASV
    # + try the rev. complement
    matched <- vapply(ref_seqs, function(ref) {
      rc_ref <- rc(ref)
      grepl(asv, ref, fixed = TRUE) ||
        grepl(ref, asv, fixed = TRUE) ||
        grepl(asv, rc_ref, fixed = TRUE) ||
        grepl(ref, rc_ref, fixed = TRUE)
    }, logical(1))
    
    if (any(matched)) {
      matched_descriptions <- ref_descriptions[matched]
      if (length(matched_descriptions) > 1){
        cat(paste0("Found more than one match for ASV ", asv_name, "; not adding it to the results\n"))
      } else {
        temp_df <- data.frame(ASV = asv,
                              ASV_index = i,
                              ref_match = matched_descriptions,
                              stringsAsFactors = FALSE)
        results <- rbind(results, temp_df)
      }
    }
  }
  rownames(results) <- results$ASV
  return(results[ ,-1])
}

if (db2 != ""){
  cat(paste0("Using database ", db2, " to assign species\n"))
  ASV_taxonomy2 <- addSpecies_custom(seqs = ASV_samples_table_noChim2, refFasta = db2)
  ASV_taxonomy2 <- ASV_taxonomy2 %>% 
    separate(ref_match, into = c("Species_addSp", "Strain", "Cluster"), sep = "-", remove = T) %>% 
    separate(Species_addSp, into = c("Genus_addSp", "Species_addSp"), sep = " ") %>% 
    select(-"ASV_index")
  
  # merge with taxonomy from assignTaxonomy()
  ASV_taxonomy3 <- transform(merge(ASV_taxonomy,ASV_taxonomy2,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  
  # combine taxonomy from both functions
  ASV_taxonomy3 <- ASV_taxonomy3 %>% 
    mutate(Genus_clean = str_remove(Genus, "^g__")) %>% 
    mutate(Species_clean = str_remove(Species, "^s__")) %>% 
    mutate(Genus_clean = if_else(!is.na(Genus_addSp),Genus_addSp,Genus_clean)) %>% # overwriting genus with genus from addSpecies if there is a hit
    mutate(Species_clean = case_when(
      !is.na(Genus_addSp) ~ Species_addSp,
      is.na(Genus_addSp) & !is.na(Genus_clean) ~ Species_clean,
      is.na(Genus_addSp) & is.na(Genus_clean) ~ Species_clean,
    )) %>% # overwriting species with species from addSpecies if there is a hit
    mutate(Species_full = case_when(
      is.na(Genus_clean) ~ "s__",
      !is.na(Genus_clean) & is.na(Species_clean) ~ paste("s__", Genus_clean, " sp.", sep = ""),
      !is.na(Genus_clean) & !is.na(Species_clean) ~ paste("s__", Genus_clean, " ", Species_clean, sep = "")
    )) %>% # rewrite full species name
    dplyr::select(-c("Species", "Species_clean", "Genus", "Genus_addSp", "Species_addSp")) %>%
    dplyr::rename(Genus = Genus_clean) %>% 
    dplyr::rename(Species = Species_full) %>% 
    dplyr::relocate("Species", .after = "Genus") %>%
    mutate(inferred_from = if_else(is.na(Cluster), "assignTaxonomy", "addSpecies_custom")) %>%
    arrange(Species)
  
  cat(paste0("Matched ", nrow(ASV_taxonomy3[ASV_taxonomy3$inferred_from == "addSpecies_custom", ]), " ASV(s) to reference sequences\n"))
  
} else {
  ASV_taxonomy3 <- as.data.frame(ASV_taxonomy)
  if (!("Species" %in% colnames(ASV_taxonomy3))){
    ASV_taxonomy3$Species <- NA
  }
  ASV_taxonomy3 <- ASV_taxonomy3 %>% 
    mutate(Strain = NA) %>%
    mutate(Cluster = NA) %>%
    mutate(inferred_from = "assignTaxonomy")
}

ASV_taxonomy4 <- ASV_taxonomy3 %>% 
  mutate(ASV = rownames(ASV_taxonomy3), .before = "Kingdom")

cat("Finished assigning taxonomy\n")


### Reformat species name ###

cat("Reformatting species names\n")
cat("For ASVs with no assigned genus, the lowest assigned taxonomic rank is used as species name.\n")

ASV_taxonomy3 <- ASV_taxonomy3 %>% 
  mutate(lowest_assign_taxon = case_when(
    is.na(Phylum) ~ Kingdom,
    is.na(Class) ~ Phylum,
    is.na(Order) ~ Class,
    is.na(Family) ~ Order,
    is.na(Genus) ~ Family
  )) %>%
  mutate(Species = if_else(Species == "s__", paste(lowest_assign_taxon, "sp."),Species)) %>% 
  dplyr::select(-"lowest_assign_taxon")


### Create phyloseq object ###

cat("Creating phyloseq object\n")

ps <- phyloseq(otu_table(ASV_samples_table_noChim2, taxa_are_rows=F), 
               sample_data(meta),
               tax_table(as.matrix(ASV_taxonomy3)))


### Remove non-bacterial ASVs if present ###

cat("Filtering out non-bacterial, chloroplast and mitochondrial ASVs\n")

#check for non-bacterial ASVs
cat(paste0("Represented domains: ", paste0(unique(ASV_taxonomy3$Kingdom), collapse = ", "), "\n"))
cat(paste0("Represented classes: ", paste0(sort(unique(ASV_taxonomy3$Class)), collapse = ", "), "\n"))

# remove non-bacterial ASVs
ps <- subset_taxa(ps, Kingdom == "d__Bacteria" | is.na(Kingdom))
ps <- subset_taxa(ps, Class != "c__Chloroplast" | is.na(Class))
ps <- subset_taxa(ps, Family != "f__Mitochondria" | is.na(Family))

cat(paste0("Found and removed ", length(which(ASV_taxonomy[,3] == "c__Chloroplast")), " chloroplast ASV(s)\n"))
cat(paste0("Found and removed ", length(which(ASV_taxonomy[,5] == "f__Mitochondria")), " mitochondrial ASV(s)\n"))
cat(paste0("Found and removed ", length(which(ASV_taxonomy[,1] != "d__Bacteria")), " other non-bacterial ASV(s)\n"))

# save sequences and give new names to ASVs
dna <- DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

names(dna) <- taxa_names(ps)

# write as fasta
writeXStringSet(dna, file.path(out.tax, "ASVs_dada2.fasta"))
# write as df with taxonomy
dna_df <- data.frame(
  seq = dna
) %>% 
  merge(tax_table(ps), by=0) %>% 
  dplyr::rename(ASV = Row.names)
write.table(dna_df, file.path(out.tax, "ASV_name_tax.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

saveRDS(object = ps, file = file.path(out.tax, "phyloseq_object.RDS"))
# ps <- readRDS(file.path(out.tax, "phyloseq_object.RDS"))

cat("Finished creating and filtering phyloseq object\n")

### Create a single abundance table ###

tab <- otu_table(ps) # samples are rows and ASV are columns
class(tab) <- "matrix" # warning but it's fine
tax <- as.data.frame(tax_table(ps))
tax <- tax %>% 
  mutate(ASV = rownames(tax), .before = "Kingdom")

cols <- colnames(tab)
tab2 <- as.data.frame(tab) %>% 
  mutate(SampleID = rownames(tab)) %>%
  pivot_longer(all_of(cols), names_to = "ASV", values_to = "count") %>% 
  left_join(tax, by = "ASV") %>% 
  left_join(meta, by = "SampleID")

write.table(tab2, file.path(out.tax, "sample_ASV_table_long.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)

### Check for single- and doubletons ###

cat("Checking for singletons and doubletons - note that singletons may be introduced by rarefaction\n")

tot <- sort(apply(X=tab, MARGIN=2, sum)) # total abundance of each ASV
# show ASVs with total abundance of 1 or 2
sdtons <- names(tot[tot < 3])
cat("The following ASVs have a total abundance of 1 or 2:\n")
cat(tax_table(ps)[sdtons, c("Species")])
# show only samples in which they are present
extr <- tab[,sdtons]
if (is.null(dim(extr))){
  extr_sum <- extr
  cat("\nSample(s) in which this ASV was found:\n")
  cat(names(extr_sum[extr_sum != 0]))
} else {
  extr_sum <- apply(extr, 1, sum) 
  cat("\nSamples in which these ASVs were found:\n")
  print(tab[names(extr_sum[extr_sum != 0]),sdtons])
}

# shows ASVs with low prevalence
tab_bin <- tab
tab_bin[tab_bin > 0] <- 1
occ <- sort(apply(X=tab_bin, MARGIN=2, sum))
low_occ <- names(occ[occ == 1])
cat("\nThe following ASVs appear in only 1 sample:\n")
cat(tax_table(ps)[low_occ, c("Species")])

### Plotting most abundant taxa ###

cat("\nPlotting most abundant taxa\n")

top_nested <- nested_top_taxa(
  ps,
  top_tax_level = "Genus", # most abundant order
  nested_tax_level = "Species", # most abundant genera within those orders
  n_top_taxa = 8, # top 8 most abundant genera
  n_nested_taxa = 2 # top 2 most abundant species
)

tax_plot <- plot_nested_bar(ps_obj = top_nested$ps_obj, 
                            top_level = "Genus", nested_level = "Species",
                            legend_title = "Taxonomy") + 
  theme_nested(theme_classic) + ylab("Relative abundance") + 
  theme(
    #strip.text.x = element_text(angle=90),
    strip.background=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.size = unit(0.4, "cm")
  ) + 
  guides(fill=guide_legend(ncol =1))

if (facet_var[1] != ""){
  tax_plot <- tax_plot + facet_wrap(formula(paste0("~ ", paste(facet_var, collapse = "+"))), drop=T, scales = "free")
}

ggsave(file.path(out.plots, "05_taxonomy_barplot.pdf"), tax_plot, device="pdf", width = 8, height = 6)

cat("Taxonomy assignment done\n")
