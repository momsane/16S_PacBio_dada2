# Load required libraries

if(!require(dplyr)){
  install.packages(pkgs = 'dplyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(dplyr)
}

if(!require(tidyr)){
  install.packages(pkgs = 'tidyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(tidyr)
}

if(!require(ggplot2)){
  install.packages(pkgs = 'ggplot2', repos = 'https://stat.ethz.ch/CRAN/')
  library(ggplot2)
}

if(!require(phyloseq)){
  if(!requireNamespace("BiocManager")){
    install.packages("BiocManager")
  }
  BiocManager::install("phyloseq")
  library(phyloseq)
}

# Collect arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6){
  stop(" Usage: 06_quantify_strains.R <cd-hit_clusters_tax_full> <metadata_table.tsv> <facet_var> <quant_results_dir> <plots_dir>", call.=FALSE)
} else {
  input.ps <- args[1] # phyloseq object resulting from 05_assign_taxonomy
  input.clusters <- args[2] # table of ASVs with their assigned cd-hit cluster and user-input taxonomy
  input.metadata <- args[3] # sample metadata table, tab-separated, first column is the the sample name
  facet_var <- args[4] # one column in the metadata table to facet the taxonomy plot, put "" if not needed
  out.quant <- args[5] # folder to write quantification results
  out.plots <- args[6] # folder to write plots
}

# root <- "/Volumes/D2c/mgarcia/20240708_mgarcia_syncom_invivo/exp01_inoculation_methods/pacbio_analysis"
# # root <- "/work/FAC/FBM/DMF/pengel/general_data/syncom_pacbio_analysis"
# input.ps <- file.path(root, "results", "assign_taxonomy", "phyloseq_object.RDS")
# input.clusters <- file.path(root, "workflow", "config", "all_16S_cd-hit_clusters_tax_full.tsv")
# input.metadata <- file.path(root, "workflow", "config", "metadata.tsv")
# facet_var <- "SampleType"
# out.quant <- file.path(root, "results", "quantify_strains")
# out.plots <- file.path(root, "plots")

### Create outdirs ###

cat("\nCreating directories\n")

dir.create(out.quant, recursive = TRUE, showWarnings = FALSE)
dir.create(out.plots, recursive = TRUE, showWarnings = FALSE)

### Inputs ###

cat("Reading inputs and extracting info\n")

ps <- readRDS(input.ps)
clusters <- read.table(input.clusters, sep = "\t", header = T)
meta <- read.table(input.metadata, sep = "\t", header = T)

### Get taxonomy table and OTU table ###

tab <- otu_table(ps) # samples are rows and ASV are columns
class(tab) <- "matrix" # warning but it's fine
tax <- as.data.frame(tax_table(ps))
tax <- tax %>% 
  mutate(ASV = rownames(tax), .before = "Kingdom")

# Quantification of known strains

cat("Inferring strain abundance\n")
cat("Note: only ASVs matching the provided custom database will be used to infer strain abundance.\n")

# A: ASV by sample matrix (r,c) = abundance of each ASV in each sample
# B: ASV by strain matrix (r,c) = number of copies of each ASV in each strain
# C: strain by sample matrix (r,c) = abundance of each strain in each sample <- this is what we want!
# We know A and B, and that A=B.C

# keep only ASVs matching a syncom ASV
cls <- sort(unique(tax$ASV[tax$inferred_from == "addSpecies_custom"]))
tab2 <- t(tab[ ,cls]) # = matrix A

# rename ASVs by cluster instead
rownames(tab2) <- tax$Cluster[match(rownames(tab2), tax$ASV)]

# get number of ASV copy per strain
clusters2 <- clusters %>% 
  filter(cluster %in% rownames(tab2)) %>%
  group_by(cluster, strain) %>% 
  summarize(n = n()) %>% 
  pivot_wider(names_from = strain, values_from = n)
# convert to matrix = B
clusters3 <- as.matrix(clusters2)
rownames(clusters3) <- clusters3[ ,1]
clusters3 <- clusters3[ ,-1]
clusters3[is.na(clusters3)] <- 0

# make sure ASV are ordered the same way in both matrices
length(symdiff(rownames(tab2), rownames(clusters3))) # should be 0
clusters3 <- clusters3[rownames(tab2), ]
match <- rownames(clusters3) == rownames(tab2)
length(match[match==FALSE]) # should be 0

# solve the strain by sample matrix 
C = round(qr.solve(clusters3,tab2))

# compute relative abundance
C_rel <- 100*scale(C, center = FALSE, scale = colSums(C))

# save matrices
write.table(C, file.path(out.quant, "strain_quantification_matrix_count.tsv"), sep = "\t", quote = F, col.names = T, row.names = T)
write.table(C_rel, file.path(out.quant, "strain_quantification_matrix_relative.tsv"), sep = "\t", quote = F, col.names = T, row.names = T)

# save as long table
df <- as.data.frame(C) %>% 
  mutate(Strain = rownames(C), .before = 1)

cols <- c("Genus", "Species", "Strain")
tax_unique <- unique(tax[ ,cols] %>% filter(!is.na(Strain))) %>% arrange(Strain)

df <- df %>% 
  pivot_longer(colnames(C), names_to = "SampleID", values_to = "count") %>% 
  group_by(SampleID) %>% 
  mutate(rel_abun = 100*count/sum(count)) %>% 
  left_join(tax_unique, by = "Strain") %>%
  arrange(Species)

write.table(df, file.path(out.quant, "strain_quantification_long_table.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)

cat("Generating plots\n")

# generate barplot
df <- df %>%
  left_join(meta, by = "SampleID")
  
# order strains by their species
df$Strain <- factor(df$Strain, levels = unique(df$Strain), ordered = T)

p1 <- ggplot(
  df,
  aes(
    x = SampleID,
    y = count,
    fill = Strain
  )
) +
  geom_col() +
  labs(
    y = "Count"
  ) +
  facet_wrap(vars(facet_var)) +
  theme_bw() +
  theme(
    #strip.text.x = element_text(angle=90),
    strip.background=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  guides(fill=guide_legend(ncol =1))

p2 <- ggplot(
  df,
  aes(
    x = SampleID,
    y = rel_abun,
    fill = Strain
  )
) +
  geom_col() +
  labs(
    y = "Relative abundance (%)"
  ) +
  facet_wrap(vars(facet_var)) +
  theme_bw() +
  theme(
    #strip.text.x = element_text(angle=90),
    strip.background=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  guides(fill=guide_legend(ncol =1))

if (facet_var[1] != ""){
  p1 <- p1 + facet_wrap(formula(paste0("~ ", paste(facet_var, collapse = "+"))), drop=T, scales = "free_x")
  p2 <- p2 + facet_wrap(formula(paste0("~ ", paste(facet_var, collapse = "+"))), drop=T, scales = "free_x")
}

ggsave(file.path(out.plots, "06_strain_barplot_count.pdf"), p1, device="pdf", width = 8, height = 6)
ggsave(file.path(out.plots, "06_strain_barplot_relative.pdf"), p2, device="pdf", width = 8, height = 6)


cat("Strain quantification done\n")