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

if(!require(rlang)){
  install.packages(pkgs = 'rlang', repos = 'https://stat.ethz.ch/CRAN/')
  library(rlang)
}

if(!require(phyloseq)){
  if(!requireNamespace("BiocManager")){
    install.packages("BiocManager")
  }
  BiocManager::install("phyloseq")
  library(phyloseq)
}

if(!require(RColorBrewer)){
  install.packages(pkgs = 'RColorBrewer', repos = 'https://stat.ethz.ch/CRAN/')
  library(RColorBrewer)
}

if(!require(grDevices)){
  install.packages(pkgs = 'grDevices', repos = 'https://stat.ethz.ch/CRAN/')
  library(grDevices)
}

if(!require(iNEXT)){
  install.packages(pkgs = 'iNEXT', repos = 'https://stat.ethz.ch/CRAN/')
  library(iNEXT)
}

# Collect arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7){
  stop(" Usage: 06_quantify_strains.R <cd-hit_clusters_tax_full> <metadata_table.tsv> <facet_var> <max_cells_raref> <quant_results_dir> <plots_dir>", call.=FALSE)
} else {
  input.ps <- args[1] # phyloseq object resulting from 05_assign_taxonomy
  input.clusters <- args[2] # table of ASVs with their assigned cd-hit cluster and user-input taxonomy
  input.metadata <- args[3] # sample metadata table, tab-separated, first column is the the sample name
  facet_var <- args[4] # one column in the metadata table to facet the taxonomy plot, put "" if not needed
  maxraref <- args[5] # maximum number of 'cells' to extrapolate rarefaction curves
  out.quant <- args[6] # folder to write quantification results
  out.plots <- args[7] # folder to write plots
}

# root <- "/Volumes/D2c/mgarcia/20240708_mgarcia_syncom_assembly/pacbio_analysis/run1_bees"
# input.ps <- file.path(root, "results", "assign_taxonomy", "phyloseq_object_filtered.RDS")
# input.clusters <- file.path(root, "workflow", "config", "all_16S_cd-hit_clusters_tax_full.tsv")
# input.metadata <- file.path(root, "workflow", "config", "metadata.tsv")
# facet_var <- "SampleType"
# maxraref <- 2500
# out.quant <- file.path(root, "results", "quantify_strains")
# out.plots <- file.path(root, "plots")

if (facet_var %in% c("Kingdom", "Phylum", "Class", "Family", "Order", "Genus", "Species")){
  cat("Error: the provided facet_var is conflicting with taxonomic rank names! Please change the name of this variable before continuing.\n")
  quit(save="no")
}

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

# remove ASV clusters that the user does not want to use for quantification
cls_use <- clusters$cluster[clusters$use_to_quantify == TRUE]

# keep only ASVs matching a syncom ASV and used for quantification
cls <- sort(unique(tax$ASV[tax$inferred_from == "addSpecies_custom" & tax$Cluster %in% cls_use]))
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

set.seed(42)
# solve the strain by sample matrix 
C = round(qr.solve(clusters3,tab2))

# compute relative abundance (there will be NA values for blanks or MD bees)
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

# generate barplots
df <- df %>%
  left_join(meta, by = "SampleID")

# custom label with species and strain name
df <- df %>% 
  mutate(label = paste(substr(Genus, 1, 1),". ", word(Species, 2), " ", Strain, sep = ""))

# order strains by their species
df$Strain <- factor(df$Strain, levels = unique(df$Strain), ordered = T)
df$label <- factor(df$label, levels = unique(df$label), ordered = T)

# palette
generate_palette <- function(df) {
  
  # get genera
  genera <- unique(df[["Genus"]])
  
  # assign color to genera
  n_colors <- length(genera)
  if (n_colors <= 8) {
    base_colors <- brewer.pal(n_colors, "Dark2")
  } else {
    base_colors <- rainbow(n_colors)
  }
  
  base_colors <- setNames(base_colors, genera)
  
  # Prepare result vector
  colors <- character(nrow(df))
  
  # Loop over each genus and assign alpha variations per species
  for (g in genera) {
    species_in_genus <- df$label[df$Genus == g]
    n_species <- length(species_in_genus)
    
    # Generate transparencies between 0.4 and 1 (evenly spaced)
    alphas <- seq(0.4, 1, length.out = n_species)
    
    # Make transparent versions of the base color
    genus_colors <- scales::alpha(base_colors[g], alphas)
    
    # Assign these transparencies to the corresponding rows
    for (i in seq_along(species_in_genus)) {
      rows <- which(df[["Genus"]] == g & df[["label"]] == species_in_genus[i])
      colors[rows] <- genus_colors[i]
    }
  }
  
  df$color <- colors
  return(df)
}


pal <- generate_palette(unique(df[ ,c("Genus", "label")]))

p1 <- ggplot(
  df,
  aes(
    x = SampleID,
    y = count,
    fill = label
  )
) +
  geom_col() +
  labs(
    y = "Count",
    fill = ""
  ) +
  scale_fill_manual(values = setNames(pal$color,pal$label)) +
  theme_bw() +
  theme(
    strip.background=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    # legend
    legend.position = "right",
    legend.key.size = unit(0.4, 'cm'),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  guides(fill=guide_legend(ncol =1))

p2 <- ggplot(
  df,
  aes(
    x = SampleID,
    y = rel_abun,
    fill = label
  )
) +
  geom_col() +
  labs(
    y = "Relative abundance (%)",
    fill = ""
  ) +
  scale_fill_manual(values = setNames(pal$color,pal$label)) +
  theme_bw() +
  theme(
    strip.background=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    # legend
    legend.position = "right",
    legend.key.size = unit(0.4, 'cm'),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
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

### Rarefaction curves ###

cat("Generating rarefaction curves\n")

# filter out samples with total abundance of 0
strain_ab <- apply(C, 2, sum)
samples_keep <- names(which(strain_ab != 0))
C2 <- C[,samples_keep]

# using iNEXT so I can also estimate the sampling coverage
dt <- iNEXT(
  C2, # samples must be as columns
  q = 0,
  datatype = "abundance",
  endpoint = maxraref,
  knots = 50,
  se = TRUE,
  conf = 0.95,
  nboot = 20
)

inextqd <- dt$iNextEst$size_based %>%
  dplyr::rename(SampleID = Assemblage)

sc.plot <- ggplot(
  inextqd[inextqd$Method != "Extrapolation", ],
  aes(
    x = m,
    y = SC,
    group = SampleID
  )
) +
  geom_vline(aes(xintercept = min(inextqd$m[inextqd$Method == "Observed"]), color = "low"), linetype = "dashed") + # sample with lowest number of reads
  geom_vline(aes(xintercept = max(inextqd$m[inextqd$Method == "Observed"]), color = "high"), linetype = "dashed") + # sample with highest number of reads
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey80") +
  scale_color_manual(name = "", values = c(low = "#669bbc", high = "#e76f51"), labels = c(low = "Lowest depth", high = "Highest depth")) +
  geom_line(alpha = 0.4) +
  geom_ribbon(
    aes(
      ymin = SC.LCL,
      ymax = SC.UCL
    ), alpha = 0.2
  ) +
  theme_bw() +
  labs(
    x = "# of cells",
    y = "Sampling coverage (iNEXT)"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8,0.1),
    legend.background = element_rect(fill=alpha('white', 0.4))
  )

qd.plot <- ggplot(
  inextqd[inextqd$Method != "Extrapolation", ],
  aes(
    x = m,
    y = qD,
    group = SampleID
  )
) +
  geom_vline(aes(xintercept = min(inextqd$m[inextqd$Method == "Observed"]), color = "low"), linetype = "dashed") + # sample with lowest number of reads
  geom_vline(aes(xintercept = max(inextqd$m[inextqd$Method == "Observed"]), color = "high"), linetype = "dashed") + # sample with highest number of reads
  scale_color_manual(name = "", values = c(low = "#669bbc", high = "#e76f51"), labels = c(low = "Lowest depth", high = "Highest depth")) +
  geom_line(alpha = 0.4) +
  geom_ribbon(
    aes(
      ymin = qD.LCL,
      ymax = qD.UCL
    ), alpha = 0.1
  ) +
  geom_label(
    data = inextqd[inextqd$Method == "Observed", ],
    aes(
      x = m,
      y = qD,
      label = SampleID
    ), size = 3, nudge_x = 70
  ) +
  scale_y_continuous(breaks = seq(0,max(inextqd$qD[inextqd$Method != "Extrapolation"])+5,5)) +
  theme_bw() +
  labs(
    x = "# of cells",
    y = "# of strains"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1,0.9),
    legend.background = element_rect(fill=alpha('white', 0.4))
  )

ggsave(file.path(out.plots, "06_sampling_coverage_strains.pdf"), sc.plot, device="pdf", width = 8, height = 6)
ggsave(file.path(out.plots, "06_rarefaction_curves_strains.pdf"), qd.plot, device="pdf", width = 10, height = 8)

write.table(
  inextqd,
  file = file.path(out.quant, "inext_data.tsv"),
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

cat("Strain quantification done\n")