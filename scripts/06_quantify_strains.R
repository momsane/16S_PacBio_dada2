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

if (length(args) != 8){
  stop(" Usage: 06_quantify_strains.R <phyloseq> <cd-hit_clusters_tax_full> <qpcr> <abundance_col> <facet_var> <max_cells_raref> <quant_results_dir> <plots_dir>", call.=FALSE)
} else {
  input.ps <- args[1] # phyloseq object resulting from 05_assign_taxonomy
  input.clusters <- args[2] # table of ASVs with their assigned cd-hit cluster and user-input taxonomy
  input.qpcr <- args[3] # qpcr data, if available
  abundance_col <- args[4] # name of the column in qpcr containing the abundance in the sample
  facet_var <- args[5] # one column in the metadata table to facet the taxonomy plot, put "" if not needed
  maxraref <- args[6] # maximum number of 'cells' to extrapolate rarefaction curves
  out.quant <- args[7] # folder to write quantification results
  out.plots <- args[8] # folder to write plots
}

# root <- "/Volumes/RECHERCHE/FAC/FBM/DMF/pengel/general_data/D2c/mgarcia/20240708_mgarcia_syncom_assembly/pacbio_analysis/run1_MD_bees"
# input.ps <- file.path(root, "results", "assign_taxonomy", "phyloseq_object_filtered_nonrarefied.RDS")
# input.clusters <- file.path(root, "workflow", "config", "all_16S_cd-hit_clusters_tax_full.tsv")
# input.qpcr <- "/Volumes/RECHERCHE/FAC/FBM/DMF/pengel/general_data/D2c/mgarcia/20240708_mgarcia_syncom_assembly/absolute_quantification_results/qPCR_results_analyzed.tsv"
# abundance_col <- "normalized_16S_copies_gut"
# facet_var <- "SampleType"
# maxraref <- -1
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

maxraref <- as.numeric(maxraref)

### Get tables ###

tab <- otu_table(ps, taxa_are_rows=F) %>% as("matrix") # samples are rows and ASV are columns

## remove samples with zero reads
zeroes <- names(which(apply(tab,1,sum)==0))
# ps <- subset_samples(ps, !(SampleID %in% zeroes))

tax <- as.data.frame(tax_table(ps))
tax <- tax %>% 
  mutate(ASV = rownames(tax), .before = "Kingdom")
cols <- c("Genus", "Species", "Strain")
tax_unique <- unique(tax[ ,cols] %>% filter(!is.na(Strain))) %>% arrange(Strain)

meta <- sample_data(ps)

# Add qPCR data if available

# function to convert to df to matrix
df_to_matrix <- function(df = data.frame(), rownames_col = character(), names_col = character(), values_col = character()){ 
  # pivot wider
  df2 <- df[,c(rownames_col, names_col, values_col)] %>% 
    pivot_wider(names_from = !!sym(names_col), values_from = !!sym(values_col)) %>% 
    as.data.frame()
  # set rownames
  rownames(df2) <- df2 %>% select(!!sym(rownames_col)) %>% unlist() %>% unname()
  mat <- as.matrix(df2 %>% select(-!!sym(rownames_col)))
  # convert to numeric
  mat2 <- mat
  mat2 <- apply(mat2,2,as.numeric)
  rownames(mat2) <- rownames(mat)
  return(mat2)
}

if ((input.qpcr != "") & (abundance_col != "")){
  cat("Adding qPCR data\n")
  qpcr <- read.table(input.qpcr, sep = "\t", header = T)
  # compute ASV absolute abundance
  long <- psmelt(ps) %>% 
    select(-"Sample") %>% 
    dplyr::rename(
      ASV = OTU,
      read_count = Abundance
    ) %>%
    arrange(SampleID, ASV, Cluster) %>% 
    group_by(SampleID) %>% 
    mutate(
      rel_abun = read_count/sum(read_count),
      rel_abun = if_else(!is.na(rel_abun),rel_abun,0), # replace NA (samples with zero reads)
      .after="read_count"
    ) %>% 
    left_join(qpcr[ ,c("SampleID", abundance_col)], by = "SampleID") %>% 
    relocate(all_of(abundance_col), .after="rel_abun") %>%
    mutate(ASV_abs_abun = round(rel_abun*.data[[abundance_col]]), .after="rel_abun")
  
  # extract samples with no absolute abundance
  samples_noqpcr <- unique(long$SampleID[is.na(long$normalized_16S_copies_gut)])
  
  if (length(samples_noqpcr) != 0){
    cat("The samples below do not have any qPCR data. Strain counts will be computed separately and based on read counts.\n")
    cat(paste0(samples_noqpcr))
    cat("\n")
    
    # re-shape data as matrix like tab
    tab2 <- df_to_matrix(
      long %>% filter(!(SampleID %in% samples_noqpcr)),
      rownames_col = "SampleID",
      names_col = "ASV",
      values_col = "ASV_abs_abun"
    )
    
    tab3 <- df_to_matrix(
      long %>% filter(SampleID %in% samples_noqpcr),
      rownames_col = "SampleID",
      names_col = "ASV",
      values_col = "read_count"
    )
    
  } else {
    
    tab2 <- df_to_matrix(
      long,
      rownames_col = "SampleID",
      names_col = "ASV",
      values_col = "ASV_abs_abun"
    )
    
  }
  
  ## save long table
  write.table(
    long,
    file.path(out.quant, "sample_ASV_table_long_qpcr.tsv"),
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F
  )
  
  ## create new phyloseq object using abs abun
  ps_abs <- phyloseq(
    otu_table(tab2, taxa_are_rows = F),
    tax_table(ps),
    sample_data(ps)
  )
  saveRDS(ps_abs, file.path(out.quant, "phyloseq_object_filtered_absabun.RDS"))
  
} else {
  tab2 <- tab
}

# Quantification of known strains

cat("Inferring strain abundance\n")
cat("Note: only ASVs matching the provided custom database will be used to infer strain abundance.\n")

# A: ASV by sample matrix (r,c) = abundance of each ASV in each sample
# B: ASV by strain matrix (r,c) = number of copies of each ASV in each strain
# C: strain by sample matrix (r,c) = abundance of each strain in each sample <- this is what we want!
# We know A and B, and that A=B.C

# make matrix B

## remove ASV clusters that the user does not want to use for quantification
cls_use <- clusters$cluster[clusters$use_to_quantify == TRUE]
## keep only ASVs matching a syncom ASV and used for quantification
asv_cls <- sort(unique(tax$ASV[tax$inferred_from == "addSpecies_custom" & tax$Cluster %in% cls_use]))

## get number of ASV copy per strain
clusters2 <- clusters %>% 
  filter(cluster %in% cls_use) %>%
  group_by(cluster, strain) %>% 
  summarize(n = n()) %>% 
  pivot_wider(names_from = strain, values_from = n)
## convert to matrix = B
clusters3 <- as.matrix(clusters2)
rownames(clusters3) <- clusters3[ ,1]
clusters3 <- clusters3[ ,-1]
clusters3[is.na(clusters3)] <- 0
## convert to numeric
clusters4 <- apply(clusters3,2,as.numeric)
rownames(clusters4) <- rownames(clusters3)

# make matrix A

make_matA <- function(matrix){
  matA <- t(matrix[ ,asv_cls])
  ## rename ASVs by cluster instead
  rownames(matA) <- tax$Cluster[match(rownames(matA), tax$ASV)]
  return(matA)
}

tab4 <- make_matA(tab2)

if (length(samples_noqpcr) != 0){
  tab5 <- make_matA(tab3)
}

# solve the strain by sample matrix = C

solve_matC <- function(matA, matB, name){
  ## make sure ASVs are exactly the same and in order in both matrices
  cat("Checking that matrices are ordered in the same way: ")
  matB <- matB[rownames(matA), ] # reorder manually
  match <- rownames(matB) == rownames(matA)
  cat(all(match[match==FALSE])) # should be TRUE
  cat("\n")
  
  # remove strains with no ASVs
  asv_presence <- apply(matB,2,sum)
  strains_keep <- names(asv_presence[asv_presence != 0])
  matB <- matB[, strains_keep]
  
  set.seed(42)
  matC = round(qr.solve(matB,matA))
  ## compute relative abundance
  matC_rel <- 100*scale(matC, center = FALSE, scale = colSums(matC))
  
  ## save matrices
  write.table(matC, file.path(out.quant, paste0("strain_quant_matrix_count_", name, ".tsv")), sep = "\t", quote = F, col.names = T, row.names = T)
  write.table(matC_rel, file.path(out.quant, paste0("strain_quant_matrix_relative_", name, ".tsv")), sep = "\t", quote = F, col.names = T, row.names = T)
  
  ## make long table
  df <- as.data.frame(matC) %>% 
    mutate(Strain = rownames(matC), .before = 1)
  
  return(list(matC, df))
}

# solve and reformat table

if ((input.qpcr != "") & (abundance_col != "")){
  
  # solve for samples with absolute abundance
  
  solved_qpcr <- solve_matC(matA = tab4, matB = clusters4, "qpcr")
  
  C <- solved_qpcr[[1]]
  df <- solved_qpcr[[2]]
  
  df <- df %>% 
    pivot_longer(colnames(C), names_to = "SampleID", values_to = "abs_abun") %>% 
    group_by(SampleID) %>% 
    mutate(
      rel_abun = 100*abs_abun/sum(abs_abun),
      rel_abun = if_else(!is.na(rel_abun),rel_abun,0),
      .after="abs_abun"
    ) %>% 
    left_join(tax_unique, by = "Strain") %>%
    left_join(meta, by = "SampleID") %>%
    arrange(Species)
  
  write.table(df, file.path(out.quant, "strain_quant_long_table_qpcr.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
  
  # solve for samples without absolute abundance
  
  if (length(samples_noqpcr) != 0){
    solved_noqpcr <- solve_matC(matA = tab5, matB = clusters4, "samples_noqpcr")
    C_noqpcr <- solved_noqpcr[[1]]
    df_noqpcr <- solved_noqpcr[[2]]
    
    df_noqpcr <- df_noqpcr %>% 
      pivot_longer(colnames(C_noqpcr), names_to = "SampleID", values_to = "count") %>% 
      group_by(SampleID) %>% 
      mutate(
        rel_abun = 100*count/sum(count),
        rel_abun = if_else(!is.na(rel_abun),rel_abun,0),
        .after="count"
      ) %>% 
      left_join(tax_unique, by = "Strain") %>%
      left_join(meta, by = "SampleID") %>%
      arrange(Species)
    
    write.table(df_noqpcr, file.path(out.quant, "strain_quant_long_table_samples_noqpcr.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
  }
  
} else {
  
  solved_noqpcr <- solve_matC(matA = tab4, matB = clusters4, "noqpcr")
  
  C <- solved_noqpcr[[1]]
  df <- solved_noqpcr[[2]]
  
  df <- df %>% 
    pivot_longer(colnames(C), names_to = "SampleID", values_to = "count") %>% 
    group_by(SampleID) %>% 
    mutate(rel_abun = 100*count/sum(count)) %>% 
    left_join(tax_unique, by = "Strain") %>%
    left_join(meta, by = "SampleID") %>%
    arrange(Species)
  
  write.table(df, file.path(out.quant, "strain_quant_long_table_noqpcr.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
}


cat("Generating plots\n")

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

pal <- generate_palette(unique(df[ ,c("Genus", "label")])) %>% 
  left_join(unique(df[ ,c("Strain", "label")]), by = "label")

# save palette
write.table(
  pal,
  file.path(out.quant, "color_palette.tsv"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

if ((input.qpcr != "") & (abundance_col != "")){
  logmax <- ceiling(max(log10(df$abs_abun)))-1
  # absolute abundance if qpcr provided
  p_abs <- ggplot(
    df %>% filter(abs_abun!=0),
    aes(
      x = SampleID,
      y = abs_abun/10^logmax,
      fill = label
    )
  ) +
    geom_col() +
    labs(
      y = bquote("Genome equivalents/gut (" * 10^.(logmax) * ")"),
      x = "",
      fill = ""
    ) +
    scale_fill_manual(values = setNames(pal$color,pal$label)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    theme(
      strip.background=element_blank(),
      strip.text.x = element_text(angle=90, hjust=0, vjust=0.5),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.position = "right",
      legend.key.size = unit(0.4, 'cm'),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) + 
    guides(fill=guide_legend(ncol =1))
  
  if (facet_var[1] != ""){
    p_abs <- p_abs + facet_wrap(formula(paste0("~ ", paste(facet_var, collapse = "+"))), drop=T, scales = "free_x", space = "free_x")
  }
  ggsave(file.path(out.plots, "06_strain_barplot_absolute.pdf"), p_abs, device="pdf", width = 8, height = 6)
}

# relative abundance
if (length(samples_noqpcr) != 0){
  df_noqpcr <- df_noqpcr %>% 
    mutate(label = paste(substr(Genus, 1, 1),". ", word(Species, 2), " ", Strain, sep = ""))
  df_rel <- rbind(df, df_noqpcr)
} else {
  df_rel <- df 
}

p_rel <- ggplot(
  df_rel %>% filter(rel_abun!=0),
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
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(
    strip.background=element_blank(),
    strip.text.x = element_text(angle=90, hjust=0, vjust=0.5),
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
  guides(fill=guide_legend(ncol=1))
# p_rel

if (facet_var[1] != ""){
  p_rel <- p_rel + facet_wrap(formula(paste0("~ ", paste(facet_var, collapse = "+"))), drop=T, scales = "free_x", space = "free_x")
}

ggsave(file.path(out.plots, "06_strain_barplot_relative.pdf"), p_rel, device="pdf", width = 8, height = 6)

### Rarefaction curves ###

if (maxraref <= 0){
  cat("Skipping rarefaction curves\n")
} else {
  cat("Generating rarefaction curves\n")
  
  # filter out samples with total abundance of 0
  strain_ab <- apply(C, 2, sum)
  samples_keep <- names(which(strain_ab != 0))
  C2 <- C[,samples_keep]
  
  # using iNEXT so I can also estimate the sampling coverage
  dt <- iNEXT(
    C2, # samples must be as columns
    q = c(0,1),
    datatype = "abundance",
    endpoint = maxraref,
    knots = 50,
    se = TRUE,
    conf = 0.95,
    nboot = 20
  )
  
  inextqd <- dt$iNextEst$size_based %>%
    dplyr::rename(SampleID = Assemblage)
  
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
    ) +
    facet_wrap( ~ Order.q, scales = "free_y")
  
  ggsave(file.path(out.plots, "06_rarefaction_curves_strains.pdf"), qd.plot, device="pdf", width = 12, height = 8)
  
  write.table(
    inextqd,
    file = file.path(out.quant, "inext_data.tsv"),
    sep = "\t",
    quote = F,
    row.names = F,
    col.names = T
  )
}

cat("Strain quantification done\n")