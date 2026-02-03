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

if(!require(microDecon)){
  if(!require(devtools)){
    install.packages(pkgs = 'devtools', repos = 'https://stat.ethz.ch/CRAN/')
  }
  devtools::install_github("donaldtmcknight/microDecon") # installing through GitHub to get latest updates
  library(microDecon)
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

if (length(args) != 6){
  stop(" Usage: decontaminate.R <phyloseq> <group_var> <blank_name> <rarefy_to> <results_dir> <plots_dir>", call.=FALSE)
} else {
  input.ps <- args[1] # phyloseq object resulting from 05_assign_taxonomy
  group_var <- args[2] # one column in the metadata table indicating sample grouping, including blanks
  blank_name <- args[3] # how blanks are named in the metadata table
  rarefy_to <- args[4] # number of reads to rarefy to; if <=0, no rarefaction
  out.decont <- args[5] # folder to write decontamination results
  out.plots <- args[6] # folder to write plots
}

root <- "/Volumes/RECHERCHE/FAC/FBM/DMF/pengel/general_data/D2c/mgarcia/20240708_mgarcia_syncom_assembly/pacbio_analysis/run3_colonized_bees"
input.ps <- file.path(root, "results", "assign_taxonomy", "phyloseq_object_filtered.RDS")
group_var <- "SampleType"
blank_name <- "blank"
rarefy_to <- 18000
out.decont <- file.path(root, "results", "decontaminate")
out.plots <- file.path(root, "plots")

### Create outdirs ###

cat("\nCreating directories\n")

dir.create(out.decont, recursive = TRUE, showWarnings = FALSE)
dir.create(out.plots, recursive = TRUE, showWarnings = FALSE)

### Inputs ###

cat("Reading inputs and extracting info\n")

ps <- readRDS(input.ps)
rarefy_to <- as.numeric(rarefy_to)

### Get tables ###

tab <- otu_table(ps,taxa_are_rows=F) %>% as("matrix") # samples are rows and ASV are columns
meta <- sample_data(ps) %>% as("data.frame")
tax <- tax_table(ps) %>% as("matrix")
abundance_table_long <- psmelt(ps) %>% 
  select(-"Sample") %>% 
  dplyr::rename(
    ASV = OTU,
    read_count = Abundance
  ) %>%
  group_by(SampleID) %>% 
  mutate(rel_abun_original = 100*read_count/sum(read_count)) %>%
  ungroup() %>%
  arrange(SampleID, ASV)

### Format data for microDecon ### 

cat("Preparing data for microDecon\n")

# goal: a df of ASVs (rows) by samples (columns), with the first column as ASV ID, then blank samples, then other samples grouped
df <- as.data.frame(t(tab)) %>% 
  mutate(ASV=colnames(tab), .before=1)

# check that no sample has only zeroes
total_reads <- apply(tab,1,sum)
if (length(total_reads[total_reads == 0]) != 0){
  cat(paste0("Warning: ", length(total_reads[total_reads == 0]), " samples have zero reads and will be removed from the analysis\n"))
  samples_zero <- names(total_reads[total_reads == 0])
  df <- df %>%
    select(-samples_zero)
}

group_column <- meta[[group_var]]
groups <- unique(group_column)
groups_ordered <- c(blank_name, groups[groups != blank_name])

# order the metadata table
meta_ordered <- meta %>%
  mutate(
    !!group_var := factor(.data[[group_var]], levels = groups_ordered)
  ) %>%
  arrange(.data[[group_var]])

# order the ASV table
sample_order <- meta_ordered$SampleID
df_ordered <- df %>%
  select(ASV, all_of(sample_order))

# get number of samples in each group
group_column <- meta_ordered[[group_var]]
group_counts <- table(group_column[group_column != blank_name])
group_counts <- unlist(group_counts[group_counts != 0])

### Run microDecon ### 

cat("Running microDecon\n")

n_blanks <- length(group_column[group_column == blank_name])

decontam <- decon(
  df_ordered,
  numb.blanks = n_blanks,
  numb.ind = group_counts,
  taxa = F,
  runs = 1,
  thresh = 1, # disable filtering
  prop.thresh = 0 # disable filtering
)

# save object
saveRDS(decontam, file.path(out.decont, "microDecon_object.RDS"))


### Analyze microDecon results ###

cat("Analyzing microDecon results\n")

## ASVs removed

### rename groups
ASV_rm <- decontam$OTUs.removed
old_names <- paste0("Group", seq_along(group_counts))
new_names <- names(group_counts)
### add taxonomy
ASV_rm <- ASV_rm %>%
  dplyr::rename(!!!setNames(old_names, new_names)) %>%
  merge(tax, by=0)
ASV_rm <- ASV_rm[ ,-1]
ASV_rm[ASV_rm == "-"]<-NA

# save table
write.table(ASV_rm, file.path(out.decont, "ASV_removed.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)

## Reads before and after microdecon
tax2 <- cbind(tax, data.frame(ASV=rownames(tax)))
reads_rm <- decontam$reads.removed
reads_rm <- reads_rm %>%
  pivot_longer(colnames(reads_rm)[-1], names_to = "SampleID", values_to = "reads_microdecon") %>% 
  left_join(tax2, by="ASV") %>%
  left_join(abundance_table_long[ ,c("SampleID", "ASV", "read_count", "rel_abun_original")], by = c("SampleID", "ASV")) %>% 
  replace_na(list(reads_microdecon = 0, read_count = 0)) %>%
  mutate(
    perc_reads_removed = if_else((SampleID != "Mean.blank") & (read_count != 0), 100*reads_microdecon/read_count,NA),
    reads_left = read_count-reads_microdecon
    )

### plot reads removed
plot_rm <- ggplot(
  reads_rm %>% filter(!is.na(perc_reads_removed)),
  aes(
    x = read_count,
    y = perc_reads_removed
  )
) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_continuous(breaks=seq(0,100,20)) +
  labs(
    x = "Initial # of reads",
    y = "% reads removed"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )
plot_rm

# save table
write.table(reads_rm, file.path(out.decont, "reads_removed.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
# save plot
ggsave(file.path(out.plots, "reads_removed_hist.pdf"), plot_rm, device="pdf", width = 6, height = 6)


### Create new phyloseq object with decontaminated table ###

cat("Creating final phyloseq object")

tab2 <- decontam$decon.table[ ,-2] # remove Mean.blank
# add blanks from old phyloseq
tab3 <- cbind(tab2, df_ordered[ ,colnames(df_ordered)[2:(n_blanks+1)]])
rownames(tab3) <- tab3[,1]
tab3 <- tab3[,-1]

ps.decont <- phyloseq(
  otu_table(t(tab3), taxa_are_rows=F),
  sample_data(ps),
  tax_table(ps)
)

abundance_table_long_decontam <- psmelt(ps.decont) %>% 
  select(-"Sample") %>% 
  dplyr::rename(
    ASV = OTU,
    read_count = Abundance
  ) %>%
  arrange(SampleID, ASV)

# save results
saveRDS(ps.decont, file.path(out.decont, "ps_filtered_decontam.RDS"))
write.table(abundance_table_long_decontam, file.path(out.decont, "sample_ASV_table_long_decontam.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)

### Rarefy reads if needed ###

if (rarefy_to <= 0){
  cat("No rarefaction\n")
} else {
  # rarefy reads
  cat(paste0("Rarefaction to ", rarefy_to, " reads; samples with fewer reads will be removed; ASVs with a subsequent total abundance of 0 will be removed\n"))
  ps.final <- rarefy_even_depth(
    ps.decont,
    sample.size = rarefy_to,
    rngseed = 42,
    replace = FALSE,
    trimOTUs = TRUE,
    verbose = FALSE
  )
  # list samples removed
  samples_rm <- setdiff(rownames(otu_table(ps.decont)),rownames(otu_table(ps.final)))
  cat(paste0(length(samples_rm), " sample(s) were removed because they contained fewer than ", rarefy_to, " reads:\n"))
  print(samples_rm)
  # show ASVs removed
  ASVs_rm <- setdiff(colnames(otu_table(ps.decont,taxa_are_rows=T)),colnames(otu_table(ps.final,taxa_are_rows=T)))
  cat(paste0(length(ASVs_rm), " ASV(s) were removed because they were no longer present in any sample:\n"))
  print(tax[ASVs_rm,])
  
  saveRDS(ps.final, file.path(out.decont, "ps_filtered_decontam_rarefied.RDS"))
}

### Make new barplots ###

cat("Generating new barplots")

top_nested <- nested_top_taxa(
  ps.final,
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
  guides(fill=guide_legend(ncol =1)) +
  facet_wrap(group_var, scales = "free")
#tax_plot

ggsave(file.path(out.plots, "taxonomy_barplot_decontam.pdf"), tax_plot, device="pdf", width = 8, height = 6)

cat("Decontamination done\n")
