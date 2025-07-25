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
  if(!require(devtools)){
    install.packages("devtools")
  }
  devtools::install_github("benjjneb/dada2") # installing through GitHub to get latest updates
  library(dada2)
}

if(!require(iNEXT)){
  install.packages(pkgs = 'iNEXT', repos = 'https://stat.ethz.ch/CRAN/')
  library(iNEXT)
}

# Collect arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9){
  stop(" Usage: 04_denoising.R <processed_reads_dir> <read_count_table> <max_reads_derep> <error_model> <max_bases_errormodel> <removeSingletons_T_F> <max_reads_raref> <denoise_results_dir> <plots_dir>", call.=FALSE)
} else {
  input.reads <- args[1] # folder with all pre-processed reads
  input.readcounts <- args[2] # table reporting the number of reads at each step
  maxReads <- args[3] # max number of reads to load at once for dereplication
  errModel <- args[4] # dada2-provided function to estimate the error model
  maxBases <- args[5] # max number of bases to use for error model inference
  removeSingletons <- args[6] # "T" or "F", whether to remove singletons during denoising or not
  maxraref <- args[7] # maximum number of reads to extrapolate rarefaction curves
  out.denois <- args[8] # folder to write denoising results
  out.plots <- args[9] # folder to write plots
}

# input.reads <- "/Volumes/D2c/mgarcia/20240708_mgarcia_syncom_invivo/exp01_inoculation_methods/pacbio_analysis/results/preprocessing/trimmed_filtered_reads"
# input.readcounts <- "/Volumes/D2c/mgarcia/20240708_mgarcia_syncom_invivo/exp01_inoculation_methods/pacbio_analysis/results/preprocessing/read_count_before_after.tsv"
# maxReads <- 1E6
# errModel <- "binnedQualErrfun" # or binnedQualErrfun
# maxBases <- 1E10
# removeSingletons <- "F"
# maxraref <- 5000
# out.denois <- "/Volumes/D2c/mgarcia/20240708_mgarcia_syncom_invivo/exp01_inoculation_methods/pacbio_analysis/results/denoising"
# out.plots <- "/Volumes/D2c/mgarcia/20240708_mgarcia_syncom_invivo/exp01_inoculation_methods/pacbio_analysis/plots"

### Create outdirs ###

print("Creating directories")

dir.create(out.plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out.denois, recursive = TRUE, showWarnings = FALSE)

### Inputs ###

maxReads <- as.numeric(maxReads)
maxBases <- as.numeric(maxBases)
maxraref <- as.numeric(maxraref)

filtered_trimmed_reads_paths <- list.files(input.reads, full.names=TRUE)
reads_df <- read.table(input.readcounts, sep = "\t", header = T)

### Dereplicate sequences ###

print("Dereplicating sequences")

dereps <- derepFastq(filtered_trimmed_reads_paths, verbose=TRUE, n = maxReads)

print("Dereplication done")



### Build error model ###

print("Building error model")

set.seed(42)

if (errModel == "binnedQualErrfun"){
  print("Building error model with bins [3, 10, 17, 22, 27, 35, 40]")
  binnedQs <- c(3, 10, 17, 22, 27, 35, 40)
  binnedQualErrfun <- makeBinnedQualErrfun(binnedQs)
  error_model <- learnErrors(
    dereps,
    errorEstimationFunction = binnedQualErrfun,
    nbases = maxBases,
    randomize = T,
    BAND_SIZE = 32,
    multithread = T,
    verbose = T
  )
} else if (errModel == "PacBioErrfun") {
  print("Building error model with standard PacBio model")
  error_model <- learnErrors(
    dereps,
    errorEstimationFunction = PacBioErrfun,
    nbases = maxBases,
    randomize = T,
    BAND_SIZE = 32,
    multithread = T,
    verbose = T
  )
}

err.plot <- plotErrors(error_model, nominalQ=TRUE)
ggsave(file.path(out.plots, paste0("04_error_plot_", errModel, ".pdf")), err.plot, device="pdf", width = 10, height = 8)

saveRDS(error_model, file.path(out.denois, "dada2_error_model.RDS"))

print("Error model built")


### ASV inference ###

print("Denoising into ASVs")

dds <- list()

if (removeSingletons == "F"){
  print("Singleton detection OFF")
  for(i in seq_along(dereps)) {
    file = names(dereps)[i]
    print(paste0("Processing:", file))
    dds[[file]] <- dada(dereps[i], err=error_model, DETECT_SINGLETONS=FALSE, multithread=TRUE, verbose = T) # may leave singletons after merging
  }
}

if (removeSingletons == "T"){
  print("Singleton detection ON")
  for(i in seq_along(dereps)) {
    file = names(dereps)[i]
    print(paste0("Processing:", file))
    dds[[file]] <- dada(dereps[i], err=error_model, DETECT_SINGLETONS=TRUE, multithread=TRUE, verbose = T) # removes any singletons
  }
}

saveRDS(dds, file.path(out.denois, "denoised_seqs.rds"))

print("Denoising done")


### Generate ASV table ###

print("Generating ASV table")

#dds <- readRDS(file = file.path(out.denois, "denoised_seqs.rds"))
ASV_samples_table <- makeSequenceTable(dds)

print(paste0("Found ", ncol(ASV_samples_table), " ASVs across the ", nrow(ASV_samples_table), " samples"))

### Remove chimera ###

print("Removing chimera")

ASV_samples_table_noChim <- removeBimeraDenovo(ASV_samples_table, verbose = T, multithread = T)

#clean sample names
names <- basename(rownames(ASV_samples_table_noChim))
names <- gsub(".fastq.gz", "", names)
rownames(ASV_samples_table_noChim) <- names

saveRDS(ASV_samples_table_noChim, file.path(out.denois, "ASV_samples_table_noChim.rds"))

print(paste0("Number of reads removed: ", sum(ASV_samples_table)-sum(ASV_samples_table_noChim)))
print(paste0("Proportion of reads removed: ", round(1-sum(ASV_samples_table_noChim)/sum(ASV_samples_table),5)))

print("Chimera removal done")


### Rarefaction curves ###

print("Generating rarefaction curves")

# using iNEXT so I can also estimate the sampling coverage
dt <- iNEXT(
  t(ASV_samples_table_noChim),
  q = 0,
  datatype = "abundance",
  endpoint = maxraref,
  knots = 50,
  se = TRUE,
  conf = 0.95,
  nboot = 100
)

inextqd <- dt$iNextEst$size_based %>%
  dplyr::rename(sample = Assemblage)

sc.plot <- ggplot(
  inextqd[inextqd$Method != "Extrapolation", ],
  aes(
    x = m,
    y = SC,
    group = sample
  )
) +
  geom_vline(aes(xintercept = min(inextqd$m[inextqd$Method == "Observed"]), color = "low"), linetype = "dashed") + # sample with lowest number of reads
  geom_vline(aes(xintercept = max(inextqd$m[inextqd$Method == "Observed"]), color = "high"), linetype = "dashed") + # sample with highest number of reads
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
    x = "# of reads",
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
    group = sample
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
      label = sample
    ), size = 4, nudge_x = 70
  ) +
  theme_bw() +
  labs(
    x = "# of reads",
    y = "# of ASVs"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1,0.9),
    legend.background = element_rect(fill=alpha('white', 0.4))
  )

ggsave(file.path(out.plots, "04_sampling_coverage.pdf"), sc.plot, device="pdf", width = 8, height = 6)
ggsave(file.path(out.plots, "04_rarefaction_curves.pdf"), qd.plot, device="pdf", width = 10, height = 8)

write.table(
  inextqd,
  file = file.path(out.denois, "inext_data.tsv"),
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

### Compute general statistics ###

print("Plotting global statistics")

parts <- str_count(pattern = "/", rownames(ASV_samples_table_noChim)[[1]])

track <- tibble(sample = rownames(ASV_samples_table_noChim),
                Denoising = rowSums(ASV_samples_table), 
                `Chimera removal` = rowSums(ASV_samples_table_noChim)) %>% 
  mutate(basename = paste(sample, ".fastq.gz", sep = "")) %>%
  dplyr::select(-"sample") %>%
  pivot_longer(c("Denoising", "Chimera removal"), names_to = "stage", values_to = "n_reads")

reads_df2 <- rbind(reads_df, track) %>% arrange(basename) %>% 
  separate(basename, into = c("sample", NA), sep = ".fastq.gz", remove = F)

reads_df2$stage <- factor(reads_df2$stage, levels = c("Input", "Primers removed", "Post processing", "Denoising", "Chimera removal"), ordered = T)

track_medians <- reads_df2 %>% 
  group_by(stage) %>% 
  summarize(median = median(n_reads))

# plots
reads.plot1 <- ggplot(
  reads_df2,
  aes(
    y=n_reads,
    x=stage,
    group=basename
  )
) +
  geom_line() +
  geom_label(
    data = reads_df2[reads_df2$stage == "Chimera removal", ],
    aes(
      y = n_reads,
      x = stage,
      label = sample
    ), size = 3, nudge_x = 0.1
  ) +
  geom_point(
    data = track_medians,
    aes(
      y=median,
      x=stage,
      color="median"
    ), inherit.aes = F, size = 3
  ) +
  scale_color_manual(values = c(median = "red"), label = c(median = "Median"), name = "") +
  theme_bw() +
  labs(
    x = "Step",
    y = "# of sequences"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.1,0.1),
    legend.background = element_rect(fill=alpha('white', 0.4))
  )

reads.plot2 <- ggplot(
  reads_df2,
  aes(
    x=n_reads
  )
) +
  geom_histogram(bins = 20, boundary = 0, closed = "left") +
  theme_bw() +
  labs(
    x = "# of sequences",
    y = "# of samples"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(vars(stage), nrow = length(unique(reads_df2$stage)), scales = "free_y")

ggsave(file.path(out.plots, "04_read_number_lines.pdf"), reads.plot1, device="pdf", width = 10, height = 8)
ggsave(file.path(out.plots, "04_read_number_hist.pdf"), reads.plot2, device="pdf", width = 4, height = 8)


print("Finished plotting global statistics")


print("Denoising done")
