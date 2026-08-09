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

if(!require(readr)){
  install.packages(pkgs = 'stringr', repos = 'https://stat.ethz.ch/CRAN/')
  library(readr)
}

if(!require(ggplot2)){
  install.packages(pkgs = 'ggplot2', repos = 'https://stat.ethz.ch/CRAN/')
  library(ggplot2)
}

if(!require(scales)){
  install.packages(pkgs = 'scales', repos = 'https://stat.ethz.ch/CRAN/')
  library(scales)
}

if(!require(dada2)){
  if(!require(devtools)){
    install.packages(pkgs = 'devtools', repos = 'https://stat.ethz.ch/CRAN/')
  }
  devtools::install_github("benjjneb/dada2") # installing through GitHub to get latest updates
  library(dada2)
}

if(!require(foreach)){
  install.packages(pkgs = 'foreach', repos = 'https://stat.ethz.ch/CRAN/')
  library(foreach)
}

if(!require(doParallel)){
  install.packages(pkgs = 'doParallel', repos = 'https://stat.ethz.ch/CRAN/')
  library(doParallel)
}

if(!require(iNEXT)){
  install.packages(pkgs = 'iNEXT', repos = 'https://stat.ethz.ch/CRAN/')
  library(iNEXT)
}

# Collect arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 10){
  stop(" Usage: 04_denoising.R <processed_reads_dir> <read_count_table> <max_reads_derep> <error_model> <max_bases_errormodel> <priors_db> <pool_T_F> <max_reads_raref> <denoise_results_dir> <plots_dir>", call.=FALSE)
} else {
  input.reads <- args[1] # folder with all pre-processed reads
  input.readcounts <- args[2] # table reporting the number of reads at each step
  maxReads <- args[3] # max number of reads to load at once for dereplication
  errModel <- args[4] # dada2-provided function to estimate the error model
  maxBases <- args[5] # max number of bases to use for error model inference
  db2 <- args[6] # database of ASVs expected in the samples
  pool <- args[7] # "T" or "pseudo" or "F", whether to pool samples for ASV inference
  maxraref <- args[8] # maximum number of reads to extrapolate rarefaction curves
  out.denois <- args[9] # folder to write denoising results
  out.plots <- args[10] # folder to write plots
}

# root <- file.path("/Volumes", "D2c", "mgarcia", "20240708_mgarcia_syncom_assembly", "pacbio_analysis", "run1_bees")
# input.reads <- file.path(root, "results", "preprocessing", "trimmed_filtered_reads")
# input.readcounts <- file.path(root, "results", "preprocessing", "read_count_before_after.tsv")
# maxReads <- 1E6
# errModel <- "binnedQualErrfun"
# maxBases <- 1E10
# db2 <- file.path(root, "data", "databases", "amplicon_based_db/syncom_custom_db_addSpecies.fa")
# pool <- "F"
# maxraref <- 8000
# out.denois <- file.path(root, "results", "denoising")
# out.plots <- file.path(root, "plots")

### Create outdirs ###

cat("\nCreating directories\n")

dir.create(out.plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out.denois, recursive = TRUE, showWarnings = FALSE)

### Inputs ###

maxReads <- as.numeric(maxReads)
maxBases <- as.numeric(maxBases)
maxraref <- as.numeric(maxraref)

filtered_trimmed_reads_paths <- list.files(input.reads, full.names=TRUE)
reads_df <- read.table(input.readcounts, sep = "\t", header = T)

n_samples <- length(unique(reads_df$basename))

### Build error model ###

# utility
checkConvergence <- function(dadaO) {
  sapply(dadaO$err_in, function(x) sum(abs(dadaO$err_out-x)))
}

# sample quality scores from 5 random samples to check whether they are binned or not
quals <- c()
for (file in sample(filtered_trimmed_reads_paths, 5, replace = F)){
  quals <- unique(append(quals, unique(plotQualityProfile(file)$data$Score)))
}

cat("Quality scores detected in the trimmed and filtered reads (five randomly chosen samples):\n")
cat(sort(quals))
cat("\n")

# if its exists, load pre-existing error model
if (file.exists(file.path(out.denois, "dada2_error_model.RDS"))){
  cat("Loading pre-existing error model\n")
  error_model <- readRDS(file.path(out.denois, "dada2_error_model.RDS"))
} else {
  cat("Building error model\n")
  
  set.seed(42)
  
  if (errModel == "binnedQualErrfun"){
    cat("Building error model with bins [3, 10, 17, 22, 27, 35, 40]\n")
    binnedQs <- c(3, 10, 17, 22, 27, 35, 40)
    binnedQualErrfun <- makeBinnedQualErrfun(binnedQs)
    error_model <- learnErrors(
      filtered_trimmed_reads_paths,
      errorEstimationFunction = binnedQualErrfun,
      nbases = maxBases,
      randomize = T,
      BAND_SIZE = 32,
      multithread = T,
      verbose = T
    )
  } else if (errModel == "PacBioErrfun") {
    cat("Building error model with standard PacBio model\n")
    error_model <- learnErrors(
      filtered_trimmed_reads_paths,
      errorEstimationFunction = PacBioErrfun,
      nbases = maxBases,
      randomize = T,
      BAND_SIZE = 32,
      multithread = T,
      verbose = T
    )
  }
  saveRDS(error_model, file.path(out.denois, "dada2_error_model.RDS"))
  cat("Error model built\n")
}

cat("Error model residuals:\n")
print(checkConvergence(error_model))

err.plot <- plotErrors(error_model, nominalQ=TRUE)
ggsave(file.path(out.plots, paste0("04_error_plot_", errModel, ".pdf")), err.plot, device="pdf", width = 10, height = 8)

### ASV inference ###

cat("Denoising into ASVs\n")

# options:
# pooling T/pseudo/F
# priors yes/no

if (file.exists(file.path(out.denois, "denoised_seqs.rds"))){
  cat("Loading pre-existing denoised samples\n")
  dds <- readRDS(file.path(out.denois, "denoised_seqs.rds"))
} else {
  if (db2 == ""){
    cat("No priors given\n")
    if (pool=="T"){
      dds <- dada(filtered_trimmed_reads_paths, err=error_model, pool=TRUE, multithread=TRUE, verbose = T)
    } else if (pool=="pseudo") {
      dds <- dada(filtered_trimmed_reads_paths, err=error_model, pool="pseudo", multithread=TRUE, verbose = T)
    } else if (pool=="F") {
      dds <- dada(filtered_trimmed_reads_paths, err=error_model, pool=FALSE, multithread=TRUE, verbose = T)
    } else {
      stop("Incorrect pooling option provided\n")
    }
  } else {
    cat(paste0("Using ", db2, " as priors\n"))
    my_priors <- getSequences(db2)
    if (pool=="T"){
      dds <- dada(filtered_trimmed_reads_paths, err=error_model, pool=TRUE, priors=my_priors, multithread=TRUE, verbose = T)
    } else if (pool=="pseudo") {
      dds <- dada(filtered_trimmed_reads_paths, err=error_model, pool="pseudo", priors=my_priors, multithread=TRUE, verbose = T)
    } else if (pool=="F") {
      dds <- dada(filtered_trimmed_reads_paths, err=error_model, pool=FALSE, priors=my_priors, multithread=TRUE, verbose = T)
    } else {
      stop("Incorrect pooling option provided\n")
    }
  }
  saveRDS(dds, file.path(out.denois, "denoised_seqs.rds"))
}

cat("Denoising done\n")


### Generate ASV table and remove chimera ###

if (file.exists(file.path(out.denois, "ASV_samples_table_noChim.rds"))){
  
  cat("Loading pre-existing ASV tables\n")
  ASV_samples_table <- readRDS(file.path(out.denois, "ASV_samples_table.rds"))
  ASV_samples_table_noChim <- readRDS(file.path(out.denois, "ASV_samples_table_noChim.rds"))
  
} else {
  
  cat("Generating ASV table and removing chimera\n")
  
  ASV_samples_table <- makeSequenceTable(dds)
  cat(paste0("Found ", ncol(ASV_samples_table), " ASVs across the ", nrow(ASV_samples_table), " samples", "\n"))
  saveRDS(ASV_samples_table, file.path(out.denois, "ASV_samples_table.rds"))
  
  if (pool %in% c("pseudo", "F")){
    ASV_samples_table_noChim <- removeBimeraDenovo(ASV_samples_table, verbose = T, multithread = T)
  } else {
    cat("Removing chimera in pooled mode")
    ASV_samples_table_noChim <- removeBimeraDenovo(ASV_samples_table, method = "pooled", verbose = T, multithread = T)
  }
  
  #clean sample names
  names <- basename(rownames(ASV_samples_table_noChim))
  names <- gsub(".fastq.gz", "", names)
  rownames(ASV_samples_table_noChim) <- names
  
  saveRDS(ASV_samples_table_noChim, file.path(out.denois, "ASV_samples_table_noChim.rds"))
  
  cat(paste0("Number of reads removed: ", sum(ASV_samples_table)-sum(ASV_samples_table_noChim), "\n"))
  cat(paste0("Proportion of reads removed: ", round(1-sum(ASV_samples_table_noChim)/sum(ASV_samples_table),5), "\n"))
  
  cat("Chimera removal done\n")
}

# cleanup
rm(dds)
invisible(gc())

### Compute general statistics ###

cat("Plotting global statistics\n")

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

if (n_samples >= 70){
  reads.plot1 <- ggplot(
    reads_df2,
    aes(
      y=n_reads,
      x=stage,
      group=basename
    )
  ) +
    geom_line(alpha=0.4) +
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
} else {
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
      ), size = 2, nudge_x = 0.1, alpha=0.4
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
}

reads.plot2 <- ggplot(
  reads_df2,
  aes(
    x=n_reads
  )
) +
  geom_histogram(bins = 20, boundary = 0, closed = "left") +
  scale_x_log10(label = label_log(), breaks=10^c(0:6)) +
  annotation_logticks(
    sides="b",
    outside=T,
    linewidth=0.3,
    short = unit(.5,"mm"),
    mid = unit(1,"mm"),
    long = unit(1.5,"mm")
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  labs(
    x = "# of sequences",
    y = "# of samples"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  facet_wrap(vars(stage), nrow = length(unique(reads_df2$stage)), scales = "free_y")

ggsave(file.path(out.plots, "04_read_number_lines.pdf"), reads.plot1, device="pdf", width = 10, height = 8)
ggsave(file.path(out.plots, "04_read_number_hist.pdf"), reads.plot2, device="pdf", width = 4, height = 8)

write.table(
  reads_df2,
  file = file.path(out.denois, "read_counts_steps.tsv"),
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

cat("Finished plotting global statistics\n")

### Rarefaction curves ###

## The range of values is split into numCores chunks to run in parallel

run_inext <- function(knots = 50, maxraref = numeric(), input, n = numCores){
  
  # define output file and log paths
  inext_file <- file.path(out.denois, "inext_data.tsv")
  log_file = file.path(out.denois, "inext.log")
  
  # clear output and log files
  if (file.exists(inext_file)) {
    cat("Deleting previous results file\n")
    file.remove(inext_file)
  }
  if (file.exists(log_file)) {
    cat("Deleting previous log file\n")
    file.remove(log_file)
  }
  
  cat("Preparing inputs\n")
  # filter out samples with total abundance of 0 in the input
  depth <- apply(input, 2, sum)
  samples_keep <- names(which(depth != 0))
  input2 <- input[,samples_keep]
  
  # create list of sizes for x knots between 1 and maxraref
  sizes <- round(c(1,c(2:knots)*maxraref/knots))
  
  # cut list into chunks for parallel processing
  if (n == 1){
    intervals = rep(1,knots)
  } else {
    if (knots %% n != 0){
      intervals <- c(cut(c(1:(knots - (knots %% n))), breaks=(n-1), label = FALSE), rep(n, knots %% n))
    } else {
      intervals <- cut(c(1:knots), breaks=n, label = FALSE)
    }
  }
  
  # run inext in parallel workers based on numCores
  
  cat("Running iNEXT in parallel\n")
  
  inext <- foreach(i = c(1:n), .packages = c("iNEXT", "dplyr")) %dopar% {
    # define sizes for which to compute diversity
    sizes_sub <- sizes[which(intervals == i)]
    
    # logging
    con <- file(log_file, open = "a")
    writeLines(paste0("Processing chunk ", i, "/", n, ": [", paste0(sizes_sub, collapse = ", "), "]"), con)
    close(con)
    
    # run iNEXT
    dt <- iNEXT(
      input2, # samples must be as columns
      q = c(0,1),
      datatype = "abundance",
      size = sizes_sub,
      se = TRUE,
      conf = 0.95,
      nboot = 10
    )
    
    # extract relevant data
    inextqd <- dt$iNextEst$size_based %>%
      dplyr::rename(SampleID = Assemblage)
    
    # save results to file (1 per chunk)
    write.table(inextqd, 
                file = paste0(sub(".tsv", "", inext_file), "_", i),
                sep = "\t", 
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE
    )
    
    # logging
    con <- file(log_file, open = "a")
    writeLines(paste0("Finished chunk ", i, "/", n), con)
    close(con)
    
    # return NULL so the master process stores nothing
    NULL
    
  }
  
  cat("Reading and plotting results\n")
  
  # read data
  inextqd_chunks <- lapply(
    paste(sub(".tsv", "", inext_file), "_", c(1:n), sep = ""),
    read_tsv,
    show_col_types = F
  )
  
  inextqd <- do.call(rbind, inextqd_chunks)
  
  # remove duplicated rows generated by iNEXT (due to sampling and rounding problems)
  inextqd <- inextqd %>%
    group_by(SampleID, m, Order.q) %>%
    mutate(dupl = row_number()) %>% 
    filter(dupl == 1) %>% 
    ungroup() %>%
    select(-"dupl") %>%
    arrange(SampleID, Order.q, m)
  
  # save results to file
  write.table(inextqd,
              file = inext_file,
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE
  )
  
  # remove old files
  lapply(
    paste(sub(".tsv", "", inext_file), "_", c(1:n), sep = ""),
    file.remove
  )
  
  # create plot
  if (length(unique(inextqd$SampleID)) >= 70){
    qd.plot <- ggplot(
      inextqd %>% filter(Method != "Extrapolation"),
      aes(
        x = m,
        y = qD,
        group = SampleID
      )
    ) +
      geom_vline(aes(xintercept = min(inextqd$m[inextqd$Method == "Observed"]), color = "low"), linetype = "dashed") + # sample with lowest number of reads
      geom_vline(aes(xintercept = max(inextqd$m[inextqd$Method == "Observed"]), color = "high"), linetype = "dashed") + # sample with highest number of reads
      geom_line(alpha = 0.4) +
      geom_ribbon(
        aes(
          ymin = qD.LCL,
          ymax = qD.UCL
        ), alpha = 0.1
      ) +
      scale_color_manual(name = "", values = c(low = "#669bbc", high = "#e76f51"), labels = c(low = "Lowest depth", high = "Highest depth")) +
      scale_y_continuous(breaks = seq(0,max(inextqd$qD[inextqd$Method != "Extrapolation"])+5,5)) +
      scale_x_log10(label = label_log(), breaks=10^c(0:6)) +
      annotation_logticks(sides = "b") +
      theme_bw() +
      labs(
        x = "# of reads",
        y = "# of ASVs"
      ) +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.1,0.9),
        legend.background = element_rect(fill=alpha('white', 0.4))
      ) +
      facet_wrap( ~ Order.q, scales = "free_y")
    
  } else {
    
    qd.plot <- ggplot(
      inextqd %>% filter(Method != "Extrapolation"),
      aes(
        x = m,
        y = qD,
        group = SampleID
      )
    ) +
      geom_vline(aes(xintercept = min(inextqd$m[inextqd$Method == "Observed"]), color = "low"), linetype = "dashed") + # sample with lowest number of reads
      geom_vline(aes(xintercept = max(inextqd$m[inextqd$Method == "Observed"]), color = "high"), linetype = "dashed") + # sample with highest number of reads
      geom_line(alpha = 0.4) +
      geom_ribbon(
        aes(
          ymin = qD.LCL,
          ymax = qD.UCL
        ), alpha = 0.1
      ) +
      geom_label(
        data = inextqd %>% filter(Method != "Extrapolation"),
        aes(
          x = m,
          y = qD,
          label = SampleID
        ), size = 3, nudge_x = log10(maxraref*1e-3)
      ) +
      scale_color_manual(name = "", values = c(low = "#669bbc", high = "#e76f51"), labels = c(low = "Lowest depth", high = "Highest depth")) +
      scale_y_continuous(breaks = seq(0,max(inextqd$qD[inextqd$Method != "Extrapolation"])+5,5)) +
      scale_x_log10(label = label_log(), breaks=10^c(0:6)) +
      annotation_logticks(sides = "b") +
      theme_bw() +
      labs(
        x = "# of reads",
        y = "# of ASVs"
      ) +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.1,0.9),
        legend.background = element_rect(fill=alpha('white', 0.4))
      ) +
      facet_wrap( ~ Order.q, scales = "free_y")
  }
  
  # save plot
  ggsave(file.path(out.plots, paste0("04_rarefaction_curves_ASVs.pdf")), qd.plot, device="pdf", width = 12, height = 8)
  
}

if (maxraref <= 0){
  cat("Skipping rarefaction curves\n")
} else {
  cat("Generating rarefaction curves\n")
  
  ### Set up parallel backend ###
  
  #numCores <- max(1, detectCores() - 2)
  numCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
  cat(paste0("Number of detected cores: ", numCores, "\n"))
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  run_inext(knots = 50, maxraref = maxraref, input = t(ASV_samples_table_noChim))
  
  ### Stop cluster ###
  stopCluster(cl)
}

cat("Denoising done\n")
