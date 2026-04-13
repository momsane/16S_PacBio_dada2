# Load required libraries

if(!require(foreach)){
  install.packages(pkgs = 'foreach', repos = 'https://stat.ethz.ch/CRAN/')
  library(foreach)
}

if(!require(doParallel)){
  install.packages(pkgs = 'doParallel', repos = 'https://stat.ethz.ch/CRAN/')
  library(doParallel)
}

if(!require(dplyr)){
  install.packages(pkgs = 'dplyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(dplyr)
}

if(!require(tidyr)){
  install.packages(pkgs = 'tidyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(tidyr)
}

if(!require(readr)){
  install.packages(pkgs = 'readr', repos = 'https://stat.ethz.ch/CRAN/')
  library(readr)
}

if(!require(stringr)){
  install.packages(pkgs = 'stringr', repos = 'https://stat.ethz.ch/CRAN/')
  library(stringr)
}

if(!require(ggplot2)){
  install.packages(pkgs = 'ggplot2', repos = 'https://stat.ethz.ch/CRAN/')
  library(ggplot2)
}

if(!require(gridExtra)){
  install.packages(pkgs = 'gridExtra', repos = 'https://stat.ethz.ch/CRAN/')
  library(gridExtra)
}

if(!require(data.table)){
  install.packages(pkgs = 'data.table', repos = 'https://stat.ethz.ch/CRAN/')
  library(data.table)
}

if(!require(ShortRead)){
  install.packages(pkgs = 'ShortRead', repos = 'https://stat.ethz.ch/CRAN/')
  library(ShortRead)
}

if(!require(dada2)){
  if(!require(devtools)){
    install.packages(pkgs = 'devtools', repos = 'https://stat.ethz.ch/CRAN/')
  }
  devtools::install_github("benjjneb/dada2") # installing through GitHub to get latest updates
  library(dada2)
}

# Collect arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9){
  stop(" Usage: 02_preprocessing.R <raw_reads_dir> <fwd_primer_sequence> <rev_primer_sequence> <min_read_length> <max_read_length> <maxEE> <overwrite> <preproc_results_dir> <plots_dir>", call.=FALSE)
} else {
  input.raw <- args[1] # folder with all raw read files (if you have pre-rarefied the reads, give the folder with the pre-rarefied fastq files)
  fwd.primer <- args[2] # forward primer sequence
  rev.primer <- args[3] # reverse primer sequence
  minLen <- args[4] # lower bound for read length, for filterAndTrim
  maxLen <- args[5] # upper bound for read length, for filterAndTrim
  maxEE <- args[6] # max number of errors per read
  overwrite <- args[7] # whether to overwrite existing files or not
  out.preproc <- args[8] # folder to write pre-processing results
  out.plots <- args[9] # folder to write plots (can be the same folder throughout the pipeline)
}

# root <- file.path("/Volumes", "D2c", "mgarcia", "20240708_mgarcia_syncom_assembly", "pacbio_analysis", "run1_bees")
# input.raw <- file.path(root, "results", "prerarefied_reads")
# fwd.primer <- "AGRGTTYGATYMTGGCTCAG"
# rev.primer <- "RGYTACCTTGTTACGACTT"
# minLen <- 1200
# maxLen <- 1700
# maxEE <- 3
# out.preproc <- file.path(root, "results", "preprocessing")
# out.plots <- file.path(root, "plots")

### Create outdirs ###

cat("\nCreating directories\n")

dir.create(out.plots, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out.preproc, "primerfree_reads"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out.preproc, "trimmed_filtered_reads"), recursive = TRUE, showWarnings = FALSE)

### Inputs ###

minLen <- as.numeric(minLen)
maxLen <- as.numeric(maxLen)
maxEE <- as.numeric(maxEE)

raw_reads_paths <- list.files(input.raw, full.names=TRUE)



### Set up parallel backend ###

# numCores <- detectCores()
numCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
cat(paste0("Number of detected cores: ", numCores, "\n"))
cl <- makeCluster(numCores)
registerDoParallel(cl)



### Orient sequences and trim primers ###

cat("Orienting sequences and trimming primers\n")

# Define paths for primers-free reads
trimmed_reads_paths <- file.path(out.preproc, "primerfree_reads", basename(raw_reads_paths))

cat(paste0("Primer-free reads written to: ", file.path(out.preproc, "primerfree_reads"), "\n"))
primer_removal_summary_path <- file.path(out.preproc, "primer_removal_summary.tsv")

# Check if some files already exist
test_exist_trimmed <- unlist(lapply(trimmed_reads_paths, file.exists))
existing_trimmed_samples <- raw_reads_paths[which(test_exist_trimmed == TRUE)]

if (length(existing_trimmed_samples)!=0){
  if (overwrite == "F"){
    cat(paste0(length(existing_trimmed_samples), " sample(s) already processed, they will not be overwritten\n"))
    samples_totrim <- raw_reads_paths[which(test_exist_trimmed == FALSE)]
    trimmed_reads_paths2 <- trimmed_reads_paths[which(test_exist_trimmed == FALSE)]
  } else if (overwrite == "T"){
    samples_totrim <- raw_reads_paths
    trimmed_reads_paths2 <- trimmed_reads_paths
    cat(paste0(length(existing_trimmed_samples), " already existing file(s) will be overwritten\n"))
    file.remove(primer_removal_summary_path)
  }
} else {
  samples_totrim <- raw_reads_paths
  trimmed_reads_paths2 <- trimmed_reads_paths
}

trim <- foreach(i = seq_along(samples_totrim), .packages = c("dada2")) %dopar% {
  # logging
  con <- file(file.path(out.preproc, "primer_removal.log"), open = "a")
  writeLines(paste0(i, " - ", basename(samples_totrim[i])), con)
  close(con)
  
  # rm primers
  res <- removePrimers(fn = samples_totrim[i],
                       fout = trimmed_reads_paths2[i], 
                       primer.fwd = fwd.primer,
                       primer.rev = rc(rev.primer),
                       orient=TRUE,
                       verbose=TRUE,
                       max.mismatch = 4
  ) # cannot multithread with dada2
  
  res <- as.data.frame(res)
  res$basename = rownames(res)
  
  # save stats to table
  write.table(res, 
              file = primer_removal_summary_path,
              sep = "\t", 
              append = TRUE, 
              col.names = !file.exists(primer_removal_summary_path),
              row.names = FALSE,
              quote = FALSE
  )
  
  # return NULL so the master process stores nothing
  NULL
}

primer_removal_summary <- read_tsv(file.path(out.preproc, "primer_removal_summary.tsv"), show_col_types = F)

cat(paste0("Mean proportion of reads retained: ", round(mean(primer_removal_summary$reads.out/primer_removal_summary$reads.in),2),"\n"))

primer_removal_summary <- primer_removal_summary %>%
  pivot_longer(c("reads.in", "reads.out"), values_to = "n_reads", names_to = "stage") %>%
  mutate(stage = ifelse(stage == "reads.in", "Input", "Primers removed"))

cat("Orientation and primer removal done\n")


### Trim and filter reads based on quality and length ###

cat("Trimming and filtering reads\n")

# Define paths for trimmed and filtered reads
filtered_trimmed_reads_paths <- file.path(out.preproc, "trimmed_filtered_reads", basename(trimmed_reads_paths))

cat(paste0("Filtered and primer-free reads written to: ", file.path(out.preproc, "trimmed_filtered_reads"), "\n"))
filtering_summary_path <- file.path(out.preproc, "filtering_summary.tsv")

# Check if some files already exist
test_exist_filtered <- unlist(lapply(filtered_trimmed_reads_paths, file.exists))
existing_filtered_samples <- trimmed_reads_paths[which(test_exist_filtered == TRUE)]

if (length(existing_filtered_samples)!=0){
  if (overwrite == "F"){
    cat(paste0(length(existing_filtered_samples), " sample(s) already processed, they will not be overwritten\n"))
    samples_tofilt <- trimmed_reads_paths[which(test_exist_filtered == FALSE)]
    filtered_trimmed_reads_paths2 <- filtered_trimmed_reads_paths[which(test_exist_filtered == FALSE)]
  } else if (overwrite == "T"){
    samples_tofilt <- trimmed_reads_paths
    filtered_trimmed_reads_paths2 <- filtered_trimmed_reads_paths
    cat(paste0(length(existing_filtered_samples), " already existing file(s) will be overwritten\n"))
    file.remove(filtering_summary_path)  
  }
} else {
  samples_tofilt <- trimmed_reads_paths
  filtered_trimmed_reads_paths2 <- filtered_trimmed_reads_paths
}

filt <- foreach(i = seq_along(samples_tofilt), .packages = c("dada2")) %dopar% {
  # logging
  con <- file(file.path(out.preproc, "filtering.log"), open = "a")
  writeLines(paste0(i, " - ", basename(samples_tofilt[i])), con)
  close(con)
  
  # rm primers
  res <- filterAndTrim(fwd = samples_tofilt[i],
                       filt = filtered_trimmed_reads_paths2[i],
                       minLen=minLen,
                       maxLen=maxLen,
                       maxN=0, # no ambiguous bases
                       rm.phix=FALSE,
                       maxEE=maxEE,
                       multithread = FALSE,
                       verbose = T) 
  
  res <- as.data.frame(res)
  res$basename = rownames(res)
  
  # save stats to table
  write.table(res, 
              file = filtering_summary_path,
              sep = "\t", 
              append = TRUE, 
              col.names = !file.exists(filtering_summary_path),
              row.names = FALSE,
              quote = FALSE
  )
  
  # return NULL so the master process stores nothing
  NULL
}

filtering_summary <- read_tsv(file.path(out.preproc, "filtering_summary.tsv"), show_col_types = F)

cat(paste0("Mean proportion of reads retained: ", round(mean(filtering_summary$reads.out/filtering_summary$reads.in),2), "\n"))

filtering_summary <- filtering_summary %>%
  pivot_longer(c("reads.in", "reads.out"), values_to = "n_reads", names_to = "stage") %>%
  mutate(stage = ifelse(stage == "reads.in", "Primers removed", "Post processing"))

cat("Trimming and filtering reads done\n")

### Summary stats on number of reads and read length ###

cat("Computing read statistics\n")

# Read length

compute_all_readL <- function(paths, stage){
  list <- foreach(i = seq_along(paths), .packages = c("ShortRead"), .inorder = FALSE) %dopar% {
    
    fq <- readFastq(paths[i])
    lens <- width(fq)
    rm(fq) # rm to cleanup memory
    
    data.frame(
      basename = rep(basename(paths[i]), length(lens)),
      read = seq_along(lens),
      length = lens,
      stage = stage,
      stringsAsFactors = FALSE
    )
  }
  df <- rbindlist(list)
  return(df)
}

lens_raw_df <- compute_all_readL(raw_reads_paths, "Input")
lens_bf_df <- compute_all_readL(trimmed_reads_paths, "Primers removed")
lens_af_df <- compute_all_readL(filtered_trimmed_reads_paths, "Post processing")

lens_df <- rbind(lens_raw_df, lens_bf_df, lens_af_df) %>% arrange(basename) # combine before and after data
lens_df$stage <- factor(lens_df$stage, levels = c("Input", "Primers removed", "Post processing"), ordered = T)

# plot read length before and after pre-processing
length.plot <- ggplot(
  lens_df,
  aes(
    x = length
  )
) +
  geom_histogram(bins = 50, boundary = 0, closed = "left") +
  scale_y_log10() +
  theme_bw() +
  labs(
    x = "Length (bp)",
    y = "# of sequences"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(vars(stage), nrow = 3)


# Read count

reads_df <- rbind(primer_removal_summary, filtering_summary[filtering_summary$stage == "Post processing", ]) %>% arrange(basename)
reads_df$stage <- factor(reads_df$stage, levels = c("Input", "Primers removed", "Post processing"), ordered = T)

ns.plot <- ggplot(
  reads_df,
  aes(
    x = n_reads
  )
) +
  geom_histogram(bins = 20, boundary = 0, closed = "left") +
  theme_bw()+
  labs(
    x = "# of sequences",
    y = "# of samples"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(vars(stage), nrow = 3)

qc_afbf <- grid.arrange(length.plot, ns.plot, ncol=2)

cat("Saving plots and tables\n")

# save plot
ggsave(file.path(out.plots, "02_read_stats_before_after.pdf"), qc_afbf, device="pdf", width = 8, height = 8)

# save stats
write.table(lens_df, file.path(out.preproc, "read_length_before_after.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(reads_df, file.path(out.preproc, "read_count_before_after.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)


cat("All pre-processing done\n")

### Stop cluster ###
stopCluster(cl)
