
library(dada2); packageVersion("dada2")
library(here)
# File parsing
path_raw <- here("data/raw/SS-11666_16s") # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
filtpath <-here("data/processed/filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
fastqFs <- sort(list.files(path_raw, pattern="R1_001.fastq.gz"))
fastqRs <- sort(list.files(path_raw, pattern="R2_001.fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(path_raw, fastqFs), filt=file.path(filtpath, fastqFs),
              rev=file.path(path_raw, fastqRs), filt.rev=file.path(filtpath, fastqRs),
              truncLen=c(240,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

