library(dada2); packageVersion("dada2")
library(here)
# File parsing
path_raw <- here("data/raw/SS-11666_16s") # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
filtpath <-here("data/processed/filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
fastqFs <- sort(list.files(path_raw, pattern="R1_001.fastq.gz"))
fastqRs <- sort(list.files(path_raw, pattern="R2_001.fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

filtFs = file.path(filtpath, paste0(fastqFs, "R1_001_filt.fastq.gz"))
filtRs = file.path(filtpath, paste0(fastqRs, "R2_001_filt.fastq.gz"))

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "results/seqtab/seqtab.rds") # CHANGE ME to where you want sequence table saved

# Merge multiple runs (if necessary)
st1 <- readRDS("results/seqtab/seqtab.rds")

# Remove chimeras
seqtab <- removeBimeraDenovo(st1, method="consensus", multithread=TRUE)
# Assign taxonomy
tax <- assignTaxonomy(seqtab, "workflow/resources/databases/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
# Write to disk
saveRDS(seqtab, "results/seqtab/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax, "results/tax_final.rds") # CHANGE ME ...
