#DADA2 for BIG DATA


#setup
library(dada2); packageVersion("dada2")
library(here)
# Filename parsing

path <- here("data/processed") # CHANGE ME to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
fns <- list.files(path2, pattern="fastq.gz") # CHANGE if different file extensions
# Filtering
filterAndTrim(file.path(path,fns), file.path(filtpath,fns), 
              truncLen=240, maxEE=1, truncQ=11, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

#Infer sequence variants



# File parsing
filtpathF <- "data/processed/filtered/FWD" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "data/processed/filtered/REV" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
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

#merge runs, remove chimeras, assign taxonomy


# Merge multiple runs (if necessary)
st1 <- readRDS("results/seqtab/seqtab.rds")
#st2 <- readRDS("path/to/run2/output/seqtab.rds")
#st3 <- readRDS("path/to/run3/output/seqtab.rds")
#st.all <- mergeSequenceTables(st1, st2, st3)
# Remove chimeras
seqtab <- removeBimeraDenovo(st.1, method="consensus", multithread=TRUE)
# Assign taxonomy
tax <- assignTaxonomy(seqtab, "C:/Users/crazy/OneDrive/Documents/FinalProjectSGS/Data/Seq/Unzipped/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
# Write to disk
saveRDS(seqtab, "results/seqtab/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax, "results/seqtab/tax_final.rds") # CHANGE ME ...
