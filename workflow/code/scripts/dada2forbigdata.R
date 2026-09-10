#DADA2 for BIG DATA

#library
install.packages("devtools")
install.packages("rlang")
install.packages("bioconductor")

library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16")

install.packages(c('Biostrings', 'ShortRead', 'RcppParallel', 'IRanges', 'XVector', 'BiocGenerics'))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead", version = "3.23") 
# Or you can use earlier versions of Biocondcutor, e.g. if you have a pre-3.6 version of R


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ShortRead")

#filter
library(dada2); packageVersion("dada2")
# Filename parsing
path <- "~/SeqData_SGS/data/processed" # CHANGE ME to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
fns <- list.files(path2, pattern="fastq.gz") # CHANGE if different file extensions
# Filtering
filterAndTrim(file.path(path,fns), file.path(filtpath,fns), 
              truncLen=240, maxEE=1, truncQ=11, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

#Infer sequence variants

library(dada2); packageVersion("dada2")
# File parsing
filtpath <- "~/SeqData_SGS/data/processed/filtered" # CHANGE ME to the directory containing your filtered fastq files
filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE) # CHANGE if different file extensions
sample.names <- sapply(strsplit(basename(filts), "_S"), `[`, 1) # Assumes filename = sample_XXX.fastq.gz
names(filts) <- sample.names
# Learn error rates
set.seed(100)
err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE)
# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
}
# Construct sequence table and write to disk
seqtab <- makeSequenceTable(dds)
saveRDS(seqtab, "C:/Users/crazy/OneDrive/Documents/FinalProjectSGS/Data/Seq/Unzipped/seqtab/seqtab.rds") # CHANGE ME to where you want sequence table saved

foo <- readRDS("C:/Users/crazy/OneDrive/Documents/FinalProjectSGS/Data/Seq/Unzipped/seqtab/seqtab.rds")

#merge runs, remove chimeras, assign taxonomy

library(dada2); packageVersion("dada2")
# Merge multiple runs (if necessary)
st1 <- readRDS("~/SeqData-SGS/results/seqtab/seqtab.rds")
#st2 <- readRDS("path/to/run2/output/seqtab.rds")
#st3 <- readRDS("path/to/run3/output/seqtab.rds")
#st.all <- mergeSequenceTables(st1, st2, st3)
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
tax <- assignTaxonomy(seqtab, "C:/Users/crazy/OneDrive/Documents/FinalProjectSGS/Data/Seq/Unzipped/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
# Write to disk
saveRDS(seqtab, "~/SeqData-SGS/results/seqtab/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax, "~/SeqData-SGS/results/seqtab/tax_final.rds") # CHANGE ME ...
