library(dada2); packageVersion("dada2")
# File parsing for the first round of sequencing
pathF_Reutigen  <- "/Users/megaptera/dada2/Reutigen/fastq/pathF" # CHANGE ME to the directory containing the fastq files before unzipping.
pathR_Reutigen <- "/Users/megaptera/dada2/Reutigen/fastq/pathR"
filtpathF_Reutigen <- file.path(pathF_Reutigen, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR_Reutigen <- file.path(pathR_Reutigen, "filtered") # ...
fastqFs_Reutigen <- sort(list.files(pathF_Reutigen, pattern="fq.gz"))
fastqRs_Reutigen <- sort(list.files(pathR_Reutigen, pattern="fq.gz"))
if(length(fastqFs_Reutigen) != length(fastqRs_Reutigen)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS; truncQ=2 for Guillaume.
filterAndTrim(fwd=file.path(pathF_Reutigen, fastqFs_Reutigen), filt=file.path(filtpathF_Reutigen, fastqFs_Reutigen),
              rev=file.path(pathR_Reutigen, fastqRs_Reutigen), filt.rev=file.path(filtpathR_Reutigen, fastqRs_Reutigen),
              truncLen=c(280,200), maxEE=c(2,5), maxN=0, trimLeft = c(28,27),
              compress=TRUE, verbose=TRUE, multithread=TRUE)


# Now quality filtering of the reads:
filtpathF_Reutigen_filtered <- "/Users/megaptera/dada2/Reutigen/fastq/pathF/filtered"
filtpathR_Reutigen_filtered <- "/Users/megaptera/dada2/Reutigen/fastq/pathR/filtered"
filtFs_Reutigen_filtered <- list.files(filtpathF_Reutigen_filtered, pattern="fq.gz", full.names = TRUE)
filtRs_Reutigen_filtered <- list.files(filtpathR_Reutigen_filtered, pattern="fq.gz", full.names = TRUE)
sample.names_Reutigen_filtered <- sapply(strsplit(basename(filtFs_Reutigen_filtered), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR_Reutigen_filtered <- sapply(strsplit(basename(filtRs_Reutigen_filtered), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names_Reutigen_filtered, sample.namesR_Reutigen_filtered)) stop("Forward and reverse files do not match.")
names(filtFs_Reutigen_filtered) <- sample.names_Reutigen_filtered
names(filtRs_Reutigen_filtered) <- sample.names_Reutigen_filtered
set.seed(20180524)
# Learn forward error rates
errF_Reutigen_filtered <- learnErrors(filtFs_Reutigen_filtered, nread=1e6, multithread=TRUE)
# Learn reverse error rates
errR_Reutigen_filtered <- learnErrors(filtRs_Reutigen_filtered, nread=1e6, multithread=TRUE)
# Sample inference and merger of paired-end reads

png(filename="/Users/megaptera/dada2/Reutigen/errF_Reutigen_filtered.png")
plotErrors(errF_Reutigen_filtered, nominalQ=TRUE)
dev.off()

png(filename="/Users/megaptera/dada2/Reutigen/errR_Reutigen_filtered.png")
plotErrors(errR_Reutigen_filtered, nominalQ=TRUE)
dev.off()

# Then merge forward and reverse:
mergers_Reutigen_filtered <- vector("list", length(sample.names_Reutigen_filtered))
names(mergers_Reutigen_filtered) <- sample.names_Reutigen_filtered
for(sam in sample.names_Reutigen_filtered) {
  cat("Processing:", sam, "\n")
  derepF_Reutigen_filtered <- derepFastq(filtFs_Reutigen_filtered[[sam]])
  ddF_Reutigen_filtered <- dada(derepF_Reutigen_filtered, err=errF_Reutigen_filtered, multithread=TRUE)
  derepR_Reutigen_filtered <- derepFastq(filtRs_Reutigen_filtered[[sam]])
  ddR_Reutigen_filtered <- dada(derepR_Reutigen_filtered, err=errR_Reutigen_filtered, multithread=TRUE)
  merger_Reutigen_filtered <- mergePairs(ddF_Reutigen_filtered, derepF_Reutigen_filtered, ddR_Reutigen_filtered, derepR_Reutigen_filtered)
  mergers_Reutigen_filtered[[sam]] <- merger_Reutigen_filtered
}
rm(derepF_Reutigen_filtered); rm(derepR_Reutigen_filtered)

# Construct sequence table
seqtab_Reutigen_filtered <- makeSequenceTable(mergers_Reutigen_filtered)
write.table(seqtab_Reutigen_filtered, file="/Users/megaptera/dada2/Reutigen/seqtab_Reutigen.txt")
saveRDS(seqtab_Reutigen_filtered, "/Users/megaptera/dada2/Reutigen/seqtab_Reutigen.rds") # CHANGE ME to where you want sequence table saved

# remove chimeras:
seqtab.chim.Reutigen <- removeBimeraDenovo(seqtab_Reutigen_filtered, method="consensus", multithread=TRUE)

# keep track of how many reads we kept after each round of filtering:
out.Reutigen <- filterAndTrim(fwd=file.path(pathF_Reutigen, fastqFs_Reutigen), filt=file.path(filtpathF_Reutigen, fastqFs_Reutigen),
                              rev=file.path(pathR_Reutigen, fastqRs_Reutigen), filt.rev=file.path(filtpathR_Reutigen, fastqRs_Reutigen),
                              truncLen=c(280,200), maxEE=c(2,5), maxN=0, trimLeft = c(28,27),
                              compress=TRUE, verbose=TRUE, multithread=TRUE)

track <- cbind(out.Reutigen, rowSums(seqtab_Reutigen_filtered), rowSums(seqtab.chim.Reutigen))
colnames(track) <- c("input","filtered","tabled", "nonchim")
rownames(track) <- sample.names_Reutigen_filtered
head(track)

write.table(track, file="/Users/megaptera/dada2/Reutigen/Reutigen_read_stats.txt")

# Assign taxonomy:
tax_Reutigen <- assignTaxonomy(seqtab.chim.Reutigen, "/Users/megaptera/dada2/silva_nr_v132_train_set.fa", multithread=TRUE) # done
tax_Reutigen_green <- assignTaxonomy(seqtab.chim.Reutigen, "/Users/megaptera/dada2/gg_13_8_train_set_97.fa", multithread=TRUE)
#tax_Reutigen_species <- addSpecies(tax_Reutigen, "/Users/megaptera/dada2/silva_species_assignment_v132.fa", multithread=TRUE)

saveRDS(seqtab.chim.Reutigen, "/Users/megaptera/dada2/Reutigen/seqtab_chim_Reutigen.rds") 
saveRDS(tax_Reutigen, "/Users/megaptera/dada2/Reutigen/tax_Reutigen.rds")
saveRDS(tax_Reutigen_green, "/Users/megaptera/dada2/Reutigen/tax_Reutigen_green.rds")
#saveRDS(tax_Reutigen_species, "/Users/megaptera/dada2/Reutigen/tax_Reutigen_species.rds")

write.table(tax_Reutigen_green, file="/Users/megaptera/dada2/Reutigen/tax_Reutigen_green.txt")
write.table(tax_Reutigen, file="/Users/megaptera/dada2/Reutigen/tax_Reutigen_silva.txt")

write.table(seqtab.chim.Reutigen, file="/Users/megaptera/dada2/Reutigen/seqtab_chim_Reutigen.txt")


# Now build a tree:

# First we need to get the right libraries.
library("knitr")
source("https://bioconductor.org/biocLite.R")
biocLite("BiocStyle")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install
  
  
  .packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}


# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

library(dada2)
library(phyloseq)
library(DECIPHER)
library(phangorn)

seqs.Reutigen <- getSequences(seqtab.chim.Reutigen)
names(seqs.Reutigen) <- seqs.Reutigen # This propagates to the tip labels of the tree
alignment.Reutigen <- AlignSeqs(DNAStringSet(seqs.Reutigen), anchor=NA,verbose=FALSE)
saveRDS(alignment.Reutigen, "/Users/megaptera/dada2/Reutigen/alignment_Reutigen.rds")

phangAlign.Reutigen <- phyDat(as(alignment.Reutigen, "matrix"), type="DNA")
dm1 <- dist.ml(phangAlign.Reutigen)
treeNJ1 <- NJ(dm1) # Note, tip order != sequence order
fit.Reutigen = pml(treeNJ1, data=phangAlign.Reutigen)
fitGTR1 <- update(fit.Reutigen, k=4, inv=0.2)
fitGTR1 <- optim.pml(fitGTR1, model="GTR", optInv=TRUE, optGamma=TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))

saveRDS(fitGTR1, "/Users/megaptera/dada2/Reutigen/fitGTR_Reutigen.rds")


library(ape)

write.tree(fitGTR1$tree, file="/Users/megaptera/dada2/Reutigen/Reutigen_tree.tre")


