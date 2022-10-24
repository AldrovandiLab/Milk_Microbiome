#!/usr/bin/Rscript

library(ggplot2)
library(plyr)
library(reshape2)
library(dada2)
library(stringi)

outdir <- "/Lab_Share/PROMISE/nwcs610"

## Quality profiles
ncores <- 16
pdf(sprintf("%s/phyloseq/qualityProfiles.pdf", outdir))
fnFs <- list.files(sprintf("%s/fastq/split.R1/", outdir), pattern=".fastq.gz", full.names=T, recursive=T)
fnRs <- list.files(sprintf("%s/fastq/split.R2/", outdir), pattern=".fastq.gz", full.names=T, recursive=T)
samples <- sapply(fnFs, function(x) gsub("_R1.fastq.gz", "", basename(x)))
# plot quality profiles
p <- plotQualityProfile(fnFs, aggregate=T) + ggtitle(sprintf("Quality profile: (R1)"))
print(p)
p <- plotQualityProfile(fnRs, aggregate=T) + ggtitle(sprintf("Quality profile: (R2)"))
print(p)
dev.off()

## set truncation length parameters - these can vary as long as the same amplicon can be reconstructed
## DO NOT use trimLeft
truncation_lengths <- c(225,225)
ncores <- 16

## DADA2 inference
# make list of FASTQ files
fnFs <- list.files(sprintf("%s/fastq/split.R1/", outdir), pattern=".fastq", full.names=T, recursive=T)
fnRs <- list.files(sprintf("%s/fastq/split.R2/", outdir), pattern=".fastq", full.names=T, recursive=T)
samples <- sapply(fnFs, function(x) gsub("_R1.fastq.gz", "", basename(x)))
# filter FASTQs
filt_path <- file.path(sprintf("%s/fastq/DADA2_filtered", outdir))
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
res <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncation_lengths, trimLeft=c(0,0), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=ncores, matchIDs=TRUE, verbose=TRUE)
res <- as.data.frame(res)
# filter out empty samples
exists <- file.exists(filtFs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]
samples <- samples[exists]
# learn error rates from subset of data
errF <- learnErrors(filtFs, multithread=ncores)
errR <- learnErrors(filtRs, multithread=ncores)
# sample inference with learned error rates
dadaFs <- dada(filtFs, err=errF, selfConsist = FALSE, multithread=ncores)
dadaRs <- dada(filtRs, err=errR, selfConsist = FALSE, multithread=ncores)
# plot error rates
pdf(sprintf("%s/phyloseq/errorRates.pdf", outdir))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()
# merge read pairs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# make OTU table
seqtab <- makeSequenceTable(mergers)
rownames(seqtab) <- gsub("_F_filt.fastq.gz", "", rownames(seqtab))
counts.summary <- data.frame(SampleID=rownames(seqtab), input=res$reads.in, filtered=res$reads.out, merged=rowSums(seqtab))
saveRDS(seqtab, sprintf("%s/fastq/merged_seqtab.rds", outdir))
save.image(file=sprintf("%s/fastq/DADA2_workspace.RData", outdir))
# remove singletons (sequences that only appear in one sample)
#	seqtab <- seqtab[, colSums(seqtab>0)>1]
# in-silico size restriction
table(nchar(getSequences(seqtab)))
#seqtab <- seqtab[, nchar(colnames(seqtab)) %in% seq(from=252,to=254)] ### EDIT MANUALLY
seqtab <- seqtab[, nchar(colnames(seqtab)) %in% seq(from=251,to=255)] ### EDIT MANUALLY
counts.summary$size.selected <- rowSums(seqtab)
saveRDS(seqtab, sprintf("%s/fastq/merged_seqtab.size_restricted.rds", outdir))
# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE, multithread=ncores)
counts.summary$chimera.filtered <- rowSums(seqtab.nochim)
saveRDS(seqtab.nochim, sprintf("%s/fastq/merged_seqtab.nochim.rds", outdir))
write.table(counts.summary, file=sprintf("%s/phyloseq/DADA2_count_summary.txt", outdir), quote=F, row.names=F, col.names=T, sep="\t")

# assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/share/DADA2/rdp_train_set_18.fa.gz", multithread=ncores)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
taxa <- addSpecies(taxa, "/share/DADA2/rdp_species_assignment_18.fa.gz", verbose=TRUE)
unname(head(taxa))
save(seqtab.nochim, taxa, file=sprintf("%s/fastq/DADA2.RData", outdir))

# save for BLAST taxonomy assignment
out <- data.frame(id=sprintf(">seq%06d", 1:ncol(seqtab.nochim)), seq=colnames(seqtab.nochim))
write.table(out, file=sprintf("%s/fastq/DADA2_sequences.fasta", outdir), quote=F, sep="\n", row.names=F, col.names=F)





