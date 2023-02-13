#!/usr/bin/Rscript

library(ggplot2)
library(ape)
library(plyr)
library(reshape2)
library(cluster)
library(RColorBrewer)
library(phyloseq)
library(decontam)
library(grid)
library(gridExtra)
library(gplots)
library(vegan)
library(irr)
library(parallel)
library(useful)
library(pscl)
library(MASS)
library(boot)
library(igraph)
library(ggfortify)
library(ggforce)
library(stringi)
library(stringr)
library(car)
library(randomForest)
library(ranger)
library(ROCR)
library(Hmisc)
library(missForest)
library(dada2)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(DESeq2)
library(NbClust)
library(tableone)
library(xgboost)
library(shapr)
library(Ckmeans.1d.dp)

source("utils.R")
source("mcc.R")
source("GMPR.R")

distance_metrics <- c("bray", "jaccard", "jsd")
alpha_metrics <- c("Chao1", "Shannon", "Simpson", "Observed")
icc_levels <- c(5,6,7)
names(icc_levels) <- c("Family", "Genus", "Species")
ncores <- 16
cols.sig <- c("black", "red", "grey"); names(cols.sig) <- c("ns", "sig", "NA")
cols.mvar <- list("Transmitter"=c("Control"="grey", "Transmitter"="red"), "Regimen"=c("Infant NVP"="blue", "Maternal triple ARV"="green"))
siglevel <- 0.05

## thresholds for association/permutation tests
nsamps_threshold <- 0.01 # fraction of relabund to call a sample positive (for Genus)
filt_threshold <- 0.1 # fraction of samples that need to be positive to keep an OTU for association testing
nsamps_threshold.species <- 0.001 # fraction of relabund to call a sample positive (for Species)
filt_threshold.species <- 0.01
thresholds <- list(nsamps=list(Genus=0.01, Species=0.001), filt=list(Genus=0.1, Species=0.01))

nperm <- 100000
siglevel <- 0.05

#################################################################
## handoff to phyloseq
dada2_fn <- "/Lab_Share/PROMISE/nwcs610/fastq/DADA2.RData"
mapping_fn <- "/Lab_Share/PROMISE/nwcs610/combined_nwcs610_ZEBS_Mapping.with_metadata.080922.txt"
output_dir <- "/Lab_Share/PROMISE/nwcs610/phyloseq"
out_txt <- sprintf("%s/phyloseq_output.%s.txt", output_dir, format(Sys.Date(), "%m%d%y"))
out_pdf <- sprintf("%s/phyloseq_output.%s.pdf", output_dir, format(Sys.Date(), "%m%d%y"))

## load mapping and get read counts DADA2 summary
mapping <- read.table(mapping_fn, header=T, sep="\t", comment.char="", quote="", as.is=T)
rownames(mapping) <- mapping$Sample.ID
read_counts <- read.table("/Lab_Share/PROMISE/nwcs610/phyloseq/DADA2_count_summary.txt", header=T, as.is=T, sep="\t", row.names=1)
sel <- intersect(rownames(mapping), rownames(read_counts))
mapping <- mapping[sel,]
mapping <- merge(mapping, read_counts, by="row.names"); rownames(mapping) <- mapping$Row.names; mapping <- mapping[,-1]


#################################################################
## load all samples run to first examine controls
load(dada2_fn)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
sel <- intersect(rownames(mapping), rownames(seqtab.nochim))
mapping <- mapping[sel,]
seqtab.nochim <- seqtab.nochim[sel,]
# fix sample names
rownames(mapping) <- make.names(rownames(mapping))
rownames(seqtab.nochim) <- make.names(rownames(seqtab.nochim))
mapping$NumReadsOTUTable <- rowSums(seqtab.nochim)[rownames(mapping)]
mapping$SampleIDstr <- sprintf("%s (%d)", rownames(mapping), rowSums(seqtab.nochim)[rownames(mapping)])
ps <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows=TRUE), sample_data(mapping), tax_table(taxa))
set.seed(prod(dim(seqtab.nochim)))

color_table <- read.table("/Lab_Share/fanli/code/Core.16S/taxa_coloring.Genus.050818.txt", header=T, as.is=T, sep="\t", comment.char="")
coloring <- color_table$Color
names(coloring) <- color_table$Genus
ordering <- rev(names(coloring))
cols <- colorRampPalette(c("white", "red"), space = "rgb")
missing_colors <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"), brewer.pal(8, "Pastel1"), brewer.pal(8, "Pastel2"))

ordering.genus <- color_table
ordering.genus$Phylum <- factor(ordering.genus$Phylum)
ordering.genus$Class <- factor(ordering.genus$Class)
inds=order(ordering.genus$Phylum, ordering.genus$Class, ordering.genus$Genus)
ordering.genus <- ordering.genus[inds,]
ordering.genus$Genus <- factor(ordering.genus$Genus, levels=unique(ordering.genus$Genus))

dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")

blank_types <- c("DNAfreewater", "BlankSwab", "EmptyWell")
pdf(out_pdf, width=12)

##################################################################################
## remove empty samples
empty_samples <- names(which(sample_sums(ps)==0))
ps <- prune_samples(setdiff(sample_names(ps), empty_samples), ps)
mapping <- as(sample_data(ps), "data.frame")

### decontam
sample_data(ps)$Sample_or_Control <- ifelse(sample_data(ps)$SampleType %in% c("BMK", "MockDNA"), "Sample", "Control") #label Sample (samples + mock) or Control (Buffer + PCR)

#inspect library sizes
mapping.sel <- as.data.frame(sample_data(ps))
mapping.sel$LibrarySize <- sample_sums(ps)
mapping.sel <- mapping.sel[order(mapping.sel$LibrarySize),]
mapping.sel$Index <- seq(nrow(mapping.sel))
ggplot(data=mapping.sel, aes(x=Index, y=LibrarySize, color=SampleType)) + geom_point()

#make contam object
set.seed(100)
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contam <- isContaminant(ps, method="auto", neg="is.neg")
ps.clean <- prune_taxa(!contam$contaminant, ps)
ps <- ps.clean

##################################################################################
## store metadata variables
metadata_variables <- read.table("/Lab_Share/PROMISE/nwcs610/metadata_variables.072722.txt", header=T, as.is=T, sep="\t", row.names=1)
sel <- intersect(rownames(metadata_variables), colnames(mapping))
metadata_variables <- metadata_variables[sel,, drop=F]
mapping.sel <- mapping[rownames(sample_data(ps)), sel]
# fix column types
for (mvar in rownames(metadata_variables)) {
  if (metadata_variables[mvar, "type"] == "factor") {
    mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
    if (!(is.na(metadata_variables[mvar, "baseline"])) && metadata_variables[mvar, "baseline"] != "") {
    	if (grepl(",", metadata_variables[mvar, "baseline"])) {
    		lvls <- unlist(strsplit(metadata_variables[mvar, "baseline"], ","))
    		mapping.sel[,mvar] <- factor(as.character(mapping.sel[,mvar]), levels=lvls)
    	} else {
      	mapping.sel[,mvar] <- relevel(mapping.sel[,mvar], metadata_variables[mvar, "baseline"])
      }
    }
  } else if (metadata_variables[mvar, "type"] == "numeric") {
    mapping.sel[,mvar] <- as.numeric(as.character(mapping.sel[,mvar]))
  } else if (metadata_variables[mvar, "type"] == "date") {
    mapping.sel[,mvar] <- as.Date(mapping.sel[,mvar], format="%Y%m%d")
    mapping.sel[,mvar] <- factor(as.character(mapping.sel[,mvar]), levels=as.character(unique(sort(mapping.sel[,mvar]))))
  }
}
sample_data(ps) <- mapping.sel


##################################################################################
## un-rarefied/filtered data for overall analysis including blanks
## remove empty samples again (post-decontam)
empty_samples <- names(which(sample_sums(ps)==0))
ps <- prune_samples(setdiff(sample_names(ps), empty_samples), ps)
mapping <- as(sample_data(ps), "data.frame")
ps.relative <- transform_sample_counts(ps, function(x) x / sum(x) )

## read counts (prior to contaminant filtering and rarefaction)
df <- melt(sample_sums(ps)); df$SampleID <- rownames(df); df <- df[order(df$value),]; df$SampleID <- factor(df$SampleID, levels=df$SampleID)
df$Group <- mapping[as.character(df$SampleID), "SampleType"]
agg <- aggregate(value ~ Group, df, median); agg <- agg[order(agg$value),]; df$Group <- factor(df$Group, levels=agg$Group)
p <- ggplot(df, aes(x=Group, y=value)) + geom_boxplot() + theme_classic() + ggtitle("Read counts before contaminant SV removal")
print(p)

## order by PC1 (Bray-Curtis)
ordi <- ordinate(ps.relative, method = "PCoA", distance = "bray")
ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))

### PCoA
## ordination - all samples
#for (distance_metric in distance_metrics) {
#	ordi <- ordinate(ps.relative, method = "PCoA", distance = distance_metric)
#	p <- plot_ordination(ps.relative, ordi, "samples", color = "SampleType") + theme_classic() + ggtitle(sprintf("%s (%s, %s)", "unfiltered", "SampleType", distance_metric))
#	print(p)
#}

## overall taxa barplots
ps.sel <- ps.relative
otu.filt <- as.data.frame(otu_table(ps.sel))
otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level="Genus")
otu.filt$Genus[which(is.na(otu.filt$Genus))] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="")] <- "Unknown"
agg <- aggregate(. ~ Genus, otu.filt, sum)
genera <- agg$Genus
agg <- agg[,-1]
agg <- sweep(agg, 2, colSums(agg), "/")
inds_to_keep <- which(rowSums(agg >= nsamps_threshold) >= ceiling(ncol(agg)*0.01)) # manually replace filt_threshold with 0.01 (to include blank taxa)
genera[setdiff(1:nrow(agg), inds_to_keep)] <- "Other"
agg$Genus <- genera
df <- melt(agg, variable.name="SampleID")
agg <- aggregate(value~Genus+SampleID, df, sum)
agg$SampleID <- as.character(agg$SampleID)
agg$SampleIDfactor <- factor(agg$SampleID, levels=ordering.pc1)
agg$SampleType <- mapping.sel[agg$SampleID, "SampleType"]
genera_to_add <- setdiff(agg$Genus, names(coloring))
coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
coloring.full <- c(coloring, coloring_to_add)
coloring.full.sel <- intersect(names(coloring.full), agg$Genus); coloring.full <- coloring.full[coloring.full.sel]
p <- ggplot(agg, aes(x=SampleIDfactor, y=value, fill=Genus, order=Genus)) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=2)) + ggtitle(sprintf("%s taxa summary (L5, ordered by PC1)", "Overall")) + scale_fill_manual(values=coloring.full, drop=T) + ylim(c(-.1, 1.01))
print(p)
for (st in levels(mapping.sel$SampleType)) {
	agg.sel <- subset(agg, SampleType==st)
	coloring.sel <- coloring.full[intersect(names(coloring.full), agg.sel$Genus)]
	p <- ggplot(agg.sel, aes(x=SampleIDfactor, y=value, fill=Genus, order=Genus)) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=2)) + ggtitle(sprintf("%s taxa summary (L5, ordered by PC1)", st)) + scale_fill_manual(values=coloring.sel, drop=T) + ylim(c(-.1, 1.01))
	print(p)
}

### overall taxa barplots by SampleType
#for (st in levels(mapping.sel$SampleType)) {
#	ps.sel <- subset_samples(ps.relative, SampleType==st)
#	otu.filt <- as.data.frame(otu_table(ps.sel))
#	otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level="Genus")
#	otu.filt$Genus[which(is.na(otu.filt$Genus))] <- "unknown"
#	otu.filt$Genus[which(otu.filt$Genus=="")] <- "unknown"
#	agg <- aggregate(. ~ Genus, otu.filt, sum)
#	genera <- agg$Genus
#	agg <- agg[,-1]
#	agg <- sweep(agg, 2, colSums(agg), "/")
#	inds_to_keep <- which(rowSums(agg >= nsamps_threshold) >= ceiling(ncol(agg)*filt_threshold)) # manually replace filt_threshold with 0.01 (to include blank taxa)
#	genera[setdiff(1:nrow(agg), inds_to_keep)] <- "Other"
#	agg$Genus <- genera
#	df <- melt(agg, variable.name="SampleID")
#	agg <- aggregate(value~Genus+SampleID, df, sum)
#	agg$SampleID <- as.character(agg$SampleID)
#	agg$SampleIDfactor <- factor(agg$SampleID, levels=ordering.pc1)
#	agg$SampleType <- mapping.sel[agg$SampleID, "SampleType"]
#	genera_to_add <- setdiff(agg$Genus, names(coloring))
#	coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
#	coloring.full <- c(coloring, coloring_to_add)
#	coloring.sel <- coloring.full[intersect(names(coloring.full), agg.sel$Genus)]
#	p <- ggplot(agg, aes(x=SampleIDfactor, y=value, fill=Genus, order=Genus)) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=2)) + ggtitle(sprintf("%s taxa summary (L5, ordered by PC1)", st)) + scale_fill_manual(values=coloring.sel, drop=T) + ylim(c(-.1, 1.01))
#	print(p)
#}

#########################################################################################
## trim to true samples and remove manually defined contaminant taxa (Prauserella, Rubrobacter, Unknown)
ps <- subset_samples(ps, SampleType=="BMK")
contam_list <- c("Prauserella", "Rubrobacter", "Unknown")
otu.filt <- as.data.frame(otu_table(ps))
genera <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps), level="Genus")
genera[which(is.na(genera))] <- "Unknown"
genera[which(genera=="")] <- "Unknown"
ps <- prune_taxa(!(genera %in% contam_list), ps)
ps.relative <- transform_sample_counts(ps, function(x) x / sum(x) )
ps.rarefied <- rarefy_even_depth(ps, sample.size=1001, rngseed=nsamples(ps))
ps <- prune_samples(sample_names(ps.rarefied), ps)
ps.relative <- prune_samples(sample_names(ps.rarefied), ps.relative)
mapping <- as(sample_data(ps), "data.frame")

## distance matrices
dm <- list()
for (distance_metric in distance_metrics) {
  dm[[length(dm)+1]] <- as.matrix(phyloseq::distance(ps.relative, method=distance_metric))
}
names(dm) <- distance_metrics

## GMPR normalization
sizeFactors <- GMPR(as.matrix(as.data.frame(otu_table(ps))), intersect.no=5)

##########################################################################################
## IPTW
#tmp <- mapping.sel; tmp$CaseControl <- as.numeric(tmp$CaseControl)
#psmod <- ps(CaseControl ~ child_age_final + child_sex + pool + bcg0 + pentayn + pneumoyn + measles + yf + ipv + routineoraldoses + siadoses + malnutrition, stop.method = c("es.mean"), data=tmp, estimand = "ATE", verbose=FALSE, n.trees=10000)
## diagnostic plots
#plot(psmod,plots=1)
#plot(psmod,plots=2)
#plot(psmod,plots=3)
#plot(psmod,plots=3,pairwiseMax = FALSE)
#plot(psmod,plots=4)   
## get weights
#dfx <- mapping.sel
#dfx$wt <- get.weights(psmod, stop.method = "es.mean")
#dfwt <- svydesign(ids = ~1, data=dfx, weights= ~wt)
#psgbm <- data.frame(SampleID = rownames(psmod$data), ps=psmod$ps[,1])
#psgbm$SampleID <- as.character(psgbm$SampleID); rownames(psgbm) <- psgbm$SampleID
#psgbm$sabin_any <- mapping.sel[psgbm$SampleID, "sabin_any"]
#t <- prop.table(table(psgbm$sabin_any))
#psgbm$ipw <- as.numeric(t[as.character(psgbm$sabin_any)]) / psgbm$ps
#mapping.sel$ipw <- psgbm[rownames(mapping.sel), "ipw"]

## store weights for ranger's case.weights argument
#tab <- table(mapping.sel$CaseControl)
#wt <- 1-tab/sum(tab)
#mapping.sel$ipw <- wt[as.character(mapping.sel$CaseControl)]
#sample_data(ps) <- mapping.sel
#sample_data(ps.relative) <- mapping.sel
#sample_data(ps.rarefied) <- mapping.sel[sample_names(ps.rarefied),]


##################################################################################
### barplots for each participant to check QC prior to rarefaction

## order by PC1 (Bray-Curtis)
ordi <- ordinate(ps.relative, method = "PCoA", distance = "bray")
ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))
mapping.sel <- as(sample_data(ps.relative), "data.frame")
mapping.sel$read_count <- sample_sums(ps)[rownames(mapping.sel)]

otu.filt <- as.data.frame(otu_table(ps.relative))
otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.relative), level="Genus")
otu.filt$Genus[which(is.na(otu.filt$Genus))] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="")] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="uncultured")] <- "Unknown"
otu.filt$Genus <- gsub("_\\d$", "", otu.filt$Genus)
agg <- aggregate(. ~ Genus, otu.filt, sum)
genera <- agg$Genus
agg <- agg[,-1]
agg <- sweep(agg, 2, colSums(agg), "/")
# for each participant (order by min read count)
tmp <- aggregate(read_count ~ Patient.ID, mapping.sel, min)
tmp <- tmp[order(tmp$read_count),]
for (pid in as.character(tmp$Patient.ID)) {
	sel <- rownames(subset(mapping.sel, Patient.ID==pid))
	agg.sel <- agg[, sel, drop=F]
	genera.sel <- genera
#	genera[which(rowMeans(agg.sel)<0.01)] <- "Other"
	inds_to_keep <- which(rowSums(agg.sel >= nsamps_threshold) >= ceiling(ncol(agg.sel)*filt_threshold))
	genera.sel[setdiff(1:nrow(agg.sel), inds_to_keep)] <- "Other"
	agg.sel$Genus <- genera.sel
	df <- melt(agg.sel, variable.name="SampleID")
	df2 <- aggregate(as.formula("value ~ Genus + SampleID"), df, sum)
	df2$SampleID <- as.character(df2$SampleID)
	df2$SampleID2 <- sprintf("%s (%d)", df2$SampleID, mapping.sel[df2$SampleID, "read_count"])
	ord <- order(as.Date(as.character(mapping.sel[df2$SampleID, "TimePoint.Date"]), format="%Y%m%d"))
	df2 <- df2[ord,]
	df2$SampleID2 <- factor(df2$SampleID2, levels=unique(df2$SampleID2))
	genera_to_add <- setdiff(df2$Genus, names(coloring))
	coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
	coloring.full <- c(coloring, coloring_to_add)
	coloring.full.sel <- intersect(names(coloring.full), df2$Genus); coloring.full <- coloring.full[coloring.full.sel]; coloring.full <- coloring.full[order(names(coloring.full))]
	p <- ggplot(df2, aes_string(x="SampleID2", y="value", fill="Genus", order="Genus")) + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.0)) + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=0, vjust=0.5, hjust=0.5, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots (%s)", pid)) + guides(col = guide_legend(ncol = 3))
	print(p)
}


## percent classified
taxa_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
classified <- !is.na(tax_table(ps.relative))
pct_classified <- {}
for (sid in rownames(mapping)) {
	tmp <- as.data.frame(otu_table(ps.relative)[,sid]); colnames(tmp) <- c("value")
	for (taxa_level in taxa_levels) {
		tmp[, taxa_level] <- classified[rownames(tmp), taxa_level]
	}
	df <- as.data.frame(do.call(rbind, lapply(taxa_levels, function(taxa_level) {
		inds <- which(tmp[, taxa_level])
		c(taxa_level, sum(tmp[inds, "value"]), 1-sum(tmp[inds, "value"]))
	})))
	colnames(df) <- c("Level", "Classified", "Unclassified")
	df$Classified <- 100*as.numeric(as.character(df$Classified)); df$Unclassified <- 100*as.numeric(as.character(df$Unclassified))
	df$SampleID <- sid
	pct_classified <- rbind(pct_classified, df)
}
p <- ggplot(pct_classified, aes(x=Level, y=Classified)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("%% classified by taxa level"))
print(p)

## Species percent classified by Genus
tt <- as.data.frame(tax_table(ps.relative)@.Data)
classified <- !is.na(tt)
res <- {}
for (g in levels(tt$Genus)) {
	inds <- which(tt$Genus==g)
	tmp <- rowSums(as(otu_table(ps.relative)[inds,], "matrix"))
	a <- ifelse(any(classified[inds, 7]), tmp[classified[inds, 7]] / sum(tmp), 0)
	meanabund <- mean(colSums(as(otu_table(ps.relative)[inds,], "matrix")))
	res <- rbind(res, c(g, a, meanabund))
}
res <- as.data.frame(res); colnames(res) <- c("Genus", "Classified_to_Species", "Mean_Abundance")
res$Classified_to_Species <- 100*as.numeric(as.character(res$Classified_to_Species))
write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs610/phyloseq/Classified_to_Species.%s.txt", "by_Genus"), quote=F, sep="\t", row.names=F, col.names=T)


##################################################################################
##################################################################################
### Aim 1
##################################################################################
##################################################################################
psaim <- subset_samples(ps, Aim1=="Yes")
psaim.relative <- subset_samples(ps.relative, Aim1=="Yes")
psaim.rarefied <- subset_samples(ps.rarefied, Aim1=="Yes")
aim <- "Aim1"
mvars.aim <- c("Transmitter", "Timepoint")
mapping.sel <- as(sample_data(psaim.relative), "data.frame")

## TableOne
demo_vars <- c("country", "instn", "cd4bl", "log10vlrna")
demo <- unique(mapping.sel[,c("Patient.ID", "Transmitter", demo_vars)])
write.table(print(CreateTableOne(vars=demo_vars, strata=c("Transmitter"), data=demo, smd=T)), file=sprintf("%s/Table_1.%s.by_Subject.txt", output_dir, aim), quote=F, sep="\t", row.names=T, col.names=T)
sample_vars <- c("age")
write.table(print(CreateTableOne(vars=sample_vars, strata=c("Timepoint", "Transmitter"), data=mapping.sel, smd=T)), file=sprintf("%s/Table_1.%s.by_Sample.txt", output_dir, aim), quote=F, sep="\t", row.names=T, col.names=T)
for (tp in levels(mapping.sel$Timepoint)) {
	write.table(print(CreateTableOne(vars=sample_vars, strata=c("Transmitter"), data=subset(mapping.sel, Timepoint==tp), smd=T)), file=sprintf("%s/Table_1.%s.by_Sample.%s.txt", output_dir, aim, tp), quote=F, sep="\t", row.names=T, col.names=T)
}

##################################################################################
### barplots by var

## order by PC1 (Bray-Curtis)
ordi <- ordinate(psaim.relative, method = "PCoA", distance = "bray")
ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))

for (mvar in rownames(subset(metadata_variables, useForPERMANOVA=="yes"))){
  otu.filt <- as.data.frame(otu_table(psaim.relative))
  otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.relative), level="Genus")
  otu.filt$Genus[which(is.na(otu.filt$Genus))] <- "Unknown"
  otu.filt$Genus[which(otu.filt$Genus=="")] <- "Unknown"
  otu.filt$Genus[which(otu.filt$Genus=="uncultured")] <- "Unknown"
  otu.filt$Genus <- gsub("_\\d$", "", otu.filt$Genus)
  agg <- aggregate(. ~ Genus, otu.filt, sum)
  genera <- agg$Genus
  agg <- agg[,-1]
  agg <- sweep(agg, 2, colSums(agg), "/")
  inds_to_keep <- which(rowSums(agg >= nsamps_threshold) >= ceiling(ncol(agg)*filt_threshold))
  genera[setdiff(1:nrow(agg), inds_to_keep)] <- "Other"
  agg$Genus <- genera
  df <- melt(agg, variable.name="SampleID")
  df2 <- aggregate(as.formula("value ~ Genus + SampleID"), df, sum)
  df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
  df2[[mvar]] <- mapping.sel[df2$SampleID, mvar]
  genera_to_add <- setdiff(df2$Genus, names(coloring))
	coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
	coloring.full <- c(coloring, coloring_to_add)
	coloring.full.sel <- intersect(names(coloring.full), df2$Genus); coloring.full <- coloring.full[coloring.full.sel]; coloring.full <- coloring.full[order(names(coloring.full))]
  p <- ggplot(df2, aes_string(x="SampleIDfactor", y="value", fill="Genus", order="Genus")) + facet_wrap(as.formula(sprintf("~%s", mvar)), scales="free_x") + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.0)) + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots - %s (Aim 1)", mvar)) + guides(col = guide_legend(ncol = 3))
  print(p)
}
# faceted by Transmitter+Timepoint
df2 <- aggregate(as.formula("value ~ Genus + SampleID"), df, sum)
df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
for (mvar in mvars.aim) {
	df2[[mvar]] <- mapping.sel[df2$SampleID, mvar]
}
genera_to_add <- setdiff(df2$Genus, names(coloring))
coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
coloring.full <- c(coloring, coloring_to_add)
coloring.full.sel <- intersect(names(coloring.full), df2$Genus); coloring.full <- coloring.full[coloring.full.sel]; coloring.full <- coloring.full[order(names(coloring.full))]
p <- ggplot(df2, aes_string(x="SampleIDfactor", y="value", fill="Genus", order="Genus")) + facet_wrap(as.formula(sprintf("~%s+%s", "Transmitter", "Timepoint")), scales="free_x") + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.0)) + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots - %s (Aim 1)", "Transmitter+Timepoint")) + guides(col = guide_legend(ncol = 3))
print(p)

# aggregated by Transmitter+Timepoint
otu.filt <- as.data.frame(otu_table(psaim.relative))
otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.relative), level="Genus")
otu.filt$Genus[which(is.na(otu.filt$Genus))] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="")] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="uncultured")] <- "Unknown"
otu.filt$Genus <- gsub("_\\d$", "", otu.filt$Genus)
agg <- aggregate(. ~ Genus, otu.filt, sum)
genera <- agg$Genus
agg <- agg[,-1]
agg <- sweep(agg, 2, colSums(agg), "/")
inds_to_keep <- which(rowSums(agg >= nsamps_threshold) >= ceiling(ncol(agg)*filt_threshold))
genera[setdiff(1:nrow(agg), inds_to_keep)] <- "Other"
agg$Genus <- genera
agg <- aggregate(. ~ Genus, agg, sum)
df <- melt(agg, variable.name="SampleID")
for (mvar in mvars.aim) {
	df[, mvar] <- mapping.sel[as.character(df$SampleID), mvar]
}
df2 <- aggregate(as.formula(sprintf("value ~ Genus + Transmitter + Timepoint")), df, mean)
genera_to_add <- setdiff(df2$Genus, names(coloring))
coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
coloring.full <- c(coloring, coloring_to_add)
coloring.full.sel <- intersect(names(coloring.full), df2$Genus); coloring.full <- coloring.full[coloring.full.sel]; coloring.full <- coloring.full[order(names(coloring.full))]
counts <- table(mapping.sel[, c(mvars.aim)])
df2[, "mvar"] <- sprintf("%s-%s\n(n=%d)", df2[, "Transmitter"], df2[, "Timepoint"], counts[cbind(as.character(df2[,"Transmitter"]), as.character(df2[,"Timepoint"]))])
p <- ggplot(df2, aes(x=Timepoint, y=value, fill=Genus)) + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.05)) + facet_wrap(~Transmitter, scales="free") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots (Aim 1)")) + guides(col = guide_legend(ncol = 3))
print(p)
p <- ggplot(df2, aes_string(x="mvar", y="value", fill="Genus", order="Genus")) + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.01)) + facet_wrap(~Transmitter, scales="free") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots (Aim 1)")) + guides(col = guide_legend(ncol = 3))
print(p)

# per-taxon violin plots over time
# aggregated by Transmitter+Timepoint
otu.filt <- as.data.frame(otu_table(psaim.relative))
otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.relative), level="Genus")
otu.filt$Genus[which(is.na(otu.filt$Genus))] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="")] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="uncultured")] <- "Unknown"
otu.filt$Genus <- gsub("_\\d$", "", otu.filt$Genus)
agg <- aggregate(. ~ Genus, otu.filt, sum)
genera <- agg$Genus
agg <- agg[,-1]
agg <- sweep(agg, 2, colSums(agg), "/")
inds_to_keep <- which(rowSums(agg >= nsamps_threshold) >= ceiling(ncol(agg)*filt_threshold))
genera[setdiff(1:nrow(agg), inds_to_keep)] <- "Other"
agg$Genus <- genera
agg <- aggregate(. ~ Genus, agg, sum)
df <- melt(agg, variable.name="SampleID")
for (mvar in mvars.aim) {
	df[, mvar] <- mapping.sel[as.character(df$SampleID), mvar]
}
for (genus in unique(df$Genus)) {
	p <- ggplot(subset(df, Genus==genus), aes(x=Timepoint, y=value, color=Transmitter)) + geom_violin() + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_color_manual(values=cols.mvar[["Transmitter"]]) + ggtitle(sprintf("%s over Timepoint (Aim 1)", genus))
	print(p)
}


##################################################################################
### PCoA + PERMANOVA
psaim.sel <- psaim.relative
mapping.sel <- as(sample_data(psaim.sel), "data.frame")
sel <- rownames(mapping.sel)

## PCoA
for (mvar in rownames(subset(metadata_variables, useForPERMANOVA=="yes"))){
  for (distance_metric in distance_metrics) {
    ordi <- ordinate(psaim.sel, method = "PCoA", distance = distance_metric)
		if (metadata_variables[mvar, "type"] == "factor") {
			bd <- betadisper(as.dist(dm[[distance_metric]][sel,sel]), mapping.sel[,mvar])
			bd$centroids <- rbind(bd$centroids, rep(NA, ncol(bd$centroids))); rownames(bd$centroids)[nrow(bd$centroids)] <- "NA" # pad with NA row
			df.vectors <- as.data.frame(ordi$vectors[,1:2]); df.vectors[,mvar] <- mapping.sel[rownames(df.vectors),mvar]
			df.vectors$PC1origin <- bd$centroids[str_replace_na(as.character(df.vectors[, mvar])), "PCoA1"]
			df.vectors$PC2origin <- bd$centroids[str_replace_na(as.character(df.vectors[, mvar])), "PCoA2"]
			df.centroids <- as.data.frame(bd$centroids[1:(nrow(bd$centroids)-1), 1:2, drop=F]); df.centroids[,mvar] <- rownames(df.centroids) # drop the NA row
		  p <- plot_ordination(psaim.sel, ordi, "samples", color = mvar) + theme_classic() + ggtitle(distance_metric) + theme_classic() + stat_ellipse(geom="polygon", alpha=0.1, type="t", aes_string(fill=mvar)) + geom_segment(data=df.vectors, aes(x=PC1origin, y=PC2origin, xend=Axis.1, yend=Axis.2), inherit.aes=F, color="grey", arrow = arrow(length = unit(0.05, "inches")), size=0.5) + geom_point(data=df.centroids, aes_string(x="PCoA1", y="PCoA2", color=mvar), size=4, inherit.aes=F)
		  print(p)
		 } else {
		 	p <- plot_ordination(psaim.sel, ordi, "samples", color = mvar) + theme_classic() + ggtitle(distance_metric) + theme_classic()
		  print(p)
		 }
  }
}
# by Transmitter+Timepoint
p <- plot_ordination(psaim.sel, ordi, "samples", color = "Timepoint", shape = "Transmitter") + theme_classic() + ggtitle(distance_metric) + theme_classic() + stat_ellipse(geom="polygon", alpha=0.05, type="t", aes_string(fill="Timepoint", lty="Transmitter")) + scale_shape_manual(values=c(20,3))
print(p)


## PERMANOVA
out_txt <- sprintf("%s/PERMANOVA.%s.%s.txt", output_dir, aim, format(Sys.Date(), "%m%d%y"))
sink(out_txt, append=F)
print(sprintf("PERMANOVA (%s)", aim))
mvars <- c("Batch", "Timepoint", "Transmitter", "cd4bl", "log10vlrna")
for (distance_metric in distance_metrics) {
  print(distance_metric)
  mapping.sel <- subset(as(sample_data(psaim.sel), "data.frame"), Aim1=="Yes")
  ids_to_remove <- names(which(unlist(apply(mapping.sel[, mvars, drop=F], 1, function(x) any(is.na(x))))))
  nonna <- setdiff(rownames(mapping.sel), ids_to_remove)
  dm.sel <- dm[[distance_metric]][nonna, nonna]
  mapping.nonna <- mapping.sel[nonna,]
  form <- as.formula(sprintf("as.dist(dm.sel) ~ %s", paste(mvars, collapse="+")))
  res <- adonis2(form , data=mapping.nonna, permutations=999, by="margin")
  res$R2 <- res$SumOfSqs / sum(res$SumOfSqs)
  print(res)
}
sink()



##################################################################################
## alpha diversity (using emmeans)
psaim.sel <- psaim.rarefied
mapping.sel <- as(sample_data(psaim.sel), "data.frame")
sel <- rownames(mapping.sel)
adiv <- estimate_richness(psaim.sel, measures=alpha_metrics)
#rownames(adiv) <- gsub("^X", "", rownames(adiv))
adiv$SampleID <- rownames(adiv)
adiv <- merge(adiv, mapping.sel, by="row.names"); rownames(adiv) <- adiv$SampleID
res <- {}
mvars <- c("Transmitter")
for (mvar in mvars) {
	adiv.sel <- adiv[!is.na(adiv[,mvar]),] # remove NAs
	# boxplots/violin plots
	df <- melt(adiv.sel[,c(alpha_metrics, mvar)])
	p <- ggplot(df, aes_string(x=mvar, y="value", color=mvar)) + geom_boxplot() + facet_wrap(~variable, scales="free") + theme_classic() + ggtitle(sprintf("alpha diversity by %s (%s)", mvar, aim))
	print(p)
	for (alpha_metric in alpha_metrics) {
		mod <- lm(as.formula(sprintf("%s ~ %s*Timepoint", alpha_metric, mvar)), data=adiv.sel)
		emm.adiv <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s | Timepoint", mvar)), adjust="none")
		tmp <- as.data.frame(emm.adiv$contrasts)
		tmp$alpha_metric <- alpha_metric; tmp$metadata_variable <- mvar
		res <- rbind(res, tmp)
	}
}
res$padj <- p.adjust(res$p.value, method="fdr")
res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
res <- res[order(res$estimate, decreasing=T),]
for (mvar in mvars) {
	df <- subset(res, metadata_variable==mvar)
	contrs <- levels(droplevels(df$contrast))
	# single level contrast
	if (length(contrs)==1) {
		lims <- max(abs(df$estimate) + abs(df$SE))*1.0
		pd <- position_dodge(0.4)
		comp <- contrs[1]
		p <- ggplot(df, aes(x=Timepoint, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=Timepoint, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_hline(yintercept=0) + facet_wrap(~alpha_metric, scales="free") + theme_classic() + ggtitle(sprintf("%s (%s %s, IPTW)", "alpha diversity", mvar, comp)) + coord_flip() + scale_color_manual(values=dircolors)
		print(p)
	}
}
write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs610/phyloseq/adiv.%s.txt", aim), quote=F, sep="\t", row.names=F, col.names=T)


##################################################################################
### differential abundance - relative abundance with lmer
for (level in c("Genus", "Species")) {
	otu.filt <- as.data.frame(otu_table(psaim.relative))
	otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.relative), level=level)
	# rename Prevotella_6, etc -> Prevotella
	otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
	agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
	if (level %in% c("Species", "Genus")){
		agg <- agg[-1,]
	}
	lvl <- agg[[level]]
	agg <- agg[,-1]
	rownames(agg) <- lvl
	ftk <- names(which(rowSums(agg >= thresholds[["nsamps"]][[level]]) >= ceiling(ncol(agg)*thresholds[["filt"]][[level]])))
	agg <- agg[ftk,]
	agg[[level]] <- rownames(agg)
	res <- {}
	out <- mclapply(agg[[level]], function(f) {
		df <- melt(agg[f,]); colnames(df) <- c(level, "SampleID", "value"); df$SampleID <- as.character(df$SampleID)
		tmp <- {}
		print(f)
		nsamps_detected <- length(which(df$value>=thresholds[["nsamps"]][[level]]))
		for (mvar in c("Transmitter")) {
		  #print(sprintf("%s %s", f, mvar))
		  df2 <- df
		  for (m in c(mvar, "Timepoint", "Patient.ID")) {
		  	df2[, m] <- mapping.sel[df2$SampleID, m]
		  }
		  df2 <- subset(df2, !is.na(df2[,mvar,drop=F]))
			df2[, mvar] <- factor(as.character(df2[,mvar]), levels=rev(levels(df2[,mvar]))) # reverse levels to get more intuitive contrast
		  mod <- lmer(as.formula(sprintf("%s ~ %s*Timepoint + (1 | Patient.ID)", "value", mvar)), data=df2); modelstr <- "LMEM"
		  # contrast Transmitter, stratify by Timepoint
			emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s | Timepoint", mvar)), adjust="none")
			tmp2 <- as.data.frame(emm$contrasts); colnames(tmp2)[2] <- "Strata"
			tmp2[, level] <- f; tmp2$metadata_variable <- mvar; tmp2$model <- modelstr
			tmp <- rbind(tmp, tmp2)
			# contrast Timepoint, stratify by Transmitter
			emm <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ Timepoint | %s", mvar)), adjust="none")
			tmp2 <- as.data.frame(emm$contrasts); colnames(tmp2)[2] <- "Strata"
			tmp2[, level] <- f; tmp2$metadata_variable <- "Timepoint"; tmp2$model <- modelstr
			tmp <- rbind(tmp, tmp2)
		}
		print(sprintf("finished %s", f))
		tmp
	}, mc.cores=16)
	res <- as.data.frame(do.call(rbind, out))
	res <- res[,c(level, setdiff(colnames(res), level))]
	colnames(res) <- c(level, "contrast", "Strata", "Estimate", "SE", "df", "t.ratio", "p.value", "metadata_variable", "model")
	res$padj <- p.adjust(res$p.value, method="fdr")
	res <- res[order(res$p.value, decreasing=F),]
	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$Estimate)==1, "up", "down"), "NS")
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs610/phyloseq/emmeans.%s.%s.txt", aim, level), quote=F, sep="\t", row.names=F, col.names=T)
	# forest plot by Transmitter
	mvar <- "Transmitter"
	df <- subset(res, metadata_variable==mvar)
	df <- df[order(df$Estimate, decreasing=T),]
	df[, level] <- factor(as.character(df[, level]), levels=unique(as.character(df[, level])))
	df$contrast <- droplevels(df$contrast)
	lims <- max(abs(as.numeric(as.character(df$Estimate))) + abs(as.numeric(as.character(df$SE))))*1.0
	pd <- position_dodge(0.8)
	p <- ggplot(df, aes_string(x=level, y="Estimate", color="dir", group="Strata")) + geom_point(position=pd) + geom_errorbar(aes_string(x=level, ymin="Estimate-SE", max="Estimate+SE"), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Strata), position=pd, hjust=1, color="black", size=1.5) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LMEM: %s (%s, %s)", levels(df$contrast), mvar, aim)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
	print(p)
	# forest plot by Timepoint
	mvar <- "Timepoint"
	df <- subset(res, metadata_variable==mvar)
	df <- df[order(df$Estimate, decreasing=T),]
	df[, level] <- factor(as.character(df[, level]), levels=unique(as.character(df[, level])))
	df$contrast <- droplevels(df$contrast)
	lims <- max(abs(as.numeric(as.character(df$Estimate))) + abs(as.numeric(as.character(df$SE))))*1.0
	pd <- position_dodge(0.8)
	p <- ggplot(df, aes_string(x=level, y="Estimate", color="dir", group="Strata")) + geom_point(position=pd) + geom_errorbar(aes_string(x=level, ymin="Estimate-SE", max="Estimate+SE"), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Strata), position=pd, hjust=1, color="black", size=1.5) + facet_wrap(~contrast, scales="free") + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LMEM: %s (%s, %s)", levels(df$contrast), mvar, aim)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
	print(p)
	# violin plots
	agg.melt <- melt(agg)
	colnames(agg.melt) <- c(level, "SampleID", "value"); agg.melt$SampleID <- as.character(agg.melt$SampleID)
	for (mvar in mvars.aim) {
		agg.melt[, mvar] <- mapping.sel[agg.melt$SampleID, mvar]
	}
	mvar <- "Transmitter"
	p <- ggplot(agg.melt, aes_string(x="Timepoint", y="value", color=mvar)) + geom_violin(position=pd) + geom_point(position=pd) + facet_wrap(as.formula(sprintf("~%s", level)), scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("relabund (%s, %s)", mvar, aim)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
	print(p)
	p <- ggplot(agg.melt, aes_string(x="Timepoint", y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(as.formula(sprintf("~%s", level)), scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("relabund (%s, %s)", mvar, aim)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
	npages <- n_pages(p)
	for (ip in 1:npages) {
		p <- ggplot(agg.melt, aes_string(x="Timepoint", y="value", color=mvar)) + geom_violin(position=pd) + geom_point(position=pd) + facet_wrap_paginate(as.formula(sprintf("~%s", level)), scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("relabund (%s, %s)", mvar, aim)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
	}
}


##################################################################################
### random forest (Transmitter, separately for each timepoint)
psaim <- psaim.relative
set.seed(nsamples(psaim))
mvar <- "Transmitter"
res.mean <- list()
res.sd <- list()
for (level in c("Genus", "Species")) {
	for (tp in levels(mapping$Timepoint)) {
		psaim.sel <- subset_samples(psaim, Timepoint==tp)
		mapping.sel <- as(sample_data(psaim.sel), "data.frame")
		tab <- table(mapping.sel[,mvar])
		wt <- 1-tab/sum(tab)
		mapping.sel$ipw <- wt[as.character(mapping.sel[,mvar])]
		
		otu.filt <- as.data.frame(otu_table(psaim.sel))
		otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.sel), level=level)
		# rename Prevotella_6, etc -> Prevotella
		otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
		otu.filt <- otu.filt[which(otu.filt[[level]]!=""),]
		agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
		lvl <- agg[[level]]
		agg <- agg[,-1]
		rownames(agg) <- lvl
		agg <- agg[which(rowSums(agg >= thresholds[["nsamps"]][[level]]) >= ceiling(ncol(agg)*thresholds[["filt"]][[level]])),]
		
		sel <- colnames(agg)
		data.sel <- as.data.frame(t(agg[,sel]))
		data.sel <- as.matrix(data.sel)
		response <- droplevels(mapping.sel[sel, mvar]); names(response) <- sel
		# subset to non-NA
		response <- subset(response, !is.na(response))
		data.sel <- data.sel[names(response),]
		agg.melt.stored <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "taxa", "value")

#		## after running for the first time, COMMENT OUT THIS BLOCK ##
#		num_iter <- 1000
#		ncores <- 20
#		ranger.seeds <- sample(1:num_iter, num_iter, replace=T) # set up a vector of seeds for the ranger mclapply
#		out <- mclapply(1:num_iter, function (dummy) {
#				importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", case.weights=mapping.sel[rownames(data.sel), "ipw"], seed=ranger.seeds[dummy], num.threads=1))
#		}, mc.cores=ncores )
#	#	out <- mclapply(1:num_iter, function (dummy) {
#	#			importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", sample.fraction=c(1, 0.25), seed=ranger.seeds[dummy], num.threads=1))
#	#	}, mc.cores=ncores )	
#		collated.importance <- do.call(cbind, out)
#		out <- mclapply(1:num_iter, function (dummy) {
#				rgcv2(trainx=data.sel, trainy=response, cv.fold=10, step=0.9, num.threads=1)$error.cv
#			}, mc.cores=ncores )
#		collated.cv <- do.call(cbind, out)

#		write.table(collated.importance, file=sprintf("%s/ranger.%s.%s.%s.%s.importance.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)
#		write.table(collated.cv, file=sprintf("%s/ranger.%s.%s.%s.%s.cv.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)
#		## END BLOCK TO COMMENT ##

		collated.importance <- read.table(sprintf("%s/ranger.%s.%s.%s.%s.importance.txt", output_dir, aim, mvar, tp, level), header=F, as.is=T, sep="\t", row.names=1, quote="")
		collated.cv <- read.table(sprintf("%s/ranger.%s.%s.%s.%s.cv.txt", output_dir, aim, mvar, tp, level), header=F, as.is=T, sep="\t", row.names=1)
		importance.mean <- rowMeans(collated.importance)
		importance.sd <- unlist(apply(collated.importance, 1, sd))
		cv.mean <- rowMeans(collated.cv)
		cv.sd <- unlist(apply(collated.cv, 1, sd))
		inds <- order(importance.mean, decreasing=T)
		inds <- inds[1:min(20, as.numeric(names(which.min(cv.mean))))] # edit as appropriate
		write.table(melt(importance.mean[inds]), file=sprintf("%s/ranger.%s.%s.%s.%s.features.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)

		## after running for the first time, COMMENT OUT THIS BLOCK ##
		# using a sparse model with N predictors
#		sparseRanger <- ranger(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, num.trees=10000, importance="permutation", case.weights=mapping.sel[rownames(data.sel), "ipw"], seed=sample(1:num_iter,1))
#		save(sparseRanger, file=sprintf("%s/ranger.%s.%s.%s.%s.model", output_dir, aim, mvar, tp, level))
		load(sprintf("%s/ranger.%s.%s.%s.%s.model", output_dir, aim, mvar, tp, level))
		# accuracy of final sparseRF model
		pred <- predictions(sparseRanger)
		pred_df <- data.frame(SampleID=names(response), predicted=pred, true=response, stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
		pred_df_out <- merge(pred_df, data.sel, by="row.names")
		write.table(pred_df_out, file=sprintf("%s/ranger.%s.%s.%s.%s.predictions.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=F, col.names=T)
		confusion_matrix <- table(pred_df[, c("true", "predicted")])
		class_errors <- unlist(lapply(levels(mapping.sel$sabin_any), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$sabin_any)
		accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
		vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
		mccvalue <- mcc(vec.pred, vec.true)
		df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", aim, mvar, tp, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
		print(p)

		write.table(confusion_matrix, file=sprintf("%s/ranger.%s.%s.%s.%s.confusion_matrix.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=T)
		## END BLOCK TO COMMENT ##

		# plotting - per-group sparse model
		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s, %s, %s, %s", aim, mvar, tp, level)))
		# plotting - per-group variables
		df <- data.frame(Taxa=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
		# load effect sizes from linear regression
	#	contr <- "Case - Control"
		lmres <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, aim, level), header=T, as.is=T, sep="\t", quote=""); lmres <- subset(lmres, metadata_variable==mvar & Strata==tp); colnames(lmres)[1] <- "Taxa"
		df <- merge(lmres[,c("Taxa", "contrast", "Estimate", "SE", "padj", "dir")], df, by="Taxa")
		df$Taxa <- factor(df$Taxa, levels=rev(names(importance.mean)[inds]))
		p <- ggplot(df, aes(x=Taxa, y=importance, label=Taxa)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=Taxa, y=0, label=Taxa), size=3, hjust=0) + ggtitle(sprintf("%s: %s explanatory %s", aim, mvar, tp)) + theme(axis.text.y=element_blank())
		print(p)	
		lims <- max(abs(df$Estimate) + abs(df$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
		p <- ggplot(df, aes(x=Taxa, y=Estimate, color=dir, group=contrast)) + geom_point(position=pd) + geom_errorbar(aes(x=Taxa, ymin=Estimate-SE, max=Estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=contrast), position=pd, hjust=1, color="black", size=1.5) + geom_tile(aes(x=Taxa, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s %s explanatory %s", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
		print(p)
		
		# shading rectangles of importance values
		df.rect <- df
		df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
		p <- ggplot(df.rect, aes(x=x, y=Taxa, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s %s explanatory %s", aim, mvar, tp)) + scale_fill_gradient(low="white", high="black")
		print(p)
		# violin plots of relabund values
		agg.melt <- agg.melt.stored
		agg.melt[,mvar] <- mapping.sel[agg.melt$SampleID, mvar]
		agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$sabin_any))
		agg.melt <- subset(agg.melt, taxa %in% levels(df$Taxa))
		agg.melt$taxa <- factor(agg.melt$taxa, levels=levels(df$Taxa))
		p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap(~taxa, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
		p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~taxa, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		npages <- n_pages(p)
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~taxa, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
			print(p)
		}
		p <- ggplot(agg.melt, aes_string(x="taxa", y="value", color=mvar)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
	}
}


##################################################################################
### random forest (Transmitter, combined over all timepoints)
psaim <- psaim.relative
set.seed(nsamples(psaim))
mvar <- "Transmitter"
tp <- "Combined"
res.mean <- list()
res.sd <- list()
for (level in c("Genus", "Species")) {
	psaim.sel <- psaim
	mapping.sel <- as(sample_data(psaim.sel), "data.frame")
	tab <- table(mapping.sel[,mvar])
	wt <- 1-tab/sum(tab)
	mapping.sel$ipw <- wt[as.character(mapping.sel[,mvar])]
	
	otu.filt <- as.data.frame(otu_table(psaim.sel))
	otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.sel), level=level)
	# rename Prevotella_6, etc -> Prevotella
	otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
	otu.filt <- otu.filt[which(otu.filt[[level]]!=""),]
	agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
	lvl <- agg[[level]]
	agg <- agg[,-1]
	rownames(agg) <- lvl
	agg <- agg[which(rowSums(agg >= thresholds[["nsamps"]][[level]]) >= ceiling(ncol(agg)*thresholds[["filt"]][[level]])),]
	
	sel <- colnames(agg)
	data.sel <- as.data.frame(t(agg[,sel]))
	data.sel <- as.matrix(data.sel)
	response <- droplevels(mapping.sel[sel, mvar]); names(response) <- sel
	# subset to non-NA
	response <- subset(response, !is.na(response))
	data.sel <- data.sel[names(response),]
	
	# cast into per-subject data (using mean abundance)
	data.sel <- as.data.frame(data.sel)
	data.sel$Patient.ID <- mapping.sel[rownames(data.sel), "Patient.ID"]
	data.sel$Timepoint <- mapping.sel[rownames(data.sel), "Timepoint"]
	agg <- aggregate(. ~ Patient.ID + Timepoint, data.sel, mean)
	agg2 <- melt(agg)
	agg3 <- dcast(agg2, Patient.ID ~ Timepoint+variable, value.var="value")
	data.sel <- agg3; rownames(data.sel) <- data.sel$Patient.ID; data.sel <- data.sel[,-1]
	tmp <- unique(mapping.sel[, c("Patient.ID", "Transmitter")])
	rownames(tmp) <- tmp$Patient.ID
	response <- tmp[rownames(data.sel), "Transmitter"]; names(response) <- rownames(tmp)
	# subset to subjects with data in all three timepoins to avoid NAs
	sel <- names(which(apply(data.sel, 1, function(x) !any(is.na(x)))))
	data.sel <- data.sel[sel,]
	response <- response[sel]	
	agg.melt.stored <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt.stored) <- c("Patient.ID", "taxa", "value")

#	## after running for the first time, COMMENT OUT THIS BLOCK ##
#	num_iter <- 1000
#	ncores <- 20
#	ranger.seeds <- sample(1:num_iter, num_iter, replace=T) # set up a vector of seeds for the ranger mclapply
#	out <- mclapply(1:num_iter, function (dummy) {
#			importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", case.weights=mapping.sel[rownames(data.sel), "ipw"], seed=ranger.seeds[dummy], num.threads=1))
#	}, mc.cores=ncores )
##	out <- mclapply(1:num_iter, function (dummy) {
##			importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", sample.fraction=c(1, 0.25), seed=ranger.seeds[dummy], num.threads=1))
##	}, mc.cores=ncores )	
#	collated.importance <- do.call(cbind, out)
#	out <- mclapply(1:num_iter, function (dummy) {
#			rgcv2(trainx=data.sel, trainy=response, cv.fold=10, step=0.9, num.threads=1)$error.cv
#		}, mc.cores=ncores )
#	collated.cv <- do.call(cbind, out)

#	write.table(collated.importance, file=sprintf("%s/ranger.%s.%s.%s.%s.importance.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)
#	write.table(collated.cv, file=sprintf("%s/ranger.%s.%s.%s.%s.cv.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)
#	## END BLOCK TO COMMENT ##

	collated.importance <- read.table(sprintf("%s/ranger.%s.%s.%s.%s.importance.txt", output_dir, aim, mvar, tp, level), header=F, as.is=T, sep="\t", row.names=1, quote="")
	collated.cv <- read.table(sprintf("%s/ranger.%s.%s.%s.%s.cv.txt", output_dir, aim, mvar, tp, level), header=F, as.is=T, sep="\t", row.names=1)
	importance.mean <- rowMeans(collated.importance)
	importance.sd <- unlist(apply(collated.importance, 1, sd))
	cv.mean <- rowMeans(collated.cv)
	cv.sd <- unlist(apply(collated.cv, 1, sd))
	inds <- order(importance.mean, decreasing=T)
	inds <- inds[1:min(20, as.numeric(names(which.min(cv.mean))))] # edit as appropriate
	write.table(melt(importance.mean[inds]), file=sprintf("%s/ranger.%s.%s.%s.%s.features.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)

	## after running for the first time, COMMENT OUT THIS BLOCK ##
	# using a sparse model with N predictors
#	sparseRanger <- ranger(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, num.trees=10000, importance="permutation", case.weights=mapping.sel[rownames(data.sel), "ipw"], seed=sample(1:num_iter,1))
#	save(sparseRanger, file=sprintf("%s/ranger.%s.%s.%s.%s.model", output_dir, aim, mvar, tp, level))
	load(sprintf("%s/ranger.%s.%s.%s.%s.model", output_dir, aim, mvar, tp, level))
	# accuracy of final sparseRF model
	pred <- predictions(sparseRanger)
	pred_df <- data.frame(SampleID=names(response), predicted=pred, true=response, stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
	pred_df_out <- merge(pred_df, data.sel, by="row.names")
	write.table(pred_df_out, file=sprintf("%s/ranger.%s.%s.%s.%s.predictions.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=F, col.names=T)
	confusion_matrix <- table(pred_df[, c("true", "predicted")])
	class_errors <- unlist(lapply(levels(mapping.sel$sabin_any), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$sabin_any)
	accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
	vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
	mccvalue <- mcc(vec.pred, vec.true)
	df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", aim, mvar, tp, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
	print(p)

	write.table(confusion_matrix, file=sprintf("%s/ranger.%s.%s.%s.%s.confusion_matrix.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=T)
	## END BLOCK TO COMMENT ##

	# plotting - per-group sparse model
	df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
	colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
	print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s, %s, %s, %s", aim, mvar, tp, level)))
	# plotting - per-group variables
	df <- data.frame(feature=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
	df$Timepoint <- unlist(lapply(as.character(df$feature), function(x) unlist(strsplit(x, "_"))[1]))
	df$Taxa <- unlist(lapply(as.character(df$feature), function(x) unlist(strsplit(x, "_"))[2]))
	# load effect sizes from linear regression
#	contr <- "Case - Control"
	lmres <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, aim, level), header=T, as.is=T, sep="\t", quote=""); lmres <- subset(lmres, metadata_variable==mvar); colnames(lmres)[1] <- "Taxa"; lmres$feature <- paste(lmres$Strata, lmres$Taxa, sep="_")
	df <- merge(lmres[,c("feature", "contrast", "Estimate", "SE", "padj", "dir")], df, by="feature")
	df <- df[order(df$importance, decreasing=T),]
	df$feature <- factor(df$feature, levels=rev(df$feature))
	p <- ggplot(df, aes(x=feature, y=importance, label=feature)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=feature, y=0, label=feature), size=3, hjust=0) + ggtitle(sprintf("%s: %s explanatory %s", aim, mvar, tp)) + theme(axis.text.y=element_blank())
	print(p)	
	lims <- max(abs(df$Estimate) + abs(df$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
	p <- ggplot(df, aes(x=feature, y=Estimate, color=dir, group=contrast)) + geom_point(position=pd) + geom_errorbar(aes(x=feature, ymin=Estimate-SE, max=Estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=contrast), position=pd, hjust=1, color="black", size=1.5) + geom_tile(aes(x=feature, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s %s explanatory %s", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
	print(p)	
	
	# shading rectangles of importance values
	df.rect <- df
	df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
	p <- ggplot(df.rect, aes(x=x, y=feature, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s %s explanatory %s", aim, mvar, tp)) + scale_fill_gradient(low="white", high="black")
	print(p)
	# violin plots of relabund values
	agg.melt <- agg.melt.stored; colnames(agg.melt)[2] <- "feature"
	agg.melt[,mvar] <- response[agg.melt$Patient.ID]
	agg.melt <- subset(agg.melt, feature %in% levels(df$feature))
	agg.melt$feature <- factor(agg.melt$feature, levels=levels(df$feature))
	p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap(~feature, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
	print(p)
	p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~feature, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
	npages <- n_pages(p)
	for (ip in 1:npages) {
		p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~feature, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
	}
	p <- ggplot(agg.melt, aes_string(x="feature", y="value", color=mvar)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
	print(p)
}


##################################################################################
### XGBoost (Transmitter, separately for each timepoint)
psaim <- psaim.relative
set.seed(nsamples(psaim))
mvar <- "Transmitter"
res.mean <- list()
res.sd <- list()
for (level in c("Genus", "Species")) {
	for (tp in levels(mapping$Timepoint)) {
		psaim.sel <- subset_samples(psaim, Timepoint==tp)
		mapping.sel <- as(sample_data(psaim.sel), "data.frame")
		tab <- table(mapping.sel[,mvar])
		wt <- 1-tab/sum(tab)
		mapping.sel$ipw <- wt[as.character(mapping.sel[,mvar])]
		
		otu.filt <- as.data.frame(otu_table(psaim.sel))
		otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.sel), level=level)
		# rename Prevotella_6, etc -> Prevotella
		otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
		otu.filt <- otu.filt[which(otu.filt[[level]]!=""),]
		agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
		lvl <- agg[[level]]
		agg <- agg[,-1]
		rownames(agg) <- lvl
		agg <- agg[which(rowSums(agg >= thresholds[["nsamps"]][[level]]) >= ceiling(ncol(agg)*thresholds[["filt"]][[level]])),]
		
		sel <- colnames(agg)
		data.sel <- as.data.frame(t(agg[,sel]))
		data.sel <- as.matrix(data.sel)
		response <- droplevels(mapping.sel[sel, mvar]); names(response) <- sel
		# subset to non-NA
		response <- subset(response, !is.na(response))
		data.sel <- data.sel[names(response),]
		agg.melt.stored <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "taxa", "value")
		
		
		# run XGBoost + SHAP
		labels <- as.numeric(response)-1
		negative_cases <- sum(labels == FALSE)
		positive_cases <- sum(labels == TRUE)
		bst <- xgboost(data=data.sel, label=labels, max.depth = 4, eta = 1, nthread = 2, nrounds = 10, objective = "binary:logistic", scale_pos_weight = negative_cases/postive_cases)
		pred <- predict(bst, data.sel)
		prediction <- as.numeric(pred > 0.5)

		importance_matrix <- xgb.importance(model = bst)
		xgb.ggplot.importance(importance_matrix = importance_matrix) + ggtitle(sprintf("xgb feature importance (%s, %s, %s)", mvar, level, tp))
		
		explainer <- shapr(data.sel, bst)
		p <- mean(labels)
		explanation <- explain(data.sel, approach = "empirical", explainer = explainer, prediction_zero = p)
		shp <- shapviz(explanation)
		
		for (i in seq(from=1, to=nrow(shp), by=2)) {
			j <- min(i+1, nrow(shp))
			plist <- list()
			for (k in i:j) {
				plist[[length(plist)+1]] <- sv_waterfall(shp, row_id=k) + ggtitle(sprintf("%s", names(response)[k]))
				plist[[length(plist)+1]] <- sv_force(shp, row_id=k) + ggtitle(sprintf("%s", names(response)[k]))
			}			
			do.call("grid.arrange", c(plist, ncol=2))
		}
		sv_importance(shp) + ggtitle(sprintf("shapr feature importance (%s, %s, %s)", mvar, level, tp))
		sv_importance(shp, kind="beeswarm") + ggtitle(sprintf("shapr feature importance (%s, %s, %s)", mvar, level, tp))
		sv_importance(shp, kind="both", show_numbers=TRUE, bee_width=0.2) + ggtitle(sprintf("shapr feature importance (%s, %s, %s)", mvar, level, tp))

#		## after running for the first time, COMMENT OUT THIS BLOCK ##
#		num_iter <- 1000
#		ncores <- 20
#		ranger.seeds <- sample(1:num_iter, num_iter, replace=T) # set up a vector of seeds for the ranger mclapply
#		out <- mclapply(1:num_iter, function (dummy) {
#				importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", case.weights=mapping.sel[rownames(data.sel), "ipw"], seed=ranger.seeds[dummy], num.threads=1))
#		}, mc.cores=ncores )
#	#	out <- mclapply(1:num_iter, function (dummy) {
#	#			importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", sample.fraction=c(1, 0.25), seed=ranger.seeds[dummy], num.threads=1))
#	#	}, mc.cores=ncores )	
#		collated.importance <- do.call(cbind, out)
#		out <- mclapply(1:num_iter, function (dummy) {
#				rgcv2(trainx=data.sel, trainy=response, cv.fold=10, step=0.9, num.threads=1)$error.cv
#			}, mc.cores=ncores )
#		collated.cv <- do.call(cbind, out)

#		write.table(collated.importance, file=sprintf("%s/ranger.%s.%s.%s.%s.importance.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)
#		write.table(collated.cv, file=sprintf("%s/ranger.%s.%s.%s.%s.cv.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)
#		## END BLOCK TO COMMENT ##

		collated.importance <- read.table(sprintf("%s/ranger.%s.%s.%s.%s.importance.txt", output_dir, aim, mvar, tp, level), header=F, as.is=T, sep="\t", row.names=1, quote="")
		collated.cv <- read.table(sprintf("%s/ranger.%s.%s.%s.%s.cv.txt", output_dir, aim, mvar, tp, level), header=F, as.is=T, sep="\t", row.names=1)
		importance.mean <- rowMeans(collated.importance)
		importance.sd <- unlist(apply(collated.importance, 1, sd))
		cv.mean <- rowMeans(collated.cv)
		cv.sd <- unlist(apply(collated.cv, 1, sd))
		inds <- order(importance.mean, decreasing=T)
		inds <- inds[1:min(20, as.numeric(names(which.min(cv.mean))))] # edit as appropriate
		write.table(melt(importance.mean[inds]), file=sprintf("%s/ranger.%s.%s.%s.%s.features.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)

		## after running for the first time, COMMENT OUT THIS BLOCK ##
		# using a sparse model with N predictors
#		sparseRanger <- ranger(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, num.trees=10000, importance="permutation", case.weights=mapping.sel[rownames(data.sel), "ipw"], seed=sample(1:num_iter,1))
#		save(sparseRanger, file=sprintf("%s/ranger.%s.%s.%s.%s.model", output_dir, aim, mvar, tp, level))
		load(sprintf("%s/ranger.%s.%s.%s.%s.model", output_dir, aim, mvar, tp, level))
		# accuracy of final sparseRF model
		pred <- predictions(sparseRanger)
		pred_df <- data.frame(SampleID=names(response), predicted=pred, true=response, stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
		pred_df_out <- merge(pred_df, data.sel, by="row.names")
		write.table(pred_df_out, file=sprintf("%s/ranger.%s.%s.%s.%s.predictions.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=F, col.names=T)
		confusion_matrix <- table(pred_df[, c("true", "predicted")])
		class_errors <- unlist(lapply(levels(mapping.sel$sabin_any), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$sabin_any)
		accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
		vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
		mccvalue <- mcc(vec.pred, vec.true)
		df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", aim, mvar, tp, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
		print(p)

		write.table(confusion_matrix, file=sprintf("%s/ranger.%s.%s.%s.%s.confusion_matrix.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=T)
		## END BLOCK TO COMMENT ##

		# plotting - per-group sparse model
		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s, %s, %s, %s", aim, mvar, tp, level)))
		# plotting - per-group variables
		df <- data.frame(Taxa=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
		# load effect sizes from linear regression
	#	contr <- "Case - Control"
		lmres <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, aim, level), header=T, as.is=T, sep="\t", quote=""); lmres <- subset(lmres, metadata_variable==mvar & Strata==tp); colnames(lmres)[1] <- "Taxa"
		df <- merge(lmres[,c("Taxa", "contrast", "Estimate", "SE", "padj", "dir")], df, by="Taxa")
		df$Taxa <- factor(df$Taxa, levels=rev(names(importance.mean)[inds]))
		p <- ggplot(df, aes(x=Taxa, y=importance, label=Taxa)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=Taxa, y=0, label=Taxa), size=3, hjust=0) + ggtitle(sprintf("%s: %s explanatory %s", aim, mvar, tp)) + theme(axis.text.y=element_blank())
		print(p)	
		lims <- max(abs(df$Estimate) + abs(df$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
		p <- ggplot(df, aes(x=Taxa, y=Estimate, color=dir, group=contrast)) + geom_point(position=pd) + geom_errorbar(aes(x=Taxa, ymin=Estimate-SE, max=Estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=contrast), position=pd, hjust=1, color="black", size=1.5) + geom_tile(aes(x=Taxa, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s %s explanatory %s", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
		print(p)
		
		# shading rectangles of importance values
		df.rect <- df
		df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
		p <- ggplot(df.rect, aes(x=x, y=Taxa, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s %s explanatory %s", aim, mvar, tp)) + scale_fill_gradient(low="white", high="black")
		print(p)
		# violin plots of relabund values
		agg.melt <- agg.melt.stored
		agg.melt[,mvar] <- mapping.sel[agg.melt$SampleID, mvar]
		agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$sabin_any))
		agg.melt <- subset(agg.melt, taxa %in% levels(df$Taxa))
		agg.melt$taxa <- factor(agg.melt$taxa, levels=levels(df$Taxa))
		p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap(~taxa, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
		p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~taxa, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		npages <- n_pages(p)
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~taxa, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
			print(p)
		}
		p <- ggplot(agg.melt, aes_string(x="taxa", y="value", color=mvar)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
	}
}



##################################################################################
##################################################################################
### Aim 1 - ZEBS only
##################################################################################
##################################################################################
psaim <- subset_samples(ps, Aim1=="Yes" & Study=="ZEBS")
psaim.relative <- subset_samples(ps.relative, Aim1=="Yes" & Study=="ZEBS")
psaim.rarefied <- subset_samples(ps.rarefied, Aim1=="Yes" & Study=="ZEBS")
aim <- "Aim1-ZEBS"
mvars.aim <- c("Transmitter", "Timepoint")
mapping.sel <- as(sample_data(psaim.relative), "data.frame")

## TableOne
demo_vars <- c("country", "instn", "cd4bl", "log10vlrna")
demo <- unique(mapping.sel[,c("Patient.ID", "Transmitter", demo_vars)])
write.table(print(CreateTableOne(vars=demo_vars, strata=c("Transmitter"), data=demo, smd=T)), file=sprintf("%s/Table_1.%s.by_Subject.txt", output_dir, aim), quote=F, sep="\t", row.names=T, col.names=T)
sample_vars <- c("age")
write.table(print(CreateTableOne(vars=sample_vars, strata=c("Timepoint", "Transmitter"), data=mapping.sel, smd=T)), file=sprintf("%s/Table_1.%s.by_Sample.txt", output_dir, aim), quote=F, sep="\t", row.names=T, col.names=T)
for (tp in levels(mapping.sel$Timepoint)) {
	write.table(print(CreateTableOne(vars=sample_vars, strata=c("Transmitter"), data=subset(mapping.sel, Timepoint==tp), smd=T)), file=sprintf("%s/Table_1.%s.by_Sample.%s.txt", output_dir, aim, tp), quote=F, sep="\t", row.names=T, col.names=T)
}

##################################################################################
### barplots by var
## order by PC1 (Bray-Curtis)
ordi <- ordinate(psaim.relative, method = "PCoA", distance = "bray")
ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))

for (mvar in rownames(subset(metadata_variables, useForPERMANOVA=="yes"))){
  otu.filt <- as.data.frame(otu_table(psaim.relative))
  otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.relative), level="Genus")
  otu.filt$Genus[which(is.na(otu.filt$Genus))] <- "Unknown"
  otu.filt$Genus[which(otu.filt$Genus=="")] <- "Unknown"
  otu.filt$Genus[which(otu.filt$Genus=="uncultured")] <- "Unknown"
  otu.filt$Genus <- gsub("_\\d$", "", otu.filt$Genus)
  agg <- aggregate(. ~ Genus, otu.filt, sum)
  genera <- agg$Genus
  agg <- agg[,-1]
  agg <- sweep(agg, 2, colSums(agg), "/")
  inds_to_keep <- which(rowSums(agg >= nsamps_threshold) >= ceiling(ncol(agg)*filt_threshold))
  genera[setdiff(1:nrow(agg), inds_to_keep)] <- "Other"
  agg$Genus <- genera
  df <- melt(agg, variable.name="SampleID")
  df2 <- aggregate(as.formula("value ~ Genus + SampleID"), df, sum)
  df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
  df2[[mvar]] <- mapping.sel[df2$SampleID, mvar]
  genera_to_add <- setdiff(df2$Genus, names(coloring))
	coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
	coloring.full <- c(coloring, coloring_to_add)
	coloring.full.sel <- intersect(names(coloring.full), df2$Genus); coloring.full <- coloring.full[coloring.full.sel]; coloring.full <- coloring.full[order(names(coloring.full))]
  p <- ggplot(df2, aes_string(x="SampleIDfactor", y="value", fill="Genus", order="Genus")) + facet_wrap(as.formula(sprintf("~%s", mvar)), scales="free_x") + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.0)) + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots - %s (%s)", mvar, aim)) + guides(col = guide_legend(ncol = 3))
  print(p)
}
# faceted by Transmitter+Timepoint
df2 <- aggregate(as.formula("value ~ Genus + SampleID"), df, sum)
df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
for (mvar in c(mvars.aim, "Batch")) {
	df2[[mvar]] <- mapping.sel[df2$SampleID, mvar]
}
genera_to_add <- setdiff(df2$Genus, names(coloring))
coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
coloring.full <- c(coloring, coloring_to_add)
coloring.full.sel <- intersect(names(coloring.full), df2$Genus); coloring.full <- coloring.full[coloring.full.sel]; coloring.full <- coloring.full[order(names(coloring.full))]
p <- ggplot(df2, aes_string(x="SampleIDfactor", y="value", fill="Genus", order="Genus")) + facet_wrap(as.formula(sprintf("~%s+%s", "Transmitter", "Timepoint")), scales="free_x") + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.0)) + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots - %s (%s)", "Transmitter+Timepoint", aim)) + guides(col = guide_legend(ncol = 3))
print(p)
for (b in levels(df2$Batch)) {
	p <- ggplot(subset(df2, Batch==b), aes_string(x="SampleIDfactor", y="value", fill="Genus", order="Genus")) + facet_wrap(as.formula(sprintf("~%s+%s", "Transmitter", "Timepoint")), scales="free_x") + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.0)) + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots - Batch %s %s (%s)", b, "Transmitter+Timepoint", aim)) + guides(col = guide_legend(ncol = 3))
	print(p)
}

# aggregated by Transmitter+Timepoint
otu.filt <- as.data.frame(otu_table(psaim.relative))
otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.relative), level="Genus")
otu.filt$Genus[which(is.na(otu.filt$Genus))] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="")] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="uncultured")] <- "Unknown"
otu.filt$Genus <- gsub("_\\d$", "", otu.filt$Genus)
agg <- aggregate(. ~ Genus, otu.filt, sum)
genera <- agg$Genus
agg <- agg[,-1]
agg <- sweep(agg, 2, colSums(agg), "/")
inds_to_keep <- which(rowSums(agg >= nsamps_threshold) >= ceiling(ncol(agg)*filt_threshold))
genera[setdiff(1:nrow(agg), inds_to_keep)] <- "Other"
agg$Genus <- genera
agg <- aggregate(. ~ Genus, agg, sum)
df <- melt(agg, variable.name="SampleID")
for (mvar in c(mvars.aim, "Batch")) {
	df[, mvar] <- mapping.sel[as.character(df$SampleID), mvar]
}
df2 <- aggregate(as.formula(sprintf("value ~ Genus + Transmitter + Timepoint")), df, mean)
genera_to_add <- setdiff(df2$Genus, names(coloring))
coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
coloring.full <- c(coloring, coloring_to_add)
coloring.full.sel <- intersect(names(coloring.full), df2$Genus); coloring.full <- coloring.full[coloring.full.sel]; coloring.full <- coloring.full[order(names(coloring.full))]
counts <- table(mapping.sel[, c(mvars.aim)])
df2[, "mvar"] <- sprintf("%s-%s\n(n=%d)", df2[, "Transmitter"], df2[, "Timepoint"], counts[cbind(as.character(df2[,"Transmitter"]), as.character(df2[,"Timepoint"]))])
p <- ggplot(df2, aes(x=Timepoint, y=value, fill=Genus)) + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.05)) + facet_wrap(~Transmitter, scales="free") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots (%s)", aim)) + guides(col = guide_legend(ncol = 3))
print(p)
p <- ggplot(df2, aes_string(x="mvar", y="value", fill="Genus", order="Genus")) + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.01)) + facet_wrap(~Transmitter, scales="free") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots (%s)", aim)) + guides(col = guide_legend(ncol = 3))
print(p)
for (b in levels(df$Batch)) {
	df2 <- aggregate(as.formula(sprintf("value ~ Genus + Transmitter + Timepoint")), subset(df, Batch==b), mean)
	genera_to_add <- setdiff(df2$Genus, names(coloring))
	coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
	coloring.full <- c(coloring, coloring_to_add)
	coloring.full.sel <- intersect(names(coloring.full), df2$Genus); coloring.full <- coloring.full[coloring.full.sel]; coloring.full <- coloring.full[order(names(coloring.full))]
	counts <- table(mapping.sel[, c(mvars.aim)])
	df2[, "mvar"] <- sprintf("%s-%s\n(n=%d)", df2[, "Transmitter"], df2[, "Timepoint"], counts[cbind(as.character(df2[,"Transmitter"]), as.character(df2[,"Timepoint"]))])

	p <- ggplot(df2, aes(x=Timepoint, y=value, fill=Genus)) + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.05)) + facet_wrap(~Transmitter, scales="free") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots (Batch %s, %s)", b, aim)) + guides(col = guide_legend(ncol = 3))
	print(p)
}



# per-taxon violin plots over time
# aggregated by Transmitter+Timepoint
otu.filt <- as.data.frame(otu_table(psaim.relative))
otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.relative), level="Genus")
otu.filt$Genus[which(is.na(otu.filt$Genus))] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="")] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="uncultured")] <- "Unknown"
otu.filt$Genus <- gsub("_\\d$", "", otu.filt$Genus)
agg <- aggregate(. ~ Genus, otu.filt, sum)
genera <- agg$Genus
agg <- agg[,-1]
agg <- sweep(agg, 2, colSums(agg), "/")
inds_to_keep <- which(rowSums(agg >= nsamps_threshold) >= ceiling(ncol(agg)*filt_threshold))
genera[setdiff(1:nrow(agg), inds_to_keep)] <- "Other"
agg$Genus <- genera
agg <- aggregate(. ~ Genus, agg, sum)
df <- melt(agg, variable.name="SampleID")
for (mvar in mvars.aim) {
	df[, mvar] <- mapping.sel[as.character(df$SampleID), mvar]
}
for (genus in unique(df$Genus)) {
	p <- ggplot(subset(df, Genus==genus), aes(x=Timepoint, y=value, color=Transmitter)) + geom_violin() + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_color_manual(values=cols.mvar[["Transmitter"]]) + ggtitle(sprintf("%s over Timepoint (%s)", genus, aim))
	print(p)
}


##################################################################################
### PCoA + PERMANOVA
psaim.sel <- psaim.relative
mapping.sel <- as(sample_data(psaim.sel), "data.frame")
sel <- rownames(mapping.sel)

## PCoA
for (mvar in rownames(subset(metadata_variables, useForPERMANOVA=="yes"))){
  for (distance_metric in distance_metrics) {
    ordi <- ordinate(psaim.sel, method = "PCoA", distance = distance_metric)
		if (metadata_variables[mvar, "type"] == "factor") {
			bd <- betadisper(as.dist(dm[[distance_metric]][sel,sel]), mapping.sel[,mvar])
			bd$centroids <- rbind(bd$centroids, rep(NA, ncol(bd$centroids))); rownames(bd$centroids)[nrow(bd$centroids)] <- "NA" # pad with NA row
			df.vectors <- as.data.frame(ordi$vectors[,1:2]); df.vectors[,mvar] <- mapping.sel[rownames(df.vectors),mvar]
			df.vectors$PC1origin <- bd$centroids[str_replace_na(as.character(df.vectors[, mvar])), "PCoA1"]
			df.vectors$PC2origin <- bd$centroids[str_replace_na(as.character(df.vectors[, mvar])), "PCoA2"]
			df.centroids <- as.data.frame(bd$centroids[1:(nrow(bd$centroids)-1), 1:2, drop=F]); df.centroids[,mvar] <- rownames(df.centroids) # drop the NA row
		  p <- plot_ordination(psaim.sel, ordi, "samples", color = mvar) + theme_classic() + ggtitle(distance_metric) + theme_classic() + stat_ellipse(geom="polygon", alpha=0.1, type="t", aes_string(fill=mvar)) + geom_segment(data=df.vectors, aes(x=PC1origin, y=PC2origin, xend=Axis.1, yend=Axis.2), inherit.aes=F, color="grey", arrow = arrow(length = unit(0.05, "inches")), size=0.5) + geom_point(data=df.centroids, aes_string(x="PCoA1", y="PCoA2", color=mvar), size=4, inherit.aes=F)
		  print(p)
		 } else {
		 	p <- plot_ordination(psaim.sel, ordi, "samples", color = mvar) + theme_classic() + ggtitle(distance_metric) + theme_classic()
		  print(p)
		 }
  }
}
# by Transmitter+Timepoint
p <- plot_ordination(psaim.sel, ordi, "samples", color = "Timepoint", shape = "Transmitter") + theme_classic() + ggtitle(distance_metric) + theme_classic() + stat_ellipse(geom="polygon", alpha=0.05, type="t", aes_string(fill="Timepoint", lty="Transmitter")) + scale_shape_manual(values=c(20,3))
print(p)


## PERMANOVA
out_txt <- sprintf("%s/PERMANOVA.%s.%s.txt", output_dir, aim, format(Sys.Date(), "%m%d%y"))
sink(out_txt, append=F)
print(sprintf("PERMANOVA (%s)", aim))
mvars <- c("Batch", "Timepoint", "Transmitter", "cd4bl", "log10vlrna")
for (distance_metric in distance_metrics) {
  print(distance_metric)
  mapping.sel <- subset(as(sample_data(psaim.sel), "data.frame"), Aim1=="Yes")
  ids_to_remove <- names(which(unlist(apply(mapping.sel[, mvars, drop=F], 1, function(x) any(is.na(x))))))
  nonna <- setdiff(rownames(mapping.sel), ids_to_remove)
  dm.sel <- dm[[distance_metric]][nonna, nonna]
  mapping.nonna <- mapping.sel[nonna,]
  form <- as.formula(sprintf("as.dist(dm.sel) ~ %s", paste(mvars, collapse="+")))
  res <- adonis2(form , data=mapping.nonna, permutations=999, by="margin")
  res$R2 <- res$SumOfSqs / sum(res$SumOfSqs)
  print(res)
}
sink()



##################################################################################
## alpha diversity (using emmeans)
psaim.sel <- psaim.rarefied
mapping.sel <- as(sample_data(psaim.sel), "data.frame")
sel <- rownames(mapping.sel)
adiv <- estimate_richness(psaim.sel, measures=alpha_metrics)
#rownames(adiv) <- gsub("^X", "", rownames(adiv))
adiv$SampleID <- rownames(adiv)
adiv <- merge(adiv, mapping.sel, by="row.names"); rownames(adiv) <- adiv$SampleID
res <- {}
mvars <- c("Transmitter")
for (mvar in mvars) {
	adiv.sel <- adiv[!is.na(adiv[,mvar]),] # remove NAs
	# boxplots/violin plots
	df <- melt(adiv.sel[,c(alpha_metrics, mvar)])
	p <- ggplot(df, aes_string(x=mvar, y="value", color=mvar)) + geom_boxplot() + facet_wrap(~variable, scales="free") + theme_classic() + ggtitle(sprintf("alpha diversity by %s (%s)", mvar, aim))
	print(p)
	for (alpha_metric in alpha_metrics) {
		mod <- lm(as.formula(sprintf("%s ~ %s*Timepoint", alpha_metric, mvar)), data=adiv.sel)
		emm.adiv <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s | Timepoint", mvar)), adjust="none")
		tmp <- as.data.frame(emm.adiv$contrasts)
		tmp$alpha_metric <- alpha_metric; tmp$metadata_variable <- mvar
		res <- rbind(res, tmp)
	}
}
res$padj <- p.adjust(res$p.value, method="fdr")
res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
res <- res[order(res$estimate, decreasing=T),]
for (mvar in mvars) {
	df <- subset(res, metadata_variable==mvar)
	contrs <- levels(droplevels(df$contrast))
	# single level contrast
	if (length(contrs)==1) {
		lims <- max(abs(df$estimate) + abs(df$SE))*1.0
		pd <- position_dodge(0.4)
		comp <- contrs[1]
		p <- ggplot(df, aes(x=Timepoint, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=Timepoint, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_hline(yintercept=0) + facet_wrap(~alpha_metric, scales="free") + theme_classic() + ggtitle(sprintf("%s (%s %s, IPTW)", "alpha diversity", mvar, comp)) + coord_flip() + scale_color_manual(values=dircolors)
		print(p)
	}
}
write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs610/phyloseq/adiv.%s.txt", aim), quote=F, sep="\t", row.names=F, col.names=T)


##################################################################################
### differential abundance - relative abundance with lmer
for (level in c("Genus", "Species")) {
	otu.filt <- as.data.frame(otu_table(psaim.relative))
	otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.relative), level=level)
	# rename Prevotella_6, etc -> Prevotella
	otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
	agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
	if (level %in% c("Species", "Genus")){
		agg <- agg[-1,]
	}
	lvl <- agg[[level]]
	agg <- agg[,-1]
	rownames(agg) <- lvl
	ftk <- names(which(rowSums(agg >= thresholds[["nsamps"]][[level]]) >= ceiling(ncol(agg)*thresholds[["filt"]][[level]])))
	agg <- agg[ftk,]
	agg[[level]] <- rownames(agg)
	res <- {}
	out <- mclapply(agg[[level]], function(f) {
		df <- melt(agg[f,]); colnames(df) <- c(level, "SampleID", "value"); df$SampleID <- as.character(df$SampleID)
		tmp <- {}
		print(f)
		nsamps_detected <- length(which(df$value>=thresholds[["nsamps"]][[level]]))
		for (mvar in c("Transmitter")) {
		  #print(sprintf("%s %s", f, mvar))
		  df2 <- df
		  for (m in c(mvar, "Timepoint", "Patient.ID")) {
		  	df2[, m] <- mapping.sel[df2$SampleID, m]
		  }
		  df2 <- subset(df2, !is.na(df2[,mvar,drop=F]))
			df2[, mvar] <- factor(as.character(df2[,mvar]), levels=rev(levels(df2[,mvar]))) # reverse levels to get more intuitive contrast
		  mod <- lmer(as.formula(sprintf("%s ~ %s*Timepoint + (1 | Patient.ID)", "value", mvar)), data=df2); modelstr <- "LMEM"
		  # contrast Transmitter, stratify by Timepoint
			emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s | Timepoint", mvar)), adjust="none")
			tmp2 <- as.data.frame(emm$contrasts); colnames(tmp2)[2] <- "Strata"
			tmp2[, level] <- f; tmp2$metadata_variable <- mvar; tmp2$model <- modelstr
			tmp <- rbind(tmp, tmp2)
			# contrast Timepoint, stratify by Transmitter
			emm <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ Timepoint | %s", mvar)), adjust="none")
			tmp2 <- as.data.frame(emm$contrasts); colnames(tmp2)[2] <- "Strata"
			tmp2[, level] <- f; tmp2$metadata_variable <- "Timepoint"; tmp2$model <- modelstr
			tmp <- rbind(tmp, tmp2)
		}
		print(sprintf("finished %s", f))
		tmp
	}, mc.cores=16)
	res <- as.data.frame(do.call(rbind, out))
	res <- res[,c(level, setdiff(colnames(res), level))]
	colnames(res) <- c(level, "contrast", "Strata", "Estimate", "SE", "df", "t.ratio", "p.value", "metadata_variable", "model")
	res$padj <- p.adjust(res$p.value, method="fdr")
	res <- res[order(res$p.value, decreasing=F),]
	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$Estimate)==1, "up", "down"), "NS")
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs610/phyloseq/emmeans.%s.%s.txt", aim, level), quote=F, sep="\t", row.names=F, col.names=T)
	# forest plot by Transmitter
	mvar <- "Transmitter"
	df <- subset(res, metadata_variable==mvar)
	df <- df[order(df$Estimate, decreasing=T),]
	df[, level] <- factor(as.character(df[, level]), levels=unique(as.character(df[, level])))
	df$contrast <- droplevels(df$contrast)
	lims <- max(abs(as.numeric(as.character(df$Estimate))) + abs(as.numeric(as.character(df$SE))))*1.0
	pd <- position_dodge(0.8)
	p <- ggplot(df, aes_string(x=level, y="Estimate", color="dir", group="Strata")) + geom_point(position=pd) + geom_errorbar(aes_string(x=level, ymin="Estimate-SE", max="Estimate+SE"), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Strata), position=pd, hjust=1, color="black", size=1.5) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LMEM: %s (%s, %s)", levels(df$contrast), mvar, aim)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
	print(p)
	# forest plot by Timepoint
	mvar <- "Timepoint"
	df <- subset(res, metadata_variable==mvar)
	df <- df[order(df$Estimate, decreasing=T),]
	df[, level] <- factor(as.character(df[, level]), levels=unique(as.character(df[, level])))
	df$contrast <- droplevels(df$contrast)
	lims <- max(abs(as.numeric(as.character(df$Estimate))) + abs(as.numeric(as.character(df$SE))))*1.0
	pd <- position_dodge(0.8)
	p <- ggplot(df, aes_string(x=level, y="Estimate", color="dir", group="Strata")) + geom_point(position=pd) + geom_errorbar(aes_string(x=level, ymin="Estimate-SE", max="Estimate+SE"), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Strata), position=pd, hjust=1, color="black", size=1.5) + facet_wrap(~contrast, scales="free") + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LMEM: %s (%s, %s)", levels(df$contrast), mvar, aim)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
	print(p)
	# violin plots
	agg.melt <- melt(agg)
	colnames(agg.melt) <- c(level, "SampleID", "value"); agg.melt$SampleID <- as.character(agg.melt$SampleID)
	for (mvar in mvars.aim) {
		agg.melt[, mvar] <- mapping.sel[agg.melt$SampleID, mvar]
	}
	mvar <- "Transmitter"
	p <- ggplot(agg.melt, aes_string(x="Timepoint", y="value", color=mvar)) + geom_violin(position=pd) + geom_point(position=pd) + facet_wrap(as.formula(sprintf("~%s", level)), scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("relabund (%s, %s)", mvar, aim)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
	print(p)
	p <- ggplot(agg.melt, aes_string(x="Timepoint", y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(as.formula(sprintf("~%s", level)), scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("relabund (%s, %s)", mvar, aim)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
	npages <- n_pages(p)
	for (ip in 1:npages) {
		p <- ggplot(agg.melt, aes_string(x="Timepoint", y="value", color=mvar)) + geom_violin(position=pd) + geom_point(position=pd) + facet_wrap_paginate(as.formula(sprintf("~%s", level)), scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("relabund (%s, %s)", mvar, aim)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
	}
}


##################################################################################
### random forest

## randomForest classification of Transmitter (separately for each time point)
psaim <- psaim.relative
set.seed(nsamples(psaim))
mvar <- "Transmitter"
res.mean <- list()
res.sd <- list()
for (level in c("Genus", "Species")) {
	for (tp in levels(mapping$Timepoint)) {
		psaim.sel <- subset_samples(psaim, Timepoint==tp)
		mapping.sel <- as(sample_data(psaim.sel), "data.frame")
		tab <- table(mapping.sel[,mvar])
		wt <- 1-tab/sum(tab)
		mapping.sel$ipw <- wt[as.character(mapping.sel[,mvar])]
		
		otu.filt <- as.data.frame(otu_table(psaim.sel))
		otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.sel), level=level)
		# rename Prevotella_6, etc -> Prevotella
		otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
		otu.filt <- otu.filt[which(otu.filt[[level]]!=""),]
		agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
		lvl <- agg[[level]]
		agg <- agg[,-1]
		rownames(agg) <- lvl
		agg <- agg[which(rowSums(agg >= thresholds[["nsamps"]][[level]]) >= ceiling(ncol(agg)*thresholds[["filt"]][[level]])),]
		
		sel <- colnames(agg)
		data.sel <- as.data.frame(t(agg[,sel]))
		data.sel <- as.matrix(data.sel)
		response <- droplevels(mapping.sel[sel, mvar]); names(response) <- sel
		# subset to non-NA
		response <- subset(response, !is.na(response))
		data.sel <- data.sel[names(response),]
		agg.melt.stored <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "taxa", "value")

		## after running for the first time, COMMENT OUT THIS BLOCK ##
		num_iter <- 1000
		ncores <- 20
		ranger.seeds <- sample(1:num_iter, num_iter, replace=T) # set up a vector of seeds for the ranger mclapply
		out <- mclapply(1:num_iter, function (dummy) {
				importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", case.weights=mapping.sel[rownames(data.sel), "ipw"], seed=ranger.seeds[dummy], num.threads=1))
		}, mc.cores=ncores )
	#	out <- mclapply(1:num_iter, function (dummy) {
	#			importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", sample.fraction=c(1, 0.25), seed=ranger.seeds[dummy], num.threads=1))
	#	}, mc.cores=ncores )	
		collated.importance <- do.call(cbind, out)
		out <- mclapply(1:num_iter, function (dummy) {
				rgcv2(trainx=data.sel, trainy=response, cv.fold=10, step=0.9, num.threads=1)$error.cv
			}, mc.cores=ncores )
		collated.cv <- do.call(cbind, out)

		write.table(collated.importance, file=sprintf("%s/ranger.%s.%s.%s.%s.importance.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)
		write.table(collated.cv, file=sprintf("%s/ranger.%s.%s.%s.%s.cv.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)
		## END BLOCK TO COMMENT ##

		collated.importance <- read.table(sprintf("%s/ranger.%s.%s.%s.%s.importance.txt", output_dir, aim, mvar, tp, level), header=F, as.is=T, sep="\t", row.names=1, quote="")
		collated.cv <- read.table(sprintf("%s/ranger.%s.%s.%s.%s.cv.txt", output_dir, aim, mvar, tp, level), header=F, as.is=T, sep="\t", row.names=1)
		importance.mean <- rowMeans(collated.importance)
		importance.sd <- unlist(apply(collated.importance, 1, sd))
		cv.mean <- rowMeans(collated.cv)
		cv.sd <- unlist(apply(collated.cv, 1, sd))
		inds <- order(importance.mean, decreasing=T)
		inds <- inds[1:min(20, as.numeric(names(which.min(cv.mean))))] # edit as appropriate
		write.table(melt(importance.mean[inds]), file=sprintf("%s/ranger.%s.%s.%s.%s.features.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)

		## after running for the first time, COMMENT OUT THIS BLOCK ##
		# using a sparse model with N predictors
		sparseRanger <- ranger(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, num.trees=10000, importance="permutation", case.weights=mapping.sel[rownames(data.sel), "ipw"], seed=sample(1:num_iter,1))
		save(sparseRanger, file=sprintf("%s/ranger.%s.%s.%s.%s.model", output_dir, aim, mvar, tp, level))
		load(sprintf("%s/ranger.%s.%s.%s.%s.model", output_dir, aim, mvar, tp, level))
		# accuracy of final sparseRF model
		pred <- predictions(sparseRanger)
		pred_df <- data.frame(SampleID=names(response), predicted=pred, true=response, stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
		pred_df_out <- merge(pred_df, data.sel, by="row.names")
		write.table(pred_df_out, file=sprintf("%s/ranger.%s.%s.%s.%s.predictions.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=F, col.names=T)
		confusion_matrix <- table(pred_df[, c("true", "predicted")])
		class_errors <- unlist(lapply(levels(mapping.sel$sabin_any), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$sabin_any)
		accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
		vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
		mccvalue <- mcc(vec.pred, vec.true)
		df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", aim, mvar, tp, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
		print(p)

		write.table(confusion_matrix, file=sprintf("%s/ranger.%s.%s.%s.%s.confusion_matrix.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=T)
		## END BLOCK TO COMMENT ##

		# plotting - per-group sparse model
		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s, %s, %s, %s", aim, mvar, tp, level)))
		# plotting - per-group variables
		df <- data.frame(Taxa=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
		# load effect sizes from linear regression
	#	contr <- "Case - Control"
		lmres <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, aim, level), header=T, as.is=T, sep="\t", quote=""); lmres <- subset(lmres, metadata_variable==mvar & Strata==tp); colnames(lmres)[1] <- "Taxa"
		df <- merge(lmres[,c("Taxa", "contrast", "Estimate", "SE", "padj", "dir")], df, by="Taxa")
		df$Taxa <- factor(df$Taxa, levels=rev(names(importance.mean)[inds]))
		p <- ggplot(df, aes(x=Taxa, y=importance, label=Taxa)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=Taxa, y=0, label=Taxa), size=3, hjust=0) + ggtitle(sprintf("%s: %s explanatory %s", aim, mvar, tp)) + theme(axis.text.y=element_blank())
		print(p)	
		lims <- max(abs(df$Estimate) + abs(df$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
		p <- ggplot(df, aes(x=Taxa, y=Estimate, color=dir, group=contrast)) + geom_point(position=pd) + geom_errorbar(aes(x=Taxa, ymin=Estimate-SE, max=Estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=contrast), position=pd, hjust=1, color="black", size=1.5) + geom_tile(aes(x=Taxa, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s %s explanatory %s", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
		print(p)
		
		# shading rectangles of importance values
		df.rect <- df
		df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
		p <- ggplot(df.rect, aes(x=x, y=Taxa, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s %s explanatory %s", aim, mvar, tp)) + scale_fill_gradient(low="white", high="black")
		print(p)
		# violin plots of relabund values
		agg.melt <- agg.melt.stored
		agg.melt[,mvar] <- mapping.sel[agg.melt$SampleID, mvar]
		agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$sabin_any))
		agg.melt <- subset(agg.melt, taxa %in% levels(df$Taxa))
		agg.melt$taxa <- factor(agg.melt$taxa, levels=levels(df$Taxa))
		p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap(~taxa, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
		p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~taxa, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		npages <- n_pages(p)
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~taxa, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
			print(p)
		}
		p <- ggplot(agg.melt, aes_string(x="taxa", y="value", color=mvar)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
	}
}






##################################################################################
##################################################################################
### Aim 2
##################################################################################
##################################################################################
psaim <- subset_samples(ps, Aim2=="Yes")
psaim.relative <- subset_samples(ps.relative, Aim2=="Yes")
psaim.rarefied <- subset_samples(ps.rarefied, Aim2=="Yes")
aim <- "Aim2"
mvars.aim <- c("Regimen", "Visit")
mapping.sel <- as(sample_data(psaim.relative), "data.frame")

## TableOne
demo_vars <- c("country", "instn", "cd4bl", "log10vlrna", "delgage", "parity")
demo <- unique(mapping.sel[,c("Patient.ID", "Regimen", demo_vars)])
write.table(print(CreateTableOne(vars=demo_vars, strata=c("Regimen"), data=demo, smd=T), noSpaces=T), file=sprintf("%s/Table_1.%s.by_Subject.txt", output_dir, aim), quote=F, sep="\t", row.names=T, col.names=T)
write.table(table(demo$instn), file=sprintf("%s/instn_counts.%s.by_Subject.txt", output_dir, aim), quote=F, sep="\t", row.names=T, col.names=T)


sample_vars <- c("age")
write.table(print(CreateTableOne(vars=sample_vars, strata=c("Visit", "Regimen"), data=mapping.sel, smd=T)), file=sprintf("%s/Table_1.%s.by_Sample.txt", output_dir, aim), quote=F, sep="\t", row.names=T, col.names=T)
for (v in levels(mapping.sel$Visit)) {
	write.table(print(CreateTableOne(vars=sample_vars, strata=c("Regimen"), data=subset(mapping.sel, Visit==v), smd=T)), file=sprintf("%s/Table_1.%s.by_Sample.%s.txt", output_dir, aim, v), quote=F, sep="\t", row.names=T, col.names=T)
}

##################################################################################
### barplots by var

## order by PC1 (Bray-Curtis)
ordi <- ordinate(psaim.relative, method = "PCoA", distance = "bray")
ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))

for (mvar in rownames(subset(metadata_variables, useForPERMANOVA=="yes"))){
  otu.filt <- as.data.frame(otu_table(psaim.relative))
  otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.relative), level="Genus")
  otu.filt$Genus[which(is.na(otu.filt$Genus))] <- "Unknown"
  otu.filt$Genus[which(otu.filt$Genus=="")] <- "Unknown"
  otu.filt$Genus[which(otu.filt$Genus=="uncultured")] <- "Unknown"
  otu.filt$Genus <- gsub("_\\d$", "", otu.filt$Genus)
  agg <- aggregate(. ~ Genus, otu.filt, sum)
  genera <- agg$Genus
  agg <- agg[,-1]
  agg <- sweep(agg, 2, colSums(agg), "/")
  inds_to_keep <- which(rowSums(agg >= nsamps_threshold) >= ceiling(ncol(agg)*filt_threshold))
  genera[setdiff(1:nrow(agg), inds_to_keep)] <- "Other"
  agg$Genus <- genera
  df <- melt(agg, variable.name="SampleID")
  df2 <- aggregate(as.formula("value ~ Genus + SampleID"), df, sum)
  df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
  df2[[mvar]] <- mapping.sel[df2$SampleID, mvar]
  genera_to_add <- setdiff(df2$Genus, names(coloring))
	coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
	coloring.full <- c(coloring, coloring_to_add)
	coloring.full.sel <- intersect(names(coloring.full), df2$Genus); coloring.full <- coloring.full[coloring.full.sel]
  p <- ggplot(df2, aes_string(x="SampleIDfactor", y="value", fill="Genus", order="Genus")) + facet_wrap(as.formula(sprintf("~%s", mvar)), scales="free_x") + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.0)) + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots - %s (Aim 2)", mvar)) + guides(col = guide_legend(ncol = 3))
  print(p)
}
# faceted by Regimen+Visit
df2 <- aggregate(as.formula("value ~ Genus + SampleID"), df, sum)
df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
for (mvar in mvars.aim) {
	df2[[mvar]] <- mapping.sel[df2$SampleID, mvar]
}
genera_to_add <- setdiff(df2$Genus, names(coloring))
coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
coloring.full <- c(coloring, coloring_to_add)
coloring.full.sel <- intersect(names(coloring.full), df2$Genus); coloring.full <- coloring.full[coloring.full.sel]
p <- ggplot(df2, aes_string(x="SampleIDfactor", y="value", fill="Genus", order="Genus")) + facet_wrap(as.formula(sprintf("~%s+%s", "Regimen", "Visit")), scales="free_x", ncol=4, nrow=2) + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.0)) + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots - %s (Aim 2)", "Regimen+Visit")) + guides(col = guide_legend(ncol = 3))
print(p)

# aggregated by Regimen+Visit
otu.filt <- as.data.frame(otu_table(psaim.relative))
otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.relative), level="Genus")
otu.filt$Genus[which(is.na(otu.filt$Genus))] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="")] <- "Unknown"
otu.filt$Genus[which(otu.filt$Genus=="uncultured")] <- "Unknown"
otu.filt$Genus <- gsub("_\\d$", "", otu.filt$Genus)
agg <- aggregate(. ~ Genus, otu.filt, sum)
genera <- agg$Genus
agg <- agg[,-1]
agg <- sweep(agg, 2, colSums(agg), "/")
inds_to_keep <- which(rowSums(agg >= nsamps_threshold) >= ceiling(ncol(agg)*filt_threshold))
genera[setdiff(1:nrow(agg), inds_to_keep)] <- "Other"
agg$Genus <- genera
agg <- aggregate(. ~ Genus, agg, sum)
df <- melt(agg, variable.name="SampleID")
for (mvar in mvars.aim) {
	df[, mvar] <- mapping.sel[as.character(df$SampleID), mvar]
}
df2 <- aggregate(as.formula(sprintf("value ~ Genus + Regimen + Visit")), df, mean)
genera_to_add <- setdiff(df2$Genus, names(coloring))
coloring_to_add <- missing_colors[1:length(genera_to_add)]; names(coloring_to_add) <- genera_to_add
coloring.full <- c(coloring, coloring_to_add)
coloring.full.sel <- intersect(names(coloring.full), df2$Genus); coloring.full <- coloring.full[coloring.full.sel]
counts <- table(mapping.sel[, c(mvars.aim)])
df2[, "mvar"] <- sprintf("%s-%s\n(n=%d)", df2[, "Regimen"], df2[, "Visit"], counts[cbind(as.character(df2[,"Regimen"]), as.character(df2[,"Visit"]))])
p <- ggplot(df2, aes(x=Visit, y=value, fill=Genus)) + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.05)) + facet_wrap(~Regimen, scales="free") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots (Aim 2)")) + guides(col = guide_legend(ncol = 3))
print(p)
p <- ggplot(df2, aes_string(x="mvar", y="value", fill="Genus", order="Genus")) + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.01)) + facet_wrap(~Regimen, scales="free") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring.full) + ggtitle(sprintf("Taxa barplots (Aim 2)")) + guides(col = guide_legend(ncol = 3))
print(p)



##################################################################################
### PCoA + PERMANOVA
psaim.sel <- psaim.relative
mapping.sel <- as(sample_data(psaim.sel), "data.frame")
sel <- rownames(mapping.sel)

## PCoA
for (mvar in rownames(subset(metadata_variables, useForPERMANOVA=="yes"))){
  for (distance_metric in distance_metrics) {
    ordi <- ordinate(psaim.sel, method = "PCoA", distance = distance_metric)
		if (metadata_variables[mvar, "type"] == "factor") {
			bd <- betadisper(as.dist(dm[[distance_metric]][sel,sel]), mapping.sel[,mvar])
			bd$centroids <- rbind(bd$centroids, rep(NA, ncol(bd$centroids))); rownames(bd$centroids)[nrow(bd$centroids)] <- "NA" # pad with NA row
			df.vectors <- as.data.frame(ordi$vectors[,1:2]); df.vectors[,mvar] <- mapping.sel[rownames(df.vectors),mvar]
			df.vectors$PC1origin <- bd$centroids[str_replace_na(as.character(df.vectors[, mvar])), "PCoA1"]
			df.vectors$PC2origin <- bd$centroids[str_replace_na(as.character(df.vectors[, mvar])), "PCoA2"]
			df.centroids <- as.data.frame(bd$centroids[1:(nrow(bd$centroids)-1), 1:2, drop=F]); df.centroids[,mvar] <- rownames(df.centroids) # drop the NA row
		  p <- plot_ordination(psaim.sel, ordi, "samples", color = mvar) + theme_classic() + ggtitle(distance_metric) + theme_classic() + stat_ellipse(geom="polygon", alpha=0.1, type="t", aes_string(fill=mvar)) + geom_segment(data=df.vectors, aes(x=PC1origin, y=PC2origin, xend=Axis.1, yend=Axis.2), inherit.aes=F, color="grey", arrow = arrow(length = unit(0.05, "inches")), size=0.5) + geom_point(data=df.centroids, aes_string(x="PCoA1", y="PCoA2", color=mvar), size=4, inherit.aes=F)
		  print(p)
		 } else {
		 	p <- plot_ordination(psaim.sel, ordi, "samples", color = mvar) + theme_classic() + ggtitle(distance_metric) + theme_classic()
		  print(p)
		 }
  }
}
# by Regimen+Visit
p <- plot_ordination(psaim.sel, ordi, "samples", color = "Visit", shape = "Regimen") + theme_classic() + ggtitle(distance_metric) + theme_classic() + stat_ellipse(geom="polygon", alpha=0.05, type="t", aes_string(fill="Visit", lty="Regimen")) + scale_shape_manual(values=c(20,3))
print(p)


## PERMANOVA
out_txt <- sprintf("%s/PERMANOVA.%s.%s.txt", output_dir, aim, format(Sys.Date(), "%m%d%y"))
sink(out_txt, append=F)
print(sprintf("PERMANOVA (%s)", aim))
mvars <- c("Batch", "Visit", "Regimen", "cd4bl", "log10vlrna", "instn")
for (distance_metric in distance_metrics) {
  print(distance_metric)
  mapping.sel <- subset(as(sample_data(psaim.sel), "data.frame"), Aim2=="Yes")
  ids_to_remove <- names(which(unlist(apply(mapping.sel[, mvars, drop=F], 1, function(x) any(is.na(x))))))
  nonna <- setdiff(rownames(mapping.sel), ids_to_remove)
  dm.sel <- dm[[distance_metric]][nonna, nonna]
  mapping.nonna <- mapping.sel[nonna,]
  form <- as.formula(sprintf("as.dist(dm.sel) ~ %s", paste(mvars, collapse="+")))
  res <- adonis2(form , data=mapping.nonna, permutations=999, by="margin")
  res$R2 <- res$SumOfSqs / sum(res$SumOfSqs)
  print(res)
}
sink()


##################################################################################
## alpha diversity (using emmeans)
psaim.sel <- psaim.rarefied
mapping.sel <- as(sample_data(psaim.sel), "data.frame")
sel <- rownames(mapping.sel)
adiv <- estimate_richness(psaim.sel, measures=alpha_metrics)
#rownames(adiv) <- gsub("^X", "", rownames(adiv))
adiv$SampleID <- rownames(adiv)
adiv <- merge(adiv, mapping.sel, by="row.names"); rownames(adiv) <- adiv$SampleID
res <- {}
mvars <- c("Regimen")
for (mvar in mvars) {
	adiv.sel <- adiv[!is.na(adiv[,mvar]),] # remove NAs
	# boxplots/violin plots
	df <- melt(adiv.sel[,c(alpha_metrics, mvar)])
	p <- ggplot(df, aes_string(x=mvar, y="value", color=mvar)) + geom_boxplot() + facet_wrap(~variable, scales="free") + theme_classic() + ggtitle(sprintf("alpha diversity by %s (%s)", mvar, aim))
	print(p)
	for (alpha_metric in alpha_metrics) {
		mod <- lm(as.formula(sprintf("%s ~ %s*Visit", alpha_metric, mvar)), data=adiv.sel)
		emm.adiv <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s | Visit", mvar)), adjust="none")
		tmp <- as.data.frame(emm.adiv$contrasts)
		tmp$alpha_metric <- alpha_metric; tmp$metadata_variable <- mvar
		res <- rbind(res, tmp)
	}
}
res$padj <- p.adjust(res$p.value, method="fdr")
res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
res <- res[order(res$estimate, decreasing=T),]
for (mvar in mvars) {
	df <- subset(res, metadata_variable==mvar)
	contrs <- levels(droplevels(df$contrast))
	# single level contrast
	if (length(contrs)==1) {
		lims <- max(abs(df$estimate) + abs(df$SE))*1.0
		pd <- position_dodge(0.4)
		comp <- contrs[1]
		p <- ggplot(df, aes(x=Visit, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=Visit, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_hline(yintercept=0) + facet_wrap(~alpha_metric, scales="free") + theme_classic() + ggtitle(sprintf("%s (%s %s, IPTW)", "alpha diversity", mvar, comp)) + coord_flip() + scale_color_manual(values=dircolors)
		print(p)
	}
}
write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs610/phyloseq/adiv.%s.txt", aim), quote=F, sep="\t", row.names=F, col.names=T)


##################################################################################
### differential abundance - relative abundance with lmer
mvar <- "Regimen"
for (level in c("Genus", "Species")) {
	otu.filt <- as.data.frame(otu_table(psaim.relative))
	otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.relative), level=level)
	# rename Prevotella_6, etc -> Prevotella
	otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
	agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
	if (level %in% c("Species", "Genus")){
		agg <- agg[-1,]
	}
	lvl <- agg[[level]]
	agg <- agg[,-1]
	rownames(agg) <- lvl
	ftk <- names(which(rowSums(agg >= thresholds[["nsamps"]][[level]]) >= ceiling(ncol(agg)*thresholds[["filt"]][[level]])))
	agg <- agg[ftk,]
	agg[[level]] <- rownames(agg)
	res <- {}
	out <- mclapply(agg[[level]], function(f) {
		df <- melt(agg[f,]); colnames(df) <- c(level, "SampleID", "value"); df$SampleID <- as.character(df$SampleID)
		tmp <- {}
		print(f)
		nsamps_detected <- length(which(df$value>=thresholds[["nsamps"]][[level]]))
		for (mvar in c("Regimen")) {
		  #print(sprintf("%s %s", f, mvar))
		  df2 <- df
		  for (m in c(mvar, "Visit", "Patient.ID")) {
		  	df2[, m] <- mapping.sel[df2$SampleID, m]
		  }
		  df2 <- subset(df2, !is.na(df2[,mvar,drop=F]))
			df2[, mvar] <- factor(as.character(df2[,mvar]), levels=rev(levels(df2[,mvar]))) # reverse levels to get more intuitive contrast
		  mod <- lmer(as.formula(sprintf("%s ~ %s*Visit + (1 | Patient.ID)", "value", mvar)), data=df2); modelstr <- "LMEM"
	#    mod <- lm(as.formula(sprintf("%s ~ %s", "value", mvar)), data=df2); modelstr <- "LMEM"
			emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s | Visit", mvar)), adjust="none")
			tmp2 <- as.data.frame(emm$contrasts) 
			tmp2[, level] <- f; tmp2$metadata_variable <- mvar; tmp2$model <- modelstr
			tmp <- rbind(tmp, tmp2)
		}
		print(sprintf("finished %s", f))
		tmp
	}, mc.cores=16)
	res <- as.data.frame(do.call(rbind, out))
	res <- res[,c(level, setdiff(colnames(res), level))]
	colnames(res) <- c(level, "contrast", "Visit", "Estimate", "SE", "df", "t.ratio", "p.value", "metadata_variable", "model")
	res$padj <- p.adjust(res$p.value, method="fdr")
	res <- res[order(res$p.value, decreasing=F),]
	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$Estimate)==1, "up", "down"), "NS")
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs610/phyloseq/emmeans.%s.%s.txt", aim, level), quote=F, sep="\t", row.names=F, col.names=T)
	# forest plot
	df <- subset(res, metadata_variable==mvar)
	df <- df[order(df$Estimate, decreasing=T),]
	df[, level] <- factor(as.character(df[, level]), levels=unique(as.character(df[, level])))
	df$contrast <- droplevels(df$contrast)
	lims <- max(abs(as.numeric(as.character(df$Estimate))) + abs(as.numeric(as.character(df$SE))))*1.0
	pd <- position_dodge(0.8)
	p <- ggplot(df, aes_string(x=level, y="Estimate", color="dir", group="Visit")) + geom_point(position=pd) + geom_errorbar(aes_string(x=level, ymin="Estimate-SE", max="Estimate+SE"), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Visit), position=pd, hjust=1, color="black", size=1.5) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LMEM: %s (%s, %s)", levels(df$contrast), mvar, aim)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
	print(p)
	# violin plots
	agg.melt <- melt(agg)
	colnames(agg.melt) <- c(level, "SampleID", "value"); agg.melt$SampleID <- as.character(agg.melt$SampleID)
	for (mvar in mvars.aim) {
		agg.melt[, mvar] <- mapping.sel[agg.melt$SampleID, mvar]
	}
	
	mvar <- "Regimen"
	p <- ggplot(agg.melt, aes_string(x="Visit", y="value", color=mvar)) + geom_violin(position=pd) + geom_point(position=pd) + facet_wrap(as.formula(sprintf("~%s", level)), scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("relabund (%s, %s)", mvar, aim)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
	print(p)
	p <- ggplot(agg.melt, aes_string(x="Visit", y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(as.formula(sprintf("~%s", level)), scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("relabund (%s, %s)", mvar, aim)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
	npages <- n_pages(p)
	for (ip in 1:npages) {
		p <- ggplot(agg.melt, aes_string(x="Visit", y="value", color=mvar)) + geom_violin(position=pd) + geom_point(position=pd) + facet_wrap_paginate(as.formula(sprintf("~%s", level)), scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("relabund (%s, %s)", mvar, aim)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
	}
}



##################################################################################
### random forest

## randomForest classification of Regimen (separately for each time point)
psaim <- psaim.relative
set.seed(nsamples(psaim))
mvar <- "Regimen"
res.mean <- list()
res.sd <- list()
for (level in c("Genus", "Species")) {
	for (tp in levels(mapping$Visit)) {
		psaim.sel <- subset_samples(psaim, Visit==tp)
		mapping.sel <- as(sample_data(psaim.sel), "data.frame")
		tab <- table(mapping.sel[,mvar])
		wt <- 1-tab/sum(tab)
		mapping.sel$ipw <- wt[as.character(mapping.sel[,mvar])]
		
		otu.filt <- as.data.frame(otu_table(psaim.sel))
		otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.sel), level=level)
		# rename Prevotella_6, etc -> Prevotella
		otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
		otu.filt <- otu.filt[which(otu.filt[[level]]!=""),]
		agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
		lvl <- agg[[level]]
		agg <- agg[,-1]
		rownames(agg) <- lvl
		agg <- agg[which(rowSums(agg >= thresholds[["nsamps"]][[level]]) >= ceiling(ncol(agg)*thresholds[["filt"]][[level]])),]
		
		sel <- colnames(agg)
		data.sel <- as.data.frame(t(agg[,sel]))
		data.sel <- as.matrix(data.sel)
		response <- droplevels(mapping.sel[sel, mvar]); names(response) <- sel
		# subset to non-NA
		response <- subset(response, !is.na(response))
		data.sel <- data.sel[names(response),]
		agg.melt.stored <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "taxa", "value")

#		## after running for the first time, COMMENT OUT THIS BLOCK ##
#		num_iter <- 1000
#		ncores <- 20
#		ranger.seeds <- sample(1:num_iter, num_iter, replace=T) # set up a vector of seeds for the ranger mclapply
#		out <- mclapply(1:num_iter, function (dummy) {
#				importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", case.weights=mapping.sel[rownames(data.sel), "ipw"], seed=ranger.seeds[dummy], num.threads=1))
#		}, mc.cores=ncores )
#	#	out <- mclapply(1:num_iter, function (dummy) {
#	#			importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", sample.fraction=c(1, 0.25), seed=ranger.seeds[dummy], num.threads=1))
#	#	}, mc.cores=ncores )	
#		collated.importance <- do.call(cbind, out)
#		out <- mclapply(1:num_iter, function (dummy) {
#				rgcv2(trainx=data.sel, trainy=response, cv.fold=10, step=0.9, num.threads=1)$error.cv
#			}, mc.cores=ncores )
#		collated.cv <- do.call(cbind, out)

#		write.table(collated.importance, file=sprintf("%s/ranger.%s.%s.%s.%s.importance.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)
#		write.table(collated.cv, file=sprintf("%s/ranger.%s.%s.%s.%s.cv.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)
#		## END BLOCK TO COMMENT ##

		collated.importance <- read.table(sprintf("%s/ranger.%s.%s.%s.%s.importance.txt", output_dir, aim, mvar, tp, level), header=F, as.is=T, sep="\t", row.names=1, quote="")
		collated.cv <- read.table(sprintf("%s/ranger.%s.%s.%s.%s.cv.txt", output_dir, aim, mvar, tp, level), header=F, as.is=T, sep="\t", row.names=1)
		importance.mean <- rowMeans(collated.importance)
		importance.sd <- unlist(apply(collated.importance, 1, sd))
		cv.mean <- rowMeans(collated.cv)
		cv.sd <- unlist(apply(collated.cv, 1, sd))
		inds <- order(importance.mean, decreasing=T)
		inds <- inds[1:min(20, as.numeric(names(which.min(cv.mean))))] # edit as appropriate
		write.table(melt(importance.mean[inds]), file=sprintf("%s/ranger.%s.%s.%s.%s.features.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)

		## after running for the first time, COMMENT OUT THIS BLOCK ##
		# using a sparse model with N predictors
#		sparseRanger <- ranger(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, num.trees=10000, importance="permutation", case.weights=mapping.sel[rownames(data.sel), "ipw"], seed=sample(1:num_iter,1))
#		save(sparseRanger, file=sprintf("%s/ranger.%s.%s.%s.%s.model", output_dir, aim, mvar, tp, level))
		load(sprintf("%s/ranger.%s.%s.%s.%s.model", output_dir, aim, mvar, tp, level))
		# accuracy of final sparseRF model
		pred <- predictions(sparseRanger)
		pred_df <- data.frame(SampleID=names(response), predicted=pred, true=response, stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
		pred_df_out <- merge(pred_df, data.sel, by="row.names")
		write.table(pred_df_out, file=sprintf("%s/ranger.%s.%s.%s.%s.predictions.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=F, col.names=T)
		confusion_matrix <- table(pred_df[, c("true", "predicted")])
		class_errors <- unlist(lapply(levels(mapping.sel$sabin_any), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$sabin_any)
		accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
		vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
		mccvalue <- mcc(vec.pred, vec.true)
		df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", aim, mvar, tp, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
		print(p)

		write.table(confusion_matrix, file=sprintf("%s/ranger.%s.%s.%s.%s.confusion_matrix.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=T)
		## END BLOCK TO COMMENT ##

		# plotting - per-group sparse model
		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s, %s, %s, %s", aim, mvar, tp, level)))
		# plotting - per-group variables
		df <- data.frame(Taxa=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
		# load effect sizes from linear regression
	#	contr <- "Case - Control"
		lmres <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, aim, level), header=T, as.is=T, sep="\t", quote=""); lmres <- subset(lmres, metadata_variable==mvar & Visit==tp); colnames(lmres)[1] <- "Taxa"
		df <- merge(lmres[,c("Taxa", "contrast", "Estimate", "SE", "padj", "dir")], df, by="Taxa")
		df$Taxa <- factor(df$Taxa, levels=rev(names(importance.mean)[inds]))
		p <- ggplot(df, aes(x=Taxa, y=importance, label=Taxa)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=Taxa, y=0, label=Taxa), size=3, hjust=0) + ggtitle(sprintf("%s: %s explanatory %s", aim, mvar, tp)) + theme(axis.text.y=element_blank())
		print(p)	
		lims <- max(abs(df$Estimate) + abs(df$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
		p <- ggplot(df, aes(x=Taxa, y=Estimate, color=dir, group=contrast)) + geom_point(position=pd) + geom_errorbar(aes(x=Taxa, ymin=Estimate-SE, max=Estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=contrast), position=pd, hjust=1, color="black", size=1.5) + geom_tile(aes(x=Taxa, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s %s explanatory %s", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
		print(p)
		
		# shading rectangles of importance values
		df.rect <- df
		df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
		p <- ggplot(df.rect, aes(x=x, y=Taxa, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s %s explanatory %s", aim, mvar, tp)) + scale_fill_gradient(low="white", high="black")
		print(p)
		# violin plots of relabund values
		agg.melt <- agg.melt.stored
		agg.melt[,mvar] <- mapping.sel[agg.melt$SampleID, mvar]
		agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$sabin_any))
		agg.melt <- subset(agg.melt, taxa %in% levels(df$Taxa))
		agg.melt$taxa <- factor(agg.melt$taxa, levels=levels(df$Taxa))
		p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap(~taxa, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
		p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~taxa, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		npages <- n_pages(p)
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~taxa, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
			print(p)
		}
		p <- ggplot(agg.melt, aes_string(x="taxa", y="value", color=mvar)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
	}
}


##################################################################################
### random forest (Regimen, combined over all timepoints)
psaim <- psaim.relative
set.seed(nsamples(psaim))
mvar <- "Regimen"
tp <- "Combined"
res.mean <- list()
res.sd <- list()
for (level in c("Genus", "Species")) {
	psaim.sel <- psaim
	mapping.sel <- as(sample_data(psaim.sel), "data.frame")
	tab <- table(mapping.sel[,mvar])
	wt <- 1-tab/sum(tab)
	mapping.sel$ipw <- wt[as.character(mapping.sel[,mvar])]
	
	otu.filt <- as.data.frame(otu_table(psaim.sel))
	otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psaim.sel), level=level)
	# rename Prevotella_6, etc -> Prevotella
	otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
	otu.filt <- otu.filt[which(otu.filt[[level]]!=""),]
	agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
	lvl <- agg[[level]]
	agg <- agg[,-1]
	rownames(agg) <- lvl
	agg <- agg[which(rowSums(agg >= thresholds[["nsamps"]][[level]]) >= ceiling(ncol(agg)*thresholds[["filt"]][[level]])),]
	
	sel <- colnames(agg)
	data.sel <- as.data.frame(t(agg[,sel]))
	data.sel <- as.matrix(data.sel)
	response <- droplevels(mapping.sel[sel, mvar]); names(response) <- sel
	# subset to non-NA
	response <- subset(response, !is.na(response))
	data.sel <- data.sel[names(response),]
	
	# cast into per-subject data (using mean abundance)
	data.sel <- as.data.frame(data.sel)
	data.sel$Patient.ID <- mapping.sel[rownames(data.sel), "Patient.ID"]
	data.sel$Visit <- mapping.sel[rownames(data.sel), "Visit"]
	agg <- aggregate(. ~ Patient.ID + Visit, data.sel, mean)
	j <- which(!colnames(agg) %in% c("Patient.ID", "Visit"))  # rescale each Patient.ID/Visit combo to Z scores
	for (i in 1:nrow(agg)) {
		x <- unlist(agg[i,j])
		agg[i,j] <- (x - mean(x)) / sd(x)
	}
	agg2 <- melt(agg)
	agg3 <- dcast(agg2, Patient.ID ~ Visit+variable, value.var="value")
	data.sel <- agg3; rownames(data.sel) <- data.sel$Patient.ID; data.sel <- data.sel[,-1]
	tmp <- unique(mapping.sel[, c("Patient.ID", "Regimen")])
	rownames(tmp) <- tmp$Patient.ID
	response <- tmp[rownames(data.sel), "Regimen"]; names(response) <- rownames(tmp)
	# subset to subjects with data in all four timepoints to avoid NAs
	sel <- names(which(apply(data.sel, 1, function(x) !any(is.na(x)))))
	data.sel <- data.sel[sel,]
	response <- response[sel]	
	agg.melt.stored <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt.stored) <- c("Patient.ID", "taxa", "value")

#	## after running for the first time, COMMENT OUT THIS BLOCK ##
#	num_iter <- 1000
#	ncores <- 20
#	ranger.seeds <- sample(1:num_iter, num_iter, replace=T) # set up a vector of seeds for the ranger mclapply
#	out <- mclapply(1:num_iter, function (dummy) {
#			importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", case.weights=mapping.sel[rownames(data.sel), "ipw"], seed=ranger.seeds[dummy], num.threads=1))
#	}, mc.cores=ncores )
##	out <- mclapply(1:num_iter, function (dummy) {
##			importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", seed=ranger.seeds[dummy], num.threads=1))
##	}, mc.cores=ncores )	
#	collated.importance <- do.call(cbind, out)
#	out <- mclapply(1:num_iter, function (dummy) {
#			rgcv2(trainx=data.sel, trainy=response, cv.fold=10, step=0.9, num.threads=1)$error.cv
#		}, mc.cores=ncores )
#	collated.cv <- do.call(cbind, out)

#	write.table(collated.importance, file=sprintf("%s/ranger.%s.%s.%s.%s.importance.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)
#	write.table(collated.cv, file=sprintf("%s/ranger.%s.%s.%s.%s.cv.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)
#	## END BLOCK TO COMMENT ##

	collated.importance <- read.table(sprintf("%s/ranger.%s.%s.%s.%s.importance.txt", output_dir, aim, mvar, tp, level), header=F, as.is=T, sep="\t", row.names=1, quote="")
	collated.cv <- read.table(sprintf("%s/ranger.%s.%s.%s.%s.cv.txt", output_dir, aim, mvar, tp, level), header=F, as.is=T, sep="\t", row.names=1)
	importance.mean <- rowMeans(collated.importance)
	importance.sd <- unlist(apply(collated.importance, 1, sd))
	cv.mean <- rowMeans(collated.cv)
	cv.sd <- unlist(apply(collated.cv, 1, sd))
	inds <- order(importance.mean, decreasing=T)
	inds <- inds[1:min(20, as.numeric(names(which.min(cv.mean))))] # edit as appropriate
	write.table(melt(importance.mean[inds]), file=sprintf("%s/ranger.%s.%s.%s.%s.features.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=F)

	## after running for the first time, COMMENT OUT THIS BLOCK ##
	# using a sparse model with N predictors
#	sparseRanger <- ranger(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, num.trees=10000, importance="permutation", case.weights=mapping.sel[rownames(data.sel), "ipw"], seed=sample(1:num_iter,1))
#	save(sparseRanger, file=sprintf("%s/ranger.%s.%s.%s.%s.model", output_dir, aim, mvar, tp, level))
	load(sprintf("%s/ranger.%s.%s.%s.%s.model", output_dir, aim, mvar, tp, level))
	# accuracy of final sparseRF model
	pred <- predictions(sparseRanger)
	pred_df <- data.frame(SampleID=names(response), predicted=pred, true=response, stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
	pred_df_out <- merge(pred_df, data.sel, by="row.names")
	write.table(pred_df_out, file=sprintf("%s/ranger.%s.%s.%s.%s.predictions.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=F, col.names=T)
	confusion_matrix <- table(pred_df[, c("true", "predicted")])
	class_errors <- unlist(lapply(levels(mapping.sel$sabin_any), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$sabin_any)
	accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
	vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
	mccvalue <- mcc(vec.pred, vec.true)
	df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", aim, mvar, tp, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
	print(p)

	write.table(confusion_matrix, file=sprintf("%s/ranger.%s.%s.%s.%s.confusion_matrix.txt", output_dir, aim, mvar, tp, level), quote=F, sep="\t", row.names=T, col.names=T)
	## END BLOCK TO COMMENT ##

	# plotting - per-group sparse model
	df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
	colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
	print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s, %s, %s, %s", aim, mvar, tp, level)))
	# plotting - per-group variables
	df <- data.frame(feature=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
	df$Timepoint <- unlist(lapply(as.character(df$feature), function(x) unlist(strsplit(x, "_"))[1]))
	df$Taxa <- unlist(lapply(as.character(df$feature), function(x) unlist(strsplit(x, "_"))[2]))
	# load effect sizes from linear regression
#	contr <- "Case - Control"
	lmres <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, aim, level), header=T, as.is=T, sep="\t", quote=""); lmres <- subset(lmres, metadata_variable==mvar); colnames(lmres)[1] <- "Taxa"; lmres$feature <- paste(lmres$Visit, lmres$Taxa, sep="_")
	df <- merge(lmres[,c("feature", "contrast", "Estimate", "SE", "padj", "dir")], df, by="feature")
	df <- df[order(df$importance, decreasing=T),]
	df$feature <- factor(df$feature, levels=rev(df$feature))
	p <- ggplot(df, aes(x=feature, y=importance, label=feature)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=feature, y=0, label=feature), size=3, hjust=0) + ggtitle(sprintf("%s: %s explanatory %s", aim, mvar, tp)) + theme(axis.text.y=element_blank())
	print(p)	
	lims <- max(abs(df$Estimate) + abs(df$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
	p <- ggplot(df, aes(x=feature, y=Estimate, color=dir, group=contrast)) + geom_point(position=pd) + geom_errorbar(aes(x=feature, ymin=Estimate-SE, max=Estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=contrast), position=pd, hjust=1, color="black", size=1.5) + geom_tile(aes(x=feature, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s %s explanatory %s", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
	print(p)	
	
	# shading rectangles of importance values
	df.rect <- df
	df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
	p <- ggplot(df.rect, aes(x=x, y=feature, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s %s explanatory %s", aim, mvar, tp)) + scale_fill_gradient(low="white", high="black")
	print(p)
	# violin plots of relabund values
	agg.melt <- agg.melt.stored; colnames(agg.melt)[2] <- "feature"
	agg.melt[,mvar] <- response[agg.melt$Patient.ID]
	agg.melt <- subset(agg.melt, feature %in% levels(df$feature))
	agg.melt$feature <- factor(agg.melt$feature, levels=levels(df$feature))
	p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap(~feature, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
	print(p)
	p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~feature, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
	npages <- n_pages(p)
	for (ip in 1:npages) {
		p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~feature, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
		print(p)
	}
	p <- ggplot(agg.melt, aes_string(x="feature", y="value", color=mvar)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s, %s)", aim, mvar, tp)) + coord_flip() + scale_color_manual(values=cols.mvar[[mvar]])
	print(p)
}






dev.off()




