#!/usr/bin/Rscript

library(ggplot2)
library(ape)
library(plyr)
library(reshape2)
library(cluster)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(gplots)
library(vegan)
library(irr)
library(useful)
library(pscl)
library(parallel)
library(poLCA)
library(igraph)
library(randomForest)
library(ROCR)
library(stringi)
library(mixOmics)
library(ggfortify)
library(Rtsne)
library(ggforce)
library(emmeans)
library(tableone)
library(abind)
library(limma)
library(psych)
library(glmnet)

source("/Lab_Share/fanli/code/PROMISE/utils.R")
source("/Lab_Share/fanli/code/PROMISE/mcc.R")

cols.sig <- c("black", "red", "grey"); names(cols.sig) <- c("ns", "sig", "NA")
cols.cohort <- c("#ca0020", "#f4a582", "#0571b0", "#92c5de"); names(cols.cohort) <- c("Term.zdv", "Preterm.zdv", "Term.zdvart", "Preterm.zdvart")
siglevel <- 0.05
dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")

mapping_fn <- "/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.083019.txt"
mapping <- read.table(mapping_fn, header=T, as.is=T, sep="\t")
colnames(mapping) <- c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "cpatid", "Country", "Plasma.drawdt", "DBS.drawdt", "InfantDBS.drawdt", "GestationalAgeAtCollection", "SampleID.Mom", "SampleID.Infant", "hemaval.mom", "hemaval.infant", "StudySite", "gender", "weight0week", "weight1week", "ap_onstgage", "nbclass", "instn.mom", "instn.infant", "DaysFromEntryToDelivery", "imputed.gage", "InfantAgeInDays", "InfantAgeInDaysBinned")

metadata_variables <- read.table("/Lab_Share/PROMISE/nwcs619/metadata_variables.081219.txt", header=T, as.is=T, sep="\t", row.names=1)
for (mvar in rownames(metadata_variables)) {
	if (metadata_variables[mvar, "type"] == "factor") {
		mapping[,mvar] <- factor(mapping[,mvar])
		if (metadata_variables[mvar, "baseline"] != "") {
			mapping[,mvar] <- relevel(mapping[,mvar], metadata_variables[mvar, "baseline"])
		}
	} else if (metadata_variables[mvar, "type"] == "ordered") {
		if (metadata_variables[mvar, "baseline"] == "") {
			lvls <- unique(mapping[,mvar])
		} else {
			lvls <- unlist(strsplit(metadata_variables[mvar, "baseline"], ","))
		}
		mapping[,mvar] <- ordered(mapping[,mvar], levels=lvls)
	} else if (metadata_variables[mvar, "type"] == "numeric") {
		mapping[,mvar] <- as.numeric(as.character(mapping[,mvar]))
	}
}

## remove Malawi 6 study site subjects
mapping <- subset(mapping, !(Country == "Malawi" & StudySite == "6"))

#########################################################################################################
### read in Metabolon data
metabolite_levels <- c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY")
## parse metabolon data and Z-transform
df.metabolon <- list()
metabolon_map <- data.frame()
metabolon_sortorder <- list()
for (st in c("DBS", "Plasma")) {
	df.metabolon[[st]] <- list()
	metabolon <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/ScaledImpData.%s.txt", st), header=T, as.is=T, sep="\t", quote="", comment.char="")
	colnames(metabolon) <- gsub("^X", "", colnames(metabolon))
	metabolon <- subset(metabolon, SUPER.PATHWAY!="") # remove uncharacterized molecules
	metabolon_map <- rbind(metabolon_map, metabolon[,c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY")])
	sel <- setdiff(colnames(metabolon), metabolite_levels)
	for (mlevel in metabolite_levels) {
		tmp <- metabolon[,c(mlevel, sel)]
		agg <- aggregate(as.formula(sprintf(". ~ %s", mlevel)), tmp, sum); rownames(agg) <- agg[,mlevel]; agg <- agg[,-1]
		agg <- agg[, as.character(intersect(colnames(agg), c(mapping$patid, mapping$cpatid)))] # filter to just the samples in the mapping file
		agg <- t(log(agg))
#		agg <- apply(agg, 1, function(x) (x-mean(x))/sd(x))
		agg <- agg[, setdiff(1:ncol(agg), which(is.na(apply(agg, 2, sd))))] # remove entries with zero variation
		# fix metabolite names as necessary
#		colnames(agg) <- gsub("\\*", "", colnames(agg))
		df.metabolon[[st]][[length(df.metabolon[[st]])+1]] <- agg
	}
	names(df.metabolon[[st]]) <- metabolite_levels
}
names(df.metabolon) <- c("DBS", "Plasma")
metabolon_map <- unique(metabolon_map); rownames(metabolon_map) <- metabolon_map$BIOCHEMICAL
cols.superpathway <- c(brewer.pal(length(unique(metabolon_map$SUPER.PATHWAY)), "Set1"), "#bbbbbb"); names(cols.superpathway) <- c(unique(metabolon_map$SUPER.PATHWAY), "NOT_METABOLITE")

out_pdf <- sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/metabolon_analysis_DBS_vs_plasma.%s.%s.pdf", "nwcs619", format(Sys.Date(), "%m%d%y"))
pdf(out_pdf, width=12)


#########################################################################################################
### QC data about metabolomics

## number of metabolites detected in each sample type, Venn diagrams
mlevel <- "BIOCHEMICAL"
qc <- {}; merged <- data.frame(BIOCHEMICAL=rownames(metabolon_map), detected.maternal_DBS=NA, detected.maternal_plasma=NA, detected.infant_DBS=NA, detected.infant_plasma=NA, median.maternal_DBS=NA, median.maternal_plasma=NA, median.infant_DBS=NA, median.infant_plasma=NA); rownames(merged) <- merged$BIOCHEMICAL
for (st in c("DBS", "Plasma")) {
	metabolon <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/OrigScale.%s.txt", st), header=T, as.is=T, sep="\t", quote="", comment.char="")
	colnames(metabolon) <- gsub("^X", "", colnames(metabolon))
	metabolon <- subset(metabolon, SUPER.PATHWAY!="") # remove uncharacterized molecules
	sel <- setdiff(colnames(metabolon), metabolite_levels)
	tmp <- metabolon[,sel]; rownames(tmp) <- metabolon[, mlevel]; ids <- colnames(tmp)
	tmp <- t(apply(tmp, 1, function(x) as.numeric(gsub(",", "", x)))); colnames(tmp) <- ids
	# count detectable as any non-NA value
	counts.maternal <- apply(tmp[, as.character(mapping$patid)], 1, function(x) length(which(!is.na(x))))
	if (any(as.character(mapping$cpatid) %in% colnames(tmp))) {
		counts.infant <- apply(tmp[, as.character(mapping$cpatid)], 1, function(x) length(which(!is.na(x))))
	} else {
		counts.infant <- rep(NA, length(counts.maternal))
	}
	# summary statistics
	median.maternal <- apply(log(tmp[, as.character(mapping$patid)]), 1, function(x) median(x,na.rm=T))
	if (any(as.character(mapping$cpatid) %in% colnames(tmp))) {
		median.infant <- apply(log(tmp[, as.character(mapping$cpatid)]), 1, function(x) median(x,na.rm=T))
	} else {
		median.infant <- rep(NA, length(median.maternal))
	}
	out <- data.frame(BIOCHEMICAL=rownames(tmp), subtype=st, detected.maternal=counts.maternal, detected.infant=counts.infant, median.maternal=median.maternal, median.infant=median.infant)
	if (st == "DBS") {
		merged[as.character(out$BIOCHEMICAL), c("BIOCHEMICAL", "detected.maternal_DBS", "detected.infant_DBS", "median.maternal_DBS", "median.infant_DBS")] <- out[, c("BIOCHEMICAL", "detected.maternal", "detected.infant", "median.maternal", "median.infant")]
	} else {
		merged[as.character(out$BIOCHEMICAL), c("BIOCHEMICAL", "detected.maternal_plasma", "detected.infant_plasma", "median.maternal_plasma", "median.infant_plasma")] <- out[, c("BIOCHEMICAL", "detected.maternal", "detected.infant", "median.maternal", "median.infant")]
	}
	write.table(out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/detected_counts.%s.txt", st), quote=F, sep="\t", row.names=F, col.names=T)
}
merged$FLAG.maternal_DBS <- ifelse(is.na(merged$detected.maternal_DBS), FALSE, merged$detected.maternal_DBS > 0)
merged$FLAG.infant_DBS <- ifelse(is.na(merged$detected.infant_DBS), FALSE, merged$detected.infant_DBS > 0)
merged$FLAG.maternal_plasma <- ifelse(is.na(merged$detected.maternal_plasma), FALSE, merged$detected.maternal_plasma > 0)
vc <- vennCounts(merged[, c("FLAG.maternal_DBS", "FLAG.infant_DBS", "FLAG.maternal_plasma")])
vennDiagram(vc, cex=c(1,1,1))
vc <- vennCounts(merged[, c("FLAG.maternal_DBS", "FLAG.maternal_plasma")])
vennDiagram(vc, cex=c(1,1))
write.table(merged, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/detected_counts.%s.txt", "merged"), quote=F, sep="\t", row.names=F, col.names=T)

## ICC and coefficient of variation in paired maternal plasma/DBS samples
subtype <- "maternal"
mapping.sel <- mapping; rownames(mapping.sel) <- mapping$patid
sel.metabolites <- rownames(subset(merged, FLAG.maternal_plasma & FLAG.maternal_DBS))
data <- matrix()
data.plasma <- df.metabolon[["Plasma"]][[mlevel]][as.character(mapping$patid), sel.metabolites]
data.dbs <- df.metabolon[["DBS"]][[mlevel]][as.character(mapping$patid), sel.metabolites]
data <- abind(data.plasma, data.dbs, along=0); dimnames(data)[[1]] <- c("Plasma", "DBS")
# ICC by metabolite
final <- apply(data, 3, function(x) {
	unlist(c(mean(x[,1]), mean(x[,2]), sd(x[,1]), sd(x[,2]), mean(abs(x[,1]-x[,2])), ICC(t(x))$results["Single_raters_absolute", c("ICC", "p", "lower bound", "upper bound")]))
#	icc(t(x))$value
}); final <- as.data.frame(t(final)); colnames(final) <- c("mean.plasma", "mean.DBS", "sd.plasma", "sd.DBS", "mean_difference", "ICC", "p", "lower_bound", "upper_bound")
final$padj <- p.adjust(final$p, method="fdr")
final <- final[order(final$p), ]
write.table(final, file="/Lab_Share/PROMISE/nwcs619/metabolon/ICC.metabolite.txt", quote=F, sep="\t", row.names=T, col.names=T)
final <- subset(final, mean.plasma > -5 & mean.DBS > -5) # remove a few really low abundance outliers
test <- cor.test(~ mean.plasma + ICC, final, method="spearman")
p <- ggplot(final, aes(x=mean.plasma, y=ICC)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.plasma (by metabolite, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
print(p)
test <- cor.test(~ mean.plasma + ICC, final, method="spearman")
p <- ggplot(final, aes(x=mean.DBS, y=ICC)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.DBS (by metabolite, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
print(p)
p <- ggplot(final, aes(x=ICC)) + geom_histogram(bins=50) + theme_classic() + ggtitle(sprintf("Distribution of ICC (by metabolite, mean=%.2g, median=%.2g)", mean(final$ICC), median(final$ICC))) + geom_vline(xintercept=mean(final$ICC), color="red")
print(p)
thresholds <- seq(from=0.05, to=0.95, by=0.05)
df <- sapply(thresholds, function(thresh) {
	length(which(final$ICC >= thresh))
}); df <- data.frame(threshold=thresholds, n=df); df$pct <- 100*(df$n / nrow(final))
p <- ggplot(df, aes(x=threshold, y=pct)) + geom_point() + geom_line() + geom_text(label=sprintf("%d (%.2g%%)", df$n, df$pct)) + theme_classic() + ggtitle(sprintf("Number of metabolites at ICC threshold (n=%d padj<0.05)", length(which(final$padj<0.05))))
print(p)
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Number of metabolites at ICC threshold")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)
df <- as.matrix(summary(final$ICC))
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("ICC by metabolite summary statistics")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)

# plot of mean metabolite values (axes are values in DBS and plasma), colored by ICC
# maybe also draw ellipse for each point as SD of each metabolite?
pdf(out_pdf,width=12)
p <- ggplot(final, aes(x=mean.plasma, y=mean.DBS, color=ICC)) + geom_point() + theme_classic() + ggtitle(sprintf("mean metabolite values and ICC")) + scale_color_gradient(low="black", high="green") + geom_vline(xintercept=0) + geom_hline(yintercept=0) + geom_abline(slope=1)
print(p)
p <- ggplot(final, aes(x=mean.plasma, y=mean.DBS, color=ICC)) + geom_point() + geom_errorbarh(aes(xmin=mean.plasma-sd.plasma, xmax=mean.plasma+sd.plasma), height=0.4) + geom_errorbar(aes(ymin=mean.DBS-sd.DBS, ymax=mean.DBS+sd.DBS), width=0.4) + theme_classic() + ggtitle(sprintf("mean metabolite values and ICC")) + scale_color_gradient(low="black", high="green")
print(p)
dev.off()


# ICC by subject
final <- apply(data, 2, function(x) {
	unlist(c(mean(x[,1]), mean(x[,2]), mean(abs(x[,1]-x[,2])), ICC(t(x))$results["Single_raters_absolute", c("ICC", "p", "lower bound", "upper bound")]))
#	icc(t(x))$value
}); final <- as.data.frame(t(final)); colnames(final) <- c("mean.plasma", "mean.DBS", "mean_difference", "ICC", "p", "lower_bound", "upper_bound")
final$padj <- p.adjust(final$p, method="fdr")
final <- final[order(final$p), ]
write.table(final, file="/Lab_Share/PROMISE/nwcs619/metabolon/ICC.subject.txt", quote=F, sep="\t", row.names=T, col.names=T)
test <- cor.test(~ mean.plasma + ICC, final, method="spearman")
p <- ggplot(final, aes(x=mean.plasma, y=ICC)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.plasma (by subject, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
print(p)
test <- cor.test(~ mean.plasma + ICC, final, method="spearman")
p <- ggplot(final, aes(x=mean.DBS, y=ICC)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.DBS (by subject, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
print(p)
p <- ggplot(final, aes(x=ICC)) + geom_histogram(bins=50) + theme_classic() + ggtitle(sprintf("Distribution of ICC (by subject, mean=%.2g, median=%.2g)", mean(final$ICC), median(final$ICC))) + geom_vline(xintercept=mean(final$ICC), color="red")
print(p)
thresholds <- seq(from=0.05, to=0.95, by=0.05)
df <- sapply(thresholds, function(thresh) {
	length(which(final$ICC >= thresh))
}); df <- data.frame(threshold=thresholds, n=df); df$pct <- 100*(df$n / nrow(final))
p <- ggplot(df, aes(x=threshold, y=pct)) + geom_point() + geom_line() + geom_text(label=sprintf("%d (%.2g%%)", df$n, df$pct)) + theme_classic() + ggtitle(sprintf("Number of metabolites at ICC threshold (by subject, n=%d padj<0.05)", length(which(final$padj<0.05))))
print(p)
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Number of metabolites at ICC threshold (by subject)")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)
df <- as.matrix(summary(final$ICC))
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("ICC by metabolite summary statistics (by subject)")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)
# tSNE by subject with both plasma and DBS
rownames(data.plasma) <- sprintf("%s.Plasma", rownames(data.plasma)); rownames(data.dbs) <- sprintf("%s.DBS", rownames(data.dbs))
data2 <- rbind(data.plasma, data.dbs)
tsne.out <- Rtsne(data2, perplexity=20)
df <- data.frame(PC1 = tsne.out$Y[,1], PC2 = tsne.out$Y[,2], SampleID=rownames(data2)); rownames(df) <- df$SampleID
df$SubjectID <- unlist(lapply(rownames(df), function(x) unlist(strsplit(x, "\\."))[1]))
df$SampleType <- unlist(lapply(rownames(df), function(x) unlist(strsplit(x, "\\."))[2]))
p <- ggplot(df, aes_string(x="PC1", y="PC2", colour="SampleType")) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("tSNE by subject (%s)", mvar)) + scale_color_brewer(palette="Set1")
print(p)
for (mvar in c("Group", "Delivery")) {
	df[, mvar] <- mapping.sel[df$SubjectID, mvar]
	p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("tSNE by subject (%s)", mvar)) + scale_color_brewer(palette="Set1")
	print(p)
}
# PERMANOVA
res <- adonis2(data2 ~ SampleType + Group + SubjectID, data=df, permutations=999, method='euclidean')
sink(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/metabolon_PERMANOVA_by_SampleType.txt"))
print(res)
sink()
# distance boxplots
distance_metric <- "euclidean"
dm <- as.matrix(dist(data2, method=distance_metric))
pair_df <- as.data.frame(t(combn(rownames(df), 2))); colnames(pair_df) <- c("s1", "s2")
pair_df$s1 <- as.character(pair_df$s1); pair_df$s2 <- as.character(pair_df$s2)
pair_df$id1 <- df[pair_df$s1, "SubjectID"]; pair_df$id2 <- df[pair_df$s2, "SubjectID"]
pair_df$Group <- ifelse(pair_df$id1 == pair_df$id2, "Within-Subject", "Between-Subject")
pair_df$dist <- dm[cbind(pair_df$s1, pair_df$s2)]
test <- wilcox.test(dist ~ Group, pair_df)
p <- ggplot(pair_df, aes(x=Group, y=dist)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("%s distance comparison (Wilcoxon p=%.4g)", distance_metric, test$p.value))
print(p)


dev.off()





