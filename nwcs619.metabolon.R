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
library(doParallel)
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
cols.cohort <- c("#808080", "#cccccc", "#0571b0", "#92c5de", "#ca0020", "#f4a582", "#ffa500", "#ffdb99"); names(cols.cohort) <- c("Term.untreated", "Preterm.untreated", "Term.zdv", "Preterm.zdv", "Term.PI-ART", "Preterm.PI-ART", "Term.other", "Preterm.other")
cols.regimen <- c("#808080", "#0571b0", "#ca0020"); names(cols.regimen) <- c("untreated", "zdv", "PI-ART")
cols.daystodelivery <- c("black", "green", "blue", "red"); names(cols.daystodelivery) <- c("(14,200]", "(7,14]", "(2,7]", "(0,2]")
cols.delivery2 <- c("#808080", "orange", "red"); names(cols.delivery2) <- c("Term", "Preterm", "VeryPreterm")
siglevel <- 0.05
siglevel_loose <- 0.2 # less conservative alpha for infant DBS
dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")

mapping_fn <- "/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.020320.txt"
mapping <- read.table(mapping_fn, header=T, as.is=T, sep="\t")

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

## (full) Cohort demographics and some QC data about metabolomics
write.table(print(CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "delgage", "hemaval.mom", "DaysPTDPlasma", "DaysPTDDBS"), strata=c("MaternalGroup"), data=mapping, smd=T)), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_1_full.MaternalGroup.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(print(CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "delgage", "hemaval.mom", "DaysPTDPlasma", "DaysPTDDBS"), strata=c("Delivery"), data=mapping, smd=T)), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_1_full.Delivery.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(print(CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "InfantAgeInDays", "delgage", "hemaval.infant"), strata=c("InfantGroup"), data=mapping, smd=T)), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_1_full.InfantGroup.txt", quote=F, sep="\t", row.names=T, col.names=T)
for (regi in c("untreated", "zdv", "PI-ART", "other")) {
	write.table(print(CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "GestationalAgeAtCollection", "delgage", "hemaval.mom", "hemaval.infant", "DaysPTDPlasma", "DaysPTDDBS"), strata=c("Delivery"), data=subset(mapping, MaternalRegimen==regi), smd=T)), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/Table_1_full.MaternalRegimen.%s.txt", regi), quote=F, sep="\t", row.names=T, col.names=T)
	write.table(print(CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "GestationalAgeAtCollection", "delgage", "hemaval.mom", "hemaval.infant"), strata=c("Delivery"), data=subset(mapping, InfantRegimen==regi), smd=T)), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/Table_1_full.InfantRegimen.%s.txt", regi), quote=F, sep="\t", row.names=T, col.names=T)
}

## remove Malawi 6 study site subjects
mapping.full <- mapping
mapping <- subset(mapping, !(Country == "Malawi" & StudySite == "6"))
mappinglist <- list(mapping, mapping.full)
#mappinglist <- list(mapping.full, mapping.full) # to keep all n=100 subjects for both DBS and Plasma
names(mappinglist) <- c("DBS", "Plasma")

out_pdf <- sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/metabolon_analysis.%s.%s.pdf", "nwcs619", format(Sys.Date(), "%m%d%y"))
pdf(out_pdf, width=12)


#########################################################################################################
### combined violin plot of metabolites that separate Malawi-6 samples from Malawi-3, with stability data from Metabolon
### read in data for n=100 dyads
#metabolite_levels <- c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY")
### parse metabolon data and Z-transform
#df.metabolon <- list()
#metabolon_map <- data.frame()
#metabolon_map_by_assay <- list()
#metabolon_sortorder <- list()
#for (st in c("DBS", "Plasma")) {
#	df.metabolon[[st]] <- list()
#	metabolon <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/ScaledImpData.%s.txt", st), header=T, as.is=T, sep="\t", quote="", comment.char="")
#	colnames(metabolon) <- gsub("^X", "", colnames(metabolon))
#	# store for data submission to PROMISE
#	out <- metabolon
#	exclusion <- ifelse(colnames(out) %in% c(mapping.full$patid, mapping.full$cpatid), "", "Unreliable")
#	exclusion[1:5] <- ""
#	out <- rbind(exclusion, out)
#	out[2:nrow(out), which(out[1,] == "Unreliable")] <- NA
#	write.table(out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/ScaledImputedData.for_DMC.%s.txt", st), row.names=F, col.names=T, sep="\t", quote=F)
##	metabolon <- subset(metabolon, SUPER.PATHWAY!="") # remove uncharacterized molecules
#	metabolon_map <- rbind(metabolon_map, metabolon[,c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY")])
#	tmp <- metabolon[,c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY", "COMP_ID", "PLATFORM")]
#	rownames(tmp) <- tmp$BIOCHEMICAL
#	metabolon_map_by_assay[[length(metabolon_map_by_assay)+1]] <- tmp
#	sel <- setdiff(colnames(metabolon), c(metabolite_levels, "COMP_ID", "PLATFORM"))
#	for (mlevel in metabolite_levels) {
#		tmp <- metabolon[,c(mlevel, sel)]
#		agg <- aggregate(as.formula(sprintf(". ~ %s", mlevel)), tmp, sum); rownames(agg) <- agg[,mlevel]; agg <- agg[,-1]
#		agg <- agg[, as.character(intersect(colnames(agg), c(mapping.full$patid, mapping.full$cpatid)))] # filter to just the samples in the mapping file
#		agg <- t(log(agg))
##		agg <- apply(agg, 1, function(x) (x-mean(x))/sd(x))
#		agg <- agg[, setdiff(1:ncol(agg), which(is.na(apply(agg, 2, sd))))] # remove entries with zero variation
#		# fix metabolite names as necessary
##		colnames(agg) <- gsub("\\*", "", colnames(agg))
#		df.metabolon[[st]][[length(df.metabolon[[st]])+1]] <- agg
#	}
#	names(df.metabolon[[st]]) <- metabolite_levels
#}
#names(df.metabolon) <- c("DBS", "Plasma")
#metabolon_map <- unique(metabolon_map); rownames(metabolon_map) <- metabolon_map$BIOCHEMICAL
#all_colors <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set3"))
#cols.superpathway <- c(all_colors[1:length(unique(metabolon_map$SUPER.PATHWAY))], "#bbbbbb"); names(cols.superpathway) <- c(unique(metabolon_map$SUPER.PATHWAY), "NOT_METABOLITE")
#names(metabolon_map_by_assay) <- c("DBS", "Plasma")


#subtype <- "maternal"
#stability_metabolites <- c("1-methylguanidine", "4-guanidinobutanoate", "N-acetylproline", "erythronate*", "glycerate", "arabonate/xylonate", "succinate", "thioproline")
#for (st in c("DBS")) {
#	for (mlevel in "BIOCHEMICAL") {
#		data <- df.metabolon[[st]][[mlevel]]
#		mapping.sel <- mapping.full[,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "weight0week", "StudySite")]
#		colnames(mapping.sel) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weight0week", "StudySite")
#		mapping.sel <- subset(mapping.sel, Country == "Malawi")
#		rownames(mapping.sel) <- mapping.sel$patid
#		data.sel <- as.data.frame(data[rownames(mapping.sel), stability_metabolites]) # subset to just the maternal samples from Malawi
#		agg.melt <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt) <- c("SampleID", "metabolite", "value")
#		agg.melt$StudySite <- droplevels(mapping.sel[agg.melt$SampleID, "StudySite"])
#		p <- ggplot(agg.melt, aes(x=metabolite, y=value, color=StudySite)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of stability metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#		print(p)
#	}
#}


#########################################################################################################
### read in Metabolon data
metabolite_levels <- c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY")
## parse metabolon data and Z-transform
df.metabolon <- list()
metabolon_map <- data.frame()
metabolon_map_by_assay <- list()
metabolon_sortorder <- list()
for (st in c("DBS", "Plasma")) {
	df.metabolon[[st]] <- list()
	metabolon <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/ScaledImpData.%s.txt", st), header=T, as.is=T, sep="\t", quote="", comment.char="")
	colnames(metabolon) <- gsub("^X", "", colnames(metabolon))
#	metabolon <- subset(metabolon, SUPER.PATHWAY!="") # remove uncharacterized molecules
	metabolon_map <- rbind(metabolon_map, metabolon[,c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY")])
	tmp <- metabolon[,c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY", "COMP_ID", "PLATFORM")]
	rownames(tmp) <- tmp$BIOCHEMICAL
	metabolon_map_by_assay[[length(metabolon_map_by_assay)+1]] <- tmp
	sel <- setdiff(colnames(metabolon), metabolite_levels)
	# use n=100 for plasma, n=79 for DBS (exclude Malawi-6)
	sel_subjects <- c(as.character(mappinglist[[st]]$patid), as.character(mappinglist[[st]]$cpatid))
	sel <- intersect(sel, sel_subjects)
	for (mlevel in metabolite_levels) {
		tmp <- metabolon[,c(mlevel, sel)]
		labels <- tmp[, mlevel]; tmp <- tmp[, setdiff(colnames(tmp), mlevel)]; ids <- colnames(tmp)
		tmp <- as.data.frame(t(apply(tmp, 1, function(x) {
			x <- as.numeric(gsub(",", "", x))
			to_impute <- which(is.na(x))
			x[to_impute] <- min(x, na.rm=T)
			x
		})))
		colnames(tmp) <- ids; tmp[, mlevel] <- labels
		agg <- aggregate(as.formula(sprintf(". ~ %s", mlevel)), tmp, sum); rownames(agg) <- agg[,mlevel]; agg <- agg[,-1]
		agg <- agg[, as.character(intersect(colnames(agg), sel_subjects))] # filter to just the samples in the mapping file
#		agg <- agg[, as.character(intersect(colnames(agg), mappinglist[[st]]$patid))] # filter to just the maternal samples in the mapping file
#		agg <- t(agg)
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
cols.superpathway <- c(brewer.pal(9, "Set1"), "#bbbbbb", "#111111")[1:(length(unique(metabolon_map$SUPER.PATHWAY))+1)]; names(cols.superpathway) <- c(unique(metabolon_map$SUPER.PATHWAY), "NOT_METABOLITE")

#######################
## parse raw metabolon data
df.metabolon.raw <- list()
for (st in c("DBS", "Plasma")) {
	df.metabolon.raw[[st]] <- list()
	metabolon <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/OrigScale.%s.txt", st), header=T, as.is=T, sep="\t", quote="", comment.char="")
	colnames(metabolon) <- gsub("^X", "", colnames(metabolon))
#	metabolon <- subset(metabolon, SUPER.PATHWAY!="") # remove uncharacterized molecules
#	metabolon_map <- rbind(metabolon_map, metabolon[,c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY")])
	sel <- setdiff(colnames(metabolon), metabolite_levels)
	# use n=100 for plasma, n=79 for DBS (exclude Malawi-6)
	sel_subjects <- c(as.character(mappinglist[[st]]$patid), as.character(mappinglist[[st]]$cpatid))
	sel <- intersect(sel, sel_subjects)
	for (mlevel in c("BIOCHEMICAL")) {
		tmp <- metabolon[,c(mlevel, sel)]
		labels <- tmp[, mlevel]; tmp <- tmp[, setdiff(colnames(tmp), mlevel)]; ids <- colnames(tmp)
		tmp <- as.data.frame(t(apply(tmp, 1, function(x) {
			x <- as.numeric(gsub(",", "", x))
			x
		})))
		to_remove <- which(apply(tmp, 1, function(x) all(is.na(x)))) # remove rows that are all NA (features found only in infant DBS)
		tmp <- tmp[setdiff(rownames(tmp), to_remove),,drop=F]; labels <- labels[setdiff(1:length(labels), to_remove)]
		colnames(tmp) <- ids; tmp[, mlevel] <- labels
		# fix metabolite names as necessary
#		colnames(agg) <- gsub("\\*", "", colnames(agg))
		tmp <- tmp[, as.character(intersect(colnames(tmp), sel_subjects))]
		rownames(tmp) <- labels
		df.metabolon.raw[[st]][[length(df.metabolon.raw[[st]])+1]] <- tmp
	}
	names(df.metabolon.raw[[st]]) <- "BIOCHEMICAL"
}
names(df.metabolon.raw) <- c("DBS", "Plasma")
#######################
## replace Metabolon ScaledImpData with impute/log-transform/Z-transform separately for maternal/infant
for (st in c("DBS")) {
	for (mlevel in c("BIOCHEMICAL")) {
		# maternal
		sel <- as.character(mapping$patid)
		tmp <- df.metabolon.raw[[st]][[mlevel]][, sel]
		tmp2 <- apply(tmp, 1, function(x) {
			x[which(is.na(x))] <- min(x, na.rm=T) # impute to min value
			x <- log(x) # log-transform
			x <- (x-mean(x))/sd(x) # Z-transform
			x
		})
		# infant
		sel <- as.character(mapping$cpatid)
		tmp.infant <- df.metabolon.raw[[st]][[mlevel]][, sel]
		tmp3 <- apply(tmp.infant, 1, function(x) {
			x[which(is.na(x))] <- min(x, na.rm=T) # impute to min value
			x <- log(x) # log-transform
			x <- (x-mean(x))/sd(x) # Z-transform
			x
		})
		replace <- rbind(tmp2, tmp3); replace <- replace[rownames(df.metabolon[[st]][[mlevel]]),]
		df.metabolon[[st]][[mlevel]] <- replace
	}
}
for (st in c("Plasma")) {
	for (mlevel in c("BIOCHEMICAL")) {
		# maternal
		sel <- as.character(mapping.full$patid)
		tmp <- df.metabolon.raw[[st]][[mlevel]][, sel]
		tmp2 <- apply(tmp, 1, function(x) {
			x[which(is.na(x))] <- min(x, na.rm=T) # impute to min value
			x <- log(x) # log-transform
			x <- (x-mean(x))/sd(x) # Z-transform
			x
		})
		replace <- tmp2; replace <- replace[rownames(df.metabolon[[st]][[mlevel]]),]
		df.metabolon[[st]][[mlevel]] <- replace
	}
}
metabolon_map <- unique(metabolon_map); rownames(metabolon_map) <- metabolon_map$BIOCHEMICAL
all_colors <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set3"))
cols.superpathway <- c(brewer.pal(9, "Set1"), "#bbbbbb", "#111111")[1:(length(unique(metabolon_map$SUPER.PATHWAY))+1)]; names(cols.superpathway) <- c(unique(metabolon_map$SUPER.PATHWAY), "NOT_METABOLITE")
names(metabolon_map_by_assay) <- c("DBS", "Plasma")



#########################################################################################################
### Cohort demographics and some QC data about metabolomics
# Table 1 for n=76 (DBS)
mapping.demo <- subset(mapping, MaternalRegimen %in% c("untreated", "zdv", "PI-ART"))
mapping.demo$MaternalGroup <- droplevels(mapping.demo$MaternalGroup)
mapping.demo$MaternalRegimen <- droplevels(mapping.demo$MaternalRegimen)
tab1 <- CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "delgage", "hemaval.mom", "InfantAgeInDays", "DaysPTDPlasma", "DaysPTDDBS"), strata=c("MaternalGroup"), data=mapping.demo, smd=T)
write.table(print(tab1), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_1.MaternalGroup.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(print(summary(tab1$ContTable)), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_1_detailed.MaternalGroup.txt", quote=F, sep="\t", row.names=T, col.names=T)
tab1 <- CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "delgage", "hemaval.mom", "InfantAgeInDays", "DaysPTDPlasma", "DaysPTDDBS"), strata=c("Delivery"), data=mapping.demo, smd=T)
write.table(print(tab1), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_1.Delivery.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(print(CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "InfantAgeInDays", "delgage", "hemaval.infant"), strata=c("InfantGroup"), data=mapping.demo, smd=T)), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_1.InfantGroup.txt", quote=F, sep="\t", row.names=T, col.names=T)
for (regi in c("untreated", "zdv", "PI-ART", "other")) {
	write.table(print(CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "GestationalAgeAtCollection", "delgage", "hemaval.mom", "hemaval.infant", "DaysPTDPlasma", "DaysPTDDBS"), strata=c("Delivery"), data=subset(mapping, MaternalRegimen==regi), smd=T)), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/Table_1.MaternalRegimen.%s.txt", regi), quote=F, sep="\t", row.names=T, col.names=T)
	write.table(print(CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "GestationalAgeAtCollection", "delgage", "hemaval.mom", "hemaval.infant"), strata=c("Delivery"), data=subset(mapping, InfantRegimen==regi), smd=T)), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/Table_1.InfantRegimen.%s.txt", regi), quote=F, sep="\t", row.names=T, col.names=T)
}

# Table 1 for n=97 (plasma)
mapping.demo <- subset(mappinglist[["Plasma"]], MaternalRegimen %in% c("untreated", "zdv", "PI-ART"))
mapping.demo$MaternalGroup <- droplevels(mapping.demo$MaternalGroup)
mapping.demo$MaternalRegimen <- droplevels(mapping.demo$MaternalRegimen)
tab1 <- CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "delgage", "hemaval.mom", "InfantAgeInDays", "DaysPTDPlasma", "DaysPTDDBS"), strata=c("MaternalGroup"), data=mapping.demo, smd=T)
write.table(print(tab1), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_1_DBS.MaternalGroup.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(print(summary(tab1$ContTable)), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_1_DBS_detailed.MaternalGroup.txt", quote=F, sep="\t", row.names=T, col.names=T)
tab1 <- CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "delgage", "hemaval.mom", "InfantAgeInDays", "DaysPTDPlasma", "DaysPTDDBS"), strata=c("Delivery"), data=mapping.demo, smd=T)
write.table(print(tab1), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_1_DBS.Delivery.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(print(CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "InfantAgeInDays", "delgage", "hemaval.infant"), strata=c("InfantGroup"), data=mapping.demo, smd=T)), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_1_DBS.InfantGroup.txt", quote=F, sep="\t", row.names=T, col.names=T)
for (regi in c("untreated", "zdv", "PI-ART", "other")) {
	write.table(print(CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "GestationalAgeAtCollection", "delgage", "hemaval.mom", "hemaval.infant", "DaysPTDPlasma", "DaysPTDDBS"), strata=c("Delivery"), data=subset(mapping, MaternalRegimen==regi), smd=T)), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/Table_1_DBS.MaternalRegimen.%s.txt", regi), quote=F, sep="\t", row.names=T, col.names=T)
	write.table(print(CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "GestationalAgeAtCollection", "delgage", "hemaval.mom", "hemaval.infant"), strata=c("Delivery"), data=subset(mapping, InfantRegimen==regi), smd=T)), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/Table_1_DBS.InfantRegimen.%s.txt", regi), quote=F, sep="\t", row.names=T, col.names=T)
}

### number of metabolites detected in each sample type, Venn diagrams
#mlevel <- "BIOCHEMICAL"
#qc <- {}; merged <- data.frame(BIOCHEMICAL=rownames(metabolon_map), detected.maternal_DBS=NA, detected.maternal_plasma=NA, detected.infant_DBS=NA, detected.infant_plasma=NA, median.maternal_DBS=NA, median.maternal_plasma=NA, median.infant_DBS=NA, median.infant_plasma=NA); rownames(merged) <- merged$BIOCHEMICAL
#for (st in c("DBS", "Plasma")) {
#	metabolon <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/OrigScale.%s.txt", st), header=T, as.is=T, sep="\t", quote="", comment.char="")
#	colnames(metabolon) <- gsub("^X", "", colnames(metabolon))
#	metabolon <- subset(metabolon, SUPER.PATHWAY!="") # remove uncharacterized molecules
#	sel <- setdiff(colnames(metabolon), metabolite_levels)
#	tmp <- metabolon[,sel]; rownames(tmp) <- metabolon[, mlevel]; ids <- colnames(tmp)
#	tmp <- t(apply(tmp, 1, function(x) as.numeric(gsub(",", "", x)))); colnames(tmp) <- ids
#	# count detectable as any non-NA value
#	counts.maternal <- apply(tmp[, as.character(mapping$patid)], 1, function(x) length(which(!is.na(x))))
#	if (any(as.character(mapping$cpatid) %in% colnames(tmp))) {
#		counts.infant <- apply(tmp[, as.character(mapping$cpatid)], 1, function(x) length(which(!is.na(x))))
#	} else {
#		counts.infant <- rep(NA, length(counts.maternal))
#	}
#	# summary statistics
#	median.maternal <- apply(log(tmp[, as.character(mapping$patid)]), 1, function(x) median(x,na.rm=T))
#	if (any(as.character(mapping$cpatid) %in% colnames(tmp))) {
#		median.infant <- apply(log(tmp[, as.character(mapping$cpatid)]), 1, function(x) median(x,na.rm=T))
#	} else {
#		median.infant <- rep(NA, length(median.maternal))
#	}
#	out <- data.frame(BIOCHEMICAL=rownames(tmp), subtype=st, detected.maternal=counts.maternal, detected.infant=counts.infant, median.maternal=median.maternal, median.infant=median.infant)
#	if (st == "DBS") {
#		merged[as.character(out$BIOCHEMICAL), c("BIOCHEMICAL", "detected.maternal_DBS", "detected.infant_DBS", "median.maternal_DBS", "median.infant_DBS")] <- out[, c("BIOCHEMICAL", "detected.maternal", "detected.infant", "median.maternal", "median.infant")]
#	} else {
#		merged[as.character(out$BIOCHEMICAL), c("BIOCHEMICAL", "detected.maternal_plasma", "detected.infant_plasma", "median.maternal_plasma", "median.infant_plasma")] <- out[, c("BIOCHEMICAL", "detected.maternal", "detected.infant", "median.maternal", "median.infant")]
#	}
#	write.table(out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/detected_counts.%s.txt", st), quote=F, sep="\t", row.names=F, col.names=T)
#}
#merged$FLAG.maternal_DBS <- ifelse(is.na(merged$detected.maternal_DBS), FALSE, merged$detected.maternal_DBS > 0)
#merged$FLAG.infant_DBS <- ifelse(is.na(merged$detected.infant_DBS), FALSE, merged$detected.infant_DBS > 0)
#merged$FLAG.maternal_plasma <- ifelse(is.na(merged$detected.maternal_plasma), FALSE, merged$detected.maternal_plasma > 0)
#vc <- vennCounts(merged[, c("FLAG.maternal_DBS", "FLAG.infant_DBS", "FLAG.maternal_plasma")])
#vennDiagram(vc, cex=c(1,1,1))
#vc <- vennCounts(merged[, c("FLAG.maternal_DBS", "FLAG.maternal_plasma")])
#vennDiagram(vc, cex=c(1,1))
#write.table(merged, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/detected_counts.%s.txt", "merged"), quote=F, sep="\t", row.names=F, col.names=T)

### ICC and coefficient of variation in paired maternal plasma/DBS samples
#subtype <- "maternal"
#mapping.sel <- mapping; rownames(mapping.sel) <- mapping$patid
#sel.metabolites <- rownames(subset(merged, FLAG.maternal_plasma & FLAG.maternal_DBS))
#data <- matrix()
#data.plasma <- df.metabolon[["Plasma"]][[mlevel]][as.character(mapping$patid), sel.metabolites]
#data.dbs <- df.metabolon[["DBS"]][[mlevel]][as.character(mapping$patid), sel.metabolites]
#data <- abind(data.plasma, data.dbs, along=0); dimnames(data)[[1]] <- c("Plasma", "DBS")
## ICC by metabolite
#final <- apply(data, 3, function(x) {
#	unlist(c(mean(x[,1]), mean(x[,2]), mean(abs(x[,1]-x[,2])), ICC(t(x))$results["Single_raters_absolute", c("ICC", "p", "lower bound", "upper bound")]))
##	icc(t(x))$value
#}); final <- as.data.frame(t(final)); colnames(final) <- c("mean.plasma", "mean.DBS", "mean_difference", "ICC", "p", "lower_bound", "upper_bound")
#final$padj <- p.adjust(final$p, method="fdr")
#final <- final[order(final$p), ]
#write.table(final, file="/Lab_Share/PROMISE/nwcs619/metabolon/ICC.metabolite.txt", quote=F, sep="\t", row.names=T, col.names=T)
#final <- subset(final, mean.plasma > -5 & mean.DBS > -5) # remove a few really low abundance outliers
#test <- cor.test(~ mean.plasma + ICC, final, method="spearman")
#p <- ggplot(final, aes(x=mean.plasma, y=ICC)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.plasma (by metabolite, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
#print(p)
#test <- cor.test(~ mean.plasma + ICC, final, method="spearman")
#p <- ggplot(final, aes(x=mean.DBS, y=ICC)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.DBS (by metabolite, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
#print(p)
#p <- ggplot(final, aes(x=ICC)) + geom_histogram(bins=50) + theme_classic() + ggtitle(sprintf("Distribution of ICC (by metabolite, mean=%.2g, median=%.2g)", mean(final$ICC), median(final$ICC))) + geom_vline(xintercept=mean(final$ICC), color="red")
#print(p)
#thresholds <- seq(from=0.05, to=0.95, by=0.05)
#df <- sapply(thresholds, function(thresh) {
#	length(which(final$ICC >= thresh))
#}); df <- data.frame(threshold=thresholds, n=df); df$pct <- 100*(df$n / nrow(final))
#p <- ggplot(df, aes(x=threshold, y=pct)) + geom_point() + geom_line() + geom_text(label=sprintf("%d (%.2g%%)", df$n, df$pct)) + theme_classic() + ggtitle(sprintf("Number of metabolites at ICC threshold (n=%d padj<0.05)", length(which(final$padj<0.05))))
#print(p)
#p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Number of metabolites at ICC threshold")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#print(p)
#df <- as.matrix(summary(final$ICC))
#p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("ICC by metabolite summary statistics")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#print(p)
## ICC by subject
#final <- apply(data, 2, function(x) {
#	unlist(c(mean(x[,1]), mean(x[,2]), mean(abs(x[,1]-x[,2])), ICC(t(x))$results["Single_raters_absolute", c("ICC", "p", "lower bound", "upper bound")]))
##	icc(t(x))$value
#}); final <- as.data.frame(t(final)); colnames(final) <- c("mean.plasma", "mean.DBS", "mean_difference", "ICC", "p", "lower_bound", "upper_bound")
#final$padj <- p.adjust(final$p, method="fdr")
#final <- final[order(final$p), ]
#write.table(final, file="/Lab_Share/PROMISE/nwcs619/metabolon/ICC.subject.txt", quote=F, sep="\t", row.names=T, col.names=T)
#test <- cor.test(~ mean.plasma + ICC, final, method="spearman")
#p <- ggplot(final, aes(x=mean.plasma, y=ICC)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.plasma (by subject, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
#print(p)
#test <- cor.test(~ mean.plasma + ICC, final, method="spearman")
#p <- ggplot(final, aes(x=mean.DBS, y=ICC)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.DBS (by subject, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
#print(p)
#p <- ggplot(final, aes(x=ICC)) + geom_histogram(bins=50) + theme_classic() + ggtitle(sprintf("Distribution of ICC (by subject, mean=%.2g, median=%.2g)", mean(final$ICC), median(final$ICC))) + geom_vline(xintercept=mean(final$ICC), color="red")
#print(p)
#thresholds <- seq(from=0.05, to=0.95, by=0.05)
#df <- sapply(thresholds, function(thresh) {
#	length(which(final$ICC >= thresh))
#}); df <- data.frame(threshold=thresholds, n=df); df$pct <- 100*(df$n / nrow(final))
#p <- ggplot(df, aes(x=threshold, y=pct)) + geom_point() + geom_line() + geom_text(label=sprintf("%d (%.2g%%)", df$n, df$pct)) + theme_classic() + ggtitle(sprintf("Number of metabolites at ICC threshold (by subject, n=%d padj<0.05)", length(which(final$padj<0.05))))
#print(p)
#p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Number of metabolites at ICC threshold (by subject)")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#print(p)
#df <- as.matrix(summary(final$ICC))
#p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("ICC by metabolite summary statistics (by subject)")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#print(p)
## tSNE by subject with both plasma and DBS
#rownames(data.plasma) <- sprintf("%s.Plasma", rownames(data.plasma)); rownames(data.dbs) <- sprintf("%s.DBS", rownames(data.dbs))
#data2 <- rbind(data.plasma, data.dbs)
#tsne.out <- Rtsne(data2, perplexity=20)
#df <- data.frame(PC1 = tsne.out$Y[,1], PC2 = tsne.out$Y[,2], SampleID=rownames(data2)); rownames(df) <- df$SampleID
#df$SubjectID <- unlist(lapply(rownames(df), function(x) unlist(strsplit(x, "\\."))[1]))
#df$SampleType <- unlist(lapply(rownames(df), function(x) unlist(strsplit(x, "\\."))[2]))
#p <- ggplot(df, aes_string(x="PC1", y="PC2", colour="SampleType")) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("tSNE by subject (%s)", mvar)) + scale_color_brewer(palette="Set1")
#print(p)
#for (mvar in c("Group", "Delivery")) {
#	df[, mvar] <- mapping.sel[df$SubjectID, mvar]
#	p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("tSNE by subject (%s)", mvar)) + scale_color_brewer(palette="Set1")
#	print(p)
#}
## PERMANOVA
#res <- adonis2(data2 ~ SampleType + Group + SubjectID, data=df, permutations=999, method='euclidean')
#sink(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/metabolon_PERMANOVA_by_SampleType.txt"))
#print(res)
#sink()

#########################################################################################################
### maternal metabolites (DBS, Plasma separately)
## ordination (t-SNE, PCA) and PERMANOVA
set.seed(nrow(mapping))
subtype <- "maternal"
for (st in c("DBS", "Plasma")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mappinglist[[st]][,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
	colnames(mapping.sel) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
	rownames(mapping.sel) <- mapping.sel$patid
	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("untreated", "zdv", "PI-ART"))
	data <- data[rownames(mapping.sel),] # subset to just the maternal samples
	to_remove <- names(which(apply(data, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
	data <- data[, setdiff(colnames(data), to_remove)]
	
	# PCA
	pca <- prcomp(data, center=F, scale=F)
	eigs <- pca$sdev^2
	pvar <- 100*(eigs / sum(eigs))
	df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], SampleID=rownames(pca$x))
	for (mvar in intersect(rownames(subset(metadata_variables, useForPERMANOVA=="yes")), colnames(mapping.sel))) {
		df[, mvar] <- mapping.sel[rownames(df), mvar]
		if (metadata_variables[mvar, "type"] == "factor") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s %s %s PCoA (Euclidean distance)", subtype, st, mlevel, mvar)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_brewer(palette="Set1") + stat_ellipse(type="t")
		} else if (metadata_variables[mvar, "type"] == "numeric") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s %s %s PCoA (Euclidean distance)", subtype, st, mlevel, mvar)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_gradient(low="grey", high="red")
		}
		else {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s %s %s PCoA (Euclidean distance)", subtype, st, mlevel, mvar)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + stat_ellipse(type="t")
		}
		print(p)
	}
	p <- ggplot(df, aes(x=PC1, y=PC2, colour=MaternalRegimen, shape=Delivery, group=MaternalRegimen)) + geom_point(size=2) + theme_classic() + ggtitle(sprintf("%s %s %s PCoA (Euclidean distance)", subtype, st, mlevel)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_brewer(palette="Set1") + stat_ellipse() + scale_shape_manual(values=c(1,3))
	print(p)
	
	# PERMANOVA
	res <- adonis2(data ~ Delivery + MaternalRegimen + Country + GestationalAgeAtCollection, data=mapping.sel, permutations=999, method='euclidean')
	sink(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/metabolon_PERMANOVA.%s.%s.%s.txt", subtype, st, mlevel))
	print(res)
	sink()
	# t-SNE
	tsne.out <- Rtsne(data, perplexity=20)
	df <- data.frame(PC1 = tsne.out$Y[,1], PC2 = tsne.out$Y[,2], SampleID=rownames(data)); rownames(df) <- df$SampleID
	for (mvar in intersect(rownames(subset(metadata_variables, useForPERMANOVA=="yes")), colnames(mapping.sel))) {
		df[, mvar] <- mapping.sel[rownames(df), mvar]
		if (metadata_variables[mvar, "type"] == "factor") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s %s %s tSNE", subtype, st, mlevel, mvar)) + scale_color_brewer(palette="Set1") + stat_ellipse(type="t")
		} else if (metadata_variables[mvar, "type"] == "numeric") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s %s %s tSNE", subtype, st, mlevel, mvar)) + scale_color_gradient(low="grey", high="red")
		}
		else {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s %s %s tSNE", subtype, st, mlevel, mvar)) + stat_ellipse(type="t")
		}
		print(p)
	}
}

## ordination (t-SNE, PCA) and PERMANOVA, with Plasma+DBS
## NOTE: to be fixed, rows should have Plasma and DBS separately, then need to have columns be intersection of metabolites found in both Plasma and DBS
set.seed(nrow(mapping))
subtype <- "maternal"
mlevel <- "BIOCHEMICAL"
mappinglist.sel <- lapply(mappinglist, function(x) subset(x, MaternalRegimen %in% c("untreated", "zdv", "PI-ART")))
mapping.dbs <- mappinglist.sel[["DBS"]]; mapping.dbs$SampleType <- "DBS"; rownames(mapping.dbs) <- sprintf("%s.DBS", mapping.dbs$patid)
mapping.plasma <- mappinglist.sel[["Plasma"]]; mapping.plasma$SampleType <- "Plasma"; rownames(mapping.plasma) <- sprintf("%s.Plasma", mapping.plasma$patid)
mapping.combined <- rbind(mapping.dbs, mapping.plasma)
data <- {}
sel.metabolites <- intersect(rownames(metabolon_map_by_assay[["DBS"]]), rownames(metabolon_map_by_assay[["Plasma"]]))
for (st in c("DBS", "Plasma")) {
	tmp <- df.metabolon[[st]][[mlevel]]
	if (st=="DBS") {
		sel_subjects <- mapping.dbs$patid
	} else {
		sel_subjects <- mapping.plasma$patid
	}
	sel_subjects <- as.character(sel_subjects)
	tmp <- tmp[sel_subjects,]
	tmp <- tmp[, intersect(colnames(tmp), sel.metabolites)] # only use metabolites found in both DBS and Plasma
	rownames(tmp) <- sprintf("%s.%s", rownames(tmp), st)
	data <- rbind(data, tmp)
}
to_remove <- names(which(apply(data, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
data <- data[, setdiff(colnames(data), to_remove)]
# PCA
pca <- prcomp(data, center=F, scale=F)
eigs <- pca$sdev^2
pvar <- 100*(eigs / sum(eigs))
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
df <- merge(df, mapping.combined, by="row.names")
df$SampleType2 <- paste(df$SampleType, df$Delivery, sep="")
p <- ggplot(df, aes(x=PC1, y=PC2, colour=Delivery, shape=SampleType, group=Delivery)) + geom_point(size=2) + theme_classic() + ggtitle(sprintf("%s %s %s PCoA (Euclidean distance)", subtype, "Plasma+DBS", mlevel)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_brewer(palette="Set1") + stat_ellipse() + scale_shape_manual(values=c(1,16))
print(p)
p <- ggplot(df, aes(x=PC1, y=PC2, colour=SampleType, shape=Delivery, group=SampleType)) + geom_point(size=2) + theme_classic() + ggtitle(sprintf("%s %s %s PCoA (Euclidean distance)", subtype, "Plasma+DBS", mlevel)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_brewer(palette="Set1") + stat_ellipse() + scale_shape_manual(values=c(1,3))
print(p)
p <- ggplot(df, aes(x=PC1, y=PC2, colour=MaternalRegimen, shape=SampleType, group=MaternalRegimen)) + geom_point(size=2) + theme_classic() + ggtitle(sprintf("%s %s %s PCoA (Euclidean distance)", subtype, "Plasma+DBS", mlevel)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_brewer(palette="Set1") + stat_ellipse() + scale_shape_manual(values=c(1,16))
print(p)
p <- ggplot(df, aes(x=PC1, y=PC2, colour=MaternalRegimen, shape=SampleType2, group=MaternalRegimen)) + geom_point(size=2) + theme_classic() + ggtitle(sprintf("%s %s %s PCoA (Euclidean distance)", subtype, "Plasma+DBS", mlevel)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_brewer(palette="Set1") + stat_ellipse() + scale_shape_manual(values=c(1,16,2,17))
print(p)
p <- ggplot(df, aes(x=PC1, y=PC2, colour=MaternalRegimen, shape=Delivery, group=MaternalRegimen)) + geom_point(size=2) + facet_wrap(~SampleType) + theme_classic() + ggtitle(sprintf("%s %s %s PCoA (Euclidean distance)", subtype, "Plasma+DBS", mlevel)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_brewer(palette="Set1") + stat_ellipse() + scale_shape_manual(values=c(1,16))
print(p)

# PERMANOVA
res <- adonis2(data ~ Delivery + MaternalRegimen + Country + GestationalAgeAtCollection + SampleType, data=mapping.combined, permutations=999, method='euclidean')
sink(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/metabolon_PERMANOVA.%s.%s.%s.txt", subtype, "PlasmaDBS", mlevel))
print(res)
sink()

### calculation of effect size, outcome Delivery, stratified by Regimen
#subtype <- "maternal"; mvar <- "Delivery"
#for (st in c("DBS", "Plasma")) {
#	mlevel <- "BIOCHEMICAL"
#	data <- df.metabolon[[st]][[mlevel]]
#	mapping2 <- mapping[,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
#	colnames(mapping2) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
#	rownames(mapping2) <- mapping2$patid
#	for (regi in setdiff(levels(mapping2$MaternalRegimen), "other")) {
#		mapping.sel <- subset(mapping2, MaternalRegimen == regi); mapping.sel$MaternalRegimen <- droplevels(mapping.sel$MaternalRegimen) # exclude 'other' regimen because too few samples
#		data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
#		name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#		sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
#		df <- merge(data.sel, mapping.sel, by="row.names"); 
#		df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
#		res <- {}
#		for (metabolite in colnames(data.sel)) {
#			tmp <- df[, c(metabolite, mvar)]
#			m1 <- mean(subset(tmp, Delivery=="Preterm")[, metabolite])
#			m2 <- mean(subset(tmp, Delivery=="Term")[, metabolite])
#			sd1 <- sd(subset(tmp, Delivery=="Preterm")[, metabolite])
#			sd2 <- sd(subset(tmp, Delivery=="Term")[, metabolite])
#			d <- abs(m1 - m2) / sqrt((sd1^2 + sd2^2)/2)
#			tmp <- data.frame(metabolite=name_map[metabolite, "original"], d=d)
#			res <- rbind(res, tmp)
#		}
#		df <- as.matrix(summary(res$d))
#		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Estimated effect size (%s %s %s, %s)", subtype, mlevel, st, regi)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#		print(p)
#		df <- melt(quantile(res$d, probs=seq(from=0,to=1,by=0.05)))
#		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Estimated effect size (%s %s %s, %s)", subtype, mlevel, st, regi)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#		print(p)
#	}
#}

### calculation of effect size, outcome Delivery, averaged across Regimen
#subtype <- "maternal"; mvar <- "Delivery"
#for (st in c("DBS", "Plasma")) {
#	mlevel <- "BIOCHEMICAL"
#	data <- df.metabolon[[st]][[mlevel]]
#	mapping.sel <- mapping[,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
#	colnames(mapping.sel) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
#	rownames(mapping.sel) <- mapping.sel$patid
#	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("untreated", "zdv", "PI-ART")); mapping.sel$MaternalRegimen <- droplevels(mapping.sel$MaternalRegimen) # exclude 'other' regimen because too few samples
#	data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
#	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
#	df <- merge(data.sel, mapping.sel, by="row.names"); 
#	df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
#	res <- {}
#	for (metabolite in colnames(data.sel)) {
#		tmp <- df[, c(metabolite, mvar)]
#		m1 <- mean(subset(tmp, Delivery=="Preterm")[, metabolite])
#		m2 <- mean(subset(tmp, Delivery=="Term")[, metabolite])
#		sd1 <- sd(subset(tmp, Delivery=="Preterm")[, metabolite])
#		sd2 <- sd(subset(tmp, Delivery=="Term")[, metabolite])
#		d <- abs(m1 - m2) / sqrt((sd1^2 + sd2^2)/2)
#		tmp <- data.frame(metabolite=name_map[metabolite, "original"], d=d)
#		res <- rbind(res, tmp)
#	}
#	df <- as.matrix(summary(res$d))
#	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Estimated effect size (%s %s %s, all regimens)", subtype, mlevel, st)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#	print(p)
#	df <- melt(quantile(res$d, probs=seq(from=0,to=1,by=0.05)))
#	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Estimated effect size (%s %s %s, all regimens)", subtype, mlevel, st)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#	print(p)
#}

## linear model with emmeans, outcome Delivery, stratified by Regimen
subtype <- "maternal"; mvar <- "Delivery"
for (st in c("DBS", "Plasma")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mappinglist[[st]][,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "weights.maternal")]
	colnames(mapping.sel) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weights.maternal")
	rownames(mapping.sel) <- mapping.sel$patid
	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("untreated", "zdv", "PI-ART")); mapping.sel$MaternalRegimen <- droplevels(mapping.sel$MaternalRegimen) # exclude 'other' regimen because too few samples
	data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
	to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
	data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
	df <- merge(data.sel, mapping.sel, by="row.names"); 
	df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
	res <- {}
	for (metabolite in colnames(data.sel)) {
		mod <- lm(as.formula(sprintf("%s ~ %s*MaternalRegimen", metabolite, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s | MaternalRegimen", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
		res <- rbind(res, tmp)
	}
	res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
	res$padj <- p.adjust(res$p.value, method="fdr")
	res <- res[order(res$p.value, decreasing=F),]
	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
	res$exp_estimate <- exp(res$estimate)
	for (addvar in c("COMP_ID", "PLATFORM")) {
		res[, addvar] <- metabolon_map_by_assay[[st]][as.character(res$metabolite), addvar]
	}
#	res$dir2 = ifelse(res$estimate < 0, "down", "up"); tab <- table(res[,c("Regimen","dir2")])
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
	# forest plot
	resSig <- subset(res, padj < siglevel)
	df <- subset(res, metabolite %in% as.character(resSig$metabolite))
	df <- df[order(df$estimate, decreasing=T),]
	df$metabolite <- factor(as.character(df$metabolite), levels=as.character(unique(df$metabolite)))
	lims <- max(abs(df$estimate) + abs(df$SE))*1.0
	pd <- position_dodge(0.8)
	p <- ggplot(df, aes(x=metabolite, y=estimate, color=dir, group=MaternalRegimen)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=MaternalRegimen), position=pd, hjust=1, color="black", size=2) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("Delivery stratified by MaternalRegimen (%s, %s, %s)", subtype, st, mlevel)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
	print(p)
}
# plot DBS and plasma together
res <- {}
subtype <- "maternal"; mvar <- "Delivery"; mlevel <- "BIOCHEMICAL"
for (st in c("DBS", "Plasma")) {
	tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), sep="\t", header=T, as.is=T, quote="", comment.char="")
	tmp$SampleType <- st
	tmp$metabolite <- sprintf("[%s] %s", tmp$SampleType, tmp$metabolite)
	res <- rbind(res, tmp)
}
resSig <- subset(res, padj < siglevel)
df <- subset(res, metabolite %in% as.character(resSig$metabolite))
df <- df[order(df$estimate, decreasing=T),]
df$metabolite <- factor(as.character(df$metabolite), levels=as.character(unique(df$metabolite)))
lims <- max(abs(df$estimate) + abs(df$SE))*1.0
pd <- position_dodge(0.8)
p <- ggplot(df, aes(x=metabolite, y=estimate, color=dir, group=MaternalRegimen)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=MaternalRegimen), position=pd, hjust=1, color="black", size=2) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("Delivery stratified by MaternalRegimen (%s, %s, %s)", subtype, "DBS+Plasma", mlevel)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
print(p)


## linear model with emmeans, outcome Delivery2, stratified by Regimen
subtype <- "maternal"; mvar <- "Delivery2"
for (st in c("DBS", "Plasma")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mappinglist[[st]][,c("Delivery2", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
	colnames(mapping.sel) <- c("Delivery2", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
	rownames(mapping.sel) <- mapping.sel$patid
	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("untreated", "zdv", "PI-ART")); mapping.sel$MaternalRegimen <- droplevels(mapping.sel$MaternalRegimen) # exclude 'other' regimen because too few samples
	data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
	to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
	data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
	df <- merge(data.sel, mapping.sel, by="row.names"); 
#	df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
	res <- {}
	for (metabolite in colnames(data.sel)) {
		mod <- lm(as.formula(sprintf("%s ~ %s*MaternalRegimen", metabolite, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s | MaternalRegimen", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
		res <- rbind(res, tmp)
	}
	res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
	res$padj <- p.adjust(res$p.value, method="fdr")
	res <- res[order(res$p.value, decreasing=F),]
	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
	res$exp_estimate <- exp(res$estimate)
	for (addvar in c("COMP_ID", "PLATFORM")) {
		res[, addvar] <- metabolon_map_by_assay[[st]][as.character(res$metabolite), addvar]
	}
#	res$dir2 = ifelse(res$estimate < 0, "down", "up"); tab <- table(res[,c("Regimen","dir2")])
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_Delivery2_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
	# forest plot
	resSig <- subset(res, padj < siglevel)
	df <- subset(res, metabolite %in% as.character(resSig$metabolite))
	if (nrow(df)>0) {
		df <- df[order(df$estimate, decreasing=T),]
		df$metabolite <- factor(as.character(df$metabolite), levels=as.character(unique(df$metabolite)))
		lims <- max(abs(df$estimate) + abs(df$SE))*1.0
		pd <- position_dodge(0.8)
		p <- ggplot(df, aes(x=metabolite, y=estimate, color=dir, group=MaternalRegimen)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims*0.9, label=MaternalRegimen), position=pd, hjust=1, color="black", size=2) + facet_wrap(~contrast) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("Delivery2 stratified by MaternalRegimen (%s, %s, %s)", subtype, st, mlevel)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
		print(p)
	}
}
# plot DBS and plasma together
res <- {}
subtype <- "maternal"; mvar <- "Delivery2"; mlevel <- "BIOCHEMICAL"
for (st in c("DBS", "Plasma")) {
	tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_Delivery2_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), sep="\t", header=T, as.is=T, quote="", comment.char="")
	tmp$SampleType <- st
	tmp$metabolite <- sprintf("[%s] %s", tmp$SampleType, tmp$metabolite)
	res <- rbind(res, tmp)
}
resSig <- subset(res, padj < siglevel)
df <- subset(res, metabolite %in% as.character(resSig$metabolite))
df <- df[order(df$estimate, decreasing=T),]
df$metabolite <- factor(as.character(df$metabolite), levels=as.character(unique(df$metabolite)))
lims <- max(abs(df$estimate) + abs(df$SE))*1.0
pd <- position_dodge(0.8)
p <- ggplot(df, aes(x=metabolite, y=estimate, color=dir, group=MaternalRegimen)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims*0.9, label=MaternalRegimen), position=pd, hjust=1, color="black", size=2) + facet_wrap(~contrast) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("Delivery2 stratified by MaternalRegimen (%s, %s, %s)", subtype, "DBS+Plasma", mlevel)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
print(p)


## linear model with emmeans, outcome Delivery, averaged across Regimen
subtype <- "maternal"; mvar <- "Delivery"
for (st in c("DBS", "Plasma")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mappinglist[[st]][,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
	colnames(mapping.sel) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
	rownames(mapping.sel) <- mapping.sel$patid
	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("untreated", "zdv", "PI-ART")); mapping.sel$MaternalRegimen <- droplevels(mapping.sel$MaternalRegimen) # exclude 'other' regimen because too few samples
	data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
	to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
	data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
	df <- merge(data.sel, mapping.sel, by="row.names"); 
	df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
	res <- {}
	for (metabolite in colnames(data.sel)) {
		mod <- lm(as.formula(sprintf("%s ~ %s+MaternalRegimen", metabolite, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
		res <- rbind(res, tmp)
	}
	res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
	res$padj <- p.adjust(res$p.value, method="fdr")
	res <- res[order(res$p.value, decreasing=F),]
	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
	res$exp_estimate <- exp(res$estimate)
	for (addvar in c("COMP_ID", "PLATFORM")) {
		res[, addvar] <- metabolon_map_by_assay[[st]][as.character(res$metabolite), addvar]
	}
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_averaged_across_Regimen.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
	# forest plot
	resSig <- subset(res, padj < siglevel)
	df <- subset(res, metabolite %in% as.character(resSig$metabolite))
	df <- df[order(df$estimate, decreasing=T),]
	df$metabolite <- factor(as.character(df$metabolite), levels=as.character(unique(df$metabolite)))
	lims <- max(abs(df$estimate) + abs(df$SE))*1.0
	pd <- position_dodge(0.8)
	p <- ggplot(df, aes(x=metabolite, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("Delivery averaged across MaternalRegimen (%s, %s, %s)", subtype, st, mlevel)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
	print(p)
}

## linear model with emmeans, outcome Delivery, stratified by Regimen+Country
subtype <- "maternal"; mvar <- "Delivery"
for (st in c("DBS", "Plasma")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mappinglist[[st]][,c("Delivery", "MaternalRegimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
	colnames(mapping.sel) <- c("Delivery", "MaternalRegimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
	rownames(mapping.sel) <- mapping.sel$patid
	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("untreated", "zdv", "PI-ART")); mapping.sel$MaternalRegimen <- droplevels(mapping.sel$MaternalRegimen) # exclude 'other' regimen because too few samples
	data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
	to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
	data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
	df <- merge(data.sel, mapping.sel, by="row.names"); 
	df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
	res <- {}
	for (metabolite in colnames(data.sel)) {
		mod <- lm(as.formula(sprintf("%s ~ %s*MaternalRegimen + %s*Country", metabolite, mvar, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s | MaternalRegimen + Country", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
		res <- rbind(res, tmp)
	}
	res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
	res$padj <- p.adjust(res$p.value, method="fdr")
	res <- res[order(res$p.value, decreasing=F),]
	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
#	res$dir2 = ifelse(res$estimate < 0, "down", "up"); tab <- table(res[,c("Regimen","dir2")])
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_RegimenCountry.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
}


### linear model with emmeans, outcome Regimen, stratified by Delivery
#subtype <- "maternal"; mvar <- "Regimen"
#for (st in c("DBS", "Plasma")) {
#	mlevel <- "BIOCHEMICAL"
#	data <- df.metabolon[[st]][[mlevel]]
#	mapping.sel <- mappinglist[[st]][,c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
#	colnames(mapping.sel) <- c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
#	rownames(mapping.sel) <- mapping.sel$patid
#	data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
#	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
#	df <- merge(data.sel, mapping.sel, by="row.names"); 
#	df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
#	res <- {}
#	for (metabolite in colnames(data.sel)) {
#		mod <- lm(as.formula(sprintf("%s ~ %s*Delivery", metabolite, mvar)), data=df); modelstr <- "LM"
#		emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s | Delivery", mvar)), adjust="none")
#		tmp <- as.data.frame(emm$contrasts) 
#		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
#		res <- rbind(res, tmp)
#	}
#	res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
#	res$padj <- p.adjust(res$p.value, method="fdr")
#	res <- res[order(res$p.value, decreasing=F),]
#	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
##	res$dir2 = ifelse(res$estimate < 0, "down", "up"); tab <- table(res[,c("Delivery","dir2")])
#	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Delivery.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
#	# forest plot
#	resSig <- subset(res, padj < siglevel)
#	df <- subset(res, metabolite %in% as.character(resSig$metabolite))
#	df <- df[order(df$estimate, decreasing=T),]
#	df$metabolite <- factor(as.character(df$metabolite), levels=as.character(unique(df$metabolite)))
#	lims <- max(abs(df$estimate) + abs(df$SE))*1.0
#	pd <- position_dodge(0.8)
#	p <- ggplot(df, aes(x=metabolite, y=estimate, color=dir, group=Delivery)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Delivery), position=pd, hjust=1, color="black", size=2) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("Regimen stratified by Delivery (%s, %s, %s)", subtype, st, mlevel)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
#	print(p)
#}

### linear model with emmeans, outcome Regimen, averaged across Delivery
#subtype <- "maternal"; mvar <- "Regimen"
#for (st in c("DBS", "Plasma")) {
#	mlevel <- "BIOCHEMICAL"
#	data <- df.metabolon[[st]][[mlevel]]
#	mapping.sel <- mappinglist[[st]][,c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
#	colnames(mapping.sel) <- c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
#	rownames(mapping.sel) <- mapping.sel$patid
#	data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
#	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
#	df <- merge(data.sel, mapping.sel, by="row.names"); 
#	df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
#	res <- {}
#	for (metabolite in colnames(data.sel)) {
#		mod <- lm(as.formula(sprintf("%s ~ %s+Delivery", metabolite, mvar)), data=df); modelstr <- "LM"
#		emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s", mvar)), adjust="none")
#		tmp <- as.data.frame(emm$contrasts) 
#		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
#		res <- rbind(res, tmp)
#	}
#	res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
#	res$padj <- p.adjust(res$p.value, method="fdr")
#	res <- res[order(res$p.value, decreasing=F),]
#	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
#	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_averaged_across_Delivery.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
#	# forest plot
#	resSig <- subset(res, padj < siglevel)
#	df <- subset(res, metabolite %in% as.character(resSig$metabolite))
#	df <- df[order(df$estimate, decreasing=T),]
#	df$metabolite <- factor(as.character(df$metabolite), levels=as.character(unique(df$metabolite)))
#	lims <- max(abs(df$estimate) + abs(df$SE))*1.0
#	pd <- position_dodge(0.8)
#	p <- ggplot(df, aes(x=metabolite, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("Regimen averaged across Delivery (%s, %s, %s)", subtype, st, mlevel)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
#	print(p)
#}


## mediation analysis (IV=ART, DV=PTB, mediator=microbiome; ZDV vs PI-ART)
subtype <- "maternal"; iv <- "MaternalRegimen"; dv <- "Delivery"
set.seed(nrow(mapping))
for (st in c("DBS", "Plasma")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mappinglist[[st]][,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "weights.maternal")]
	colnames(mapping.sel) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weights.maternal")
	rownames(mapping.sel) <- mapping.sel$patid
	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("untreated", "zdv", "PI-ART")); mapping.sel$MaternalRegimen <- droplevels(mapping.sel$MaternalRegimen) # exclude 'other' regimen because too few samples
	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("zdv", "PI-ART")); mapping.sel$MaternalRegimen <- droplevels(mapping.sel$MaternalRegimen) # subset to just ZDV/PI-ART comparison
	data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
	to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
	data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
	df <- merge(data.sel, mapping.sel, by="row.names"); 
	res <- {}
	
	for (iv_level in c("PI-ART")) {
#	for (iv_level in setdiff(levels(mapping.sel[, iv]), "untreated")) {
		for (metabolite in colnames(data.sel)) {
			fit.mediator <- lm(as.formula(sprintf("%s ~ %s", metabolite, iv)), data=df, weights=df$weights.maternal)
			fit.dv <- glm(as.formula(sprintf("%s ~ %s + %s", dv, iv, metabolite)), family="binomial", data=df, weights=df$weights.maternal)
#			med2 <- try(mediation::mediate(fit.mediator, fit.dv, treat=iv, mediator=metabolite, boot=T, control.value="untreated", treat.value=iv_level, parallel="multicore", ncpus=16), silent=T)
			med2 <- try(mediation::mediate(fit.mediator, fit.dv, treat=iv, mediator=metabolite, boot=T, control.value="zdv", treat.value=iv_level, parallel="multicore", ncpus=16), silent=T)
			if (class(med2)=="mediate") {
				summ <- summary(med2)
				res <- rbind(res, c(metabolite, iv_level, unlist(summ[c("d0", "d0.ci", "d0.p", "d1", "d1.ci", "d1.p", "z0", "z0.ci", "z0.p", "z1", "z1.ci", "z1.p", "d.avg", "d.avg.ci", "d.avg.p", "z.avg", "z.avg.ci", "z.avg.p", "tau.coef", "tau.ci", "tau.p", "n0", "n0.ci", "n0.p", "n1", "n1.ci", "n1.p", "n.avg", "n.avg.ci", "n.avg.p")]) ))
			} else {
				res <- rbind(res, c(metabolite, iv_level, rep(NA, 40)))
			}
		}
	}
	res <- as.data.frame(res)
	colnames(res) <- c("metabolite", "Regimen", "ACME.control", "ACME.control.ci.lower", "ACME.control.ci.upper", "ACME.control.p", "ACME.treated", "ACME.treated.ci.lower", "ACME.treated.ci.upper", "ACME.treated.p", "ADE.control", "ADE.control.ci.lower", "ADE.control.ci.upper", "ADE.control.p", "ADE.treated", "ADE.treated.ci.lower", "ADE.treated.ci.upper", "ADE.treated.p", "ACME.average", "ACME.average.ci.lower", "ACME.average.ci.upper", "ACME.average.p", "ADE.average", "ADE.average.ci.lower", "ADE.average.ci.upper", "ADE.average.ci.p", "totaleffect", "totaleffect.ci.lower", "totaleffect.ci.upper", "totaleffect.p", "propmediated.control.avg", "propmediated.control.avg.ci.lower", "propmediated.control.avg.ci.upper", "propmediated.control.p", "propmediated.treatment.avg", "propmediated.treatment.avg.ci.lower", "propmediated.treatment.avg.ci.upper", "propmediated.treatment.p", "propmediated.avg", "propmediated.avg.ci.lower", "propmediated.avg.ci.upper", "propmediated.p")
	res$metabolite <- as.character(name_map[as.character(res$metabolite), "original"])
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/mediation_analysis.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
}


## mediation analysis (IV=ART, DV=PTB, mediator=microbiome; untreated baseline)
subtype <- "maternal"; iv <- "MaternalRegimen"; dv <- "Delivery"
set.seed(nrow(mapping))
for (st in c("DBS", "Plasma")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mappinglist[[st]][,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "weights.maternal")]
	colnames(mapping.sel) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weights.maternal")
	rownames(mapping.sel) <- mapping.sel$patid
	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("untreated", "zdv", "PI-ART")); mapping.sel$MaternalRegimen <- droplevels(mapping.sel$MaternalRegimen) # exclude 'other' regimen because too few samples
#	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("zdv", "PI-ART")); mapping.sel$MaternalRegimen <- droplevels(mapping.sel$MaternalRegimen) # subset to just ZDV/PI-ART comparison
	data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
	to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
	data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
	df <- merge(data.sel, mapping.sel, by="row.names"); 
	res <- {}
	
#	for (iv_level in c("PI-ART")) {
	for (iv_level in setdiff(levels(mapping.sel[, iv]), "untreated")) {
		for (metabolite in colnames(data.sel)) {
			fit.mediator <- lm(as.formula(sprintf("%s ~ %s", metabolite, iv)), data=df, weights=df$weights.maternal)
			fit.dv <- glm(as.formula(sprintf("%s ~ %s + %s", dv, iv, metabolite)), family="binomial", data=df, weights=df$weights.maternal)
			med2 <- try(mediation::mediate(fit.mediator, fit.dv, treat=iv, mediator=metabolite, boot=T, control.value="untreated", treat.value=iv_level, parallel="multicore", ncpus=16), silent=T)
#			med2 <- try(mediation::mediate(fit.mediator, fit.dv, treat=iv, mediator=metabolite, boot=T, control.value="zdv", treat.value=iv_level, parallel="multicore", ncpus=16), silent=T)
			if (class(med2)=="mediate") {
				summ <- summary(med2)
				res <- rbind(res, c(metabolite, iv_level, unlist(summ[c("d0", "d0.ci", "d0.p", "d1", "d1.ci", "d1.p", "z0", "z0.ci", "z0.p", "z1", "z1.ci", "z1.p", "d.avg", "d.avg.ci", "d.avg.p", "z.avg", "z.avg.ci", "z.avg.p", "tau.coef", "tau.ci", "tau.p", "n0", "n0.ci", "n0.p", "n1", "n1.ci", "n1.p", "n.avg", "n.avg.ci", "n.avg.p")]) ))
			} else {
				res <- rbind(res, c(metabolite, iv_level, rep(NA, 40)))
			}
		}
	}
	res <- as.data.frame(res)
	colnames(res) <- c("metabolite", "Regimen", "ACME.control", "ACME.control.ci.lower", "ACME.control.ci.upper", "ACME.control.p", "ACME.treated", "ACME.treated.ci.lower", "ACME.treated.ci.upper", "ACME.treated.p", "ADE.control", "ADE.control.ci.lower", "ADE.control.ci.upper", "ADE.control.p", "ADE.treated", "ADE.treated.ci.lower", "ADE.treated.ci.upper", "ADE.treated.p", "ACME.average", "ACME.average.ci.lower", "ACME.average.ci.upper", "ACME.average.p", "ADE.average", "ADE.average.ci.lower", "ADE.average.ci.upper", "ADE.average.ci.p", "totaleffect", "totaleffect.ci.lower", "totaleffect.ci.upper", "totaleffect.p", "propmediated.control.avg", "propmediated.control.avg.ci.lower", "propmediated.control.avg.ci.upper", "propmediated.control.p", "propmediated.treatment.avg", "propmediated.treatment.avg.ci.lower", "propmediated.treatment.avg.ci.upper", "propmediated.treatment.p", "propmediated.avg", "propmediated.avg.ci.lower", "propmediated.avg.ci.upper", "propmediated.p")
	res$metabolite <- as.character(name_map[as.character(res$metabolite), "original"])
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/mediation_analysis_untreated_baseline.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
}


### mediation analysis (IV=ART, DV=PTB, mediator=microbiome; use mediation package)
#subtype <- "maternal"; iv <- "MaternalRegimen"; dv <- "Delivery"
#set.seed(nrow(mapping))
#for (st in c("DBS", "Plasma")) {
#	mlevel <- "BIOCHEMICAL"
#	data <- df.metabolon[[st]][[mlevel]]
#	mapping.sel <- mappinglist[[st]][,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
#	colnames(mapping.sel) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
#	rownames(mapping.sel) <- mapping.sel$patid
#	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("untreated", "zdv", "PI-ART")); mapping.sel$MaternalRegimen <- droplevels(mapping.sel$MaternalRegimen) # exclude 'other' regimen because too few samples
##	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("zdv", "PI-ART")); mapping.sel$MaternalRegimen <- droplevels(mapping.sel$MaternalRegimen) # subset to just ZDV/PI-ART comparison
#	data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
#	to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
#	data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
#	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
##	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
#	pca <- prcomp(data.sel, center=T, scale=T)
#	data.sel <- pca$x[,1:10]
#	df <- merge(data.sel, mapping.sel, by="row.names")
#	res <- {}
#	
##	for (iv_level in c("PI-ART")) {
#	for (iv_level in setdiff(levels(mapping.sel[, iv]), "untreated")) {
#		for (metabolite in colnames(data.sel)) {
#			fit.mediator <- lm(as.formula(sprintf("%s ~ %s", metabolite, iv)), data=df)
#			fit.dv <- glm(as.formula(sprintf("%s ~ %s + %s", dv, iv, metabolite)), family="binomial", data=df)
#			med2 <- try(mediation::mediate(fit.mediator, fit.dv, treat=iv, mediator=metabolite, boot=T, control.value="untreated", treat.value=iv_level, parallel="multicore", ncpus=16), silent=T)
##			med2 <- try(mediation::mediate(fit.mediator, fit.dv, treat=iv, mediator=metabolite, boot=T, control.value="zdv", treat.value=iv_level, parallel="multicore", ncpus=16), silent=T)
#			if (class(med2)=="mediate") {
#				summ <- summary(med2)
#				res <- rbind(res, c(metabolite, iv_level, unlist(summ[c("d0", "d0.ci", "d0.p", "d1", "d1.ci", "d1.p", "z0", "z0.ci", "z0.p", "z1", "z1.ci", "z1.p", "d.avg", "d.avg.ci", "d.avg.p", "z.avg", "z.avg.ci", "z.avg.p", "tau.coef", "tau.ci", "tau.p", "n0", "n0.ci", "n0.p", "n1", "n1.ci", "n1.p", "n.avg", "n.avg.ci", "n.avg.p")]) ))
#			} else {
#				res <- rbind(res, c(metabolite, iv_level, rep(NA, 40)))
#			}
#		}
#	}
#	res <- as.data.frame(res)
#	colnames(res) <- c("metabolite", "Regimen", "ACME.control", "ACME.control.ci.lower", "ACME.control.ci.upper", "ACME.control.p", "ACME.treated", "ACME.treated.ci.lower", "ACME.treated.ci.upper", "ACME.treated.p", "ADE.control", "ADE.control.ci.lower", "ADE.control.ci.upper", "ADE.control.p", "ADE.treated", "ADE.treated.ci.lower", "ADE.treated.ci.upper", "ADE.treated.p", "ACME.average", "ACME.average.ci.lower", "ACME.average.ci.upper", "ACME.average.p", "ADE.average", "ADE.average.ci.lower", "ADE.average.ci.upper", "ADE.average.ci.p", "totaleffect", "totaleffect.ci.lower", "totaleffect.ci.upper", "totaleffect.p", "propmediated.control.avg", "propmediated.control.avg.ci.lower", "propmediated.control.avg.ci.upper", "propmediated.control.p", "propmediated.treatment.avg", "propmediated.treatment.avg.ci.lower", "propmediated.treatment.avg.ci.upper", "propmediated.treatment.p", "propmediated.avg", "propmediated.avg.ci.lower", "propmediated.avg.ci.upper", "propmediated.p")
#	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/mediation_analysis_PCs.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
#}


## LASSO regression of MaternalGroup (separately by Regimen); using METABOLITE data [DBS, Plasma]
ncvreps <- 100
set.seed(nrow(mapping))
subtype <- "maternal"; mvar <- "Delivery"
for (st in c("DBS", "Plasma")) {
	for (mlevel in "BIOCHEMICAL") {
		data <- df.metabolon[[st]][[mlevel]]
		mapping2 <- mappinglist[[st]][,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "weight0week", "DaysPTDPlasma2")]
		colnames(mapping2) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weight0week", "DaysPTDPlasma2")
		rownames(mapping2) <- mapping2$patid
		res.mean <- {}; res.sd <- {}
		for (regi in setdiff(levels(mapping2$MaternalRegimen), "other")) {
			mapping.sel <- subset(mapping2, MaternalRegimen==regi); mapping.sel$MaternalGroup <- droplevels(mapping.sel$MaternalGroup)
			data.sel <- data[rownames(mapping.sel),] # subset to just the desired maternal samples from regimen
			to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
			data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
			response <- droplevels(mapping.sel$MaternalGroup); names(response) <- rownames(mapping.sel)
			# add Country as covariates
#			data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
			agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
#			fit <- glmnet(x = data.sel, y = response, alpha = 1, family="binomial")
#			cvfit <- cv.glmnet(data.sel, response, family="binomial", type.measure="class")
			
			# do CV 100 reps because it seems to have high variability
			out <- mclapply(1:ncvreps, function(dummy) {
				fit <- cv.glmnet(data.sel, response, family="binomial", type.measure="class")
				errors = data.frame(fit$lambda,fit$cvm)
				errors
			}, mc.cores=16)
			lambdas <- do.call(rbind, out)
			
			# take mean cvm for each lambda
			df <- aggregate(fit.cvm ~ fit.lambda, lambdas, mean); colnames(df) <- c("lambda", "cvm")
			df$cvsd <- aggregate(fit.cvm ~ fit.lambda, lambdas, sd)[,2]
			df$loglambda <- log(df$lambda)
			# select the best one
			bestindex = which.min(df$cvm); bestlambda = df[bestindex, "lambda"]
			fit <- glmnet(data.sel, response, alpha=1, family="binomial")
			lasso_coef = predict(fit, type="coefficients", s=bestlambda, exact=T)
			# plot CV fit (log-lambda vs CV error)
			p <- ggplot(df, aes(x=loglambda, y=cvm)) + geom_point(color="red") + geom_errorbar(aes(x=loglambda, ymin=cvm-cvsd, ymax=cvm+cvsd)) + geom_vline(xintercept=log(bestlambda), linetype="dotted") + theme_classic() + ggtitle(sprintf("CV fit for LASSO (%s %s %s %s)", regi, subtype, st, mlevel))
			print(p)
			
#			res <- as.data.frame(as.matrix(coef(cvfit, s="lambda.min"))); colnames(res) <- c("value")
			res <- as.data.frame(as.matrix(lasso_coef)); colnames(res) <- c("value")
			res <- subset(res, value != 0); res$taxa <- rownames(res)
			res <- subset(res, taxa != "(Intercept)")
			if (nrow(res)>0) {
				res <- res[order(res$value, decreasing=T),]; res$taxa <- factor(res$taxa, levels=res$taxa)
				res.mean <- rbind(res.mean, res)
				p <- ggplot(res, aes(x=taxa, y=value)) + geom_bar(stat="identity") + theme_classic() + ggtitle(sprintf("[LASSO] %s ~ selected features (%s %s %s, lambda=%.4g)", regi, subtype, st, mlevel, bestlambda)) + theme(plot.title=element_text(size=8), axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5))
				print(p)
				
				# violin plots of metabolite values
				agg.melt <- agg.melt.stored
				agg.melt$MaternalGroup <- mapping.sel[agg.melt$SampleID, "MaternalGroup"]
				agg.melt$delgage <- mapping.sel[agg.melt$SampleID, "delgage"]
				agg.melt$weight0week <- mapping.sel[agg.melt$SampleID, "weight0week"]
				agg.melt$DaysPTDPlasma2 <- mapping.sel[agg.melt$SampleID, "DaysPTDPlasma2"]
				agg.melt <- subset(agg.melt, metabolite %in% rownames(res))
				agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(res))
				p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of LASSO metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
				p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of LASSO metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				npages <- n_pages(p)
				for (ip in 1:npages) {
					p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of LASSO metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
					print(p)
				}
				for (ip in 1:npages) {
					p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=delgage)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of LASSO metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
					print(p)
				}
				for (ip in 1:npages) {
					p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=weight0week)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of LASSO metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
					print(p)
				}
				for (ip in 1:npages) {
					p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=DaysPTDPlasma2)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of LASSO metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.daystodelivery)
					print(p)
				}
			}
			write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/LASSO.%s.%s.%s.%s.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
		}
	}
}


## elastic net regression of MaternalGroup (separately by Regimen); using METABOLITE data [DBS, Plasma]
num_cores <- 16
cluster <- makePSOCKcluster(num_cores)
registerDoParallel(cluster)

ncvreps <- 100
ctrl1 <- trainControl(method="repeatedcv", number=10, repeats=ncvreps, returnResamp="all", classProbs=TRUE)

set.seed(nrow(mapping))
subtype <- "maternal"; mvar <- "Delivery"
for (st in c("DBS", "Plasma")) {
	for (mlevel in "BIOCHEMICAL") {
		data <- df.metabolon[[st]][[mlevel]]
		mapping2 <- mappinglist[[st]][,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "weight0week", "DaysPTDPlasma2")]
		colnames(mapping2) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weight0week", "DaysPTDPlasma2")
		rownames(mapping2) <- mapping2$patid
		res.mean <- {}; res.sd <- {}
		for (regi in setdiff(levels(mapping2$MaternalRegimen), "other")) {
			mapping.sel <- subset(mapping2, MaternalRegimen==regi); mapping.sel$MaternalGroup <- droplevels(mapping.sel$MaternalGroup)
			data.sel <- data[rownames(mapping.sel),] # subset to just the desired maternal samples from regimen
			to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
			data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
			response <- droplevels(mapping.sel$MaternalGroup); names(response) <- rownames(mapping.sel)
			tmp <- factor(make.names(response), levels=make.names(levels(response))); names(tmp) <- names(response); response <- tmp
			# add Country as covariates
#			data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
			agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
			
			grid.glmnet <- expand.grid(alpha=seq(from=0, to=1, by=0.05), lambda=seq(from=0, to=0.3, by=0.01))
			caret_res <- train(data.sel, response, method = "glmnet", trControl = ctrl1, tuneGrid=grid.glmnet)
			
			# plot caret results
			p <- plot(caret_res)
			print(p)
			
			# build best model
			fit <- glmnet(data.sel, response, alpha=caret_res$bestTune$alpha, family="binomial") # same as caret_res$finalModel (?)
			enet_coef = predict(fit, x=data.sel, y=response, type="coefficients", s=caret_res$bestTune$lambda, exact=T)
			
#			res <- as.data.frame(as.matrix(coef(cvfit, s="lambda.min"))); colnames(res) <- c("value")
			res <- as.data.frame(as.matrix(enet_coef)); colnames(res) <- c("value")
			res <- subset(res, value != 0); res$taxa <- rownames(res)
			res <- subset(res, taxa != "(Intercept)")
			if (nrow(res)>0) {
				res <- res[order(res$value, decreasing=T),]; res$taxa <- factor(res$taxa, levels=res$taxa)
				res.mean <- rbind(res.mean, res)
				p <- ggplot(res, aes(x=taxa, y=value)) + geom_bar(stat="identity") + theme_classic() + ggtitle(sprintf("[Elastic net] %s ~ selected features (%s %s %s, alpha=%.4g lambda=%.4g)", regi, subtype, st, mlevel, caret_res$bestTune$alpha, caret_res$bestTune$lambda)) + theme(plot.title=element_text(size=8), axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5))
				print(p)
				
				# violin plots of metabolite values
				agg.melt <- agg.melt.stored
				agg.melt$MaternalGroup <- mapping.sel[agg.melt$SampleID, "MaternalGroup"]
				agg.melt$delgage <- mapping.sel[agg.melt$SampleID, "delgage"]
				agg.melt$weight0week <- mapping.sel[agg.melt$SampleID, "weight0week"]
				agg.melt$DaysPTDPlasma2 <- mapping.sel[agg.melt$SampleID, "DaysPTDPlasma2"]
				agg.melt <- subset(agg.melt, metabolite %in% rownames(res))
				agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(res))
				cols.cohort.sel <- cols.cohort[intersect(names(cols.cohort), levels(agg.melt$MaternalGroup))]
				p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of Elastic net metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort.sel)
				print(p)
				p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of Elastic net metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort.sel)
				npages <- n_pages(p)
				for (ip in 1:npages) {
					p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of Elastic net metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort.sel)
					print(p)
				}
				for (ip in 1:npages) {
					p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=delgage)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of Elastic net metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
					print(p)
				}
				for (ip in 1:npages) {
					p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=weight0week)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of Elastic net metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
					print(p)
				}
				for (ip in 1:npages) {
					p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=DaysPTDPlasma2)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of Elastic net metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.daystodelivery)
					print(p)
				}
			}
			write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/elastic_net.%s.%s.%s.%s.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
		}
	}
}


## randomForest classification of MaternalGroup (separately by Regimen); using METABOLITE data [DBS, Plasma]
set.seed(nrow(mapping))
# MAIN LOOP for random forest (through metabolite levels)
subtype <- "maternal"
results.performance <- {}
for (st in c("DBS", "Plasma")) {
	for (mlevel in "BIOCHEMICAL") {
		data <- df.metabolon[[st]][[mlevel]]
		mapping2 <- mappinglist[[st]][,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "weight0week", "DaysPTDPlasma2")]
		colnames(mapping2) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weight0week", "DaysPTDPlasma2")
		rownames(mapping2) <- mapping2$patid
		for (regi in setdiff(levels(mapping2$MaternalRegimen), "other")) {
			mapping.sel <- subset(mapping2, MaternalRegimen==regi); mapping.sel$MaternalGroup <- droplevels(mapping.sel$MaternalGroup)
			data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the desired maternal samples from regimen
			to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
			data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
			response <- droplevels(mapping.sel$MaternalGroup); names(response) <- rownames(mapping.sel)
			# add Country as covariates
			data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
			agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
			
#			## after running for the first time, COMMENT OUT THIS BLOCK ##
#			num_iter <- 100
#			ncores <- 20
#			out <- mclapply(1:num_iter, function (dummy) {
#					importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
#				}, mc.cores=ncores )
#			collated.importance <- do.call(cbind, out)
#			out <- mclapply(1:num_iter, function (dummy) {
#					rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.9)$error.cv
#				}, mc.cores=ncores )
#			collated.cv <- do.call(cbind, out)

#			write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
#			write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
#			## END BLOCK TO COMMENT ##

			collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1, quote="")
			collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
			importance.mean <- rowMeans(collated.importance)
			importance.sd <- unlist(apply(collated.importance, 1, sd))
			cv.mean <- rowMeans(collated.cv)
			cv.sd <- unlist(apply(collated.cv, 1, sd))
			inds <- order(importance.mean, decreasing=T)
#			inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate
			inds <- inds[1:min(20, as.numeric(names(cv.mean)[which.min(cv.mean)]))] # select number of features based on cross-validation
			write.table(melt(importance.mean[inds]), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

			## after running for the first time, COMMENT OUT THIS BLOCK ##
			# using a sparse model with N predictors
#			sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
#			save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", regi, subtype, mlevel, st))
			load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", regi, subtype, mlevel, st))
			# accuracy of final sparseRF model
			pred <- predict(sparseRF, type="prob")
			pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
			pred_df_out <- merge(pred_df, data.sel, by="row.names")
			write.table(pred_df_out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.predictions.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
			confusion_matrix <- table(pred_df[, c("true", "predicted")])
			class_errors <- unlist(lapply(levels(mapping.sel$MaternalGroup), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$MaternalGroup)
			accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
			vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
			mccvalue <- mcc(vec.pred, vec.true)
			cf <- confusionMatrix(confusion_matrix, positive=sprintf("Preterm.%s", regi))
			cflabel <- sprintf("Positive: %s Sens: %.4g  Spec: %.4g\n PPV: %.4g  NPV: %.4g", cf$positive, cf$byClass[["Sensitivity"]], cf$byClass[["Specificity"]], cf$byClass[["Pos Pred Value"]], cf$byClass[["Neg Pred Value"]])
			df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
			p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", regi, subtype, mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + annotate(geom="text", x=2.5, y=10, label=cflabel)
			print(p)
			pred2 <- prediction(pred[,2], ordered(response))
			perf <- performance(pred2, "tpr", "fpr")
			perf.auc <- performance(pred2, "auc")
			plot(perf, main=sprintf("ROC %s %s %s %s (sparseRF final model)", regi, subtype, mlevel, st)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", unlist(perf.auc@y.values)))
			to_store <- data.frame(fpr=perf@x.values[[1]], tpr=perf@y.values[[1]], alpha=perf@alpha.values[[1]], SampleType=st, Regimen=regi, Subtype=subtype, level=mlevel)
			results.performance <- rbind(results.performance, to_store)
			
			write.table(confusion_matrix, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.confusion_matrix.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=T)
			## END BLOCK TO COMMENT ##
			
			# plotting - per-group sparse model
			df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
			colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
			print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s %s", regi, subtype, mlevel, st)))
			# plotting - per-group variables
			df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
			df$metabolite_name <- as.character(df$OTU)
			if (mlevel == "BIOCHEMICAL") {
				df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
				df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
				df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
			} else if (mlevel == "SUB.PATHWAY") {
				df$subpathway <- df$metabolite_name
				df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
				df$OTU_string <- sprintf("%s", df$OTU)
			} else {
				df$superpathway <- df$metabolite_name
				df$OTU_string <- sprintf("%s", df$OTU)
			}
			df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
			# load effect sizes from linear regression
			lmres <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote=""); lmres <- subset(lmres, MaternalRegimen==regi); rownames(lmres) <- lmres$metabolite
			for (lmvar in c("estimate", "SE", "padj", "dir")) {
				df[,lmvar] <- lmres[df$metabolite_name, lmvar]
			}
			p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("%s - %s %s explanatory %s", regi, subtype, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
			print(p)
			lims <- max(abs(df$estimate) + abs(df$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
			p <- ggplot(df, aes(x=OTU, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=OTU, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_tile(aes(x=OTU, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s - %s %s explanatory %s", regi, subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
			print(p)
			# stratified by country
			lmres.bycountry <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_RegimenCountry.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote="")
			tmp <- levels(df$OTU); lmres.bycountry <- subset(lmres.bycountry, metabolite %in% tmp & MaternalRegimen==regi)
			lmres.bycountry$metabolite <- factor(lmres.bycountry$metabolite, levels=levels(df$OTU))
			lims <- max(abs(lmres.bycountry$estimate) + abs(lmres.bycountry$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
			dircolors2 <- c("blue", "red", "black"); names(dircolors2) <- c("down", "up", "NS")
			lmres.bycountry[, "importance"] <- df[match(lmres.bycountry$metabolite, df$OTU), "importance"]
			p <- ggplot(lmres.bycountry, aes(x=metabolite, y=estimate, color=dir, group=Country)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Country), position=pd, hjust=1, color="black", size=1.5) + geom_tile(mapping=aes(x=metabolite, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates - %s %s explanatory %s", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
			print(p)
			p <- ggplot(lmres.bycountry, aes(x=metabolite, y=estimate, color=dir, group=Country)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Country), position=pd, hjust=1, color="black", size=1.5) + geom_tile(mapping=aes(x=metabolite, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + geom_vline(xintercept=as.numeric(df$OTU)+0.5, color="grey") + theme_classic() + ggtitle(sprintf("LM estimates - %s %s explanatory %s", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors2) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
			print(p)
			# shading rectangles of importance values
			df.rect <- df
			df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
	#		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
			p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s %s explanatory %s", regi, subtype, mlevel, st)) + scale_fill_gradient(low="white", high="black")
			print(p)
			# violin plots of metabolite values
			agg.melt <- agg.melt.stored
			agg.melt$MaternalGroup <- mapping.sel[agg.melt$SampleID, "MaternalGroup"]
			agg.melt$delgage <- mapping.sel[agg.melt$SampleID, "delgage"]
			agg.melt$weight0week <- mapping.sel[agg.melt$SampleID, "weight0week"]
			agg.melt$DaysPTDPlasma2 <- mapping.sel[agg.melt$SampleID, "DaysPTDPlasma2"]
			agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$MaternalGroup))
			agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
			agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			npages <- n_pages(p)
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
			}
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=delgage)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
				print(p)
			}
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=weight0week)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
				print(p)
			}
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=DaysPTDPlasma2)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.daystodelivery)
				print(p)
			}
			
			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=MaternalGroup)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
			agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
			p <- ggplot(agg.melt2, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=MaternalGroup)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
#			# violin plots stratified by prediction/truth
#			for (met in levels(agg.melt$metabolite)) {
#				tmp <- subset(agg.melt, metabolite==met)
#				# color scheme - manual
#				p <- ggplot(tmp, aes(x=MaternalGroup, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_fill_manual(values=cols.cohort)
#				print(p)
#				p <- ggplot(tmp, aes(x=MaternalGroup, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_color_manual(values=cols.cohort)
#				print(p)
#			}
		}
	}
}
## save ROC data and plot combined ROC
write.table(results.performance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.performance_data.txt"), quote=F, sep="\t", row.names=F, col.names=T)
for (st in c("DBS", "Plasma")) {
	p <- ggplot(subset(results.performance, SampleType==st), aes(x=fpr, y=tpr, color=Regimen)) + geom_point() + geom_line() + theme_classic() + ggtitle(sprintf("Combined ROC plot (maternal %s)", st)) + scale_color_manual(values=cols.regimen)
	print(p)
}

## combined heatmap of stratified Group classification models + ROC plots; combined across DBS and plasma
featurelist <- {}
subtype <- "maternal"
for (st in c("DBS", "Plasma")) {
	for (regi in c("untreated", "zdv", "PI-ART")) {
		tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", quote=""); colnames(tmp) <- c("feature", "importance")
		featurelist <- c(featurelist, tmp$feature)
	}
}
featurelist <- unique(featurelist)
df <- matrix(0, ncol=6, nrow=length(featurelist)); rownames(df) <- featurelist; colnames(df) <- c("Plasma - untreated", "DBS - untreated", "Plasma - zdv", "DBS - zdv", "Plasma - PI-ART", "DBS - PI-ART")
for (st in c("DBS", "Plasma")) {
	for (regi in c("untreated", "zdv", "PI-ART")) {
		tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", quote=""); colnames(tmp) <- c("feature", "importance")
		df[tmp$feature, sprintf("%s - %s", st, regi)] <- tmp$importance
	}
}
heatmap.2(df, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.3, main=sprintf("RF importance values (%s, %s)", subtype, mlevel))
heatmap.2(df, Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.3, main=sprintf("RF importance values (%s, %s)", subtype, mlevel))
sel <- names(which(rowSums(df>0)>1)) # metabolites found in >1 model
df <- df>0; df <- df[sel,] # convert to binary flag
#write.table(df, file="/Lab_Share/PROMISE/nwcs619/metabolon/RF_features_in_multiple_models.txt", row.names=T, col.names=T, sep="\t", quote=F)
heatmap.2(df+1-1, Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.6, main=sprintf("Metabolites found in >1 model (%s, %s)", subtype, mlevel))
#sel <- c(names(which(df[,"DBS - untreated"] & df[,"Plasma - untreated"])), names(which(df[,"DBS - zdv"] & df[,"Plasma - zdv"])), names(which(df[,"DBS - PI-ART"] & df[,"Plasma - PI-ART"])))
#heatmap.2(df[sel,]+1-1, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.6, main=sprintf("Metabolites found consistently (%s, %s)", subtype, mlevel))
#heatmap.2(df[setdiff(rownames(df),sel),,drop=F]+1-1, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.6, main=sprintf("Metabolites found otherwise > 1 model (%s, %s)", subtype, mlevel))


## violin plots of all metabolite values
subtype <- "maternal"
for (st in c("DBS", "Plasma")) {
	for (mlevel in "BIOCHEMICAL") {
		data <- df.metabolon[[st]][[mlevel]]
		mapping.sel <- mappinglist[[st]][,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "weight0week")]
		colnames(mapping.sel) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weight0week")
		mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("untreated", "zdv", "PI-ART")) # exclude other regimen because insufficient numbers
		rownames(mapping.sel) <- mapping.sel$patid
		data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the maternal samples
		agg.melt <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt) <- c("SampleID", "metabolite", "value")
		agg.melt$MaternalGroup <- droplevels(mapping.sel[agg.melt$SampleID, "MaternalGroup"])
		agg.melt$MaternalGroup <- factor(as.character(agg.melt$MaternalGroup), levels=rev(c("Preterm.untreated", "Term.untreated", "Preterm.zdv", "Term.zdv", "Preterm.PI-ART", "Term.PI-ART")))
		# manual pagination
		metabolites <- unique(agg.melt$metabolite)
		for (i in seq(from=1,to=length(metabolites), by=12)) {
			j <- min(i+11, length(metabolites))
			agg.melt.sel <- subset(agg.melt, metabolite %in% metabolites[i:j])
			p <- ggplot(agg.melt.sel, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3, nrow=4) + theme_classic() + ggtitle(sprintf("Rel. abund. of all metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort) + theme(strip.text = element_text(size=7))
			print(p)
		}
	}
}

### randomForest classification of Group (multiclass); using METABOLITE data [DBS, Plasma]
#set.seed(nrow(mapping))	
## MAIN LOOP for random forest (through metabolite levels)
#subtype <- "maternal"
#for (st in c("DBS", "Plasma")) {
#	for (mlevel in "BIOCHEMICAL") {
#		data <- df.metabolon[[st]][[mlevel]]
#		mapping.sel <- mappinglist[[st]][,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "weight0week", "DaysPTDPlasma2")]
#		colnames(mapping.sel) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weight0week", "DaysPTDPlasma2")
#		mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("untreated", "zdv", "PI-ART")) # exclude other regimen because insufficient numbers
#		rownames(mapping.sel) <- mapping.sel$patid
#		data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the maternal samples
#		to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
#		data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
#		response <- droplevels(mapping.sel$MaternalGroup); names(response) <- rownames(mapping.sel); group_levels <- levels(response)
#		# add Country as covariates
#		data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
#		agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
#		
#		## after running for the first time, COMMENT OUT THIS BLOCK ##
#		num_iter <- 100
#		ncores <- 20
#		out <- mclapply(1:num_iter, function (dummy) {
#				importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
#			}, mc.cores=ncores )
#		collated.importance <- do.call(cbind, out)
#		out <- mclapply(1:num_iter, function (dummy) {
#				rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.9)$error.cv
#			}, mc.cores=ncores )
#		collated.cv <- do.call(cbind, out)

#		write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
#		write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
#		## END BLOCK TO COMMENT ##

#		collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", "multiclass", subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1, quote="")
#		collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", "multiclass", subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
#		importance.mean <- rowMeans(collated.importance)
#		importance.sd <- unlist(apply(collated.importance, 1, sd))
#		cv.mean <- rowMeans(collated.cv)
#		cv.sd <- unlist(apply(collated.cv, 1, sd))
#		inds <- order(importance.mean, decreasing=T)
##		inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate
#		inds <- inds[1:min(50, as.numeric(names(cv.mean)[which.min(cv.mean)]))] # select number of features based on cross-validation
#		write.table(melt(importance.mean[inds]), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

#		## after running for the first time, COMMENT OUT THIS BLOCK ##
#		# using a sparse model with N predictors
#		sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
#		save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", "multiclass", subtype, mlevel, st))
#		load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", "multiclass", subtype, mlevel, st))
#		# accuracy of final sparseRF model
#		pred <- predict(sparseRF, type="prob")
#		pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
#		pred_df_out <- merge(pred_df, data.sel, by="row.names")
#		write.table(pred_df_out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.predictions.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
#		confusion_matrix <- table(pred_df[, c("true", "predicted")])
#		class_errors <- unlist(lapply(group_levels, function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- group_levels
#		accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
#		vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
#		mccvalue <- mcc(vec.pred, vec.true)
#		df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
#		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", "multiclass", subtype, mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#		print(p)
#		
#		write.table(confusion_matrix, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.confusion_matrix.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=T)
#		## END BLOCK TO COMMENT ##
#		
#		# plotting - per-group sparse model
#		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
#		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
#		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s %s", "multiclass", subtype, mlevel, st)))
#		# plotting - per-group variables
#		df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
#		df$metabolite_name <- as.character(df$OTU)
#		if (mlevel == "BIOCHEMICAL") {
#			df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
#			df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
#			df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
#		} else if (mlevel == "SUB.PATHWAY") {
#			df$subpathway <- df$metabolite_name
#			df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
#			df$OTU_string <- sprintf("%s", df$OTU)
#		} else {
#			df$superpathway <- df$metabolite_name
#			df$OTU_string <- sprintf("%s", df$OTU)
#		}
#		df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
#		# load effect sizes from linear regression
#		lmres <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote="")
#		tmp <- levels(df$OTU); lmres <- subset(lmres, metabolite %in% tmp)
#		lmres$metabolite <- factor(lmres$metabolite, levels=levels(df$OTU))
#		lims <- max(abs(lmres$estimate) + abs(lmres$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
#		dircolors2 <- c("blue", "red", "black"); names(dircolors2) <- c("down", "up", "NS")
#		lmres[, "importance"] <- df[match(lmres$metabolite, df$OTU), "importance"]
#		p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("%s - %s %s explanatory %s", "multiclass", subtype, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
#		print(p)
#		p <- ggplot(lmres, aes(x=metabolite, y=estimate, color=dir, group=MaternalRegimen)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=MaternalRegimen), position=pd, hjust=1, color="black", size=2) + geom_tile(mapping=aes(x=metabolite, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates - %s %s explanatory %s", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
#		print(p)
#		p <- ggplot(lmres, aes(x=metabolite, y=estimate, color=dir, group=MaternalRegimen)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=MaternalRegimen), position=pd, hjust=1, color="black", size=2) + geom_tile(mapping=aes(x=metabolite, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + geom_vline(xintercept=as.numeric(df$OTU)+0.5, color="grey") + theme_classic() + ggtitle(sprintf("LM estimates - %s %s explanatory %s", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors2) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
#		print(p)
#		# stratified by country
#		lmres.bycountry <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_RegimenCountry.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote="")
#		tmp <- levels(df$OTU); lmres.bycountry <- subset(lmres.bycountry, metabolite %in% tmp)
#		lmres.bycountry$metabolite <- factor(lmres.bycountry$metabolite, levels=levels(df$OTU))
#		lims <- max(abs(lmres.bycountry$estimate) + abs(lmres.bycountry$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
#		dircolors2 <- c("blue", "red", "black"); names(dircolors2) <- c("down", "up", "NS")
#		lmres.bycountry[, "importance"] <- df[match(lmres.bycountry$metabolite, df$OTU), "importance"]
#		for (i in seq(from=1, to=nlevels(lmres.bycountry$metabolite), by=10)) {
#			sel <- rev(levels(lmres.bycountry$metabolite))[i:min(i+9, nlevels(lmres.bycountry$metabolite))]
#			lmres.bycountry.sel <- subset(lmres.bycountry, metabolite %in% sel)
#			p <- ggplot(lmres.bycountry.sel, aes(x=metabolite, y=estimate, color=dir, group=Country)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Country), position=pd, hjust=1, color="black", size=1.5) + geom_tile(mapping=aes(x=metabolite, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + geom_vline(xintercept=as.numeric(df$OTU)+0.5, color="grey") + facet_wrap(~MaternalRegimen) + theme_classic() + ggtitle(sprintf("LM estimates - %s %s explanatory %s", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors2) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
#			print(p)
#		}
#		# shading rectangles of importance values
#		df.rect <- df
#		df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
##		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
#		p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s %s explanatory %s", "multiclass", subtype, mlevel, st)) + scale_fill_gradient(low="white", high="black")
#		print(p)
#		# violin plots of metabolite values
#		agg.melt <- agg.melt.stored
#		agg.melt$MaternalGroup <- mapping.sel[agg.melt$SampleID, "MaternalGroup"]
#		agg.melt$delgage <- mapping.sel[agg.melt$SampleID, "delgage"]
#		agg.melt$weight0week <- mapping.sel[agg.melt$SampleID, "weight0week"]
#		agg.melt$DaysPTDPlasma2 <- mapping.sel[agg.melt$SampleID, "DaysPTDPlasma2"]
#		agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$MaternalGroup))
#		agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
#		agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
#		p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#		print(p)
#		p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#		npages <- n_pages(p)
#		for (ip in 1:npages) {
#			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
#		}
#		for (ip in 1:npages) {
#			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=delgage)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
#			print(p)
#		}
#		for (ip in 1:npages) {
#			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=weight0week)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
#			print(p)
#		}
#		for (ip in 1:npages) {
#			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=DaysPTDPlasma2)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.daystodelivery)
#			print(p)
#		}
#		
#		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
#		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=MaternalGroup)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#		print(p)
#		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
#		agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
#		p <- ggplot(agg.melt2, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#		print(p)
#		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=MaternalGroup)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#		print(p)
##		# violin plots stratified by prediction/truth
##		for (met in levels(agg.melt$metabolite)) {
##			tmp <- subset(agg.melt, metabolite==met)
##			# color scheme - manual
##			p <- ggplot(tmp, aes(x=MaternalGroup, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_fill_manual(values=cols.cohort)
##			print(p)
##			p <- ggplot(tmp, aes(x=MaternalGroup, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_color_manual(values=cols.cohort)
##			print(p)
##		}
#	}
#}



#########################################################################################################
### maternal multiomics (Plasma and DBS combined)
## randomForest classification of MaternalGroup (separately by Regimen); using METABOLITE data [DBS+Plasma]
set.seed(nrow(mapping))
# MAIN LOOP for random forest (through metabolite levels)
subtype <- "maternal"
results.performance <- {}
for (mlevel in "BIOCHEMICAL") {
	mapping2 <- mapping[,c("Delivery", "Delivery2", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "weight0week", "DaysPTDPlasma2")]
	colnames(mapping2) <- c("Delivery", "Delivery2", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weight0week", "DaysPTDPlasma2")
	rownames(mapping2) <- mapping2$patid
	# combine Plasma and DBS data
	data <- {}
	for (st in c("Plasma", "DBS")) {
		tmp <- df.metabolon[[st]][[mlevel]]
		colnames(tmp) <- sprintf("[%s] %s", st, colnames(tmp))
		tmp <- tmp[rownames(mapping2),]
		data <- cbind(data, tmp)
	}
	for (regi in setdiff(levels(mapping2$MaternalRegimen), "other")) {
		mapping.sel <- subset(mapping2, MaternalRegimen==regi); mapping.sel$MaternalGroup <- droplevels(mapping.sel$MaternalGroup)
		data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the desired maternal samples from regimen
		to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
		data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
		response <- droplevels(mapping.sel$MaternalGroup); names(response) <- rownames(mapping.sel)
		# add Country as covariates
		data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
		agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
		
#		## after running for the first time, COMMENT OUT THIS BLOCK ##
#		num_iter <- 100
#		ncores <- 20
#		out <- mclapply(1:num_iter, function (dummy) {
#				importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
#			}, mc.cores=ncores )
#		collated.importance <- do.call(cbind, out)
#		out <- mclapply(1:num_iter, function (dummy) {
#				rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.9)$error.cv
#			}, mc.cores=ncores )
#		collated.cv <- do.call(cbind, out)

#		write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", regi, subtype, mlevel, "multiomics"), quote=F, sep="\t", row.names=T, col.names=F)
#		write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", regi, subtype, mlevel, "multiomics"), quote=F, sep="\t", row.names=T, col.names=F)
#		## END BLOCK TO COMMENT ##

		collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", regi, subtype, mlevel, "multiomics"), header=F, as.is=T, sep="\t", row.names=1, quote="")
		collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", regi, subtype, mlevel, "multiomics"), header=F, as.is=T, sep="\t", row.names=1)
		importance.mean <- rowMeans(collated.importance)
		importance.sd <- unlist(apply(collated.importance, 1, sd))
		cv.mean <- rowMeans(collated.cv)
		cv.sd <- unlist(apply(collated.cv, 1, sd))
		inds <- order(importance.mean, decreasing=T)
#		inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate
		inds <- inds[1:min(20, as.numeric(names(cv.mean)[which.min(cv.mean)]))] # select number of features based on cross-validation
		write.table(melt(importance.mean[inds]), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, "multiomics"), quote=F, sep="\t", row.names=T, col.names=F)

		## after running for the first time, COMMENT OUT THIS BLOCK ##
		# using a sparse model with N predictors
#		sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
#		save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", regi, subtype, mlevel, "multiomics"))
		load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", regi, subtype, mlevel, "multiomics"))
		# accuracy of final sparseRF model
		pred <- predict(sparseRF, type="prob")
		pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
		pred_df_out <- merge(pred_df, data.sel, by="row.names")
		write.table(pred_df_out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.predictions.txt", regi, subtype, mlevel, "multiomics"), quote=F, sep="\t", row.names=F, col.names=T)
		confusion_matrix <- table(pred_df[, c("true", "predicted")])
		class_errors <- unlist(lapply(levels(mapping.sel$MaternalGroup), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$MaternalGroup)
		accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
		vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
		mccvalue <- mcc(vec.pred, vec.true)
		cf <- confusionMatrix(confusion_matrix, positive=sprintf("Preterm.%s", regi))
		cflabel <- sprintf("Positive: %s Sens: %.4g  Spec: %.4g\n PPV: %.4g  NPV: %.4g", cf$positive, cf$byClass[["Sensitivity"]], cf$byClass[["Specificity"]], cf$byClass[["Pos Pred Value"]], cf$byClass[["Neg Pred Value"]])
		df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", regi, subtype, mlevel, "multiomics", accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + annotate(geom="text", x=2.5, y=10, label=cflabel)
		print(p)
		
		pred2 <- prediction(pred[,2], ordered(response))
		perf <- performance(pred2, "tpr", "fpr")
		perf.auc <- performance(pred2, "auc")
		plot(perf, main=sprintf("ROC %s %s %s %s (sparseRF final model)", regi, subtype, mlevel, "multiomics")) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", unlist(perf.auc@y.values)))
		to_store <- data.frame(fpr=perf@x.values[[1]], tpr=perf@y.values[[1]], alpha=perf@alpha.values[[1]], SampleType="multiomics", Regimen=regi, Subtype=subtype, level=mlevel)
		results.performance <- rbind(results.performance, to_store)
		
		write.table(confusion_matrix, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.confusion_matrix.txt", regi, subtype, mlevel, "multiomics"), quote=F, sep="\t", row.names=T, col.names=T)
		## END BLOCK TO COMMENT ##
		
		# plotting - per-group sparse model
		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s %s", regi, subtype, mlevel, "multiomics")))
		# plotting - per-group variables
		df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
		df$metabolite_name <- unlist(lapply(as.character(df$OTU), function(metstr) {
			x <- unlist(strsplit(metstr, " "))
			paste(x[2:length(x)], collapse=" ")
		}))
		df$SampleType <- unlist(lapply(as.character(df$OTU), function(metstr) {
			gsub("\\[", "", gsub("\\]", "", unlist(strsplit(metstr, " "))[1]))
		}))
		if (mlevel == "BIOCHEMICAL") {
			df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
			df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
			df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
		} else if (mlevel == "SUB.PATHWAY") {
			df$subpathway <- df$metabolite_name
			df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
			df$OTU_string <- sprintf("%s", df$OTU)
		} else {
			df$superpathway <- df$metabolite_name
			df$OTU_string <- sprintf("%s", df$OTU)
		}
		df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
		# load effect sizes from linear regression
		lmres <- {}
		for (st in c("Plasma", "DBS")) {
			tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote=""); tmp <- subset(tmp, MaternalRegimen==regi); rownames(tmp) <- tmp$metabolite
			tmp$SampleType <- st
			lmres <- rbind(lmres, tmp)
		}
		rownames(lmres) <- sprintf("[%s] %s", lmres$SampleType, lmres$metabolite)
		for (lmvar in c("estimate", "SE", "padj", "dir")) {
			df[,lmvar] <- lmres[rownames(df), lmvar]
		}
		p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("%s - %s %s explanatory %s", regi, subtype, mlevel, "multiomics")) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
		print(p)
		lims <- max(abs(df$estimate) + abs(df$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
		p <- ggplot(df, aes(x=OTU, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=OTU, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_tile(aes(x=OTU, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s - %s %s explanatory %s", regi, subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
		print(p)
		# shading rectangles of importance values
		df.rect <- df
		df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
#		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
		p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s %s explanatory %s", regi, subtype, mlevel, "multiomics")) + scale_fill_gradient(low="white", high="black")
		print(p)
		# violin plots of metabolite values
		agg.melt <- agg.melt.stored
		agg.melt$MaternalGroup <- mapping.sel[agg.melt$SampleID, "MaternalGroup"]
		agg.melt$delgage <- mapping.sel[agg.melt$SampleID, "delgage"]
		agg.melt$weight0week <- mapping.sel[agg.melt$SampleID, "weight0week"]
		agg.melt$DaysPTDPlasma2 <- mapping.sel[agg.melt$SampleID, "DaysPTDPlasma2"]
		agg.melt$Delivery2 <- mapping.sel[agg.melt$SampleID, "Delivery2"]
		agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$MaternalGroup))
		agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
		agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
		p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.cohort)
		npages <- n_pages(p)
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
		}
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=delgage)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_gradient(low="red", high="black")
			print(p)
		}
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=weight0week)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_gradient(low="red", high="black")
			print(p)
		}
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=DaysPTDPlasma2)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.daystodelivery)
			print(p)
		}
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=Delivery2)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.delivery2)
			print(p)
		}
		
		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=MaternalGroup)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
		agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
		p <- ggplot(agg.melt2, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=MaternalGroup)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
	}
}



## elastic net regression of MaternalGroup (separately by Regimen); using METABOLITE data [DBS, Plasma]
ncvreps <- 100
ctrl1 <- trainControl(method="repeatedcv", number=10, repeats=ncvreps, returnResamp="all", classProbs=TRUE)
set.seed(nrow(mapping))
subtype <- "maternal"; mvar <- "Delivery"; st <- "multiomics"
for (mlevel in "BIOCHEMICAL") {
	mapping2 <- mapping[,c("Delivery", "Delivery2", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "weight0week", "DaysPTDPlasma2")]
	colnames(mapping2) <- c("Delivery", "Delivery2", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weight0week", "DaysPTDPlasma2")
	rownames(mapping2) <- mapping2$patid
	# combine Plasma and DBS data
	data <- {}
	for (s in c("Plasma", "DBS")) {
		tmp <- df.metabolon[[s]][[mlevel]]
		colnames(tmp) <- sprintf("[%s] %s", s, colnames(tmp))
		tmp <- tmp[rownames(mapping2),]
		data <- cbind(data, tmp)
	}
	for (regi in setdiff(levels(mapping2$MaternalRegimen), "other")) {
		mapping.sel <- subset(mapping2, MaternalRegimen==regi); mapping.sel$MaternalGroup <- droplevels(mapping.sel$MaternalGroup)
		data.sel <- data[rownames(mapping.sel),] # subset to just the desired maternal samples from regimen
		to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
		data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
		response <- droplevels(mapping.sel$MaternalGroup); names(response) <- rownames(mapping.sel)
		tmp <- factor(make.names(response), levels=make.names(levels(response))); names(tmp) <- names(response); response <- tmp
		# add Country as covariates
#			data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
		agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
		
		grid.glmnet <- expand.grid(alpha=seq(from=0, to=1, by=0.05), lambda=seq(from=0, to=0.3, by=0.01))
		caret_res <- train(data.sel, response, method = "glmnet", trControl = ctrl1, tuneGrid=grid.glmnet)
		
		# plot caret results
		p <- plot(caret_res)
		print(p)
		
		# build best model
		fit <- glmnet(data.sel, response, alpha=caret_res$bestTune$alpha, family="binomial") # same as caret_res$finalModel (?)
		enet_coef = predict(fit, x=data.sel, y=response, type="coefficients", s=caret_res$bestTune$lambda, exact=T)
		
#			res <- as.data.frame(as.matrix(coef(cvfit, s="lambda.min"))); colnames(res) <- c("value")
		res <- as.data.frame(as.matrix(enet_coef)); colnames(res) <- c("value")
		res <- subset(res, value != 0); res$taxa <- rownames(res)
		res <- subset(res, taxa != "(Intercept)")
		if (nrow(res)>0) {
			res <- res[order(res$value, decreasing=T),]; res$taxa <- factor(res$taxa, levels=res$taxa)
			res.mean <- rbind(res.mean, res)
			p <- ggplot(res, aes(x=taxa, y=value)) + geom_bar(stat="identity") + theme_classic() + ggtitle(sprintf("[Elastic net] %s ~ selected features (%s %s %s, alpha=%.4g lambda=%.4g)", regi, subtype, st, mlevel, caret_res$bestTune$alpha, caret_res$bestTune$lambda)) + theme(plot.title=element_text(size=8), axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5))
			print(p)
			
			# violin plots of metabolite values
			agg.melt <- agg.melt.stored
			agg.melt$MaternalGroup <- mapping.sel[agg.melt$SampleID, "MaternalGroup"]
			agg.melt$delgage <- mapping.sel[agg.melt$SampleID, "delgage"]
			agg.melt$weight0week <- mapping.sel[agg.melt$SampleID, "weight0week"]
			agg.melt$DaysPTDPlasma2 <- mapping.sel[agg.melt$SampleID, "DaysPTDPlasma2"]
			agg.melt <- subset(agg.melt, metabolite %in% rownames(res))
			agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(res))
			cols.cohort.sel <- cols.cohort[intersect(names(cols.cohort), levels(agg.melt$MaternalGroup))]
			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of Elastic net metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort.sel)
			print(p)
			p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of Elastic net metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort.sel)
			npages <- n_pages(p)
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of Elastic net metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort.sel)
				print(p)
			}
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=delgage)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of Elastic net metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
				print(p)
			}
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=weight0week)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of Elastic net metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
				print(p)
			}
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=DaysPTDPlasma2)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of Elastic net metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.daystodelivery)
				print(p)
			}
		}
		write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/elastic_net.%s.%s.%s.%s.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
	}
}




### randomForest classification of Group (multiclass); using METABOLITE data [DBS+Plasma]
#set.seed(nrow(mapping))	
## MAIN LOOP for random forest (through metabolite levels)
#subtype <- "maternal"
#for (mlevel in "BIOCHEMICAL") {
#	mapping.sel <- mapping[,c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "weight0week", "DaysPTDPlasma2")]
#	colnames(mapping.sel) <- c("Delivery", "MaternalRegimen", "MaternalGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weight0week", "DaysPTDPlasma2")
#	mapping.sel <- subset(mapping.sel, MaternalRegimen %in% c("untreated", "zdv", "PI-ART")) # exclude other regimen because insufficient numbers
#	rownames(mapping.sel) <- mapping.sel$patid
#	# combine Plasma and DBS data
#	data <- {}
#	for (st in c("Plasma", "DBS")) {
#		tmp <- df.metabolon[[st]][[mlevel]]
#		colnames(tmp) <- sprintf("[%s] %s", st, colnames(tmp))
#		tmp <- tmp[rownames(mapping2),]
#		data <- cbind(data, tmp)
#	}
#	data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the maternal samples
#	to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
#	data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
#	response <- droplevels(mapping.sel$MaternalGroup); names(response) <- rownames(mapping.sel); group_levels <- levels(response)
#	# add Country as covariates
#	data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
#	agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
#	
##	## after running for the first time, COMMENT OUT THIS BLOCK ##
##	num_iter <- 100
##	ncores <- 20
##	out <- mclapply(1:num_iter, function (dummy) {
##			importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
##		}, mc.cores=ncores )
##	collated.importance <- do.call(cbind, out)
##	out <- mclapply(1:num_iter, function (dummy) {
##			rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.9)$error.cv
##		}, mc.cores=ncores )
##	collated.cv <- do.call(cbind, out)

##	write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", "multiclass", subtype, mlevel, "multiomics"), quote=F, sep="\t", row.names=T, col.names=F)
##	write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", "multiclass", subtype, mlevel, "multiomics"), quote=F, sep="\t", row.names=T, col.names=F)
##	## END BLOCK TO COMMENT ##

#	collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", "multiclass", subtype, mlevel, "multiomics"), header=F, as.is=T, sep="\t", row.names=1, quote="")
#	collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", "multiclass", subtype, mlevel, "multiomics"), header=F, as.is=T, sep="\t", row.names=1)
#	importance.mean <- rowMeans(collated.importance)
#	importance.sd <- unlist(apply(collated.importance, 1, sd))
#	cv.mean <- rowMeans(collated.cv)
#	cv.sd <- unlist(apply(collated.cv, 1, sd))
#	inds <- order(importance.mean, decreasing=T)
##	inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate
#	inds <- inds[1:min(50, as.numeric(names(cv.mean)[which.min(cv.mean)]))] # select number of features based on cross-validation
#	write.table(melt(importance.mean[inds]), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", "multiclass", subtype, mlevel, "multiomics"), quote=F, sep="\t", row.names=T, col.names=F)

#	## after running for the first time, COMMENT OUT THIS BLOCK ##
#	# using a sparse model with N predictors
##	sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
##	save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", "multiclass", subtype, mlevel, "multiomics"))
#	load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", "multiclass", subtype, mlevel, "multiomics"))
#	# accuracy of final sparseRF model
#	pred <- predict(sparseRF, type="prob")
#	pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
#	pred_df_out <- merge(pred_df, data.sel, by="row.names")
#	write.table(pred_df_out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.predictions.txt", "multiclass", subtype, mlevel, "multiomics"), quote=F, sep="\t", row.names=F, col.names=T)
#	confusion_matrix <- table(pred_df[, c("true", "predicted")])
#	class_errors <- unlist(lapply(group_levels, function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- group_levels
#	accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
#	vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
#	mccvalue <- mcc(vec.pred, vec.true)
#	df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
#	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", "multiclass", subtype, mlevel, "multiomics", accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#	print(p)
#	
#	write.table(confusion_matrix, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.confusion_matrix.txt", "multiclass", subtype, mlevel, "multiomics"), quote=F, sep="\t", row.names=T, col.names=T)
#	## END BLOCK TO COMMENT ##
#	
#	# plotting - per-group sparse model
#	df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
#	colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
#	print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s %s", "multiclass", subtype, mlevel, "multiomics")))
#	# plotting - per-group variables
#	df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
#	df$metabolite_name <- unlist(lapply(as.character(df$OTU), function(metstr) {
#		x <- unlist(strsplit(metstr, " "))
#		paste(x[2:length(x)], collapse=" ")
#	}))
#	df$SampleType <- unlist(lapply(as.character(df$OTU), function(metstr) {
#		gsub("\\[", "", gsub("\\]", "", unlist(strsplit(metstr, " "))[1]))
#	}))
#	if (mlevel == "BIOCHEMICAL") {
#		df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
#		df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
#		df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
#	} else if (mlevel == "SUB.PATHWAY") {
#		df$subpathway <- df$metabolite_name
#		df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
#		df$OTU_string <- sprintf("%s", df$OTU)
#	} else {
#		df$superpathway <- df$metabolite_name
#		df$OTU_string <- sprintf("%s", df$OTU)
#	}
#	df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
#	# load effect sizes from linear regression
#	lmres <- {}
#	for (st in c("Plasma", "DBS")) {
#		tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote=""); tmp <- subset(tmp, MaternalRegimen==regi); rownames(tmp) <- tmp$metabolite
#		tmp$SampleType <- st
#		lmres <- rbind(lmres, tmp)
#	}
#	lmres$metabolite_str <- sprintf("[%s] %s", lmres$SampleType, lmres$metabolite)
#	rownames(lmres) <- lmres$metabolite_str
#	tmp <- levels(df$OTU); lmres <- subset(lmres, metabolite_str %in% tmp)
#	lmres$metabolite_str <- factor(lmres$metabolite_str, levels=levels(df$OTU))
#	lims <- max(abs(lmres$estimate) + abs(lmres$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
#	dircolors2 <- c("blue", "red", "black"); names(dircolors2) <- c("down", "up", "NS")
#	lmres[, "importance"] <- df[match(lmres$metabolite_str, df$OTU), "importance"]
#	p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("%s - %s %s explanatory %s", "multiclass", subtype, mlevel, "multiomics")) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
#	print(p)
#	p <- ggplot(lmres, aes(x=metabolite, y=estimate, color=dir, group=MaternalRegimen)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=MaternalRegimen), position=pd, hjust=1, color="black", size=2) + geom_tile(mapping=aes(x=metabolite, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates - %s %s explanatory %s", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
#	print(p)
#	p <- ggplot(lmres, aes(x=metabolite_str, y=estimate, color=dir, group=MaternalRegimen)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=MaternalRegimen), position=pd, hjust=1, color="black", size=2) + geom_tile(mapping=aes(x=metabolite, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + geom_vline(xintercept=as.numeric(df$OTU)+0.5, color="grey") + theme_classic() + ggtitle(sprintf("LM estimates - %s %s explanatory %s", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=dircolors2) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
#	print(p)
#	# shading rectangles of importance values
#	df.rect <- df
#	df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
##		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
#	p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s %s explanatory %s", "multiclass", subtype, mlevel, "multiomics")) + scale_fill_gradient(low="white", high="black")
#	print(p)
#	# violin plots of metabolite values
#	agg.melt <- agg.melt.stored
#	agg.melt$MaternalGroup <- mapping.sel[agg.melt$SampleID, "MaternalGroup"]
#	agg.melt$delgage <- mapping.sel[agg.melt$SampleID, "delgage"]
#	agg.melt$weight0week <- mapping.sel[agg.melt$SampleID, "weight0week"]
#	agg.melt$DaysPTDPlasma2 <- mapping.sel[agg.melt$SampleID, "DaysPTDPlasma2"]
#	agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$MaternalGroup))
#	agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
#	agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
#	p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.cohort)
#	print(p)
#	p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.cohort)
#	npages <- n_pages(p)
#	for (ip in 1:npages) {
#		p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.cohort)
#		print(p)
#	}
#	for (ip in 1:npages) {
#		p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=delgage)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_gradient(low="red", high="black")
#		print(p)
#	}
#	for (ip in 1:npages) {
#		p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=weight0week)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_gradient(low="red", high="black")
#		print(p)
#	}
#	for (ip in 1:npages) {
#		p <- ggplot(agg.melt, aes(x=MaternalGroup, y=value, color=DaysPTDPlasma2)) + geom_violin(aes(x=MaternalGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.daystodelivery)
#		print(p)
#	}
#	
#	agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
#	p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=MaternalGroup)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.cohort)
#	print(p)
#	agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
#	agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
#	p <- ggplot(agg.melt2, aes(x=MaternalGroup, y=value, color=MaternalGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.cohort)
#	print(p)
#	p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=MaternalGroup)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, "multiomics")) + coord_flip() + scale_color_manual(values=cols.cohort)
#	print(p)
#}



#########################################################################################################
### randomForest classification of Group (separately by Regimen); using METABOLITE data [DBS, Plasma] from Malawi StudySite 6 subgroup
#set.seed(nrow(mapping))	
## MAIN LOOP for random forest (through metabolite levels)
#subtype <- "maternal"
#for (st in c("DBS", "Plasma")) {
#	for (mlevel in "BIOCHEMICAL") {
#		data <- df.metabolon[[st]][[mlevel]]
#		mapping2 <- mapping[,c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "StudySite")]
#		colnames(mapping2) <- c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "StudySite")
#		rownames(mapping2) <- mapping2$patid
#		mapping2 <- subset(mapping2, Country=="Malawi" & StudySite==6)
#		for (regi in levels(mapping2$Regimen)) {
#			mapping.sel <- subset(mapping2, Regimen==regi); mapping.sel$Group <- droplevels(mapping.sel$Group)
#			data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the desired maternal samples from regimen
#			response <- droplevels(mapping.sel$Group); names(response) <- rownames(mapping.sel)
#			# add Country as covariates
#			data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
#			agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
#			
#			## after running for the first time, COMMENT OUT THIS BLOCK ##
#			num_iter <- 100
#			ncores <- 20
#			out <- mclapply(1:num_iter, function (dummy) {
#					importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
#				}, mc.cores=ncores )
#			collated.importance <- do.call(cbind, out)
#			out <- mclapply(1:num_iter, function (dummy) {
#					rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.9)$error.cv
#				}, mc.cores=ncores )
#			collated.cv <- do.call(cbind, out)

#			write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_Malawi6.%s.%s.%s.%s.importance.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
#			write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_Malawi6.%s.%s.%s.%s.cv.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
#			## END BLOCK TO COMMENT ##

#			collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_Malawi6.%s.%s.%s.%s.importance.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1, quote="")
#			collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_Malawi6.%s.%s.%s.%s.cv.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
#			importance.mean <- rowMeans(collated.importance)
#			importance.sd <- unlist(apply(collated.importance, 1, sd))
#			cv.mean <- rowMeans(collated.cv)
#			cv.sd <- unlist(apply(collated.cv, 1, sd))
#			inds <- order(importance.mean, decreasing=T)
#			inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate
#			write.table(melt(importance.mean[inds]), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_Malawi6.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

#			## after running for the first time, COMMENT OUT THIS BLOCK ##
#			# using a sparse model with N predictors
#			sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
#			save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_Malawi6.%s.%s.%s.%s.model", regi, subtype, mlevel, st))
#			load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_Malawi6.%s.%s.%s.%s.model", regi, subtype, mlevel, st))
#			# accuracy of final sparseRF model
#			pred <- predict(sparseRF, type="prob")
#			pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
#			pred_df_out <- merge(pred_df, data.sel, by="row.names")
#			write.table(pred_df_out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_Malawi6.%s.%s.%s.%s.predictions.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
#			confusion_matrix <- table(pred_df[, c("true", "predicted")])
#			class_errors <- unlist(lapply(levels(mapping.sel$Group), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$Group)
#			accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
#			vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
#			mccvalue <- mcc(vec.pred, vec.true)
#			df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
#			p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix Malawi6 (%s, %s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", regi, subtype, mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#			print(p)
#			
#			write.table(confusion_matrix, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_Malawi6.%s.%s.%s.%s.confusion_matrix.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=T)
#			## END BLOCK TO COMMENT ##
#			
#			# plotting - per-group sparse model
#			df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
#			colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
#			print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - Malawi6 %s %s %s %s", regi, subtype, mlevel, st)))
#			# plotting - per-group variables
#			df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
#			df$metabolite_name <- as.character(df$OTU)
#			if (mlevel == "BIOCHEMICAL") {
#				df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
#				df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
#				df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
#			} else if (mlevel == "SUB.PATHWAY") {
#				df$subpathway <- df$metabolite_name
#				df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
#				df$OTU_string <- sprintf("%s", df$OTU)
#			} else {
#				df$superpathway <- df$metabolite_name
#				df$OTU_string <- sprintf("%s", df$OTU)
#			}
#			df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
#			p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("Malawi6 %s - %s %s explanatory %s", regi, subtype, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
#			print(p) 
#			# shading rectangles of importance values
#			df.rect <- df
#			df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
#	#		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
#			p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s %s explanatory %s", regi, subtype, mlevel, st)) + scale_fill_gradient(low="white", high="black")
#			print(p)
#			# violin plots of metabolite values
#			agg.melt <- agg.melt.stored
#			agg.melt$Group <- mapping.sel[agg.melt$SampleID, "Group"]
#			agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$Group))
#			agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
#			agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
#			p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (Malawi6 %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
#			p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (Malawi6 %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			npages <- n_pages(p)
#			for (ip in 1:npages) {
#				p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#				print(p)
#			}
#			
#			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
#			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Group)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (Malawi6 %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
#			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
#			agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
#			p <- ggplot(agg.melt2, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (Malawi6 %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
#			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Group)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (Malawi6 %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
#			# violin plots stratified by prediction/truth
#			for (met in levels(agg.melt$metabolite)) {
#				tmp <- subset(agg.melt, metabolite==met)
#				# color scheme - manual
#				p <- ggplot(tmp, aes(x=Group, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (Malawi6 %s, %s, %s)", met, subtype, mlevel, st)) + scale_fill_manual(values=cols.cohort)
#				print(p)
#				p <- ggplot(tmp, aes(x=Group, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (Malawi6 %s, %s, %s)", met, subtype, mlevel, st)) + scale_color_manual(values=cols.cohort)
#				print(p)
#			}
#		}
#	}
#}


### randomForest classification of Group (separately by Regimen); using METABOLITE data [DBS, Plasma] from NON- Malawi StudySite 6 subgroup
#set.seed(nrow(mapping))	
## MAIN LOOP for random forest (through metabolite levels)
#subtype <- "maternal"
#for (st in c("DBS", "Plasma")) {
#	for (mlevel in "BIOCHEMICAL") {
#		data <- df.metabolon[[st]][[mlevel]]
#		mapping2 <- mapping[,c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "StudySite")]
#		colnames(mapping2) <- c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "StudySite")
#		rownames(mapping2) <- mapping2$patid
#		mapping2 <- subset(mapping2, !(Country=="Malawi" & StudySite==6))
#		for (regi in levels(mapping2$Regimen)) {
#			mapping.sel <- subset(mapping2, Regimen==regi); mapping.sel$Group <- droplevels(mapping.sel$Group)
#			data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the desired maternal samples from regimen
#			response <- droplevels(mapping.sel$Group); names(response) <- rownames(mapping.sel)
#			# add Country as covariates
#			data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
#			agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
#			
#			## after running for the first time, COMMENT OUT THIS BLOCK ##
#			num_iter <- 100
#			ncores <- 20
#			out <- mclapply(1:num_iter, function (dummy) {
#					importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
#				}, mc.cores=ncores )
#			collated.importance <- do.call(cbind, out)
#			out <- mclapply(1:num_iter, function (dummy) {
#					rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.9)$error.cv
#				}, mc.cores=ncores )
#			collated.cv <- do.call(cbind, out)

#			write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_nonMalawi6.%s.%s.%s.%s.importance.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
#			write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_nonMalawi6.%s.%s.%s.%s.cv.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
#			## END BLOCK TO COMMENT ##

#			collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_nonMalawi6.%s.%s.%s.%s.importance.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1, quote="")
#			collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_nonMalawi6.%s.%s.%s.%s.cv.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
#			importance.mean <- rowMeans(collated.importance)
#			importance.sd <- unlist(apply(collated.importance, 1, sd))
#			cv.mean <- rowMeans(collated.cv)
#			cv.sd <- unlist(apply(collated.cv, 1, sd))
#			inds <- order(importance.mean, decreasing=T)
#			inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate
#			write.table(melt(importance.mean[inds]), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_nonMalawi6.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

#			## after running for the first time, COMMENT OUT THIS BLOCK ##
#			# using a sparse model with N predictors
#			sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
#			save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_nonMalawi6.%s.%s.%s.%s.model", regi, subtype, mlevel, st))
#			load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_nonMalawi6.%s.%s.%s.%s.model", regi, subtype, mlevel, st))
#			# accuracy of final sparseRF model
#			pred <- predict(sparseRF, type="prob")
#			pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
#			pred_df_out <- merge(pred_df, data.sel, by="row.names")
#			write.table(pred_df_out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_nonMalawi6.%s.%s.%s.%s.predictions.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
#			confusion_matrix <- table(pred_df[, c("true", "predicted")])
#			class_errors <- unlist(lapply(levels(mapping.sel$Group), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$Group)
#			accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
#			vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
#			mccvalue <- mcc(vec.pred, vec.true)
#			df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
#			p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix nonMalawi6 (%s, %s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", regi, subtype, mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#			print(p)
#			
#			write.table(confusion_matrix, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_nonMalawi6.%s.%s.%s.%s.confusion_matrix.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=T)
#			## END BLOCK TO COMMENT ##
#			
#			# plotting - per-group sparse model
#			df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
#			colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
#			print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - nonMalawi6 %s %s %s %s", regi, subtype, mlevel, st)))
#			# plotting - per-group variables
#			df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
#			df$metabolite_name <- as.character(df$OTU)
#			if (mlevel == "BIOCHEMICAL") {
#				df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
#				df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
#				df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
#			} else if (mlevel == "SUB.PATHWAY") {
#				df$subpathway <- df$metabolite_name
#				df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
#				df$OTU_string <- sprintf("%s", df$OTU)
#			} else {
#				df$superpathway <- df$metabolite_name
#				df$OTU_string <- sprintf("%s", df$OTU)
#			}
#			df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
#			p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("nonMalawi6 %s - %s %s explanatory %s", regi, subtype, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
#			print(p) 
#			# shading rectangles of importance values
#			df.rect <- df
#			df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
#	#		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
#			p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s %s explanatory %s", regi, subtype, mlevel, st)) + scale_fill_gradient(low="white", high="black")
#			print(p)
#			# violin plots of metabolite values
#			agg.melt <- agg.melt.stored
#			agg.melt$Group <- mapping.sel[agg.melt$SampleID, "Group"]
#			agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$Group))
#			agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
#			agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
#			p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (nonMalawi6 %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
#			p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (nonMalawi6 %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			npages <- n_pages(p)
#			for (ip in 1:npages) {
#				p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (nonMalawi6 %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#				print(p)
#			}
#			
#			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
#			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Group)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (nonMalawi6 %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
#			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
#			agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
#			p <- ggplot(agg.melt2, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (nonMalawi6 %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
#			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Group)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (nonMalawi6 %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
#			# violin plots stratified by prediction/truth
#			for (met in levels(agg.melt$metabolite)) {
#				tmp <- subset(agg.melt, metabolite==met)
#				# color scheme - manual
#				p <- ggplot(tmp, aes(x=Group, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (nonMalawi6 %s, %s, %s)", met, subtype, mlevel, st)) + scale_fill_manual(values=cols.cohort)
#				print(p)
#				p <- ggplot(tmp, aes(x=Group, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (nonMalawi6 %s, %s, %s)", met, subtype, mlevel, st)) + scale_color_manual(values=cols.cohort)
#				print(p)
#			}
#		}
#	}
#}



### randomForest classification of Regimen (separately by Delivery); using METABOLITE data [DBS, Plasma]
#set.seed(nrow(mapping))	
## MAIN LOOP for random forest (through metabolite levels)
#subtype <- "maternal"
#for (st in c("DBS", "Plasma")) {
#	for (mlevel in "BIOCHEMICAL") {
#		data <- df.metabolon[[st]][[mlevel]]
#		mapping2 <- mapping[,c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
#		colnames(mapping2) <- c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
#		rownames(mapping2) <- mapping2$patid
#		for (deli in levels(mapping2$Delivery)) {
#			mapping.sel <- subset(mapping2, Delivery==deli); mapping.sel$Group <- droplevels(mapping.sel$Group)
#			data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the desired maternal samples from delivery
#			response <- droplevels(mapping.sel$Group); names(response) <- rownames(mapping.sel)
#			# add Country as covariates
#			data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
#			agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
#			
##			## after running for the first time, COMMENT OUT THIS BLOCK ##
##			num_iter <- 100
##			ncores <- 20
##			out <- mclapply(1:num_iter, function (dummy) {
##					importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
##				}, mc.cores=ncores )
##			collated.importance <- do.call(cbind, out)
##			out <- mclapply(1:num_iter, function (dummy) {
##					rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.9)$error.cv
##				}, mc.cores=ncores )
##			collated.cv <- do.call(cbind, out)

##			write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", deli, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
##			write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", deli, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
##			## END BLOCK TO COMMENT ##

#			collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", deli, subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1, quote="")
#			collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", deli, subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
#			importance.mean <- rowMeans(collated.importance)
#			importance.sd <- unlist(apply(collated.importance, 1, sd))
#			cv.mean <- rowMeans(collated.cv)
#			cv.sd <- unlist(apply(collated.cv, 1, sd))
#			inds <- order(importance.mean, decreasing=T)
#			inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate
#			write.table(melt(importance.mean[inds]), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", deli, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

#			## after running for the first time, COMMENT OUT THIS BLOCK ##
#			# using a sparse model with N predictors
##			sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
##			save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", deli, subtype, mlevel, st))
#			load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", deli, subtype, mlevel, st))
#			# accuracy of final sparseRF model
#			pred <- predict(sparseRF, type="prob")
#			pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
#			pred_df_out <- merge(pred_df, data.sel, by="row.names")
#			write.table(pred_df_out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.predictions.txt", deli, subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
#			confusion_matrix <- table(pred_df[, c("true", "predicted")])
#			class_errors <- unlist(lapply(levels(mapping.sel$Group), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$Group)
#			accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
#			vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
#			mccvalue <- mcc(vec.pred, vec.true)
#			df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
#			p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix for Regimen (%s, %s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", deli, subtype, mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#			print(p)
#			
#			write.table(confusion_matrix, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.confusion_matrix.txt", deli, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=T)
#			## END BLOCK TO COMMENT ##
#			
#			# plotting - per-group sparse model
#			df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
#			colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
#			print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection for Regimen - %s %s %s %s", deli, subtype, mlevel, st)))
#			# plotting - per-group variables
#			df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
#			df$metabolite_name <- as.character(df$OTU)
#			if (mlevel == "BIOCHEMICAL") {
#				df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
#				df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
#				df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
#			} else if (mlevel == "SUB.PATHWAY") {
#				df$subpathway <- df$metabolite_name
#				df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
#				df$OTU_string <- sprintf("%s", df$OTU)
#			} else {
#				df$superpathway <- df$metabolite_name
#				df$OTU_string <- sprintf("%s", df$OTU)
#			}
#			df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
#			# load effect sizes from linear regression
#			lmres <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Delivery.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote=""); lmres <- subset(lmres, Delivery==deli); rownames(lmres) <- lmres$metabolite
#			for (lmvar in c("estimate", "SE", "padj", "dir")) {
#				df[,lmvar] <- lmres[df$metabolite_name, lmvar]
#			}
#			p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("%s - %s %s explanatory %s", deli, subtype, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
#			print(p)
#			lims <- max(abs(df$estimate) + abs(df$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
#			p <- ggplot(df, aes(x=OTU, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=OTU, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_tile(aes(x=OTU, y=-lims*0.96, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s - %s %s explanatory %s", regi, subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
#			print(p)
#			# shading rectangles of importance values
#			df.rect <- df
#			df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
#	#		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
#			p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s %s explanatory %s", deli, subtype, mlevel, st)) + scale_fill_gradient(low="white", high="black")
#			print(p)
#			# violin plots of metabolite values
#			agg.melt <- agg.melt.stored
#			agg.melt$Group <- mapping.sel[agg.melt$SampleID, "Group"]
#			agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$Group))
#			agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
#			agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
#			p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of Regimen RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
#			p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of Regimen RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			npages <- n_pages(p)
#			for (ip in 1:npages) {
#				p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of Regimen RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#				print(p)
#			}
#			
#			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
#			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Group)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of Regimen RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
#			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
#			agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
#			p <- ggplot(agg.melt2, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of Regimen RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
#			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Group)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of Regimen RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
#			print(p)
##			# violin plots stratified by prediction/truth
##			for (met in levels(agg.melt$metabolite)) {
##				tmp <- subset(agg.melt, metabolite==met)
##				# color scheme - manual
##				p <- ggplot(tmp, aes(x=Group, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_fill_manual(values=cols.cohort)
##				print(p)
##				p <- ggplot(tmp, aes(x=Group, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of Regimen %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_color_manual(values=cols.cohort)
##				print(p)
##			}
#		}
#	}
#}

### randomForest classification of StudySite; using METABOLITE data [DBS, Plasma] from Malawi only
#set.seed(nrow(mapping))	
## MAIN LOOP for random forest (through metabolite levels)
#subtype <- "maternal"
#for (st in c("DBS", "Plasma")) {
#	for (mlevel in "BIOCHEMICAL") {
#		data <- df.metabolon[[st]][[mlevel]]
#		mapping2 <- mapping[,c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "StudySite")]
#		colnames(mapping2) <- c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "StudySite")
#		rownames(mapping2) <- mapping2$patid
#		mapping.sel <- subset(mapping2, Country=="Malawi")
#		data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the desired maternal samples
#		response <- droplevels(mapping.sel$StudySite); names(response) <- rownames(mapping.sel)
#		agg.melt.stored <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
#		
##		## after running for the first time, COMMENT OUT THIS BLOCK ##
##		num_iter <- 100
##		ncores <- 20
##		out <- mclapply(1:num_iter, function (dummy) {
##				importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
##			}, mc.cores=ncores )
##		collated.importance <- do.call(cbind, out)
##		out <- mclapply(1:num_iter, function (dummy) {
##				rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.9)$error.cv
##			}, mc.cores=ncores )
##		collated.cv <- do.call(cbind, out)

##		write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.importance.txt", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
##		write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.cv.txt", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
##		## END BLOCK TO COMMENT ##

#		collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.importance.txt", subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1, quote="")
#		collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.cv.txt", subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
#		importance.mean <- rowMeans(collated.importance)
#		importance.sd <- unlist(apply(collated.importance, 1, sd))
#		cv.mean <- rowMeans(collated.cv)
#		cv.sd <- unlist(apply(collated.cv, 1, sd))
#		inds <- order(importance.mean, decreasing=T)
#		inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate
#		write.table(melt(importance.mean[inds]), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.features.txt", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

#		## after running for the first time, COMMENT OUT THIS BLOCK ##
#		# using a sparse model with N predictors
##		sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
##		save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.model", subtype, mlevel, st))
#		load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.model", subtype, mlevel, st))
#		# accuracy of final sparseRF model
#		pred <- predict(sparseRF, type="prob")
#		pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
#		pred_df_out <- merge(pred_df, data.sel, by="row.names")
#		write.table(pred_df_out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.predictions.txt", subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
#		confusion_matrix <- table(pred_df[, c("true", "predicted")])
#		class_errors <- unlist(lapply(levels(mapping.sel$StudySite), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$StudySite)
#		accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
#		vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
#		mccvalue <- mcc(vec.pred, vec.true)
#		df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
#		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix StudySite (%s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", subtype, mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#		print(p)
#		
#		write.table(confusion_matrix, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.confusion_matrix.txt", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=T)
#		## END BLOCK TO COMMENT ##
#		
#		# plotting - per-group sparse model
#		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
#		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
#		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - StudySite %s %s %s", subtype, mlevel, st)))
#		# plotting - per-group variables
#		df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
#		df$metabolite_name <- as.character(df$OTU)
#		if (mlevel == "BIOCHEMICAL") {
#			df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
#			df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
#			df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
#		} else if (mlevel == "SUB.PATHWAY") {
#			df$subpathway <- df$metabolite_name
#			df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
#			df$OTU_string <- sprintf("%s", df$OTU)
#		} else {
#			df$superpathway <- df$metabolite_name
#			df$OTU_string <- sprintf("%s", df$OTU)
#		}
#		df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
#		p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("StudySite - %s %s explanatory %s", subtype, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
#		print(p) 
#		# shading rectangles of importance values
#		df.rect <- df
#		df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
##		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
#		p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("StudySite - %s %s explanatory %s", subtype, mlevel, st)) + scale_fill_gradient(low="white", high="black")
#		print(p)
#		# violin plots of metabolite values
#		agg.melt <- agg.melt.stored
#		agg.melt$StudySite <- mapping.sel[agg.melt$SampleID, "StudySite"]
#		agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$StudySite))
#		agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
#		agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
#		p <- ggplot(agg.melt, aes(x=StudySite, y=value, color=StudySite)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (StudySite %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#		print(p)
#		p <- ggplot(agg.melt, aes(x=StudySite, y=value, color=StudySite)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (StudySite %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#		npages <- n_pages(p)
#		for (ip in 1:npages) {
#			p <- ggplot(agg.melt, aes(x=StudySite, y=value, color=StudySite)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#			print(p)
#		}
#		
#		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
#		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=StudySite)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (StudySite %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#		print(p)
#		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
#		agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
#		p <- ggplot(agg.melt2, aes(x=StudySite, y=value, color=StudySite)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (StudySite %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#		print(p)
#		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=StudySite)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (StudySite %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#		print(p)
#		
#		# t-tests
#		agg.melt.stored$StudySite <- mapping.sel[agg.melt.stored$SampleID, "StudySite"]
#		res <- {}
#		for (met in unique(agg.melt.stored$metabolite)) {
#			tmp <- subset(agg.melt.stored, metabolite==met)
#			if (length(unique(tmp$value))==1) {
#				next
#			}
#			test <- t.test(value ~ StudySite, tmp, alternative="two.sided")
#			res <- rbind(res, c(met, test$estimate, test$statistic, test$p.value, test$method))
#		}
#		res <- as.data.frame(res); colnames(res) <- c("metabolite", colnames(res)[2:3], "t", "pval", "method")
#		res$pval <- as.numeric(as.character(res$pval)); res <- res[order(res$pval),]
#		res$padj <- p.adjust(res$pval, method="fdr")
#		write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/differences_by_StudySite.t-test.%s.%s.%s.txt", subtype, mlevel, st), quote=F, row.names=F, col.names=T, sep="\t")
#	}
#}


#########################################################################################################
### infant metabolites (DBS only) (first 2 days of life)
mapping <- subset(mappinglist[["DBS"]], InfantAgeInDays <= 3) # first 3 days of life
#mapping <- subset(mappinglist[["DBS"]], InfantAgeInDays <= 1) # first 1 day of life
#mapping <- subset(mappinglist[["DBS"]], InfantAgeInDays <= 30) # or first 30 days of life?? pretty much everyone is in first month, except a 43 and 183

## Cohort demographics and some QC data about metabolomics
mapping.demo <- subset(mapping, InfantRegimen %in% c("zdv", "PI-ART"))
mapping.demo$InfantGroup <- droplevels(mapping.demo$InfantGroup)
mapping.demo$InfantRegimen <- droplevels(mapping.demo$InfantRegimen)
tab1 <- CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "delgage", "hemaval.infant", "InfantAgeInDays", "DaysPTDDBS2"), strata=c("InfantGroup"), data=mapping.demo, smd=T)
write.table(print(tab1), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_2.InfantGroup.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(print(summary(tab1$ContTable)), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_2_detailed.InfantGroup.txt", quote=F, sep="\t", row.names=T, col.names=T)
tab1 <- CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "delgage", "hemaval.infant", "InfantAgeInDays", "DaysPTDDBS2"), strata=c("Delivery"), data=mapping.demo, smd=T)
write.table(print(tab1), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_2.Delivery.txt", quote=F, sep="\t", row.names=T, col.names=T)

mapping.demo <- subset(mapping, InfantRegimen %in% c("other", "zdv", "PI-ART"))
mapping.demo$InfantGroup <- droplevels(mapping.demo$InfantGroup)
mapping.demo$InfantRegimen <- droplevels(mapping.demo$InfantRegimen)
tab1 <- CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "delgage", "hemaval.infant", "InfantAgeInDays", "DaysPTDDBS2"), strata=c("InfantGroup"), data=mapping.demo, smd=T)
write.table(print(tab1), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_2_full.InfantGroup.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(print(summary(tab1$ContTable)), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_2_full_detailed.InfantGroup.txt", quote=F, sep="\t", row.names=T, col.names=T)
tab1 <- CreateTableOne(vars=c("gender", "weight0week", "weight1week", "ap_onstgage", "Country", "GestationalAgeAtCollection", "delgage", "hemaval.infant", "InfantAgeInDays", "DaysPTDDBS2"), strata=c("Delivery"), data=mapping.demo, smd=T)
write.table(print(tab1), file="/Lab_Share/PROMISE/nwcs619/metabolon/Table_2_full.Delivery.txt", quote=F, sep="\t", row.names=T, col.names=T)


## ordination (t-SNE, PCA) and PERMANOVA
set.seed(nrow(mapping))
subtype <- "infant"
for (st in c("DBS")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mapping[,c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Infant", "hemaval.infant", "InfantAgeInDays", "InfantAgeInDaysBinned")]
	colnames(mapping.sel) <- c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "InfantAgeInDays", "InfantAgeInDaysBinned")	
	rownames(mapping.sel) <- mapping.sel$cpatid
	data <- data[rownames(mapping.sel),] # subset to just the infant samples
	to_remove <- names(which(apply(data, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
	data <- data[, setdiff(colnames(data), to_remove)]
	
	# PCA
	pca <- prcomp(data, center=F, scale=F)
	eigs <- pca$sdev^2
	pvar <- 100*(eigs / sum(eigs))
	df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], SampleID=rownames(pca$x))
	for (mvar in intersect(rownames(subset(metadata_variables, useForPERMANOVA=="yes")), colnames(mapping.sel))) {
		df[, mvar] <- mapping.sel[rownames(df), mvar]
		if (metadata_variables[mvar, "type"] == "factor") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s %s %s PCoA (Euclidean distance)", subtype, st, mlevel, mvar)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_brewer(palette="Set1") + stat_ellipse(type="t")
		} else if (metadata_variables[mvar, "type"] == "numeric") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s %s %s PCoA (Euclidean distance)", subtype, st, mlevel, mvar)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_gradient(low="grey", high="red")
		}
		else {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s %s %s PCoA (Euclidean distance)", subtype, st, mlevel, mvar)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + stat_ellipse(type="t")
		}
		print(p)
	}
	p <- ggplot(df, aes(x=PC1, y=PC2, colour=InfantRegimen, shape=Delivery, group=InfantRegimen)) + geom_point(size=2) + theme_classic() + ggtitle(sprintf("%s %s %s PCoA (Euclidean distance)", subtype, st, mlevel)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_brewer(palette="Set1") + stat_ellipse() + scale_shape_manual(values=c(1,3))
	print(p)
	
	# PERMANOVA
	res <- adonis2(data ~ Delivery + InfantRegimen + Country + GestationalAgeAtCollection + InfantAgeInDays, data=mapping.sel, permutations=999, method='euclidean')
	sink(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/metabolon_PERMANOVA.%s.%s.%s.txt", subtype, st, mlevel))
	print(res)
	sink()
	# t-SNE
	tsne.out <- Rtsne(data, perplexity=10)
	df <- data.frame(PC1 = tsne.out$Y[,1], PC2 = tsne.out$Y[,2], SampleID=rownames(data)); rownames(df) <- df$SampleID
	for (mvar in intersect(rownames(subset(metadata_variables, useForPERMANOVA=="yes")), colnames(mapping.sel))) {
		df[, mvar] <- mapping.sel[rownames(df), mvar]
		if (metadata_variables[mvar, "type"] == "factor") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s %s %s tSNE", subtype, st, mlevel, mvar)) + scale_color_brewer(palette="Set1") + stat_ellipse(type="t")
		} else if (metadata_variables[mvar, "type"] == "numeric") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s %s %s tSNE", subtype, st, mlevel, mvar)) + scale_color_gradient(low="grey", high="red")
		}
		else {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s %s %s tSNE", subtype, st, mlevel, mvar)) + stat_ellipse(type="t")
		}
		print(p)
	}
}

## mediation analysis (IV=ART, DV=metabolome, mediator=PTB; use mediation package)
subtype <- "infant"; iv <- "InfantRegimen"; mediator <- "Delivery"
set.seed(nrow(mapping))
for (st in c("DBS")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mapping[,c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.infant", "weight0week")]
	colnames(mapping.sel) <- c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weight0week")
	mapping.sel <- subset(mapping.sel, InfantRegimen %in% c("zdv", "PI-ART")) # exclude other regimen because insufficient numbers
	mapping.sel$InfantRegimen <- droplevels(mapping.sel$InfantRegimen)
	rownames(mapping.sel) <- mapping.sel$cpatid
	data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the infant samples
	to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
	data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
	df <- merge(data.sel, mapping.sel, by="row.names"); df[, mediator] <- as.numeric(df[, mediator])-1
	res <- {}
	for (metabolite in colnames(data.sel)) {
		fit.mediator <- glm(as.formula(sprintf("%s ~ %s", mediator, iv)), family="binomial", data=df)
		fit.dv <- lm(as.formula(sprintf("%s ~ %s + %s", metabolite, iv, mediator)), data=df)
		med2 <- try(mediation::mediate(fit.mediator, fit.dv, treat=iv, mediator=mediator, boot=T, control.value="zdv", treat.value="PI-ART", parallel="multicore", ncpus=16), silent=T)
		if (class(med2)=="mediate") {
			summ <- summary(med2)
			res <- rbind(res, c(metabolite, unlist(summ[c("d0", "d0.ci", "d0.p", "d1", "d1.ci", "d1.p", "z0", "z0.ci", "z0.p", "z1", "z1.ci", "z1.p", "d.avg", "d.avg.ci", "d.avg.p", "z.avg", "z.avg.ci", "z.avg.p", "tau.coef", "tau.ci", "tau.p", "n0", "n0.ci", "n0.p", "n1", "n1.ci", "n1.p", "n.avg", "n.avg.ci", "n.avg.p")]) ))
		} else {
			res <- rbind(res, c(metabolite, rep(NA, 40)))
		}
	}
	
	res <- as.data.frame(res)
	colnames(res) <- c("metabolite", "ACME.control", "ACME.control.ci.lower", "ACME.control.ci.upper", "ACME.control.p", "ACME.treated", "ACME.treated.ci.lower", "ACME.treated.ci.upper", "ACME.treated.p", "ADE.control", "ADE.control.ci.lower", "ADE.control.ci.upper", "ADE.control.p", "ADE.treated", "ADE.treated.ci.lower", "ADE.treated.ci.upper", "ADE.treated.p", "ACME.average", "ACME.average.ci.lower", "ACME.average.ci.upper", "ACME.average.p", "ADE.average", "ADE.average.ci.lower", "ADE.average.ci.upper", "ADE.average.ci.p", "totaleffect", "totaleffect.ci.lower", "totaleffect.ci.upper", "totaleffect.p", "propmediated.control.avg", "propmediated.control.avg.ci.lower", "propmediated.control.avg.ci.upper", "propmediated.control.p", "propmediated.treatment.avg", "propmediated.treatment.avg.ci.lower", "propmediated.treatment.avg.ci.upper", "propmediated.treatment.p", "propmediated.avg", "propmediated.avg.ci.lower", "propmediated.avg.ci.upper", "propmediated.p")
	res$metabolite <- as.character(name_map[as.character(res$metabolite), "original"])
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/mediation_analysis.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)

}


## violin plots of all metabolite values
subtype <- "infant"
for (st in c("DBS")) {
	for (mlevel in "BIOCHEMICAL") {
		data <- df.metabolon[[st]][[mlevel]]
		mapping.sel <- mapping[,c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.infant", "weight0week")]
		colnames(mapping.sel) <- c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weight0week")
		mapping.sel <- subset(mapping.sel, InfantRegimen %in% c("zdv", "PI-ART")) # exclude other regimen because insufficient numbers
		rownames(mapping.sel) <- mapping.sel$cpatid
		data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the infant samples
		agg.melt <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt) <- c("SampleID", "metabolite", "value")
		agg.melt$InfantGroup <- droplevels(mapping.sel[agg.melt$SampleID, "InfantGroup"])
		agg.melt$InfantGroup <- factor(as.character(agg.melt$InfantGroup), levels=rev(c("Preterm.zdv", "Term.zdv", "Preterm.PI-ART", "Term.PI-ART")))
		# manual pagination
		metabolites <- unique(agg.melt$metabolite)
		for (i in seq(from=1,to=length(metabolites), by=12)) {
			j <- min(i+11, length(metabolites))
			agg.melt.sel <- subset(agg.melt, metabolite %in% metabolites[i:j])
			p <- ggplot(agg.melt.sel, aes(x=InfantGroup, y=value, color=InfantGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3, nrow=4) + theme_classic() + ggtitle(sprintf("Rel. abund. of all metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort) + theme(strip.text = element_text(size=7))
			print(p)
		}
	}
}

## linear model with emmeans, stratified by Regimen
subtype <- "infant"; mvar <- "Delivery"
for (st in c("DBS")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mapping[,c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Infant", "hemaval.infant", "InfantAgeInDays")]
	colnames(mapping.sel) <- c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "InfantAgeInDays")
	rownames(mapping.sel) <- mapping.sel$cpatid
	mapping.sel <- subset(mapping.sel, InfantRegimen %in% c("zdv", "PI-ART")); mapping.sel$InfantRegimen <- droplevels(mapping.sel$InfantRegimen) # exclude 'untreated and 'other' regimen because too few samples
	data.sel <- data[rownames(mapping.sel),] # subset to just the infant samples
	to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
	data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
	df <- merge(data.sel, mapping.sel, by="row.names"); 
	df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
	res <- {}
	for (metabolite in colnames(data.sel)) {
		mod <- lm(as.formula(sprintf("%s ~ %s*InfantRegimen + InfantAgeInDays", metabolite, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s | InfantRegimen", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
		res <- rbind(res, tmp)
	}
	res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
	res$padj <- p.adjust(res$p.value, method="fdr")
	res <- res[order(res$p.value, decreasing=F),]
	res$dir <- ifelse(res$padj < siglevel_loose, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
	res$dir2 = ifelse(res$estimate < 0, "down", "up")
	res$exp_estimate <- exp(res$estimate)
	for (addvar in c("COMP_ID", "PLATFORM")) {
		res[, addvar] <- metabolon_map_by_assay[[st]][as.character(res$metabolite), addvar]
	}
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
	# forest plot
	resSig <- subset(res, padj < siglevel_loose)
	df <- subset(res, metabolite %in% as.character(resSig$metabolite))
	if (nrow(df)>0) {
		df <- df[order(df$estimate, decreasing=T),]
		df$metabolite <- factor(as.character(df$metabolite), levels=as.character(unique(df$metabolite)))
		lims <- max(abs(df$estimate) + abs(df$SE))*1.0
		pd <- position_dodge(0.8)
		p <- ggplot(df, aes(x=metabolite, y=estimate, color=dir, group=InfantRegimen)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=InfantRegimen), position=pd, hjust=1, color="black", size=2) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("Delivery stratified by InfantRegimen (%s, %s, %s)", subtype, st, mlevel)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
		print(p)
	}
}

## linear model with emmeans, averaged across Regimen
subtype <- "infant"; mvar <- "Delivery"
for (st in c("DBS")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mapping[,c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Infant", "hemaval.infant", "InfantAgeInDays")]
	colnames(mapping.sel) <- c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "InfantAgeInDays")
	rownames(mapping.sel) <- mapping.sel$cpatid
	mapping.sel <- subset(mapping.sel, InfantRegimen %in% c("zdv", "PI-ART")); mapping.sel$InfantRegimen <- droplevels(mapping.sel$InfantRegimen) # exclude 'untreated and 'other' regimen because too few samples
	data.sel <- data[rownames(mapping.sel),] # subset to just the infant samples
	to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
	data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
	df <- merge(data.sel, mapping.sel, by="row.names"); 
	df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
	res <- {}
	for (metabolite in colnames(data.sel)) {
		mod <- lm(as.formula(sprintf("%s ~ %s+InfantRegimen+InfantAgeInDays", metabolite, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
		res <- rbind(res, tmp)
	}
	res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
	res$padj <- p.adjust(res$p.value, method="fdr")
	res <- res[order(res$p.value, decreasing=F),]
	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
	res$exp_estimate <- exp(res$estimate)
	for (addvar in c("COMP_ID", "PLATFORM")) {
		res[, addvar] <- metabolon_map_by_assay[[st]][as.character(res$metabolite), addvar]
	}
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_averaged_across_Regimen.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
	# forest plot
	resSig <- subset(res, padj < siglevel)
	df <- subset(res, metabolite %in% as.character(resSig$metabolite))
	if (nrow(df)>0) {
		df <- df[order(df$estimate, decreasing=T),]
		df$metabolite <- factor(as.character(df$metabolite), levels=as.character(unique(df$metabolite)))
		lims <- max(abs(df$estimate) + abs(df$SE))*1.0
		pd <- position_dodge(0.8)
		p <- ggplot(df, aes(x=metabolite, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("Delivery averaged across InfantRegimen (%s, %s, %s)", subtype, st, mlevel)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
		print(p)
	}
}

## linear model with emmeans, outcome Delivery, stratified by Regimen+Country
subtype <- "infant"; mvar <- "Delivery"
for (st in c("DBS")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mapping[,c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "InfantAgeInDays")]
	colnames(mapping.sel) <- c("Delivery", "InfantRegimen", "InfantGroup", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "InfantAgeInDays")
	rownames(mapping.sel) <- mapping.sel$patid
	mapping.sel <- subset(mapping.sel, InfantRegimen %in% c("zdv", "PI-ART")); mapping.sel$InfantRegimen <- droplevels(mapping.sel$InfantRegimen) # exclude 'untreated and 'other' regimen because too few samples
	data.sel <- data[rownames(mapping.sel),] # subset to just the infant samples
	to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
	data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
#	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
	df <- merge(data.sel, mapping.sel, by="row.names"); 
	df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
	res <- {}
	for (metabolite in colnames(data.sel)) {
		mod <- lm(as.formula(sprintf("%s ~ %s*InfantRegimen + %s*Country + InfantAgeInDays", metabolite, mvar, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s | InfantRegimen + Country", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
		res <- rbind(res, tmp)
	}
	res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
	res$padj <- p.adjust(res$p.value, method="fdr")
	res <- res[order(res$p.value, decreasing=F),]
	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
#	res$dir2 = ifelse(res$estimate < 0, "down", "up"); tab <- table(res[,c("Regimen","dir2")])
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_RegimenCountry.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
}


## LASSO regression of Group (separately by Regimen); using METABOLITE data [DBS]
set.seed(nrow(mapping))
subtype <- "infant"; mvar <- "Delivery"
for (st in c("DBS")) {
	for (mlevel in "BIOCHEMICAL") {
		data <- df.metabolon[[st]][[mlevel]]
		mapping2 <- mapping[,c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Infant", "hemaval.infant", "weight0week")]
		colnames(mapping2) <- c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "weight0week")
		rownames(mapping2) <- mapping2$cpatid
		res.mean <- {}; res.sd <- {}
		for (regi in setdiff(levels(mapping2$InfantRegimen), c("other", "untreated"))) {
			mapping.sel <- subset(mapping2, InfantRegimen==regi); mapping.sel$InfantGroup <- droplevels(mapping.sel$InfantGroup)
			data.sel <- data[rownames(mapping.sel),] # subset to just the desired infant samples from regimen
			to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
			data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
			response <- droplevels(mapping.sel$InfantGroup); names(response) <- rownames(mapping.sel)
			# add Country as covariates
#			data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
			agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
#			fit <- glmnet(x = data.sel, y = response, alpha = 1, family="binomial")
#			cvfit <- cv.glmnet(data.sel, response, family="binomial", type.measure="class")
			# do CV 100 reps because it seems to have high variability
			out <- mclapply(1:ncvreps, function(dummy) {
				fit <- cv.glmnet(data.sel, response, family="binomial", type.measure="class")
				errors = data.frame(fit$lambda,fit$cvm)
				errors
			}, mc.cores=16)
			lambdas <- do.call(rbind, out)
			# take mean cvm for each lambda
			df <- aggregate(fit.cvm ~ fit.lambda, lambdas, mean); colnames(df) <- c("lambda", "cvm")
			df$cvsd <- aggregate(fit.cvm ~ fit.lambda, lambdas, sd)[,2]
			df$loglambda <- log(df$lambda)
			# select the best one
			bestindex = which.min(df$cvm); bestlambda = df[bestindex, "lambda"]
			fit <- glmnet(data.sel, response, alpha=1, family="binomial")
			lasso_coef = predict(fit, type="coefficients", s=bestlambda, exact=T)
			# plot CV fit (log-lambda vs CV error)
			p <- ggplot(df, aes(x=loglambda, y=cvm)) + geom_point(color="red") + geom_errorbar(aes(x=loglambda, ymin=cvm-cvsd, ymax=cvm+cvsd)) + geom_vline(xintercept=log(bestlambda), linetype="dotted") + theme_classic() + ggtitle(sprintf("CV fit for LASSO (%s %s %s %s)", regi, subtype, st, mlevel))
			print(p)
#			print(plot(cvfit, main=sprintf("CV fit for LASSO regression (%s %s %s %s)", regi, subtype, st, mlevel)))
			res <- as.data.frame(as.matrix(lasso_coef)); colnames(res) <- c("value")
			res <- subset(res, value != 0); res$taxa <- rownames(res)
			res <- subset(res, taxa != "(Intercept)")
			if (nrow(res)>0) {
				res <- res[order(res$value, decreasing=T),]; res$taxa <- factor(res$taxa, levels=res$taxa)
				res.mean <- rbind(res.mean, res)
				p <- ggplot(res, aes(x=taxa, y=value)) + geom_bar(stat="identity") + theme_classic() + ggtitle(sprintf("[LASSO] %s ~ selected features (%s %s %s, lambda=%.4g)", regi, subtype, st, mlevel, bestlambda)) + theme(plot.title=element_text(size=8), axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5))
				print(p)
				
				# violin plots of metabolite values
				agg.melt <- agg.melt.stored
				agg.melt$InfantGroup <- mapping.sel[agg.melt$SampleID, "InfantGroup"]
				agg.melt$delgage <- mapping.sel[agg.melt$SampleID, "delgage"]
				agg.melt$weight0week <- mapping.sel[agg.melt$SampleID, "weight0week"]
				agg.melt <- subset(agg.melt, metabolite %in% rownames(res))
				agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(res))
				p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of LASSO metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
				p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of LASSO metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				npages <- n_pages(p)
				for (ip in 1:npages) {
					p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of LASSO metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
					print(p)
				}
				for (ip in 1:npages) {
					p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=delgage)) + geom_violin(aes(x=InfantGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of LASSO metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
					print(p)
				}
				for (ip in 1:npages) {
					p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=weight0week)) + geom_violin(aes(x=InfantGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of LASSO metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
					print(p)
				}
			}
			write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/LASSO.%s.%s.%s.%s.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
		}
	}
}


## randomForest classification of Group (separately by Regimen); using METABOLITE data [DBS]
set.seed(nrow(mapping))	
# MAIN LOOP for random forest (through metabolite levels)
subtype <- "infant"
for (st in c("DBS")) {
	for (mlevel in "BIOCHEMICAL") {
		data <- df.metabolon[[st]][[mlevel]]
		mapping2 <- mapping[,c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Infant", "hemaval.infant", "InfantAgeInDays", "InfantAgeInDaysBinned", "weight0week")]
		colnames(mapping2) <- c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "InfantAgeInDays", "InfantAgeInDaysBinned", "weight0week")
		rownames(mapping2) <- mapping2$cpatid
		for (regi in setdiff(levels(mapping2$InfantRegimen), c("other", "untreated"))) {
			mapping.sel <- subset(mapping2, InfantRegimen==regi); mapping.sel$InfantGroup <- droplevels(mapping.sel$InfantGroup)
			data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the desired infant samples from regimen
			to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
			data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
			response <- droplevels(mapping.sel$InfantGroup); names(response) <- rownames(mapping.sel)
			# add [Country,InfantAgeInDays] as covariates
			data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
			data.sel$InfantAgeInDays <- mapping.sel[rownames(data.sel), "InfantAgeInDays"]
			agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), c("Country", "InfantAgeInDays"))]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
			
			## after running for the first time, COMMENT OUT THIS BLOCK ##
			num_iter <- 100
			ncores <- 20
			out <- mclapply(1:num_iter, function (dummy) {
					importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
				}, mc.cores=ncores )
			collated.importance <- do.call(cbind, out)
			out <- mclapply(1:num_iter, function (dummy) {
					rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.9)$error.cv
				}, mc.cores=ncores )
			collated.cv <- do.call(cbind, out)

			write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
			write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
			## END BLOCK TO COMMENT ##

			collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1, quote="")
			collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
			importance.mean <- rowMeans(collated.importance)
			importance.sd <- unlist(apply(collated.importance, 1, sd))
			cv.mean <- rowMeans(collated.cv)
			cv.sd <- unlist(apply(collated.cv, 1, sd))
			inds <- order(importance.mean, decreasing=T)
			inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate
			write.table(melt(importance.mean[inds]), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

			## after running for the first time, COMMENT OUT THIS BLOCK ##
			# using a sparse model with N predictors
			sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
			save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", regi, subtype, mlevel, st))
			load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", regi, subtype, mlevel, st))
			# accuracy of final sparseRF model
			pred <- predict(sparseRF, type="prob")
			pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
			pred_df_out <- merge(pred_df, data.sel, by="row.names")
			write.table(pred_df_out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.predictions.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
			confusion_matrix <- table(pred_df[, c("true", "predicted")])
			class_errors <- unlist(lapply(levels(mapping.sel$InfantGroup), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$InfantGroup)
			accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
			vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
			mccvalue <- mcc(vec.pred, vec.true)
			df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
			cf <- confusionMatrix(confusion_matrix, positive=sprintf("Preterm.%s", regi))
			cflabel <- sprintf("Positive: %s Sens: %.4g  Spec: %.4g\n PPV: %.4g  NPV: %.4g", cf$positive, cf$byClass[["Sensitivity"]], cf$byClass[["Specificity"]], cf$byClass[["Pos Pred Value"]], cf$byClass[["Neg Pred Value"]])
			p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", regi, subtype, mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + annotate(geom="text", x=2.5, y=10, label=cflabel)
			print(p)
			
			pred2 <- prediction(pred[,2], ordered(response))
			perf <- performance(pred2, "tpr", "fpr")
			perf.auc <- performance(pred2, "auc")
			plot(perf, main=sprintf("ROC %s %s %s %s (sparseRF final model)", regi, subtype, mlevel, st)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", unlist(perf.auc@y.values)))
			to_store <- data.frame(fpr=perf@x.values[[1]], tpr=perf@y.values[[1]], alpha=perf@alpha.values[[1]], SampleType=st, Regimen=regi, Subtype=subtype, level=mlevel)
			
			write.table(confusion_matrix, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.confusion_matrix.txt", regi, subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=T)
			## END BLOCK TO COMMENT ##
			
			# plotting - per-group sparse model
			df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
			colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
			print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s %s", regi, subtype, mlevel, st)))
			# plotting - per-group variables
			df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
			df$metabolite_name <- as.character(df$OTU)
			if (mlevel == "BIOCHEMICAL") {
				df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
				df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
				df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
			} else if (mlevel == "SUB.PATHWAY") {
				df$subpathway <- df$metabolite_name
				df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
				df$OTU_string <- sprintf("%s", df$OTU)
			} else {
				df$superpathway <- df$metabolite_name
				df$OTU_string <- sprintf("%s", df$OTU)
			}
			df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
			lmres <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote=""); lmres <- subset(lmres, InfantRegimen==regi); rownames(lmres) <- lmres$metabolite
			for (lmvar in c("estimate", "SE", "padj", "dir")) {
				df[,lmvar] <- lmres[df$metabolite_name, lmvar]
			}
			p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("%s - %s %s explanatory %s", regi, subtype, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
			print(p)			
			lims <- max(abs(df$estimate) + abs(df$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
			p <- ggplot(df, aes(x=OTU, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=OTU, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_tile(aes(x=OTU, y=-lims*0.96, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s - %s %s explanatory %s", regi, subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
			print(p)
			# stratified by country
			lmres.bycountry <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_RegimenCountry.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote="")
			tmp <- levels(df$OTU); lmres.bycountry <- subset(lmres.bycountry, metabolite %in% tmp & InfantRegimen==regi)
			lmres.bycountry$metabolite <- factor(lmres.bycountry$metabolite, levels=levels(df$OTU))
			lims <- max(abs(lmres.bycountry$estimate) + abs(lmres.bycountry$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
			dircolors2 <- c("blue", "red", "black"); names(dircolors2) <- c("down", "up", "NS")
			lmres.bycountry[, "importance"] <- df[match(lmres.bycountry$metabolite, df$OTU), "importance"]
			p <- ggplot(lmres.bycountry, aes(x=metabolite, y=estimate, color=dir, group=Country)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Country), position=pd, hjust=1, color="black", size=1.5) + geom_tile(mapping=aes(x=metabolite, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s - %s %s explanatory %s", regi, subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
			print(p)
			p <- ggplot(lmres.bycountry, aes(x=metabolite, y=estimate, color=dir, group=Country)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Country), position=pd, hjust=1, color="black", size=1.5) + geom_tile(mapping=aes(x=metabolite, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + geom_vline(xintercept=as.numeric(df$OTU)+0.5, color="grey") + theme_classic() + ggtitle(sprintf("LM estimates %s - %s %s explanatory %s", regi, subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors2) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
			print(p)
			# shading rectangles of importance values
			df.rect <- df
			df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
	#		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
			p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s %s explanatory %s", regi, subtype, mlevel, st)) + scale_fill_gradient(low="white", high="black")
			print(p)
			# violin plots of metabolite values
			agg.melt <- agg.melt.stored
			agg.melt$InfantGroup <- mapping.sel[agg.melt$SampleID, "InfantGroup"]
			agg.melt$InfantAgeInDays <- mapping.sel[agg.melt$SampleID, "InfantAgeInDays"]
			agg.melt$InfantAgeInDaysBinned <- mapping.sel[agg.melt$SampleID, "InfantAgeInDaysBinned"]
			agg.melt$delgage <- mapping.sel[agg.melt$SampleID, "delgage"]
			agg.melt$weight0week <- mapping.sel[agg.melt$SampleID, "weight0week"]
			agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$InfantGroup))
			agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
			agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
			p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			npages <- n_pages(p)
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
			}
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantAgeInDaysBinned)) + geom_violin(aes(x=InfantGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
				print(p)
			}
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantAgeInDays)) + geom_violin(aes(x=InfantGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip()
				print(p)
			}
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=delgage)) + geom_violin(aes(x=InfantGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
				print(p)
			}
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=weight0week)) + geom_violin(aes(x=InfantGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_gradient(low="red", high="black")
				print(p)
			}
			
			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=InfantGroup)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
			agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
			p <- ggplot(agg.melt2, aes(x=InfantGroup, y=value, color=InfantGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=InfantGroup)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
#			# violin plots stratified by prediction/truth
#			for (met in levels(agg.melt$metabolite)) {
#				tmp <- subset(agg.melt, metabolite==met)
#				# color scheme - manual
#				p <- ggplot(tmp, aes(x=InfantGroup, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_fill_manual(values=cols.cohort)
#				print(p)
#				p <- ggplot(tmp, aes(x=InfantGroup, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_color_manual(values=cols.cohort)
#				print(p)
#			}
		}
	}
}


## randomForest classification of Group (multiclass); using METABOLITE data [DBS]
set.seed(nrow(mapping))	
# MAIN LOOP for random forest (through metabolite levels)
subtype <- "infant"
for (st in c("DBS")) {
	for (mlevel in "BIOCHEMICAL") {
		data <- df.metabolon[[st]][[mlevel]]
		mapping.sel <- mapping[,c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Infant", "hemaval.infant", "InfantAgeInDays", "InfantAgeInDaysBinned", "weight0week")]
		colnames(mapping.sel) <- c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "InfantAgeInDays", "InfantAgeInDaysBinned", "weight0week")
		rownames(mapping.sel) <- mapping.sel$cpatid
		mapping.sel <- subset(mapping.sel, InfantRegimen %in% c("zdv", "PI-ART")) # exclude untreated and other regimen because insufficient numbers
		data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the infant samples
		to_remove <- names(which(apply(data.sel, 2, function(x) any(is.na(x))))) # drop metabolites with NAs (due to no variation for Z-transform)
		data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)]
		response <- droplevels(mapping.sel$InfantGroup); names(response) <- rownames(mapping.sel); group_levels <- levels(response)
		# add [Country, InfantAgeInDays] as covariates
		data.sel$Country <- mapping.sel[rownames(data.sel), "Country"]
		data.sel$InfantAgeInDays <- mapping.sel[rownames(data.sel), "InfantAgeInDays"]
		agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), c("Country", "InfantAgeInDays"))]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
		
		## after running for the first time, COMMENT OUT THIS BLOCK ##
		num_iter <- 100
		ncores <- 20
		out <- mclapply(1:num_iter, function (dummy) {
				importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
			}, mc.cores=ncores )
		collated.importance <- do.call(cbind, out)
		out <- mclapply(1:num_iter, function (dummy) {
				rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.9)$error.cv
			}, mc.cores=ncores )
		collated.cv <- do.call(cbind, out)

		write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
		write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
		## END BLOCK TO COMMENT ##

		collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", "multiclass", subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1, quote="")
		collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", "multiclass", subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
		importance.mean <- rowMeans(collated.importance)
		importance.sd <- unlist(apply(collated.importance, 1, sd))
		cv.mean <- rowMeans(collated.cv)
		cv.sd <- unlist(apply(collated.cv, 1, sd))
		inds <- order(importance.mean, decreasing=T)
		inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate
		write.table(melt(importance.mean[inds]), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

		## after running for the first time, COMMENT OUT THIS BLOCK ##
		# using a sparse model with N predictors
		sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
		save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", "multiclass", subtype, mlevel, st))
		load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", "multiclass", subtype, mlevel, st))
		# accuracy of final sparseRF model
		pred <- predict(sparseRF, type="prob")
		pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
		pred_df_out <- merge(pred_df, data.sel, by="row.names")
		write.table(pred_df_out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.predictions.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
		confusion_matrix <- table(pred_df[, c("true", "predicted")])
		class_errors <- unlist(lapply(group_levels, function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- group_levels
		accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
		vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
		mccvalue <- mcc(vec.pred, vec.true)
		df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", "multiclass", subtype, mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
		print(p)
		
		write.table(confusion_matrix, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.confusion_matrix.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=T)
		## END BLOCK TO COMMENT ##
		
		# plotting - per-group sparse model
		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s %s", "multiclass", subtype, mlevel, st)))
		# plotting - per-group variables
		df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
		df$metabolite_name <- as.character(df$OTU)
		if (mlevel == "BIOCHEMICAL") {
			df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
			df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
			df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
		} else if (mlevel == "SUB.PATHWAY") {
			df$subpathway <- df$metabolite_name
			df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
			df$OTU_string <- sprintf("%s", df$OTU)
		} else {
			df$superpathway <- df$metabolite_name
			df$OTU_string <- sprintf("%s", df$OTU)
		}
		df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
		# load effect sizes from linear regression
		lmres <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote="")
		tmp <- levels(df$OTU); lmres <- subset(lmres, metabolite %in% tmp)
		lmres$metabolite <- factor(lmres$metabolite, levels=levels(df$OTU))
		lims <- max(abs(lmres$estimate) + abs(lmres$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
		dircolors2 <- c("blue", "red", "black"); names(dircolors2) <- c("down", "up", "NS")
		lmres[, "importance"] <- df[match(lmres$metabolite, df$OTU), "importance"]
		p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("%s - %s %s explanatory %s", "multiclass", subtype, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
		print(p)
		p <- ggplot(lmres, aes(x=metabolite, y=estimate, color=dir, group=InfantRegimen)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=InfantRegimen), position=pd, hjust=1, color="black", size=2) + geom_tile(mapping=aes(x=metabolite, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates - %s %s explanatory %s", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
		print(p)
		p <- ggplot(lmres, aes(x=metabolite, y=estimate, color=dir, group=InfantRegimen)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=InfantRegimen), position=pd, hjust=1, color="black", size=2) + geom_tile(mapping=aes(x=metabolite, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + geom_vline(xintercept=as.numeric(df$OTU)+0.5, color="grey") + theme_classic() + ggtitle(sprintf("LM estimates - %s %s explanatory %s", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors2) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
		print(p)
		# stratified by country
		lmres.bycountry <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_RegimenCountry.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote="")
		tmp <- levels(df$OTU); lmres.bycountry <- subset(lmres.bycountry, metabolite %in% tmp)
		lmres.bycountry$metabolite <- factor(lmres.bycountry$metabolite, levels=levels(df$OTU))
		lims <- max(abs(lmres.bycountry$estimate) + abs(lmres.bycountry$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
		dircolors2 <- c("blue", "red", "black"); names(dircolors2) <- c("down", "up", "NS")
		lmres.bycountry[, "importance"] <- df[match(lmres.bycountry$metabolite, df$OTU), "importance"]
		for (i in seq(from=1, to=nlevels(lmres.bycountry$metabolite), by=10)) {
			sel <- rev(levels(lmres.bycountry$metabolite))[i:min(i+9, nlevels(lmres.bycountry$metabolite))]
			lmres.bycountry.sel <- subset(lmres.bycountry, metabolite %in% sel)
			p <- ggplot(lmres.bycountry.sel, aes(x=metabolite, y=estimate, color=dir, group=Country)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Country), position=pd, hjust=1, color="black", size=1.5) + geom_tile(mapping=aes(x=metabolite, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + geom_vline(xintercept=as.numeric(df$OTU)+0.5, color="grey") + facet_wrap(~InfantRegimen) + theme_classic() + ggtitle(sprintf("LM estimates - %s %s explanatory %s", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=dircolors2) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
			print(p)
		}
		# shading rectangles of importance values
		df.rect <- df
		df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
#		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
		p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s %s explanatory %s", "multiclass", subtype, mlevel, st)) + scale_fill_gradient(low="white", high="black")
		print(p)
		# violin plots of metabolite values
		agg.melt <- agg.melt.stored
		agg.melt$InfantGroup <- mapping.sel[agg.melt$SampleID, "InfantGroup"]
		agg.melt$InfantAgeInDays <- mapping.sel[agg.melt$SampleID, "InfantAgeInDays"]
		agg.melt$InfantAgeInDaysBinned <- mapping.sel[agg.melt$SampleID, "InfantAgeInDaysBinned"]
		agg.melt$delgage <- mapping.sel[agg.melt$SampleID, "delgage"]
		agg.melt$weight0week <- mapping.sel[agg.melt$SampleID, "weight0week"]
		agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$InfantGroup))
		agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
		agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
		p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		npages <- n_pages(p)
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantGroup)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
		}
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantAgeInDaysBinned)) + geom_violin(aes(x=InfantGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
			print(p)
		}
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=InfantAgeInDays)) + geom_violin(aes(x=InfantGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip()
			print(p)
		}
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=delgage)) + geom_violin(aes(x=InfantGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip()
			print(p)
		}
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes(x=InfantGroup, y=value, color=weight0week)) + geom_violin(aes(x=InfantGroup, y=value), inherit.aes=F) + geom_jitter(width=0.2) + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip()
			print(p)
		}
	
		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=InfantGroup)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
		agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
		p <- ggplot(agg.melt2, aes(x=InfantGroup, y=value, color=InfantGroup)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=InfantGroup)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
#		# violin plots stratified by prediction/truth
#		for (met in levels(agg.melt$metabolite)) {
#			tmp <- subset(agg.melt, metabolite==met)
#			# color scheme - manual
#			p <- ggplot(tmp, aes(x=InfantGroup, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_fill_manual(values=cols.cohort)
#			print(p)
#			p <- ggplot(tmp, aes(x=InfantGroup, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_color_manual(values=cols.cohort)
#			print(p)
#			
#		}
	}
}


### randomForest classification of StudySite; using METABOLITE data [DBS] from Malawi only
#set.seed(nrow(mapping))	
## MAIN LOOP for random forest (through metabolite levels)
#subtype <- "infant"
#for (st in c("DBS")) {
#	for (mlevel in "BIOCHEMICAL") {
#		data <- df.metabolon[[st]][[mlevel]]
#		mapping2 <- mapping[,c("Delivery", "Regimen", "Group", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom", "StudySite")]
#		colnames(mapping2) <- c("Delivery", "Regimen", "Group", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "StudySite")
#		rownames(mapping2) <- mapping2$cpatid
#		mapping.sel <- subset(mapping2, Country=="Malawi")
#		data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the desired infant samples
#		response <- droplevels(mapping.sel$StudySite); names(response) <- rownames(mapping.sel)
#		agg.melt.stored <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
#		
##		## after running for the first time, COMMENT OUT THIS BLOCK ##
##		num_iter <- 100
##		ncores <- 20
##		out <- mclapply(1:num_iter, function (dummy) {
##				importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
##			}, mc.cores=ncores )
##		collated.importance <- do.call(cbind, out)
##		out <- mclapply(1:num_iter, function (dummy) {
##				rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.9)$error.cv
##			}, mc.cores=ncores )
##		collated.cv <- do.call(cbind, out)

##		write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.importance.txt", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
##		write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.cv.txt", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
##		## END BLOCK TO COMMENT ##

#		collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.importance.txt", subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1, quote="")
#		collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.cv.txt", subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
#		importance.mean <- rowMeans(collated.importance)
#		importance.sd <- unlist(apply(collated.importance, 1, sd))
#		cv.mean <- rowMeans(collated.cv)
#		cv.sd <- unlist(apply(collated.cv, 1, sd))
#		inds <- order(importance.mean, decreasing=T)
#		inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate
#		write.table(melt(importance.mean[inds]), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.features.txt", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

#		## after running for the first time, COMMENT OUT THIS BLOCK ##
#		# using a sparse model with N predictors
##		sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
##		save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.model", subtype, mlevel, st))
#		load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.model", subtype, mlevel, st))
#		# accuracy of final sparseRF model
#		pred <- predict(sparseRF, type="prob")
#		pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
#		pred_df_out <- merge(pred_df, data.sel, by="row.names")
#		write.table(pred_df_out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.predictions.txt", subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
#		confusion_matrix <- table(pred_df[, c("true", "predicted")])
#		class_errors <- unlist(lapply(levels(mapping.sel$StudySite), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$StudySite)
#		accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
#		vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
#		mccvalue <- mcc(vec.pred, vec.true)
#		df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
#		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix StudySite (%s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", subtype, mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#		print(p)
#		
#		write.table(confusion_matrix, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE_StudySite.%s.%s.%s.confusion_matrix.txt", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=T)
#		## END BLOCK TO COMMENT ##
#		
#		# plotting - per-group sparse model
#		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
#		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
#		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - StudySite %s %s %s", subtype, mlevel, st)))
#		# plotting - per-group variables
#		df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
#		df$metabolite_name <- as.character(df$OTU)
#		if (mlevel == "BIOCHEMICAL") {
#			df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
#			df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
#			df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
#		} else if (mlevel == "SUB.PATHWAY") {
#			df$subpathway <- df$metabolite_name
#			df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
#			df$OTU_string <- sprintf("%s", df$OTU)
#		} else {
#			df$superpathway <- df$metabolite_name
#			df$OTU_string <- sprintf("%s", df$OTU)
#		}
#		df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
#		p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("StudySite - %s %s explanatory %s", subtype, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
#		print(p) 
#		# shading rectangles of importance values
#		df.rect <- df
#		df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
##		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
#		p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("StudySite - %s %s explanatory %s", subtype, mlevel, st)) + scale_fill_gradient(low="white", high="black")
#		print(p)
#		# violin plots of metabolite values
#		agg.melt <- agg.melt.stored
#		agg.melt$StudySite <- mapping.sel[agg.melt$SampleID, "StudySite"]
#		agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$StudySite))
#		agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
#		agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
#		p <- ggplot(agg.melt, aes(x=StudySite, y=value, color=StudySite)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (StudySite %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#		print(p)
#		p <- ggplot(agg.melt, aes(x=StudySite, y=value, color=StudySite)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (StudySite %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#		npages <- n_pages(p)
#		for (ip in 1:npages) {
#			p <- ggplot(agg.melt, aes(x=StudySite, y=value, color=StudySite)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#			print(p)
#		}
#		
#		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
#		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=StudySite)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (StudySite %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#		print(p)
#		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
#		agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
#		p <- ggplot(agg.melt2, aes(x=StudySite, y=value, color=StudySite)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (StudySite %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#		print(p)
#		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=StudySite)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (StudySite %s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_brewer(palette="Set1")
#		print(p)
#		
#		# t-tests
#		agg.melt.stored$StudySite <- mapping.sel[agg.melt.stored$SampleID, "StudySite"]
#		res <- {}
#		for (met in unique(agg.melt.stored$metabolite)) {
#			tmp <- subset(agg.melt.stored, metabolite==met)
#			if (length(unique(tmp$value))==1) {
#				next
#			}
#			test <- t.test(value ~ StudySite, tmp, alternative="two.sided")
#			res <- rbind(res, c(met, test$estimate, test$statistic, test$p.value, test$method))
#		}
#		res <- as.data.frame(res); colnames(res) <- c("metabolite", colnames(res)[2:3], "t", "pval", "method")
#		res$pval <- as.numeric(as.character(res$pval)); res <- res[order(res$pval),]
#		res$padj <- p.adjust(res$pval, method="fdr")
#		write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/differences_by_StudySite.t-test.%s.%s.%s.txt", subtype, mlevel, st), quote=F, row.names=F, col.names=T, sep="\t")
#		
#	}
#}

### combine StudySite t-test data
#data.maternal <- read.table("/Lab_Share/PROMISE/nwcs619/metabolon/differences_by_StudySite.t-test.maternal.BIOCHEMICAL.DBS.txt", header=T, as.is=T, sep="\t", quote="")
#data.infant <- read.table("/Lab_Share/PROMISE/nwcs619/metabolon/differences_by_StudySite.t-test.infant.BIOCHEMICAL.DBS.txt", header=T, as.is=T, sep="\t", quote="")
#data.merged <- merge(data.maternal, data.infant, by="metabolite", all=T)
#write.table(data.merged, file="/Lab_Share/PROMISE/nwcs619/metabolon/differences_by_StudySite.t-test.BIOCHEMICAL.DBS.txt", quote=F, col.names=T, row.names=F, sep="\t")

### gestational age prediction
#set.seed(nrow(mapping))	
#subtype <- "infant"
#for (st in c("DBS")) {
#	for (mlevel in "BIOCHEMICAL") {
#		data <- df.metabolon[[st]][[mlevel]]
#		mapping.sel <- mapping[,c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Infant", "hemaval.infant", "gender", "weight0week", "weight1week")]
#		colnames(mapping.sel) <- c("Delivery", "InfantRegimen", "InfantGroup", "cpatid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval", "gender", "weight0week", "weight1week")
#		rownames(mapping.sel) <- mapping.sel$cpatid
#		mapping.sel <- subset(mapping.sel, InfantRegimen %in% c("zdv", "PI-ART")) # exclude other and untreated groups since too small
#		mapping.sel <- subset(mapping.sel, delgage >= 27 & delgage <= 46 & !is.na(weight0week) & !is.na(weight1week)) # remove some crazy outliers (~24wks, 53wks), anyone with missing data for now (can consider imputation later)
#		mapping.sel$InfantGroup <- droplevels(mapping.sel$InfantGroup)
#		data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the infant samples
#		# add Country, gender, weight0week as covariates
#		for (mvar in c("Country", "gender", "weight0week")) {
#			data.sel[, mvar] <- mapping.sel[rownames(data.sel), mvar]
#		}
#		# split into training/testing sets (stratified by InfantGroup)
##		training.sel <- sample(rownames(mapping.sel), size=round(nrow(mapping.sel)*0.75))
#		training.sel <- rownames(mapping.sel)[sample.stratified(mapping.sel$InfantGroup, pct=0.75)]
#		testing.sel <- setdiff(rownames(mapping.sel), training.sel)		
#		data.training <- data.sel[training.sel,]; data.testing <- data.sel[testing.sel,]
#		response <- mapping.sel[training.sel, "delgage"]; names(response) <- training.sel
#		response.testing <- mapping.sel[testing.sel, "delgage"]; names(response.testing) <- testing.sel

#		agg.melt.stored <- melt(as.matrix(data.training[, setdiff(colnames(data.training), c("Country", "gender", "weight0week"))]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
#		
##		## after running for the first time, COMMENT OUT THIS BLOCK ##
##		num_iter <- 100
##		ncores <- 20
##		out <- mclapply(1:num_iter, function (dummy) {
##				importance(randomForest(x=data.training, y=response, ntree=10000, importance=T), type=1, scale=F)
##			}, mc.cores=ncores )
##		collated.importance <- do.call(cbind, out)
##		out <- mclapply(1:num_iter, function (dummy) {
##				rfcv(trainx=data.training, trainy=response, cv.fold=10, step=0.9)$error.cv
##			}, mc.cores=ncores )
##		collated.cv <- do.call(cbind, out)

##		write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/GestationalAgePrediction.%s.%s.%s.%s.importance.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
##		write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/GestationalAgePrediction.%s.%s.%s.%s.cv.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
##		## END BLOCK TO COMMENT ##

#		collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/GestationalAgePrediction.%s.%s.%s.%s.importance.txt", "multiclass", subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1, quote="")
#		collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/GestationalAgePrediction.%s.%s.%s.%s.cv.txt", "multiclass", subtype, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
#		importance.mean <- rowMeans(collated.importance)
#		importance.sd <- unlist(apply(collated.importance, 1, sd))
#		cv.mean <- rowMeans(collated.cv)
#		cv.sd <- unlist(apply(collated.cv, 1, sd))
#		inds <- order(importance.mean, decreasing=T)
#		inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.01)))] # edit as appropriate
#		write.table(melt(importance.mean[inds]), file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/GestationalAgePrediction.%s.%s.%s.%s.features.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

#		## after running for the first time, COMMENT OUT THIS BLOCK ##
#		# using a sparse model with N predictors
##		sparseRF <- randomForest(x=data.training[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
##		save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/GestationalAgePrediction.%s.%s.%s.%s.model", "multiclass", subtype, mlevel, st))
#		load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/GestationalAgePrediction.%s.%s.%s.%s.model", "multiclass", subtype, mlevel, st))
#		# accuracy of final sparseRF model
#		pred.training <- predict(sparseRF, type="response", newdata=data.training)
#		pred.testing <- predict(sparseRF, type="response", newdata=data.testing)
#		res <- data.frame(StudyID=c(names(pred.training), names(pred.testing)), gage.true=c(response, response.testing), gage.predicted=c(pred.training, pred.testing), Set=c(rep("Training", length(pred.training)), rep("Testing", length(pred.testing))))
#		res$InfantGroup <- mapping.sel[rownames(res), "InfantGroup"]
#		write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/GestationalAgePrediction.%s.%s.%s.%s.predictions.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
#		# predicted vs actual gestational age, RMSE, residuals
#		p <- ggplot(res, aes(x=gage.true, y=gage.predicted, color=Set)) + geom_point() + stat_smooth(method="lm") + theme_bw() + ggtitle(sprintf("Predicted vs actual gestational age (%s, %s, %s, %s)", "multiclass", subtype, mlevel, st)) + scale_color_brewer(palette="Set1")
#		print(p)
#		p <- ggplot(res, aes(x=gage.true, y=gage.predicted-gage.true, color=Set)) + geom_point() + stat_smooth(method="lm") + theme_bw() + ggtitle(sprintf("Residual plot (%s, %s, %s, %s)", "multiclass", subtype, mlevel, st)) + scale_color_brewer(palette="Set1")
#		print(p)
#		res$gweeks.true <- cut(res$gage.true, breaks=c(0, seq(from=27, to=43, by=1), 60))
#		agg <- ddply(res, .(Set, gweeks.true), function(x) {
#			c(sqrt(mean(x$gage.predicted - x$gage.true)^2), mean(x$gage.predicted - x$gage.true), sd(x$gage.predicted - x$gage.true))
#		})
#		colnames(agg) <- c("Set", "gweeks.true", "RMSE", "error.mean", "error.sd")
#		p <- ggplot(agg, aes(x=gweeks.true, y=RMSE, color=Set, group=Set)) + geom_point() + geom_line() + geom_text(aes(y=RMSE+0.5, label=round(RMSE, digits=2))) + ggtitle(sprintf("RMSE by weeks gestation (%s, %s, %s, %s)", "multiclass", subtype, mlevel, st)) + theme_bw() + scale_color_brewer(palette="Set1") + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
#		print(p)
#		p <- ggplot(agg, aes(x=gweeks.true, y=error.mean, color=Set, group=Set)) + geom_point() + geom_line() + geom_errorbar(aes(ymin=error.mean-error.sd, ymax=error.mean+error.sd), width=0.5) + ggtitle(sprintf("Mean error by weeks gestation (%s, %s, %s, %s)", "multiclass", subtype, mlevel, st)) + theme_bw() + scale_color_brewer(palette="Set1") + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
#		print(p)
#		# plotting - per-group sparse model
#		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
#		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
#		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s %s", "multiclass", subtype, mlevel, st)))
#		# plotting - per-group variables
#		df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
#		df$metabolite_name <- as.character(df$OTU)
#		if (mlevel == "BIOCHEMICAL") {
#			df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
#			df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
#			df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
#		} else if (mlevel == "SUB.PATHWAY") {
#			df$subpathway <- df$metabolite_name
#			df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
#			df$OTU_string <- sprintf("%s", df$OTU)
#		} else {
#			df$superpathway <- df$metabolite_name
#			df$OTU_string <- sprintf("%s", df$OTU)
#		}
#		df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
#		p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("%s - %s %s explanatory %s", "multiclass", subtype, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
#		print(p) 
#		# shading rectangles of importance values
#		df.rect <- df
#		df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
##		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
#		p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s %s explanatory %s", "multiclass", subtype, mlevel, st)) + scale_fill_gradient(low="white", high="black")
#		print(p)
#	}
#}



## aggregate lists of random forests features (separately for maternal and infant)
mlevel <- "BIOCHEMICAL"
subtype <- "maternal"
for (st in c("DBS", "Plasma")) {
	out <- {}
	for (regi in c("untreated", "zdv", "PI-ART")) {
		# RF
		tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", quote=""); colnames(tmp) <- c("feature", "importance")
		tmp <- cbind(sprintf("RF - %s", regi), tmp)
		colnames(tmp) <- c("Analysis", "feature", "value")
		out <- rbind(out, tmp)
		# LM
		tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote="")
		tmp <- subset(tmp, MaternalRegimen == regi & dir %in% c("up", "down"))
		if (nrow(tmp) > 0) {
			tmp <- cbind(sprintf("LM - %s", regi), tmp[, c("metabolite", "estimate")])
			colnames(tmp) <- c("Analysis", "feature", "value")
			out <- rbind(out, tmp)
		}
#		# LASSO
#		tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/LASSO.%s.%s.%s.%s.txt", regi, subtype, mlevel, st), header=T, as.is=T, sep="\t", quote=""); colnames(tmp) <- c("value", "feature")
#		tmp <- cbind(sprintf("LASSO - %s", regi), tmp[, c("feature", "value")])
#		colnames(tmp) <- c("Analysis", "feature", "value")
#		out <- rbind(out, tmp)
	}
#	# RF multiclass
#	regi <- "multiclass"
#	tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", quote=""); colnames(tmp) <- c("feature", "importance")
#	tmp <- cbind(sprintf("RF - %s", regi), tmp)
#	colnames(tmp) <- c("Analysis", "feature", "value")
#	out <- rbind(out, tmp)
	# LM averaged
	regi <- "averaged"
	tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_averaged_across_Regimen.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote="")
	tmp <- subset(tmp, dir %in% c("up", "down"))
	tmp <- cbind(sprintf("LM - %s", regi), tmp[, c("metabolite", "estimate")])
	colnames(tmp) <- c("Analysis", "feature", "value")
	out <- rbind(out, tmp)
	
	# put in useful table form
	df <- dcast(out, feature ~ Analysis)
	rownames(df) <- df$feature; df <- df[,-1]
	df[which(is.na(df), arr.ind=T)] <- 0
	df <- df != 0
	write.table(df, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/summarized_features_all_analyses.%s.%s.txt", subtype, st), row.names=T, col.names=T, quote=F, sep="\t")
}

mlevel <- "BIOCHEMICAL"
subtype <- "infant"
for (st in c("DBS")) {
	out <- {}
	for (regi in c("zdv", "PI-ART")) {
		# RF
		tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", quote=""); colnames(tmp) <- c("feature", "importance")
		tmp <- cbind(sprintf("RF - %s", regi), tmp)
		colnames(tmp) <- c("Analysis", "feature", "value")
		out <- rbind(out, tmp)
		# LM
		tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote="")
		tmp <- subset(tmp, InfantRegimen == regi & dir %in% c("up", "down"))
		if (nrow(tmp) > 0) {
			tmp <- cbind(sprintf("LM - %s", regi), tmp[, c("metabolite", "estimate")])
			colnames(tmp) <- c("Analysis", "feature", "value")
			out <- rbind(out, tmp)
		}
#		# LASSO
#		tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/LASSO.%s.%s.%s.%s.txt", regi, subtype, mlevel, st), header=T, as.is=T, sep="\t", quote=""); colnames(tmp) <- c("value", "feature")
#		tmp <- cbind(sprintf("LASSO - %s", regi), tmp[, c("feature", "value")])
#		colnames(tmp) <- c("Analysis", "feature", "value")
#		out <- rbind(out, tmp)
	}
	# RF multiclass
	regi <- "multiclass"
	tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", quote=""); colnames(tmp) <- c("feature", "importance")
	tmp <- cbind(sprintf("RF - %s", regi), tmp)
	colnames(tmp) <- c("Analysis", "feature", "value")
	out <- rbind(out, tmp)
	# LM averaged
	regi <- "averaged"
	tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_averaged_across_Regimen.%s.%s.%s.txt", subtype, st, mlevel), header=T, as.is=T, sep="\t", quote="")
	tmp <- subset(tmp, dir %in% c("up", "down"))
	if (nrow(tmp)>0) {
		tmp <- cbind(sprintf("LM - %s", regi), tmp[, c("metabolite", "estimate")])
		colnames(tmp) <- c("Analysis", "feature", "value")
		out <- rbind(out, tmp)
	}
	# put in useful table form
	df <- dcast(out, feature ~ Analysis)
	rownames(df) <- df$feature; df <- df[,-1]
	df[which(is.na(df), arr.ind=T)] <- 0
	df <- df != 0
	write.table(df, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/summarized_features_all_analyses.%s.%s.txt", subtype, st), row.names=T, col.names=T, quote=F, sep="\t")
}


#pdf(out_pdf, width=12)

### PLS-DA classification of Delivery (stratified by Regimen); maternal DBS data
### NOTE: multiclass version (4 groups) does not work well due to single ARV being selected on component 1
#set.seed(nrow(mapping))
#nrepeats <- 500
#ncomps <- 10

## no looping for PLS-DA/sPLS-DA as need to manually tune parameters
## maternal DBS BIOCHEMICAL
#subtype <- "maternal"; st <- "DBS"; mlevel <- "BIOCHEMICAL"; regi <- "zdv"
#data <- df.metabolon[[st]][[mlevel]]
#mapping.sel <- subset(mapping[,c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")], Regimen==regi)
#colnames(mapping.sel) <- c("Delivery", "Regimen", "Group", "patid", "delgage", "deldtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
#rownames(mapping.sel) <- mapping.sel$patid
#data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the maternal samples
#response <- mapping.sel$Delivery; names(response) <- rownames(mapping.sel)
#agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
## PLS-DA
#plsda.res <- plsda(data.sel, response, ncomp = ncomps) # where ncomp is the number of components wanted
#perf.plsda <- perf(plsda.res, validation = "Mfold", folds = 5, progressBar = FALSE, auc = TRUE, nrepeat = nrepeats)
#plot(perf.plsda, sd = TRUE, legend.position = "horizontal")
#plotIndiv(plsda.res, comp = 1:2, group = plsda.res$Y, ind.names = FALSE, ellipse = TRUE, legend = TRUE, title=sprintf("PLS-DA (%s, %s, %s, %s)", subtype, st, mlevel, regi))
#plotIndiv(plsda.res, comp = 1:2, group = plsda.res$Y, ind.names = FALSE, ellipse = TRUE, star = TRUE, legend = TRUE, title=sprintf("PLS-DA (%s, %s, %s, %s)", subtype, st, mlevel, regi))
#background <- background.predict(plsda.res, comp.predicted=2, dist = "max.dist") 
#plotIndiv(plsda.res, comp = 1:2, group = response, ind.names = FALSE, title = sprintf("PLS-DA (%s, %s, %s, %s)", subtype, st, mlevel, regi), legend = TRUE,  background = background)
## sPLS-DA
#list.keepX <- c(1:10,  seq(from=20, to=100, by=10))
#tune.splsda.multiclass <- tune.splsda(data.sel, response, ncomp = ncomps, validation = 'Mfold', folds = 5, progressBar = TRUE, dist = 'max.dist', measure = "BER", test.keepX = list.keepX, nrepeat = nrepeats, cpus = 16)
#error <- tune.splsda.multiclass$error.rate
#ncomp <- tune.splsda.multiclass$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
#ncomp
#select.keepX <- tune.splsda.multiclass$choice.keepX[1:ncomp]  # optimal number of variables to select
#select.keepX
#plot(tune.splsda.multiclass)
#splsda.final <- splsda(data.sel, response, ncomp = ncomp, keepX = select.keepX)
#plotIndiv(splsda.final, ind.names = FALSE, legend=TRUE, ellipse = TRUE, title=sprintf("SPLS-DA final result (%s, %s, %s, %s)", subtype, st, mlevel, regi))
#plotIndiv(splsda.final, ind.names = FALSE, legend=TRUE, ellipse = TRUE, star = TRUE, title=sprintf("SPLS-DA final result (%s, %s, %s, %s)", subtype, st, mlevel, regi))
#background <- background.predict(splsda.final, comp.predicted=2, dist = "max.dist") 
#plotIndiv(splsda.final, comp = 1:2, group = response, ind.names = FALSE, title = sprintf("SPLS-DA (%s, %s, %s, %s)", subtype, st, mlevel, regi), legend = TRUE,  background = background)
#perf.final <- perf(splsda.final, validation = "Mfold", folds = 5, dist = 'max.dist', nrepeat = nrepeats, progressBar = FALSE) 
#plot(perf.final)
## compute mean+SD for performance metrics (using ncomp components as this is model with the selected number of components)
#cmlist <- lapply(1:nrepeats, function(repi) get.confusion_matrix(truth=response, predicted=perf.final$class[["max.dist"]][,repi,ncomp]))
#cmarr <- abind(cmlist, along=3)
#cm.mean <- apply(cmarr, c(1,2), mean); cm.sd <- apply(cmarr, c(1,2), sd)
#accuracy <- 100*apply(cmarr, 3, function(x) sum(diag(x))/sum(x))
#mccvalues <- {}
#for (i in 1:nrepeats) {
#	predicted <- factor(perf.final$class[["max.dist"]][,i,ncomp][names(response)], levels=levels(response))
#	mccvalue <- mcc(as.numeric(predicted)-1, as.numeric(response)-1)
#	mccvalues <- c(mccvalues, mccvalue)
#}
#metrics <- data.frame(accuracy=accuracy, mcc=mccvalues)
#cmtab <- matrix(paste(formatC(cm.mean, format="f", digits=3), formatC(cm.sd, format="f", digits=3), sep=" +/- "), nrow=nrow(cm.mean), dimnames=dimnames(cm.mean))
#p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix mean+/-SD (%s, %s, %s, %s)", subtype, st, mlevel, regi)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(cmtab), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#print(p)
#p <- ggplot(metrics, aes(accuracy)) + geom_histogram() + theme_classic() + ggtitle(sprintf("Accuracy (%s, %s, %s, %s) mean=%.4g SD=%.4g", subtype, st, mlevel, regi, mean(metrics$accuracy), sd(metrics$accuracy))) + geom_vline(xintercept=mean(metrics$accuracy), col="red")
#print(p)
#p <- ggplot(metrics, aes(mcc)) + geom_histogram() + theme_classic() + ggtitle(sprintf("MCC (%s, %s, %s, %s) mean=%.4g SD=%.4g", subtype, st, mlevel, regi, mean(metrics$mcc), sd(metrics$mcc))) + geom_vline(xintercept=mean(metrics$mcc), col="red")
#print(p)
## arrow plots, loadings, feature stability
#plotArrow(splsda.final, legend=T)
#for (i in 1:ncomp) {
#	plotLoadings(splsda.final, comp = i, title = sprintf("Loadings on comp %d (%s, %s, %s, %s)", i, subtype, st, mlevel, regi), contrib = 'max', method = 'mean')
#}
#for (i in 1:ncomp) {
#	inds <- match(selectVar(splsda.final, comp = i)$name, names(perf.final$features$stable[[i]]))
#	freq <- as.numeric(perf.final$features$stable[[i]][inds])
#	df <- data.frame(selectVar(splsda.final, comp = i)$value, freq); df <- df[order(df$freq, decreasing=F),]; df$feature <- factor(rownames(df), levels=rownames(df)) # feature stability
#	p <- ggplot(df, aes(x=feature, y=freq)) + geom_bar(stat="identity") + coord_flip() + theme_classic() + ggtitle(sprintf("Feature stability comp %d (%s, %s, %s, %s)", i, subtype, st, mlevel, regi))
#	print(p)
#}

dev.off()





