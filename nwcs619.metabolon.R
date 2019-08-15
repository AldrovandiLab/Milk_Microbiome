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

source("/Lab_Share/fanli/code/PROMISE/utils.R")
source("/Lab_Share/fanli/code/PROMISE/mcc.R")

cols.sig <- c("black", "red", "grey"); names(cols.sig) <- c("ns", "sig", "NA")
cols.cohort <- c("#ca0020", "#f4a582", "#0571b0", "#92c5de"); names(cols.cohort) <- c("Term.zdv", "Preterm.zdv", "Term.zdvart", "Preterm.zdvart")
siglevel <- 0.05
dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")

mapping_fn <- "/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.081219.txt"
mapping <- read.table(mapping_fn, header=T, as.is=T, sep="\t")
colnames(mapping) <- c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "cpatid", "Country", "Plasma.drawdt", "DBS.drawdt", "InfantDBS.drawdt", "GestationalAgeAtCollection", "SampleID.Mom", "SampleID.Infant", "hemaval.mom", "hemaval.infant")

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
		agg <- t(log(agg))
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

out_pdf <- sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/metabolon_analysis.%s.%s.pdf", "nwcs619", format(Sys.Date(), "%m%d%y"))
pdf(out_pdf, width=12)


#########################################################################################################
### Cohort demographics and some QC data about metabolomics
mlevel <- "BIOCHEMICAL"
qc <- {}
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
	out <- data.frame(BIOCHEMICAL=rownames(tmp), detected.maternal=counts.maternal, detected.infant=counts.infant)
	write.table(out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/detected_counts.%s.txt", st), quote=F, sep="\t", row.names=F, col.names=T)
}
	
#########################################################################################################
### maternal metabolites (DBS, Plasma separately)
## ordination (t-SNE, PCA) and PERMANOVA
set.seed(nrow(mapping))
subtype <- "maternal"
for (st in c("DBS", "Plasma")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mapping[,c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
	colnames(mapping.sel) <- c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
	rownames(mapping.sel) <- mapping.sel$patid
	data <- data[rownames(mapping.sel),] # subset to just the maternal samples
	
	# PCA
	pca <- prcomp(data, center=F, scale=F)
	eigs <- pca$sdev^2
	pvar <- 100*(eigs / sum(eigs))
	df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], SampleID=rownames(pca$x))
	for (mvar in rownames(subset(metadata_variables, useForPERMANOVA=="yes"))) {
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
	# PERMANOVA
	res <- adonis2(data ~ Delivery + Regimen + Country + GestationalAgeAtCollection, data=mapping.sel, permutations=999, method='euclidean')
	sink(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/metabolon_PERMANOVA.%s.%s.%s.txt", subtype, st, mlevel))
	print(res)
	sink()
	# t-SNE
	tsne.out <- Rtsne(data, perplexity=20)
	df <- data.frame(PC1 = tsne.out$Y[,1], PC2 = tsne.out$Y[,2], SampleID=rownames(data)); rownames(df) <- df$SampleID
	for (mvar in rownames(subset(metadata_variables, useForPERMANOVA=="yes"))) {
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

## linear model with emmeans, stratified by Regimen
subtype <- "maternal"; mvar <- "Delivery"
for (st in c("DBS", "Plasma")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mapping[,c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
	colnames(mapping.sel) <- c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
	rownames(mapping.sel) <- mapping.sel$patid
	data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
	df <- merge(data.sel, mapping.sel, by="row.names"); 
	df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
	res <- {}
	for (metabolite in colnames(data.sel)) {
		mod <- lm(as.formula(sprintf("%s ~ %s*Regimen", metabolite, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s | Regimen", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
		res <- rbind(res, tmp)
	}
	res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
	res$padj <- p.adjust(res$p.value, method="fdr")
	res <- res[order(res$p.value, decreasing=F),]
	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_by_Regimen.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
	# forest plot
	resSig <- subset(res, padj < siglevel)
	df <- subset(res, metabolite %in% as.character(resSig$metabolite))
	df <- df[order(df$estimate, decreasing=T),]
	df$metabolite <- factor(as.character(df$metabolite), levels=as.character(unique(df$metabolite)))
	lims <- max(abs(df$estimate) + abs(df$SE))*1.0
	pd <- position_dodge(0.8)
	p <- ggplot(df, aes(x=metabolite, y=estimate, color=dir, group=Regimen)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_text(aes(y=-lims, label=Regimen), position=pd, hjust=1, color="black", size=2) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("Delivery stratified by Regimen (%s, %s, %s)", subtype, st, mlevel)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
	print(p)
}

## linear model with emmeans, averaged across Regimen
subtype <- "maternal"; mvar <- "Delivery"
for (st in c("DBS", "Plasma")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- mapping[,c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
	colnames(mapping.sel) <- c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
	rownames(mapping.sel) <- mapping.sel$patid
	data.sel <- data[rownames(mapping.sel),] # subset to just the maternal samples
	name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
	sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds < quantile(sds, probs=0.05))); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites in the lowest 5% for variance
	df <- merge(data.sel, mapping.sel, by="row.names"); 
	df[, mvar] <- factor(as.character(df[,mvar]), levels=rev(levels(df[,mvar]))) # reverse levels to get more intuitive contrast
	res <- {}
	for (metabolite in colnames(data.sel)) {
		mod <- lm(as.formula(sprintf("%s ~ %s+Regimen", metabolite, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("pairwise ~ %s", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
		res <- rbind(res, tmp)
	}
	res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
	res$padj <- p.adjust(res$p.value, method="fdr")
	res <- res[order(res$p.value, decreasing=F),]
	res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/emmeans_averaged_across_Regimen.%s.%s.%s.txt", subtype, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
	# forest plot
	resSig <- subset(res, padj < siglevel)
	df <- subset(res, metabolite %in% as.character(resSig$metabolite))
	df <- df[order(df$estimate, decreasing=T),]
	df$metabolite <- factor(as.character(df$metabolite), levels=as.character(unique(df$metabolite)))
	lims <- max(abs(df$estimate) + abs(df$SE))*1.0
	pd <- position_dodge(0.8)
	p <- ggplot(df, aes(x=metabolite, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=metabolite, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("Delivery averaged across Regimen (%s, %s, %s)", subtype, st, mlevel)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
	print(p)
}


## randomForest classification of Group (separately by Regimen); using METABOLITE data [plasma, fecal]
set.seed(nrow(mapping))	
# MAIN LOOP for random forest (through metabolite levels)
subtype <- "maternal"
for (st in c("DBS", "Plasma")) {
	for (mlevel in "BIOCHEMICAL") {
		data <- df.metabolon[[st]][[mlevel]]
		mapping2 <- mapping[,c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
		colnames(mapping2) <- c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
		rownames(mapping2) <- mapping2$patid
		for (regi in levels(mapping2$Regimen)) {
			mapping.sel <- subset(mapping2, Regimen==regi); mapping.sel$Group <- droplevels(mapping.sel$Group)
			data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the desired maternal samples from regimen
			response <- droplevels(mapping.sel$Group); names(response) <- rownames(mapping.sel)
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
			inds <- inds[1:min(50, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate
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
			class_errors <- unlist(lapply(levels(mapping.sel$Group), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$Group)
			accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
			vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
			mccvalue <- mcc(vec.pred, vec.true)
			df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
			p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", regi, subtype, mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
			print(p)
			
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
			p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("%s - %s %s explanatory %s", regi, subtype, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
			print(p) 
			# shading rectangles of importance values
			df.rect <- df
			df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
	#		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
			p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s %s explanatory %s", regi, subtype, mlevel, st)) + scale_fill_gradient(low="white", high="black")
			print(p)
			# violin plots of metabolite values
			agg.melt <- agg.melt.stored
			agg.melt$Group <- mapping.sel[agg.melt$SampleID, "Group"]
			agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$Group))
			agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
			agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
			p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			npages <- n_pages(p)
			for (ip in 1:npages) {
				p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
			}
			
			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Group)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
			agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
			p <- ggplot(agg.melt2, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Group)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			# violin plots stratified by prediction/truth
			for (met in levels(agg.melt$metabolite)) {
				tmp <- subset(agg.melt, metabolite==met)
				# color scheme - manual
				p <- ggplot(tmp, aes(x=Group, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_fill_manual(values=cols.cohort)
				print(p)
				p <- ggplot(tmp, aes(x=Group, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_color_manual(values=cols.cohort)
				print(p)
			}
		}
	}
}


## randomForest classification of Group (multiclass); using METABOLITE data [plasma, fecal]
set.seed(nrow(mapping))	
# MAIN LOOP for random forest (through metabolite levels)
subtype <- "maternal"
for (st in c("DBS", "Plasma")) {
	for (mlevel in "BIOCHEMICAL") {
		data <- df.metabolon[[st]][[mlevel]]
		mapping.sel <- mapping[,c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
		colnames(mapping.sel) <- c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
		rownames(mapping.sel) <- mapping.sel$patid
		data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the maternal samples
		response <- mapping.sel$Group; names(response) <- rownames(mapping.sel)
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

#		write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.importance.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
#		write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.cv.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
#		## END BLOCK TO COMMENT ##

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
#		sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
#		save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", "multiclass", subtype, mlevel, st))
		load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.model", "multiclass", subtype, mlevel, st))
		# accuracy of final sparseRF model
		pred <- predict(sparseRF, type="prob")
		pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
		pred_df_out <- merge(pred_df, data.sel, by="row.names")
		write.table(pred_df_out, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s.%s.%s.%s.predictions.txt", "multiclass", subtype, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
		confusion_matrix <- table(pred_df[, c("true", "predicted")])
		class_errors <- unlist(lapply(levels(mapping.sel$Group), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$Group)
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
		p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("%s - %s %s explanatory %s", "multiclass", subtype, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
		print(p) 
		# shading rectangles of importance values
		df.rect <- df
		df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
#		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
		p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s %s explanatory %s", "multiclass", subtype, mlevel, st)) + scale_fill_gradient(low="white", high="black")
		print(p)
		# violin plots of metabolite values
		agg.melt <- agg.melt.stored
		agg.melt$Group <- mapping.sel[agg.melt$SampleID, "Group"]
		agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$Group))
		agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
		agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
		p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		npages <- n_pages(p)
		for (ip in 1:npages) {
			p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
		}
		
		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Group)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
		agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
		p <- ggplot(agg.melt2, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Group)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", subtype, mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		# violin plots stratified by prediction/truth
		for (met in levels(agg.melt$metabolite)) {
			tmp <- subset(agg.melt, metabolite==met)
			# color scheme - manual
			p <- ggplot(tmp, aes(x=Group, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_fill_manual(values=cols.cohort)
			print(p)
			p <- ggplot(tmp, aes(x=Group, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s, %s)", met, subtype, mlevel, st)) + scale_color_manual(values=cols.cohort)
			print(p)
			
		}
	}
}




#pdf(out_pdf, width=12)

### PLS-DA classification of Group (multiclass); using maternal DBS/plasma data
#set.seed(nrow(mapping))
#nrepeats <- 500
#ncomps <- 5

## no looping for PLS-DA/sPLS-DA as need to manually tune parameters
## maternal DBS BIOCHEMICAL
#subtype <- "maternal"; st <- "DBS"; mlevel <- "BIOCHEMICAL"
#data <- df.metabolon[[st]][[mlevel]]
#mapping.sel <- mapping[,c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "Country", "GestationalAgeAtCollection", "SampleID.Mom", "hemaval.mom")]
#colnames(mapping.sel) <- c("Delivery", "Regimen", "Group", "patid", "del.gage", "del.dtup", "Country", "GestationalAgeAtCollection", "SampleID", "hemaval")
#rownames(mapping.sel) <- mapping.sel$patid
#data.sel <- as.data.frame(data[rownames(mapping.sel),]) # subset to just the maternal samples
#response <- mapping.sel$Group; names(response) <- rownames(mapping.sel)
#agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Country")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
## PLS-DA
#plsda.res <- plsda(data.sel, response, ncomp = ncomps) # where ncomp is the number of components wanted
#perf.plsda <- perf(plsda.res, validation = "Mfold", folds = 5, progressBar = FALSE, auc = TRUE, nrepeat = nrepeats)
#plot(perf.plsda, sd = TRUE, legend.position = "horizontal")
#plotIndiv(plsda.res, comp = 1:2, group = plsda.res$Y, ind.names = FALSE, ellipse = TRUE, legend = TRUE, title=sprintf("PLS-DA (%s, %s, %s)", subtype, st, mlevel))
#plotIndiv(plsda.res, comp = 1:2, group = plsda.res$Y, ind.names = FALSE, ellipse = TRUE, star = TRUE, legend = TRUE, title=sprintf("PLS-DA (%s, %s, %s)", subtype, st, mlevel))
#background <- background.predict(plsda.res, comp.predicted=2, dist = "max.dist") 
#plotIndiv(plsda.res, comp = 1:2, group = response, ind.names = FALSE, title = sprintf("PLS-DA (%s, %s)", st, mlevel), legend = TRUE,  background = background)
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
#plotIndiv(splsda.final, ind.names = FALSE, legend=TRUE, ellipse = TRUE, title=sprintf("SPLS-DA final result (%s, %s)", st, mlevel))
#plotIndiv(splsda.final, ind.names = FALSE, legend=TRUE, ellipse = TRUE, star = TRUE, title=sprintf("SPLS-DA final result (%s, %s)", st, mlevel))
#background <- background.predict(splsda.final, comp.predicted=2, dist = "max.dist") 
#plotIndiv(splsda.final, comp = 1:2, group = response, ind.names = FALSE, title = sprintf("SPLS-DA (%s, %s)", st, mlevel), legend = TRUE,  background = background)
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
#p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix mean+/-SD (%s, %s)", st, mlevel)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(cmtab), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#print(p)
#p <- ggplot(metrics, aes(accuracy)) + geom_histogram() + theme_classic() + ggtitle(sprintf("Accuracy (%s, %s) mean=%.4g SD=%.4g", st, mlevel, mean(metrics$accuracy), sd(metrics$accuracy))) + geom_vline(xintercept=mean(metrics$accuracy), col="red")
#print(p)
#p <- ggplot(metrics, aes(mcc)) + geom_histogram() + theme_classic() + ggtitle(sprintf("MCC (%s, %s) mean=%.4g SD=%.4g", st, mlevel, mean(metrics$mcc), sd(metrics$mcc))) + geom_vline(xintercept=mean(metrics$mcc), col="red")
#print(p)
## arrow plots, loadings, feature stability
#plotArrow(splsda.final, legend=T)
#for (i in 1:ncomp) {
#	plotLoadings(splsda.final, comp = i, title = sprintf("Loadings on comp %d (%s, %s)", i, st, mlevel), contrib = 'max', method = 'mean')
#}
#for (i in 1:ncomp) {
#	inds <- match(selectVar(splsda.final, comp = i)$name, names(perf.final$features$stable[[i]]))
#	freq <- as.numeric(perf.final$features$stable[[i]][inds])
#	df <- data.frame(selectVar(splsda.final, comp = i)$value, freq); df <- df[order(df$freq, decreasing=F),]; df$feature <- factor(rownames(df), levels=rownames(df)) # feature stability
#	p <- ggplot(df, aes(x=feature, y=freq)) + geom_bar(stat="identity") + coord_flip() + theme_classic() + ggtitle(sprintf("Feature stability (%s, %s)", st, mlevel))
#	print(p)
#}

dev.off()


###############################################################################################################
#### METABOLITES (plasma and rectal combined)

#pdf(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/phyloseq_output.%s.%s.pdf", "Kathy_NAFLD_merged_metabolites", format(Sys.Date(), "%m%d%y")), width=12)
#cols.sig <- c("black", "red", "grey"); names(cols.sig) <- c("ns", "sig", "NA")

### randomForest classification of Group (one-vs-one setup, X vs NMLBMI baseline); using METABOLITE data [plasma, fecal]
#set.seed(122517)	
## MAIN LOOP for random forest (through metabolite levels)
#st <- "merged"
#metabolite_sources <- c("plasma", "fecal")
#for (mlevel in metabolite_levels) {
#	data <- {}
#	sel <- intersect(rownames(df.metabolon[["plasma"]][[mlevel]]), rownames(df.metabolon[["fecal"]][[mlevel]]))
#	for (mst in metabolite_sources) {
#		tmp <- df.metabolon[[mst]][[mlevel]][sel,]
#		colnames(tmp) <- sprintf("[%s] %s", mst, colnames(tmp))
#		data <- cbind(data, tmp)		
#	}
#	mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% rownames(data))
#	# random forest classification
#	res.mean <- {}; res.sd <- {}; res.sig <- list()
#	for (mvar_level in setdiff(levels(mapping.sel[, "Group"]), "NMLBMI")) {
#		sel <- intersect(rownames(data), unique(subset(sample_data(ps), Group %in% c(mvar_level, "NMLBMI"))$StudyID))
#		data.sel <- as.data.frame(data[sel,])
#		response <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Group"]); names(response) <- sel
#		# add BMI, Age, Sex as covariates
#		data.sel$BMI <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "BMI"])
#		data.sel$Age <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Age"])
#		data.sel$Sex <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Sex"])
#		agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Sex")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
#		
##		## after running for the first time, COMMENT OUT THIS BLOCK ##
##		num_iter <- 100
##		ncores <- 20
##		out <- mclapply(1:num_iter, function (dummy) {
##				importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
##			}, mc.cores=ncores )
##		collated.importance <- do.call(cbind, out)
##		out <- mclapply(1:num_iter, function (dummy) {
##				rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
##			}, mc.cores=ncores )
##		collated.cv <- do.call(cbind, out)

##		write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.importance.txt", mvar_level, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
##		write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.cv.txt", mvar_level, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
##		## END BLOCK TO COMMENT ##

#		collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.importance.txt", mvar_level, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
#		collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.cv.txt", mvar_level, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
#		importance.mean <- rowMeans(collated.importance)
#		importance.sd <- unlist(apply(collated.importance, 1, sd))
#		cv.mean <- rowMeans(collated.cv)
#		cv.sd <- unlist(apply(collated.cv, 1, sd))
#		res.mean <- rbind(res.mean, importance.mean)
#		res.sd <- rbind(res.sd, importance.sd)
#		inds <- order(importance.mean, decreasing=T)
#		inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.001)))] # edit as appropriate
#		res.sig[[length(res.sig)+1]] <- names(importance.mean[inds])

#		## after running for the first time, COMMENT OUT THIS BLOCK ##
#		# using a sparse model with N predictors
##		sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
##		save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.model", mvar_level, mlevel, st))
#		load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.model", mvar_level, mlevel, st))
#		# ROC and AUC of final sparseRF model
#		pred <- predict(sparseRF, type="prob")
#		pred2 <- prediction(pred[,2], ordered(response))
#		perf <- performance(pred2, "tpr", "fpr")
#		perf.auc <- performance(pred2, "auc")
#		pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
#		confusion_matrix <- table(pred_df[, c("true", "predicted")])
#		accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
#		vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
#		mccvalue <- mcc(vec.pred, vec.true)
#		plot(perf, main=sprintf("ROC %s %s %s (sparseRF final model)", mvar_level, mlevel, st)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g, accuracy=%.2g%%, MCC=%.4g", unlist(perf.auc@y.values), accuracy, mccvalue))			
#		## END BLOCK TO COMMENT ##
#		
#		# plotting - per-group sparse model
#		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
#		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
#		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s", mvar_level, mlevel, st)))
#		# plotting - per-group variables
#		df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
#		df$metabolite_name <- gsub("[plasma] ", "", gsub("[fecal] ", "", as.character(df$OTU), fixed=T), fixed=T)
#		df$metabolite_source <- unlist(lapply(as.character(df$OTU), function(x) gsub("]", "", gsub("[", "", unlist(strsplit(x, " "))[1], fixed=T), fixed=T)))
#		res <- {}
#		for (mst in metabolite_sources) {
#			tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/METABOLITE_linear_regression.%s.%s.txt", mst, mlevel), header=T, sep="\t", as.is=T, quote="")
#			tmp <- subset(tmp, coefficient == paste("Group", mvar_level, sep="")); rownames(tmp) <- sprintf("[%s] %s", mst, tmp$metabolite)
#			res <- rbind(res, tmp)
#		}
#		df$Z <- res[as.character(df$OTU), "Z"]; df$padj <- res[as.character(df$OTU), "padj"]
#		if (mlevel == "BIOCHEMICAL") {
#			df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
#			df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
#			df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g; %s; %s)", df$OTU, df$Z, df$padj, df$subpathway, df$superpathway)
#		} else if (mlevel == "SUB.PATHWAY") {
#			df$subpathway <- df$metabolite_name
#			df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
#			df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g)", df$OTU, df$Z, df$padj)
#		} else {
#			df$superpathway <- df$metabolite_name
#			df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g)", df$OTU, df$Z, df$padj)
#		}
#		df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
#		df$sig <- ifelse(df$padj < 0.1, "sig", "ns"); df$sig[which(is.na(df$sig))] <- "ns"
#		p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string, color=sig), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory %s", mvar_level, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
#		print(p)
#		# violin plots of relative abundance
#		agg.melt <- agg.melt.stored
#		agg.melt$Group <- mapping.sel[match(agg.melt$SampleID, mapping.sel$StudyID), "Group"]
#		agg.melt <- subset(agg.melt, metabolite %in% rownames(df) & Group %in% c(mvar_level, "NMLBMI"))
#		agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
#		p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mvar_level, st)) + coord_flip()
#		print(p)			
#	}
#	# plotting - all groups importance values
#	rownames(res.mean) <- setdiff(levels(mapping.sel[, "Group"]), "NMLBMI"); rownames(res.sd) <- rownames(res.mean)
#	df <- res.mean[, unique(unlist(res.sig))]; df[which(df<0.001, arr.ind=T)] <- 0
#	df.sig <- matrix("", nrow=nrow(df), ncol=ncol(df)); rownames(df.sig) <- rownames(df); colnames(df.sig) <- colnames(df)
#	for (i in 1:length(res.sig)) {
#		df.sig[i, res.sig[[i]]] <- "*"
#	}
#	heatmap.2(df, Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(10,10), cexCol=0.6, cexRow=0.6, cellnote=df.sig, notecex=0.8, notecol="red", main=sprintf("RF importance values (%s, %s)", st, mlevel))
#}


### randomForest classification of Group (multiclass); using METABOLITE data [plasma, fecal]
#set.seed(122617)	
## MAIN LOOP for random forest (through metabolite levels)
#st <- "merged"
#for (mlevel in metabolite_levels) {
#	data <- {}
#	sel <- intersect(rownames(df.metabolon[["plasma"]][[mlevel]]), rownames(df.metabolon[["fecal"]][[mlevel]]))
#	for (mst in c("plasma", "fecal")) {
#		tmp <- df.metabolon[[mst]][[mlevel]][sel,]
#		colnames(tmp) <- sprintf("[%s] %s", mst, colnames(tmp))
#		data <- cbind(data, tmp)		
#	}
#	mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% rownames(data))
#	data.sel <- as.data.frame(data)
#	response <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Group"]); names(response) <- rownames(data.sel)
#	# add BMI, Age, Sex as covariates
#	data.sel$BMI <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "BMI"])
#	data.sel$Age <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Age"])
#	data.sel$Sex <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Sex"])
#	agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Sex")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
#	
##	## after running for the first time, COMMENT OUT THIS BLOCK ##
##	num_iter <- 100
##	ncores <- 20
##	out <- mclapply(1:num_iter, function (dummy) {
##			importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
##		}, mc.cores=ncores )
##	collated.importance <- do.call(cbind, out)
##	out <- mclapply(1:num_iter, function (dummy) {
##			rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
##		}, mc.cores=ncores )
##	collated.cv <- do.call(cbind, out)

##	write.table(collated.importance, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.importance.txt", "multiclass", mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
##	write.table(collated.cv, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.cv.txt", "multiclass", mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
##	## END BLOCK TO COMMENT ##

#	collated.importance <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.importance.txt", "multiclass", mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
#	collated.cv <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.cv.txt", "multiclass", mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
#	importance.mean <- rowMeans(collated.importance)
#	importance.sd <- unlist(apply(collated.importance, 1, sd))
#	cv.mean <- rowMeans(collated.cv)
#	cv.sd <- unlist(apply(collated.cv, 1, sd))
#	inds <- order(importance.mean, decreasing=T)
#	inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.0005)))] # edit as appropriate

#	## after running for the first time, COMMENT OUT THIS BLOCK ##
#	# using a sparse model with N predictors
##	sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
##	save(sparseRF, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.model", "multiclass", mlevel, st))
#	load(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.model", "multiclass", mlevel, st))
#	# accuracy of final sparseRF model
#	pred <- predict(sparseRF, type="prob")
#	pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
#	confusion_matrix <- table(pred_df[, c("true", "predicted")])
#	class_errors <- unlist(lapply(levels(mapping.sel$Group), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$Group)
#	accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
#	vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
#	mccvalue <- mcc(vec.pred, vec.true)
#	df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
#	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", "multiclass", mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#	print(p)
#	
#	write.table(confusion_matrix, file=sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/randomForest_METABOLITE.%s_%s.%s.confusion_matrix.txt", "multiclass", mlevel, st), quote=F, sep="\t", row.names=T, col.names=T)
#	## END BLOCK TO COMMENT ##
#	
#	# plotting - per-group sparse model
#	df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
#	colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
#	print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s", "multiclass", mlevel, st)))
#	# plotting - per-group variables
#	df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
#	df$metabolite_name <- gsub("[plasma] ", "", gsub("[fecal] ", "", as.character(df$OTU), fixed=T), fixed=T)
#	df$metabolite_source <- unlist(lapply(as.character(df$OTU), function(x) gsub("]", "", gsub("[", "", unlist(strsplit(x, " "))[1], fixed=T), fixed=T)))
#	res <- {}
#	for (mst in metabolite_sources) {
#		tmp <- read.table(sprintf("/Lab_Share/PROMISE/nwcs619/metabolon/METABOLITE_linear_regression.%s.%s.txt", mst, mlevel), header=T, sep="\t", as.is=T, quote="")
#		tmp <- tmp[grep("Group", tmp$coefficient),]
#		tmp <- aggregate(padj ~ metabolite, tmp, min); rownames(tmp) <- sprintf("[%s] %s", mst, tmp$metabolite)
#		res <- rbind(res, tmp)
#	}
#	df$padj <- res[as.character(df$OTU), "padj"]
#	df$sig <- ifelse(df$padj < 0.1, "sig", "ns"); df$sig[which(is.na(df$sig))] <- "ns"
#	if (mlevel == "BIOCHEMICAL") {
#		df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
#		df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
#		df$OTU_string <- sprintf("%s (min padj=%.4g; %s; %s)", df$OTU, df$padj, df$subpathway, df$superpathway)
#	} else if (mlevel == "SUB.PATHWAY") {
#		df$subpathway <- df$metabolite_name
#		df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
#		df$OTU_string <- sprintf("%s (min padj=%.4g)", df$OTU, df$padj)
#	} else {
#		df$superpathway <- df$metabolite_name
#		df$OTU_string <- sprintf("%s (min padj=%.4g)", df$OTU, df$padj)
#	}
#	df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
#	p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU, color=sig), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory %s", "multiclass", mlevel, st)) + scale_fill_manual(values=cols.superpathway) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
#	print(p)
#	# violin plots of relative abundance
#	agg.melt <- agg.melt.stored
#	agg.melt$Group <- mapping.sel[match(agg.melt$SampleID, mapping.sel$StudyID), "Group"]
#	agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$Group))
#	agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
#	agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
#	p <- ggplot(agg.melt, aes(x=Group, y=value, color=Group)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mlevel, st)) + coord_flip()
#	print(p)
#	# violin plots stratified by prediction/truth
#	for (met in levels(agg.melt$metabolite)) {
#		tmp <- subset(agg.melt, metabolite==met)
#		p <- ggplot(tmp, aes(x=Group, y=value, fill=Prediction, color=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s)", met, mlevel, st))
#		print(p)
#	}
#}


#dev.off()





