
library(ggplot2)
library(reshape2)
library(plyr)
library(useful)
library(lubridate)
library(RColorBrewer)

## subject counts
#ddply(subset(data, Aim1=="Yes"), .(Transmitter), function(x) {
#	length(unique(x$Patient.ID))
#})
#ddply(subset(data, Aim2=="Yes"), .(Regimen), function(x) {
#	length(unique(x$Patient.ID))
#})

nbatches <- 40
samples_per_batch <- 36

## perform randomization for BMK_1410 samples for Metabolon
data <- read.table("/Lab_Share/PROMISE/nwcs610/BMK_1410_Master_Sample_Sheet.092022.txt", header=T, as.is=T, sep="\t", quote="")
set.seed(nrow(data))
df <- ddply(data, .(Patient.ID), function(x) {
	if (unique(x$Study)=="NWCS610") {
		x$Group <- ifelse(any(x$Aim1=="Yes"), sprintf("Aim1-%s", x$Transmitter), sprintf("Aim2-%s", x$Regimen))
	} else if (unique(x$Study=="ZEBS")) {
		x$Group <- sprintf("ZEBS-%s", x$Transmitter)
	} else if (unique(x$Study=="Haiti")) {
		x$Group <- sprintf("Haiti-%s", x$HIVStatus)
	} else {
		x$Group <- "Unknown"
	}
	x
})

# create subject-level table
res <- unique(df[,c("Patient.ID", "Group", "instn")])
rownames(res) <- res$Patient.ID
res$Batch <- NA

# assign Aim1/ZEBS (one Transmitter and two Control to each of 40 batches, extra 9 subjects to random batches)
# use matching from SELECTED.Aim1.030419.txt for NWCS610-Aim1, matched_1to2.Aim1.010821.txt for ZEBS
matching <- read.table("/Lab_Share/PROMISE/nwcs610/SELECTED.Aim1.030419.txt", header=T, as.is=T, sep="\t", quote="")
rownames(matching) <- matching$mpatid
i <- 1
for (pid in as.character(subset(res, Group=="Aim1-Transmitter")$Patient.ID)) {
	res[pid, "Batch"] <- i
	res[as.character(matching[pid, "Match1"]), "Batch"] <- i
	res[as.character(matching[pid, "Match2"]), "Batch"] <- i
	i <- i+1
}
matching <- read.table("/Lab_Share/ZEBS/matched_1to2.Aim1.010821.txt", header=T, as.is=T, sep="\t", quote="")
rownames(matching) <- matching$studyid
for (pid in as.character(subset(res, Group=="ZEBS-Transmitter")$Patient.ID)) {
	res[pid, "Batch"] <- i
	for (m in rownames(subset(matching, match==pid))) {
		res[m, "Batch"] <- i
	}
	i <- i+1
}
# reassign n=9 subjects in batches 41-43 to 1-40 that have more slots open
pids.reassign <- rownames(subset(res, Batch>40))
tmp <- table(res$Batch)
fixed <- sample(as.numeric(names(which(tmp==2))))
float <- sample(intersect(1:nbatches, as.numeric(names(which(tmp==3)))), length(pids.reassign)-length(fixed))
res[pids.reassign, "Batch"] <- c(fixed, float)

# assign Aim2 (one NVP and one maternal triple ART to each of 22 batches, add leftovers to their appropriate Aim1 matches)
matching <- read.table("/Lab_Share/PROMISE/nwcs610/SELECTED.Aim2.030519.txt", header=T, as.is=T, sep="\t", quote="")
id1 <- as.character(matching$id1); id2 <- as.character(matching$id2)
inds <- which(!is.na(res[id1, "Batch"]))
used <- unique(res[id1[inds], "Batch"])
res[id2[inds], "Batch"] <- res[id1[inds], "Batch"]
inds <- which(is.na(res[id1, "Batch"]))
slots <- setdiff(as.numeric(names(which(table(res$Batch)<4))), used)
for (i in inds) {
	x <- sample(slots, 1)
	res[id1[i], "Batch"] <- x
	res[id2[i], "Batch"] <- x
	slots <- setdiff(slots, x)
}

# assign McGrath HEU/HUU
matching <- read.table("/Lab_Share/PROMISE/nwcs610/TunzaMwanaStudy_Del_6wk.10132022.txt", header=T, as.is=T, sep="\t", quote="")
rownames(matching) <- matching$PID
matching$Group <- sprintf("McGrath-%s", matching$HIV.Status)
matching$Delivery.Visit.Date <- as.Date(matching$Delivery.Visit.Date, format="%m/%d/%Y")
matching$X6.Week.Visit.Date <- as.Date(matching$X6.Week.Visit.Date, format="%m/%d/%Y")
matching <- subset(matching, Study.Status=="Active")
matching <- subset(matching, !(is.na(Delivery.Visit.Date) & is.na(X6.Week.Visit.Date)))

# calculate open slots (6x HEU/HUU for batches with 3 assigned, 5x for batches with 4 assigned, 3x+1 for batches with 5 assigned)
tab <- table(res$Batch)
slots <- c(rep(which(tab==3), each=6), rep(which(tab==4), each=5), rep(which(tab==5), each=3))
names(slots) <- {}
sel <- rownames(subset(matching, HIV.Status=="Negative"))
matching[sel, "Batch"] <- sample(slots, size=length(sel))
sel <- rownames(subset(matching, HIV.Status=="Positive"))
matching[sample(sel, size=min(length(sel), length(slots))), "Batch"] <- sample(slots, size=min(length(sel), length(slots)))
sel <- rownames(matching)[which(is.na(matching$Batch))]
mc <- table(matching$Batch)
assign <- c(which(mc<6), sample(which(mc==6), size=length(sel)-length(which(mc<6))))
matching[sel, "Batch"] <- assign
mcgrath <- matching[, c("PID", "Group", "Study.Site", "Batch")]
colnames(mcgrath) <- colnames(res)

# convert to samples
out <- data
out$Batch <- res[as.character(out$Patient.ID), "Batch"]

# add McGrath samples
vars <- c("HIV.Status", "Study.Site", "Ethnicity", "Date.of.Birth", "Enrolment.Height", "Enrolment.Weight")
tmp <- melt(matching[,c("PID", "Delivery.Visit.Date", "X6.Week.Visit.Date")], id.vars=c("PID"))
colnames(tmp) <- c("Patient.ID", "Timepoint2", "TimePoint.Date")
tmp$SampleID <- sprintf("%s.%s", tmp$Patient.ID, format(tmp$TimePoint.Date, "%Y%m%d"))
tmp$Timepoint2 <- ifelse(tmp$Timepoint2=="Delivery.Visit.Date", "Delivery", "Week06")
tmp$SampleType <- "BMK"
tmp$HIVStatus <- matching[tmp$Patient.ID, "HIV.Status"]
tmp$Study <- "McGrath"
tmp$Arm <- mcgrath[tmp$Patient.ID, "Group"]
tmp$Batch <- mcgrath[tmp$Patient.ID, "Batch"]
tmp$instn <- mcgrath[tmp$Patient.ID, "instn"]
tmp <- subset(tmp, !is.na(TimePoint.Date))

out$TimePoint.Date <- as.character(out$TimePoint.Date)
tmp$TimePoint.Date <- as.character(tmp$TimePoint.Date, format="%Y%m%d")
final <- merge(out, tmp, all=T)
rownames(final) <- final$SampleID

# add Haiti to batches with most open slots
tab <- table(final$Batch)
slots <- which(tab<=32); needed <- max(table(subset(final, Study=="Haiti")$HIVStatus))
assign <- c(sample(slots, size=min(length(slots), needed), replace=F), sample(slots, size=needed-min(length(slots),needed), replace=F))
sel <- rownames(subset(final, Study=="Haiti" & HIVStatus=="Positive"))
final[sel, "Batch"] <- assign[1:length(sel)]
sel <- rownames(subset(final, Study=="Haiti" & HIVStatus=="Negative"))
final[sel, "Batch"] <- assign[1:length(sel)]


## fix 3 batches that have 37 samples
tab <- table(final$Batch)
inds <- which(tab>36) # Batch 7 has 37 samples
reassign <- which(tab==32) # Batch 13 only has 32 samples
# move PID 8-01-002 and 8-01-024 from Batch 13 to Batch 7
sel <- c("8-01-002", "8-01-024")
final[rownames(subset(final, Patient.ID %in% sel)), "Batch"] <- reassign

write.table(final, file="/Lab_Share/PROMISE/nwcs610/BMK_1410_Master_Sample_Sheet.with_batching.101722.txt", quote=F, sep="\t", row.names=F, col.names=T)

################################################################################################################
## check randomization
final$instn2 <- final$instn
final[which(final$instn %in% names(which(table(final$instn) <= 12))), "instn2"] <- "Other CRS"
color_list <- list()
colors <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"), brewer.pal(8, "Pastel1"), brewer.pal(8, "Pastel2"))
cols.instn2 <- colors[1:length(unique(final$instn2))]; names(cols.instn2) <- unique(final$instn2)
color_list[["instn2"]] <- cols.instn2

out_pdf <- sprintf("/Lab_Share/PROMISE/nwcs610/BMK_1410_randomization.pdf")
pdf(out_pdf, width=12)

# distribution of Study/Arm across batches
mvars <- c("Study", "Arm")
for (mvar in mvars) {
	df <- melt(table(final[,c(mvar, "Batch")]))
	p <- ggplot(df, aes_string(x="Batch", y="value", fill=mvar)) + geom_bar(stat="identity", position="stack") + theme_classic() + ggtitle(sprintf("%s across batches", mvar)) + geom_hline(yintercept=samples_per_batch, linetype="dashed") + scale_fill_brewer(palette="Set1")
	print(p)
}

# distribution of instn/instn2 across batches
mvar <- "instn"
df <- melt(table(final[,c(mvar, "Batch")]))
p <- ggplot(df, aes_string(x="Batch", y="value", fill=mvar)) + geom_bar(stat="identity", position="stack") + theme_classic() + ggtitle(sprintf("%s across batches", mvar)) + geom_hline(yintercept=samples_per_batch, linetype="dashed")
print(p)
mvar <- "instn2"
df <- melt(table(final[,c(mvar, "Batch")]))
p <- ggplot(df, aes_string(x="Batch", y="value", fill=mvar)) + geom_bar(stat="identity", position="stack") + theme_classic() + ggtitle(sprintf("%s across batches", mvar)) + geom_hline(yintercept=samples_per_batch, linetype="dashed") + scale_fill_manual(values=color_list[[mvar]])
print(p)

# distribution of instn2 stratified by Arm
df <- melt(table(final[,c("Arm", mvar, "Batch")]))
p <- ggplot(df, aes_string(x="Batch", y="value", fill=mvar)) + geom_bar(stat="identity", position="stack") + facet_wrap(~Arm) + theme_classic() + ggtitle(sprintf("%s across batches, stratified by Arm", mvar)) + scale_fill_manual(values=color_list[[mvar]])
print(p)


dev.off()



