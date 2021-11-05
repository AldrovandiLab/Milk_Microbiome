
library(ggplot2)
library(reshape2)
library(plyr)
library(useful)
library(lubridate)

## perform joint matching for NWCS 610 (BMK transmission) and ZEBS combined
data.610 <- read.table("/Lab_Share/PROMISE/nwcs610/Selected 610 BMK FINAL.txt", header=T, as.is=T, sep="\t")
metadata.610 <- read.table("/Lab_Share/PROMISE/nwcs610/SELECTED.Aim2.030519.txt", header=T, as.is=T, sep="\t")
data.zebs <- read.table("/Lab_Share/ZEBS/matched_1to2.Aim1.010821.txt", header=T, as.is=T, sep="\t")

# make a table of subjects for randomization
# nwcs610 aim 1
data.subjects <- unique(subset(data.610, Aim1=="Aim1")[,c("Patient.ID", "Category", "Aim1")])
colnames(data.subjects) <- c("Patient.ID", "Category", "Study")
data.subjects$Category <- ifelse(data.subjects$Category=="Case Subjects", "Transmitter", "Control")
data.subjects$Study <- "NWCS610-Aim1"
# nwcs610 aim 2
tmp <- data.frame(Patient.ID=c(metadata.610$id1, metadata.610$id2), Category=c(metadata.610$pp_trtl1, metadata.610$pp_trtl2), Study="NWCS610-Aim2")
data.subjects <- rbind(data.subjects, tmp)
# ZEBS
tmp <- data.frame(Patient.ID=data.zebs$studyid, Category=ifelse(data.zebs$Transmitter=="Yes", "Transmitter", "Control"), Study="ZEBS")
data.subjects <- rbind(data.subjects, tmp)
tmp1 <- aggregate(Category ~ Patient.ID, data.subjects, function(x) paste(x, collapse=","))
tmp2 <- aggregate(Study ~ Patient.ID, data.subjects, function(x) paste(x, collapse=","))
data.subjects.final <- merge(tmp1, tmp2, by="Patient.ID")
data.subjects.final$Group <- gsub("Control,", "", data.subjects.final$Category)
rownames(data.subjects.final) <- data.subjects.final$Patient.ID

# evenly assign to 9 batches (20 subjects/batch)
groups <- unique(data.subjects.final$Group)
n <- 20 # 20 subjects per batch: max of 80 samples per 96-well plate, leaving room for 16 controls
nbatches <- 9
slots_remaining <- rep(n, nbatches); names(slots_remaining) <- 1:nbatches
# random assignment in batches to ensure evenness
# e.g. if you have 5 batches, assign samples 5 at a time evenly, then the remainder assign randomly to remaining
# need to keep track of how many slots remaining in each batch and sample only from the batches with slots remaining
set.seed(nrow(data.subjects.final))
for (gr in groups) {
	tmp <- subset(data.subjects.final, Group==gr)
	open_batches <- which(slots_remaining > 0)
	vec <- rep(as.numeric(names(slots_remaining)), times=slots_remaining)
	num_requested <- floor(nrow(tmp) / length(open_batches))
	num_even_batches <- min(num_requested, slots_remaining)
	nsubj.even <- nbatches * num_even_batches # number of subjects to split evenly among batches
	nsubj.random <- nrow(tmp) - nsubj.even # leftover subjects to be randomly assigned
	# sample evenly
	vec.even <- sample(rep(open_batches, num_even_batches))
	slots_remaining <- slots_remaining - table(vec.even)
	
	# sample randomly from remaining slots
	open_batches <- which(slots_remaining > 0)
	vec <- rep(as.numeric(names(slots_remaining)), times=slots_remaining)
	vec.random <- sample(vec, size=nsubj.random)
	slots_remaining[names(table(vec.random))] = slots_remaining[names(table(vec.random))] - table(vec.random)
	
	vec <- c(vec.even, vec.random)
	data.subjects.final[rownames(tmp), "Batch"] <- vec
}
table(data.subjects.final[,c("Group", "Batch")])




#######
# map onto samples
data.610$Batch <- data.subjects.final[as.character(data.610$Patient.ID), "Batch"]
data.610$Category <- data.subjects.final[as.character(data.610$Patient.ID), "Category"]
data.610$Study <- data.subjects.final[as.character(data.610$Patient.ID), "Study"]

tmp <- melt(data.zebs[,c("studyid", "PreDate", "TransmissionDate", "Post1Date", "Post2Date")], id.vars=c("studyid"))
colnames(tmp) <- c("studyid", "Visit", "VisitDate")
tmp$studyid <- as.character(tmp$studyid)
tmp$Batch <- data.subjects.final[tmp$studyid, "Batch"]
tmp$Category <- data.subjects.final[tmp$studyid, "Category"]
tmp$Study <- data.subjects.final[tmp$studyid, "Study"]
colnames(tmp) <- c("Patient.ID", "Visit", "TimePoint.Date", "Batch", "Category", "Study")

mvars <- c("Patient.ID", "TimePoint.Date", "Category", "Study", "Batch")

# fix date strings
data.610$TimePoint.Date <- format(as.Date(data.610$TimePoint.Date, format="%m/%d/%Y"), "%Y%m%d")
tmp$TimePoint.Date <- format(as.Date(tmp$TimePoint.Date, format="%Y-%m-%d"), "%Y%m%d")

data.samples.final <- unique(rbind(data.610[,mvars], tmp[,mvars]))
data.samples.final$Sample.ID <- sprintf("P%s.%s", data.samples.final$Patient.ID, data.samples.final$TimePoint.Date)

write.table(data.samples.final, file="/Lab_Share/PROMISE/nwcs610/nwcs610_ZEBS_combined_batching.samples.052421.txt", quote=F, sep="\t", row.names=F, col.names=T)
write.table(data.subjects.final, file="/Lab_Share/PROMISE/nwcs610/nwcs610_ZEBS_combined_batching.subjects.052421.txt", quote=F, sep="\t", row.names=F, col.names=T)
tab <- table(data.samples.final[,c("Study", "Batch")])
write.table(tab, file="/Lab_Share/PROMISE/nwcs610/nwcs610_ZEBS_combined_batching.Batch_counts_by_sample.052421.txt", quote=F, sep="\t", row.names=T, col.names=T)
tab <- table(data.subjects.final[,c("Study", "Batch")])
write.table(tab, file="/Lab_Share/PROMISE/nwcs610/nwcs610_ZEBS_combined_batching.Batch_counts_by_subjects.052421.txt", quote=F, sep="\t", row.names=T, col.names=T)
tab <- table(data.samples.final[,c("Category", "Batch")])
write.table(tab, file="/Lab_Share/PROMISE/nwcs610/nwcs610_ZEBS_combined_batching.Category_counts_by_sample.052421.txt", quote=F, sep="\t", row.names=T, col.names=T)
tab <- table(data.subjects.final[,c("Category", "Batch")])
write.table(tab, file="/Lab_Share/PROMISE/nwcs610/nwcs610_ZEBS_combined_batching.Category_counts_by_subjects.052421.txt", quote=F, sep="\t", row.names=T, col.names=T)




## 09/14/2021 - merging with LDMS data for shipment
shipment <- read.table("/Lab_Share/PROMISE/nwcs610/14SEP2021_ZEBSPROMISE_PEN.txt", header=T, as.is=T, sep="\t")
batching <- read.table("/Lab_Share/PROMISE/nwcs610/nwcs610_ZEBS_combined_batching.samples.052421.txt", header=T, as.is=T, sep="\t")

shipment$TimePoint.Date <- format(as.Date(shipment$Collection.Date, format="%d-%b-%y"), "%Y%m%d")
shipment$ID2 <- ifelse(shipment$Project=="ZEBS", shipment$ID1, unlist(lapply(shipment$ID1, function(x) substr(x, 1, 7))))
shipment$Sample.ID <- sprintf("P%s.%s", shipment$ID2, shipment$TimePoint.Date)
rownames(shipment) <- shipment$Sample.ID
batching$Global.Specimen.ID <- shipment[batching$Sample.ID, "Global.Specimen.ID"]
write.table(batching, file="/Lab_Share/PROMISE/nwcs610/nwcs610_ZEBS_combined_batching_with_LDMS.samples.052421.txt", quote=F, sep="\t", row.names=F, col.names=T)




