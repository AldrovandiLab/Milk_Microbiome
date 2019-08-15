library(ggplot2)
library(reshape2)
library(plyr)
library(gtools)

## perform selection and matching for NWCS610 cohort

## Aim 1: [HIV transmitters, 1:1 matched non-transmitters]
##  4 timepoints each [pre, HIV-infection, post1, post2]
##  match by infant age at time of transmission, maternal baseline CD4 and HIV RNA, country of origin, date of breast milk sample

timing_window <- 12*7 # 12 weeks

## load mom-baby info
mpatids <- read.csv("/Lab_Share/PROMISE/nwcs610/MPATIDS.csv", row.names=1)

## load ART and transmission data
data <- read.csv("/Lab_Share/PROMISE/nwcs610/ARV_START.csv")
rownames(data) <- data$patid
tmp <- read.csv("/Lab_Share/PROMISE/nwcs610/INFECT2.csv", row.names=1)
tmp$cpatid <- rownames(tmp)
# merge into mpatids
mpatids <- merge(mpatids, data, by="row.names"); rownames(mpatids) <- mpatids[,1]; mpatids <- mpatids[, -1]
mpatids <- merge(mpatids, tmp, by="cpatid"); rownames(mpatids) <- mpatids[, "mpatid"]
mpatids$statdt2 <- as.Date(as.character(mpatids$statdt), format="%d%b%Y")
mpatids$deldtup2 <- as.Date(as.character(mpatids$deldtup), format="%d%b%Y")
# add maternal baseline characteristics
tmp <- read.csv("/Lab_Share/PROMISE/nwcs610/BASELINEANDPREGOUT.csv", row.names=1)
mpatids <- merge(mpatids, tmp[,1:4], by="row.names"); rownames(mpatids) <- mpatids[,1]; mpatids <- mpatids[, -1]
data <- mpatids
# recode some regimens
data$first_regimen[which(data$first_regimen %in% c("zdv", "zdv,sd.ftc,sd.tdf,sd.nvp", "zdv,ftc,tdf"))] <- "zdv,sd.nvp,ftc,tdf"
data$first_regimen[which(data$first_regimen %in% c("3tc,zdv,lpv/r,sd.nvp"))] <- "3tc,zdv,lpv/r"
# add country
tmp <- read.csv("/Lab_Share/PROMISE/nwcs610/OVERALL_STATUS.csv", row.names=1)
data$country1 <- tmp[rownames(data), "country1"]


## load aliquot data and subset to useful aliquots (Milk volume >= 1.8mL)
aliquot <- read.csv("/Lab_Share/PROMISE/nwcs610/VOLUME_LEFT.csv")
aliquot$drawdt2 <- as.Date(as.character(aliquot$drawdt), format="%d%b%Y")
aliquot <- subset(aliquot, aliquot_type=="Milk" & volume_available >= 1.8)

## get all transmitters and then subset for those where infection occurred after DOL 14
transmitters <- subset(data, status=="Infected")
transmitters$InfectionAgeInDays <- transmitters$statdt2 - transmitters$deldtup2
transmitters <- subset(transmitters, InfectionAgeInDays > 14) # leaves n=33 subjects

## subset for transmitters with aliquots
sel <- rownames(transmitters)[which(rownames(transmitters) %in% aliquot$patid)]
transmitters <- transmitters[sel,]
out <- transmitters

## get timepoints for each transmitter
for (id in rownames(out)) {
	aliquot.sel <- subset(aliquot, patid == id)
	tmp <- aliquot.sel$drawdt2
	out[id, "aliquots"] <- paste(tmp, collapse=",")
}
#write.table(out, file="/Lab_Share/PROMISE/nwcs610/transmitters.022119.txt", quote=F, sep="\t", row.names=T, col.names=T)

# load manually selected transmitter timepoints
manual <- read.table("/Lab_Share/PROMISE/nwcs610/manually_selected_transmitters.030419.txt", header=T, as.is=T, sep="\t")
tmp <- ddply(manual, .(mpatid), function(x) paste(x$date, collapse=",")); rownames(tmp) <- tmp$mpatid
out$Timepoints <- tmp[rownames(out), "V1"]

### get timepoints for each transmitter
#for (id in rownames(transmitters)) {
#	aliquot.sel <- subset(aliquot, patid == id)
#	tmp <- abs(aliquot.sel$drawdt2 - transmitters[id, "statdt2"])
#	i <- which.min(tmp)
#	if (tmp[i] > timing_window) {
#		transmitters[id, "TransmissionSample"] <- NA
#		transmitters[id, "PreSample"] <- NA
#		transmitters[id, "Post1Sample"] <- NA
#		transmitters[id, "Post2Sample"] <- NA
#	} else {
#		transmitters[id, "TransmissionSample"] <- aliquot.sel$drawdt2[i]
#		transmitters[id, "PreSample"] <- aliquot.sel$drawdt2[i-1]
#		transmitters[id, "Post1Sample"] <- aliquot.sel$drawdt2[i+1]
#		transmitters[id, "Post2Sample"] <- aliquot.sel$drawdt2[i+2]
#	}
#}
##transmitters$Valid <- !is.na(transmitters$TransmissionSample) & !is.na(transmitters$PreSample) & !is.na(transmitters$Post1Sample) & !is.na(transmitters$Post2Sample)
#transmitters$Valid <- !is.na(transmitters$TransmissionSample)
#transmitters <- subset(transmitters, Valid)
#transmitters$AgeAtTransmissionSample <- transmitters$TransmissionSample - transmitters$deldtup2
#transmitters$AgeAtPreSample <- transmitters$PreSample - transmitters$deldtup2
#transmitters$AgeAtPost1Sample <- transmitters$Post1Sample - transmitters$deldtup2
#transmitters$AgeAtPost2Sample <- transmitters$Post2Sample - transmitters$deldtup2

## set random seed
set.seed(nrow(data))

## match controls based on mom's ARV regiment, followed by age at infection and maternal baseline CD4/HIVRNA
controls <- subset(data, status=="Uninfected"); controls$Used <- FALSE
transmitters <- out
transmitters$Match1 <- NA; transmitters$Timepoints1 <- NA; transmitters$ClosestAge1 <- NA
for (i in 1:nrow(transmitters)) {
	print(sprintf("Matching for transmitter %d/%d", i, nrow(transmitters)))
	id.transmitter <- rownames(transmitters)[i]
	controls.sel <- subset(controls, !Used & first_regimen==transmitters[id.transmitter, "first_regimen"])
	controls.sel$Timepoints <- NA; controls.sel$Valid <- FALSE; controls.sel$ClosestAge <- NA
	# identify pre/pre/infection/post timepoints for each candidate based on transmitter timepoints
	for (j in 1:nrow(controls.sel)) {
		aliquot.sel <- subset(aliquot, patid==rownames(controls.sel)[j])
		samples.trans <- sort(as.Date(unlist(strsplit(transmitters[id.transmitter, "Timepoints"], ","))))
		if (nrow(aliquot.sel) >= length(samples.trans)) {
			aliquot.sel <- aliquot.sel[order(aliquot.sel$drawdt2),]
			perms.candidate <- permutations(n=nrow(aliquot.sel), r=length(samples.trans), v=as.character(aliquot.sel$drawdt2))
			agediffs <- apply(perms.candidate, 1, function(x) sum(abs(samples.trans - as.Date(x))))
			sel <- which.min(agediffs)
			controls.sel[j, "ClosestAge"] <- agediffs[sel]
			# add timepoints if necessary to get 4
			tp <- perms.candidate[sel,]
			additional <- data.frame(datestr=setdiff(as.character(aliquot.sel$drawdt2), tp), stringsAsFactors=F)
			additional$dist <- sapply(additional$datestr, function(x) sum(abs(as.Date(x) - as.Date(tp))))
			additional <- additional[order(additional$dist),]
			if (length(tp) < 4) {
				tp <- c(tp, additional[1:(4-length(tp)), "datestr"])
			}
			controls.sel[j, "Timepoints"] <- paste(sort(as.Date(tp), na.last=T), collapse=",")
			controls.sel[j, "Valid"] <- !grepl("NA", controls.sel[j, "Timepoints"])
		}
	}
	controls.sel <- subset(controls.sel, Valid)
	# pick best control based on ClosestAge+cd4bl+rnacopbl
	if (nrow(controls.sel) > 0) {
		controls.sel$dist.age <- as.numeric(abs(controls.sel$ClosestAge))
		controls.sel$dist.cd4bl <- as.numeric(abs(controls.sel$cd4bl - transmitters[id.transmitter, "cd4bl"]))
		controls.sel$dist.rnacopbl <- as.numeric(abs(controls.sel$rnacopbl - transmitters[id.transmitter, "rnacopbl"]))
		controls.sel$dist <- controls.sel$dist.age + controls.sel$dist.cd4bl + log10(controls.sel$dist.rnacopbl)
		selected <- rownames(controls.sel)[which.min(controls.sel$dist)]
		transmitters[id.transmitter, "Match1"] <- selected
		transmitters[id.transmitter, "Timepoints1"] <- controls.sel[selected, "Timepoints"]
		transmitters[id.transmitter, "ClosestAge1"] <- controls.sel[selected, "ClosestAge"]
		controls[selected, "Used"] <- TRUE
	}
}

# match round 2
transmitters$Match2 <- NA; transmitters$Timepoints2 <- NA; transmitters$ClosestAge2 <- NA
for (i in 1:nrow(transmitters)) {
	print(sprintf("Matching for transmitter %d/%d", i, nrow(transmitters)))
	id.transmitter <- rownames(transmitters)[i]
	controls.sel <- subset(controls, !Used & first_regimen==transmitters[id.transmitter, "first_regimen"])
	controls.sel$Timepoints <- NA; controls.sel$Valid <- FALSE; controls.sel$ClosestAge <- NA
	# identify pre/pre/infection/post timepoints for each candidate based on transmitter timepoints
	for (j in 1:nrow(controls.sel)) {
		aliquot.sel <- subset(aliquot, patid==rownames(controls.sel)[j])
		samples.trans <- sort(as.Date(unlist(strsplit(transmitters[id.transmitter, "Timepoints"], ","))))
		if (nrow(aliquot.sel) >= length(samples.trans)) {
			aliquot.sel <- aliquot.sel[order(aliquot.sel$drawdt2),]
			perms.candidate <- permutations(n=nrow(aliquot.sel), r=length(samples.trans), v=as.character(aliquot.sel$drawdt2))
			agediffs <- apply(perms.candidate, 1, function(x) sum(abs(samples.trans - as.Date(x))))
			sel <- which.min(agediffs)
			controls.sel[j, "ClosestAge"] <- agediffs[sel]
			# add timepoints if necessary to get 4
			tp <- perms.candidate[sel,]
			additional <- data.frame(datestr=setdiff(as.character(aliquot.sel$drawdt2), tp), stringsAsFactors=F)
			additional$dist <- sapply(additional$datestr, function(x) sum(abs(as.Date(x) - as.Date(tp))))
			additional <- additional[order(additional$dist),]
			if (length(tp) < 4) {
				tp <- c(tp, additional[1:(4-length(tp)), "datestr"])
			}
			controls.sel[j, "Timepoints"] <- paste(sort(as.Date(tp), na.last=T), collapse=",")
			controls.sel[j, "Valid"] <- !grepl("NA", controls.sel[j, "Timepoints"])
		}
	}
	controls.sel <- subset(controls.sel, Valid)
	# pick best control based on ClosestAge+cd4bl+rnacopbl
	if (nrow(controls.sel) > 0) {
		controls.sel$dist.age <- as.numeric(abs(controls.sel$ClosestAge))
		controls.sel$dist.cd4bl <- as.numeric(abs(controls.sel$cd4bl - transmitters[id.transmitter, "cd4bl"]))
		controls.sel$dist.rnacopbl <- as.numeric(abs(controls.sel$rnacopbl - transmitters[id.transmitter, "rnacopbl"]))
		controls.sel$dist <- controls.sel$dist.age + controls.sel$dist.cd4bl + log10(controls.sel$dist.rnacopbl)
		selected <- rownames(controls.sel)[which.min(controls.sel$dist)]
		transmitters[id.transmitter, "Match2"] <- selected
		transmitters[id.transmitter, "Timepoints2"] <- controls.sel[selected, "Timepoints"]
		transmitters[id.transmitter, "ClosestAge2"] <- controls.sel[selected, "ClosestAge"]
		controls[selected, "Used"] <- TRUE
	}
}
# pull in metadata to verify matches
transmitters$first_regimen.match1 <- data[transmitters$Match1, "first_regimen"]
transmitters$cd4bl.match1 <- data[transmitters$Match1, "cd4bl"]
transmitters$rnacopbl.match1 <- data[transmitters$Match1, "rnacopbl"]
transmitters$first_regimen.match2 <- data[transmitters$Match2, "first_regimen"]
transmitters$cd4bl.match2 <- data[transmitters$Match2, "cd4bl"]
transmitters$rnacopbl.match2 <- data[transmitters$Match2, "rnacopbl"]


write.table(transmitters, file=sprintf("/Lab_Share/PROMISE/nwcs610/SELECTED.Aim1.%s.txt", format(Sys.Date(), "%m%d%y")), quote=F, sep="\t", row.names=T, col.names=T)
countries <- melt(transmitters[,c("mpatid", "Match1", "Match2")]); colnames(countries) <- c("cohort", "id"); countries$id <- as.character(countries$id)
countries$country <- data[countries$id, "country1"]
table(countries[,c("country", "cohort")])


## Aim 2: [HIV non-transmitters, n=25 for Maternal TDF+FTC+LVP/r and Infant NVP]
##  4 timepoints each [week 6, week 26, week 50, week 74]
##  match by maternal baseline CD4/HIVRNA, country of origin, and date of randomization +/- 3 months

timing_window <- 5 # +/- 5 days from the desired timepoints
desired <- c(6, 26, 50, 74)*7
n <- 25

## load mom-baby info
mpatids <- read.csv("/Lab_Share/PROMISE/nwcs610/MPATIDS.csv", row.names=1)

## load ART and transmission data
data <- read.csv("/Lab_Share/PROMISE/nwcs610/ARV_START.csv")
rownames(data) <- data$patid
tmp <- read.csv("/Lab_Share/PROMISE/nwcs610/INFECT2.csv", row.names=1)
tmp$cpatid <- rownames(tmp)
# merge into mpatids
mpatids <- merge(mpatids, data, by="row.names"); rownames(mpatids) <- mpatids[,1]; mpatids <- mpatids[, -1]
mpatids <- merge(mpatids, tmp, by="cpatid"); rownames(mpatids) <- mpatids[, "mpatid"]
mpatids$statdt2 <- as.Date(as.character(mpatids$statdt), format="%d%b%Y")
mpatids$deldtup2 <- as.Date(as.character(mpatids$deldtup), format="%d%b%Y")
# add maternal baseline characteristics
tmp <- read.csv("/Lab_Share/PROMISE/nwcs610/BASELINEANDPREGOUT.csv", row.names=1)
mpatids <- merge(mpatids, tmp[,1:4], by="row.names"); rownames(mpatids) <- mpatids[,1]; mpatids <- mpatids[, -1]
data <- mpatids
# add postpartum therapy arm (Maternal triple ARV vs Infant NVP)
tmp <- read.csv("/Lab_Share/PROMISE/nwcs610/BP_TRT.csv", row.names=3)
data$pp_trtl <- tmp[rownames(data), "pp_trtl"]
# add country
tmp <- read.csv("/Lab_Share/PROMISE/nwcs610/OVERALL_STATUS.csv", row.names=1)
data$country1 <- tmp[rownames(data), "country1"]

## load aliquot data and subset to useful aliquots (Milk volume >= 1.8mL)
aliquot <- read.csv("/Lab_Share/PROMISE/nwcs610/VOLUME_LEFT.csv")
aliquot$drawdt2 <- as.Date(as.character(aliquot$drawdt), format="%d%b%Y")
aliquot <- subset(aliquot, aliquot_type=="Milk" & volume_available >= 1.8)

## get all non-transmitters on selected arms (ftc,tdf,lpv/r and zdv,sd.nvp,ftc,tdf)
nontransmitters <- subset(data, status=="Uninfected" & !is.na(pp_trtl))

## subset for non-transmitters with aliquots
sel <- rownames(nontransmitters)[which(rownames(nontransmitters) %in% aliquot$patid)]
nontransmitters <- nontransmitters[sel,]
nontransmitters$Valid <- FALSE

## get timepoints for each transmitter
for (id in rownames(nontransmitters)) {
	aliquot.sel <- subset(aliquot, patid == id)
	tmp <- aliquot.sel$drawdt2
	ages <- tmp - nontransmitters[id, "deldtup2"]
	valid <- all(sapply(desired, function(x) any(abs(ages-x) < timing_window)))
	selected <- sapply(desired, function(x) which.min(abs(ages-x)))
	nontransmitters[id, "aliquots"] <- paste(tmp, collapse=",")
	nontransmitters[id, "selected_aliquots"] <- paste(tmp[selected], collapse=",")
	nontransmitters[id, "ages"] <- paste(ages[selected], collapse=",")
	nontransmitters[id, "Valid"] <- valid
}
table(nontransmitters[,c("pp_trtl", "Valid")])
#nontransmitters <- subset(nontransmitters, Valid & country1 %in% c("Malawi", "South Africa"))
nontransmitters <- subset(nontransmitters, Valid)

## match by CD4 and HIVRNA, stratified by country
pairs <- as.data.frame(t(combn(rownames(nontransmitters), 2)))
colnames(pairs) <- c("id1", "id2"); pairs$id1 <- as.character(pairs$id1); pairs$id2 <- as.character(pairs$id2)
pairs$country1 <- nontransmitters[pairs$id1, "country1"]; pairs$country2 <- nontransmitters[pairs$id2, "country1"]
pairs$unif <- pairs$country1 == pairs$country2
pairs$pp_trtl1 <- nontransmitters[pairs$id1, "pp_trtl"]; pairs$pp_trtl2 <- nontransmitters[pairs$id2, "pp_trtl"]
pairs$arms <- ifelse(pairs$pp_trtl1 == pairs$pp_trtl2, "Same", "Different")
pairs <- subset(pairs, unif & arms=="Different")
pairs$cd4bl.1 <- nontransmitters[pairs$id1, "cd4bl"]; pairs$cd4bl.2 <- nontransmitters[pairs$id2, "cd4bl"]
pairs$rnacopbl.1 <- nontransmitters[pairs$id1, "rnacopbl"]; pairs$rnacopbl.2 <- nontransmitters[pairs$id2, "rnacopbl"]
pairs$deldtup2.1 <- nontransmitters[pairs$id1, "deldtup2"]; pairs$deldtup2.2 <- nontransmitters[pairs$id2, "deldtup2"]
pairs$ages.1 <- nontransmitters[pairs$id1, "ages"]; pairs$ages.2 <- nontransmitters[pairs$id2, "ages"]
pairs$selected_aliquots.1 <- nontransmitters[pairs$id1, "selected_aliquots"]
pairs$selected_aliquots.2 <- nontransmitters[pairs$id2, "selected_aliquots"]
pairs$dist <- abs(nontransmitters[pairs$id1, "cd4bl"] - nontransmitters[pairs$id2, "cd4bl"]) + log10(abs(nontransmitters[pairs$id1, "rnacopbl"] - nontransmitters[pairs$id2, "rnacopbl"]))
pairs <- pairs[order(pairs$dist),]

selected <- data.frame()
for (i in 1:n) {
	ind <- which.min(pairs$dist)
	used <- c(pairs[ind, "id1"], pairs[ind, "id2"])
	selected <- rbind(selected, pairs[ind,])
	pairs <- subset(pairs, !(id1 %in% used) & !(id2 %in% used))
}
write.table(selected, file=sprintf("/Lab_Share/PROMISE/nwcs610/SELECTED.Aim2.%s.txt", format(Sys.Date(), "%m%d%y")), quote=F, sep="\t", row.names=F, col.names=T)





