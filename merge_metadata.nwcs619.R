
library(ggplot2)
library(reshape2)
library(plyr)
library(useful)
library(lubridate)
library(dplyr)

# merge metadata for NWCS619

#mapping <- read.table("/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.081219.txt", header=T, as.is=T, sep="\t")

#infant_baseline <- read.csv("/Lab_Share/PROMISE/nwcs619/20190830/INFANT_BASELINE.csv")
#rownames(infant_baseline) <- infant_baseline$patid
#for (mvar in c("gender", "weight0week", "weight1week")) {
#	mapping[, mvar] <- infant_baseline[as.character(mapping$cpatid), mvar]
#}

#overall <- read.csv("/Lab_Share/PROMISE/nwcs619/20190830/OVERALL_STATUS.csv")
#rownames(overall) <- overall$patid
#for (mvar in c("ap_onstgage")) {
#	mapping[, mvar] <- overall[as.character(mapping$cpatid), mvar]
#}

#nbw <- read.csv("/Lab_Share/PROMISE/nwcs619/20190830/NBW0012.csv")
#rownames(nbw) <- nbw$patid
#for (mvar in c("nbclass")) {
#	mapping[, mvar] <- nbw[as.character(mapping$cpatid), mvar]
#}

#arv <- read.csv("/Lab_Share/PROMISE/nwcs619/20181207/ARV_REGIMEN_PROMISE_WITH_HOLDS.csv")
#sites <- aggregate(instn ~ patid, arv, function(x) paste(unique(x), collapse=","))
#rownames(sites) <- sites$patid
#mapping[, "instn.mom"] <- sites[as.character(mapping$patid), "instn"]
#mapping[, "instn.infant"] <- sites[as.character(mapping$cpatid), "instn"]

#pe6864t <- read.csv("/Lab_Share/PROMISE/nwcs619/20181207/PE6864T.csv")
#pe6864t.sel <- subset(pe6864t, primevnt=="Diagnosis" & dx40ap1 %in% c("MALARIA, CONFIRMED", "MALARIA, PROBABLE"))

#overall <- read.csv("/Lab_Share/PROMISE/nwcs619/20181207/OVERALL_STATUS.csv")
#rownames(overall) <- overall$patid
#overall$DaysFromEntryToDelivery <- as.Date(as.character(overall$deldtup), format="%d%B%Y") - as.Date(as.character(overall$ap_randdt1), format="%d%B%Y")
#mapping$DaysFromEntryToDelivery <- overall[as.character(mapping$patid), "DaysFromEntryToDelivery"]
#mapping$imputed.gage <- mapping$ap_onstgage + as.numeric(mapping$DaysFromEntryToDelivery / 7)
#mapping$InfantAgeInDays <- as.numeric(as.Date(mapping$InfantDBS.drawdt, format="%d%B%Y") - as.Date(mapping$deldtup, format="%d%B%Y"))
##mapping$InfantAgeInDaysBinned <- cut(mapping$InfantAgeInDays, breaks=c(-1,7,30,365))
#mapping$InfantAgeInDaysBinned <- cut(mapping$InfantAgeInDays, breaks=c(-1,1,3,30,365))

#write.table(mapping, file=sprintf("/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.083019.txt"), quote=F, sep="\t", row.names=F, col.names=T)


### fixing regimens from ARV.102319_nt-1-2.txt
#data <- read.table("/Lab_Share/PROMISE/nwcs619/20181207/ARV.102319_nt-1-2.txt", header=T, as.is=T, sep="\t", quote="", comment.char="")
#for (datevar in c("deldt1", "regimen_start", "regimen_stopdate", "Plasma.drawdt", "DBS.drawdt", "InfantDBS.drawdt")) {
#	data[, datevar] <- as.Date(data[, datevar], format="%d%B%Y")
#}
#primary_inds <- which(data$Ct == 1); names(primary_inds) <- data[primary_inds, "patid"] # primary index for each subject

## calculate maternal days on regimen
#tmp <- ddply(data[, c("patid", "Plasma.drawdt", "regimen_start", "regimen_stopdate")], .(patid), function(x) {
#	sel <- unlist(apply(x, 1, function(rowdata) {between(ymd(rowdata[2]), ymd(rowdata[3]), ymd(rowdata[4]))} ))
#	ind <- which(sel)
#	if (length(ind) == 1) {
#		x[ind, "Plasma.drawdt"] - x[ind, "regimen_start"]
#	} else {
#		-1
#	}
#})
#data[primary_inds[as.character(tmp$patid)], "Maternal.days.on.regimen"] <- tmp$V1
#data[which(is.na(data$Maternal.days.on.regimen)), "Maternal.days.on.regimen"] <- ""

## calculate infant days on regimen
#tmp <- ddply(data[, c("patid", "deldt1", "regimen_start", "regimen_stopdate")], .(patid), function(x) {
#	sel <- unlist(apply(x, 1, function(rowdata) {between(ymd(rowdata[2]), ymd(rowdata[3]), ymd(rowdata[4]))} ))
#	ind <- which(sel)
#	if (length(ind) == 1) {
#		x[ind, "deldt1"] - x[ind, "regimen_start"]
#	} else {
#		-1
#	}
#})
#data[primary_inds[as.character(tmp$patid)], "Infant.days.on.regimen"] <- tmp$V1
#data[which(is.na(data$Infant.days.on.regimen)), "Infant.days.on.regimen"] <- ""

#write.table(data, file="/Lab_Share/PROMISE/nwcs619/20181207/ARV.102319_nt-1-3.txt", quote=F, sep="\t", row.names=F, col.names=T)

## merge into mapping
#mapping <- read.table("/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.102919.txt", header=T, as.is=T, sep="\t", comment.char="", quote="")
#regimen <- data[primary_inds,]; rownames(regimen) <- regimen$patid
#rownames(mapping) <- mapping$patid
#mapping$MaternalRegimen <- regimen[rownames(mapping), "Maternal.Regimen"]; mapping$MaternalRegimen[which(mapping$MaternalRegimen %in% c("zdvart", "PI-based ART"))] <- "PI-ART"
#mapping$InfantRegimen <- regimen[rownames(mapping), "Infant.Regimen"]; mapping$InfantRegimen[which(mapping$InfantRegimen %in% c("zdvart", "PI-based ART"))] <- "PI-ART"
#mapping$MaternalGroup <- paste(mapping$delivery, mapping$MaternalRegimen, sep=".")
#mapping$InfantGroup <- paste(mapping$delivery, mapping$InfantRegimen, sep=".")

#write.table(mapping, file="/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.110419.txt", sep="\t", row.names=F, col.names=T, quote=F)


### 2/3/20: add DaysPTDPlasma, DaysPTDDBS and binned variants
#mapping <- read.table("/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.110419.txt", header=T, as.is=T, sep="\t", comment.char="", quote="")
#mapping$DaysPTDPlasma <- as.Date(mapping$deldtup, format="%d%B%Y") - as.Date(mapping$Plasma.drawdt, format="%d%B%Y")
#mapping$DaysPTDDBS <- as.Date(mapping$deldtup, format="%d%B%Y") - as.Date(mapping$DBS.drawdt, format="%d%B%Y")
#mapping$DaysPTDPlasma2 <- cut(as.numeric(mapping$DaysPTDPlasma), breaks=c(0,2,7,14,200))
#mapping$DaysPTDDBS2 <- cut(as.numeric(mapping$DaysPTDDBS), breaks=c(0,2,7,14,200))
#write.table(mapping, file="/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.020320.txt", sep="\t", row.names=F, col.names=T, quote=F)


## 2/6/2023: redefine PTB based on gesagebr variable and rerun analyses
mapping <- read.table("/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.020320.txt", header=T, as.is=T, sep="\t", comment.char="", quote="")
gesagebr <- read.csv("/Lab_Share/PROMISE/nwcs619/nwcs619/20230201/infant_gesagebr.csv", header=T, as.is=T)

inds <- match(mapping$cpatid, gesagebr$patid)
mapping$gesagebr <- gesagebr[inds, "gesagebr"]
mapping$Delivery.gesagebr <- ifelse(mapping$gesagebr <= 36, "Preterm", "Term")
mapping$MaternalGroup.gesagebr <- sprintf("%s.%s", mapping$Delivery.gesagebr, mapping$MaternalRegimen)

write.table(mapping, file="/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.020623.txt", sep="\t", row.names=F, col.names=T, quote=F)

