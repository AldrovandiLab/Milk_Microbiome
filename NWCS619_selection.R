
library(ggplot2)
library(reshape2)
library(plyr)
library(useful)

## perform selection and matching for NWCS619 cohort

## 4 groups: [preterm-ZDV, preterm-ZDVART, term-ZDV, term-ZDVART]
## n=25 for each group
## single timepoint prior to delivery
## maternal plasma and DBS

gestational_age_cutoff <- 36
timing_window <- -28
timing_window_infant <- 28

## load mom-baby info
mpatids <- read.csv("/Lab_Share/PROMISE/nwcs619/MPATIDS.csv", row.names=1)

## load gestational age and delivery data
data <- read.csv("/Lab_Share/PROMISE/nwcs619/ARV_START.csv")
rownames(data) <- data$patid
data$cpatid <- mpatids[rownames(data), "cpatid"]
## add country information
tmp <- read.csv("/Lab_Share/PROMISE/nwcs619/OVERALL_STATUS.csv", row.names=1)
data$country1 <- as.character(tmp[rownames(data), "country1"])
## remove Tanzania and Zimbabwe
data <- subset(data, !(country1 %in% c("Tanzania", "Zimbabwe")))

## load aliquot data
aliquot <- read.csv("/Lab_Share/PROMISE/nwcs619/volume_left updated 14dec18.csv")
aliquot$drawdt2 <- as.Date(as.character(aliquot$drawdt), format="%d%b%Y")

## load hematocrit data and merge into aliquot data
hema.infant <- read.csv("/Lab_Share/PROMISE/nwcs619/LBW0107T.csv")
hema.mom <- read.csv("/Lab_Share/PROMISE/nwcs619/LBW0109T.csv")
hema.infant <- subset(hema.infant, !(is.na(hemaval)))
hema.mom <- subset(hema.mom, !(is.na(hema)))
aliquot$matchstr <- sprintf("%s.%s", aliquot$patid, aliquot$drawdt)
hema.infant$matchstr <- sprintf("%s.%s", hema.infant$patid, hema.infant$visitdt)
hema.mom$matchstr <- sprintf("%s.%s", hema.mom$patid, hema.mom$visitdt)
write.table(hema.infant, file="/Lab_Share/PROMISE/nwcs619/hema.infant.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(hema.mom, file="/Lab_Share/PROMISE/nwcs619/hema.mom.txt", sep="\t", row.names=F, col.names=T, quote=F)

# hematocrit
agg.infant <- aggregate(hemaval ~ matchstr, hema.infant, mean); rownames(agg.infant) <- as.character(agg.infant$matchstr)
agg.mom <- aggregate(hema ~ matchstr, hema.mom, mean); rownames(agg.mom) <- as.character(agg.mom$matchstr); colnames(agg.mom) <- c("matchstr", "hemaval")
agg <- rbind(agg.infant, agg.mom)
# hemoglobin
agg.infant.hemo <- aggregate(hemocv ~ matchstr, hema.infant, mean); rownames(agg.infant.hemo) <- as.character(agg.infant.hemo$matchstr)
agg.mom.hemo <- aggregate(hemocv ~ matchstr, hema.mom, mean); rownames(agg.mom.hemo) <- as.character(agg.mom.hemo$matchstr); colnames(agg.mom.hemo) <- c("matchstr", "hemocv")
agg.hemo <- rbind(agg.infant.hemo, agg.mom.hemo)
#aliquot <- merge(aliquot, agg, by="matchstr", all.x=T)

## set random seed
set.seed(nrow(data))

## select preterm-ZDV moms with a delivery gestational age <=36weeks (do not enforce hematocrit, insufficient n=20)
data.sel <- subset(data, delgage <= gestational_age_cutoff & first_regimen=="zdv")
data.sel$PlasmaDateDiff <- NA
data.sel$DBSDateDiff <- NA
data.sel$Plasma.drawdt <- NA
data.sel$DBS.drawdt <- NA
data.sel$InfantDBSDateDiff <- NA
data.sel$InfantDBS.drawdt <- NA
for (i in 1:nrow(data.sel)) {
	pid <- data.sel[i, "patid"]
	cpid <- data.sel[i, "cpatid"]
	delivery_date <- as.Date(data.sel[i, "deldtup"], format="%d%b%Y")
	
	# check for plasma
	aliquot.sel <- subset(aliquot, patid==pid & aliquot_type == "Plasma" & volume_available >= 0.5)
	aliquot.sel$days_from_delivery <- aliquot.sel$drawdt2 - delivery_date
	aliquot.sel <- subset(aliquot.sel, days_from_delivery<0)
	if (nrow(aliquot.sel) > 0) {
		j <- which.max(aliquot.sel$days_from_delivery)
		data.sel[i, "PlasmaDateDiff"] <- as.numeric(aliquot.sel[j, "days_from_delivery"])
		data.sel[i, "Plasma.drawdt"] <- as.character(aliquot.sel[j, "drawdt"])
	}
	# check for DBS
	aliquot.sel <- subset(aliquot, patid==pid & aliquot_type == "DBS" & volume_available >= 2)
	aliquot.sel$days_from_delivery <- aliquot.sel$drawdt2 - delivery_date
	aliquot.sel <- subset(aliquot.sel, days_from_delivery<0)
	if (nrow(aliquot.sel) > 0) {
		j <- which.max(aliquot.sel$days_from_delivery)
		data.sel[i, "DBSDateDiff"] <- as.numeric(aliquot.sel[j, "days_from_delivery"])
		data.sel[i, "DBS.drawdt"] <- as.character(aliquot.sel[j, "drawdt"])
	}
	# check for infant DBS
	aliquot.sel <- subset(aliquot, patid==cpid & aliquot_type == "DBS" & volume_available >= 2)
	aliquot.sel$days_from_delivery <- aliquot.sel$drawdt2 - delivery_date
	aliquot.sel <- subset(aliquot.sel, days_from_delivery>0)
	if (nrow(aliquot.sel) > 0) {
		j <- which.min(aliquot.sel$days_from_delivery)
		data.sel[i, "InfantDBSDateDiff"] <- as.numeric(aliquot.sel[j, "days_from_delivery"])
		data.sel[i, "InfantDBS.drawdt"] <- as.character(aliquot.sel[j, "drawdt"])
	}
}
data.sel$Chooseable <- data.sel$PlasmaDateDiff >= timing_window & data.sel$DBSDateDiff >= timing_window & data.sel$InfantDBSDateDiff <= timing_window_infant
data.sel$GestationalAgeAtCollection <- data.sel$delgage + (data.sel$PlasmaDateDiff / 7)
data.sel.preterm.zdv <- data.sel; rownames(data.sel.preterm.zdv) <- data.sel.preterm.zdv$patid
selected.preterm.zdv <- as.character(sample(subset(data.sel, Chooseable)$patid, size=min(25, nrow(subset(data.sel, Chooseable)))))
# NOTE: only 20/26 subjects matched all criteria for preterm-ZDV, manually add 5 more samples that missed the timing cutoffs
#data.sel[setdiff(rownames(data.sel),selected.preterm.zdv), c(9,10,13,15)]
selected.preterm.zdv <- c(selected.preterm.zdv, c("3024444", "3024498", "3030374", "3035669", "3036625"))

#data.sel[,c(9,10,13,15)]

## select preterm-ZDVART moms with a delivery gestational age <=36weeks (enforce required hematocrit value, enough n=87)
data.sel <- subset(data, delgage <= gestational_age_cutoff & first_regimen=="3tc,zdv,lpv/r")
data.sel$PlasmaDateDiff <- NA
data.sel$DBSDateDiff <- NA
data.sel$Plasma.drawdt <- NA
data.sel$DBS.drawdt <- NA
data.sel$InfantDBSDateDiff <- NA
data.sel$InfantDBS.drawdt <- NA
for (i in 1:nrow(data.sel)) {
	pid <- data.sel[i, "patid"]
	cpid <- data.sel[i, "cpatid"]
	delivery_date <- as.Date(data.sel[i, "deldtup"], format="%d%b%Y")
	
	# check for plasma
	aliquot.sel <- subset(aliquot, patid==pid & aliquot_type == "Plasma" & volume_available >= 0.5 & matchstr %in% agg$matchstr)
	aliquot.sel$days_from_delivery <- aliquot.sel$drawdt2 - delivery_date
	aliquot.sel <- subset(aliquot.sel, days_from_delivery<0)
	if (nrow(aliquot.sel) > 0) {
		j <- which.max(aliquot.sel$days_from_delivery)
		data.sel[i, "PlasmaDateDiff"] <- as.numeric(aliquot.sel[j, "days_from_delivery"])
		data.sel[i, "Plasma.drawdt"] <- as.character(aliquot.sel[j, "drawdt"])
	}
	# check for DBS
	aliquot.sel <- subset(aliquot, patid==pid & aliquot_type == "DBS" & volume_available >= 2 & matchstr %in% agg$matchstr)
	aliquot.sel$days_from_delivery <- aliquot.sel$drawdt2 - delivery_date
	aliquot.sel <- subset(aliquot.sel, days_from_delivery<0)
	if (nrow(aliquot.sel) > 0) {
		j <- which.max(aliquot.sel$days_from_delivery)
		data.sel[i, "DBSDateDiff"] <- as.numeric(aliquot.sel[j, "days_from_delivery"])
		data.sel[i, "DBS.drawdt"] <- as.character(aliquot.sel[j, "drawdt"])
	}
	# check for infant DBS
	aliquot.sel <- subset(aliquot, patid==cpid & aliquot_type == "DBS" & volume_available >= 2 & matchstr %in% agg$matchstr)
	aliquot.sel$days_from_delivery <- aliquot.sel$drawdt2 - delivery_date
	aliquot.sel <- subset(aliquot.sel, days_from_delivery>0)
	if (nrow(aliquot.sel) > 0) {
		j <- which.min(aliquot.sel$days_from_delivery)
		data.sel[i, "InfantDBSDateDiff"] <- as.numeric(aliquot.sel[j, "days_from_delivery"])
		data.sel[i, "InfantDBS.drawdt"] <- as.character(aliquot.sel[j, "drawdt"])
	}
}
data.sel$Chooseable <- data.sel$PlasmaDateDiff >= timing_window & data.sel$DBSDateDiff >= timing_window & data.sel$InfantDBSDateDiff <= timing_window_infant
data.sel$GestationalAgeAtCollection <- data.sel$delgage + (data.sel$PlasmaDateDiff / 7)
data.sel.preterm.zdvart <- data.sel; rownames(data.sel.preterm.zdv) <- data.sel.preterm.zdv$patid
selected.preterm.zdvart <- as.character(sample(subset(data.sel, Chooseable)$patid, size=25))

## select term-ZDV moms with a delivery gestational age >=37weeks (insufficient infants with hematocrit to be required for infant DBS)
data.sel <- subset(data, delgage >= 37 & first_regimen=="zdv")
rownames(data.sel) <- data.sel$patid
aliquot.sel <- subset(aliquot, patid %in% data.sel$patid & ((aliquot_type == "Plasma" & volume_available >= 0.5) | aliquot_type == "DBS" & volume_available >= 2) & matchstr %in% agg$matchstr)
aliquot.sel$deldtup <- as.Date(data.sel[as.character(aliquot.sel$patid), "deldtup"], format="%d%B%Y")
aliquot.sel$delgage <- data.sel[as.character(aliquot.sel$patid), "delgage"]
aliquot.sel$GestationalAgeAtCollection <- aliquot.sel$delgage + (as.numeric(aliquot.sel$drawdt2 - aliquot.sel$deldtup) / 7)
aliquot.sel$drawdt <- as.character(aliquot.sel$drawdt); aliquot.sel$aliquot_type <- as.character(aliquot.sel$aliquot_type)
# subset to aliquot dates with sufficient Plasma+DBS
aliquot.sel <- aggregate(cbind(aliquot_type, drawdt) ~ patid + GestationalAgeAtCollection, aliquot.sel, function(x) paste(x, collapse=","))
aliquot.sel <- subset(aliquot.sel, aliquot_type == "DBS,Plasma")
aliquot.sel$patid <- as.character(aliquot.sel$patid)
aliquot.sel$country1 <- data[aliquot.sel$patid, "country1"]
# get valid infant aliquots, then filter mom aliquots to only those with good matching infant aliquots
child_samples <- subset(aliquot, patid %in% data.sel$cpatid & (aliquot_type == "DBS" & volume_available >= 2))
child_samples$deldtup <- as.Date(data.sel[match(as.character(child_samples$patid), data.sel$cpatid), "deldtup"], format="%d%B%Y")
child_samples$drawdt <- as.character(child_samples$drawdt); child_samples$aliquot_type <- as.character(child_samples$aliquot_type)
child_samples$datediff <- as.numeric(child_samples$drawdt2 - child_samples$deldtup)
child_samples <- subset(child_samples, datediff <= timing_window_infant)
child_samples <- ddply(child_samples, .(patid), function(x) {
	i <- which.min(x[, "datediff"])
	x[i,]
})
aliquot.sel <- subset(aliquot.sel, patid %in% rownames(subset(mpatids, cpatid %in% child_samples$patid)))
# greedily assign matches to each preterm-ZDV sample
match.preterm.zdv <- data.sel.preterm.zdv[selected.preterm.zdv,]
match.preterm.zdv$match.patid <- NA
match.preterm.zdv$match.cpatid <- NA
match.preterm.zdv$match.deldtup <- NA
match.preterm.zdv$match.delgage <- NA
match.preterm.zdv$match.GestationalAgeAtCollection <- NA
match.preterm.zdv$match.drawdt <- NA
match.preterm.zdv$match.cdrawdt <- NA
match.preterm.zdv$match.country1 <- NA
for (i in 1:nrow(match.preterm.zdv)) {
	aliquot.sel2 <- subset(aliquot.sel, country1==match.preterm.zdv[i, "country1"])
	sel.aliquot <- which.min(abs(aliquot.sel2$GestationalAgeAtCollection - match.preterm.zdv[i, "GestationalAgeAtCollection"]))
	sel.patid <- as.character(aliquot.sel2[sel.aliquot, "patid"])
	sel.gage <- aliquot.sel2[sel.aliquot, "GestationalAgeAtCollection"]
	match.preterm.zdv[i, "match.patid"] <- aliquot.sel2[sel.aliquot, "patid"]
	match.preterm.zdv[i, "match.cpatid"] <- data[aliquot.sel2[sel.aliquot, "patid"], "cpatid"]
	match.preterm.zdv[i, "match.deldtup"] <- as.character(data[sel.patid, "deldtup"])
	match.preterm.zdv[i, "match.delgage"] <- data[sel.patid, "delgage"]
	match.preterm.zdv[i, "match.GestationalAgeAtCollection"] <- aliquot.sel2[sel.aliquot, "GestationalAgeAtCollection"]
	match.preterm.zdv[i, "match.drawdt"] <- sapply(strsplit(aliquot.sel2[sel.aliquot, "drawdt"], ","), unique)
	aliquot.sel <- subset(aliquot.sel, !(patid %in% sel.patid))
	match.preterm.zdv[i, "match.cdrawdt"] <- subset(child_samples, patid==match.preterm.zdv[i, "match.cpatid"])$drawdt
	match.preterm.zdv[i, "match.country1"] <- as.character(data[sel.patid, "country1"])
}

## select term-ZDVART moms with a delivery gestational age >=37weeks (require both mom and baby to have hematocrit)
data.sel <- subset(data, delgage >= 37 & first_regimen=="3tc,zdv,lpv/r")
rownames(data.sel) <- data.sel$patid
aliquot.sel <- subset(aliquot, patid %in% data.sel$patid & ((aliquot_type == "Plasma" & volume_available >= 0.5) | aliquot_type == "DBS" & volume_available >= 2) & matchstr %in% agg$matchstr)
aliquot.sel$deldtup <- as.Date(data.sel[as.character(aliquot.sel$patid), "deldtup"], format="%d%B%Y")
aliquot.sel$delgage <- data.sel[as.character(aliquot.sel$patid), "delgage"]
aliquot.sel$GestationalAgeAtCollection <- aliquot.sel$delgage + (as.numeric(aliquot.sel$drawdt2 - aliquot.sel$deldtup) / 7)
aliquot.sel$drawdt <- as.character(aliquot.sel$drawdt); aliquot.sel$aliquot_type <- as.character(aliquot.sel$aliquot_type)
# subset to aliquot dates with sufficient Plasma+DBS
aliquot.sel <- aggregate(cbind(aliquot_type, drawdt) ~ patid + GestationalAgeAtCollection, aliquot.sel, function(x) paste(x, collapse=","))
aliquot.sel <- subset(aliquot.sel, aliquot_type == "DBS,Plasma")
aliquot.sel$patid <- as.character(aliquot.sel$patid)
aliquot.sel$country1 <- data[aliquot.sel$patid, "country1"]
# get valid infant aliquots, then filter mom aliquots to only those with good matching infant aliquots
child_samples <- subset(aliquot, patid %in% data.sel$cpatid & (aliquot_type == "DBS" & volume_available >= 2) & matchstr %in% agg$matchstr)
child_samples$deldtup <- as.Date(data.sel[match(as.character(child_samples$patid), data.sel$cpatid), "deldtup"], format="%d%B%Y")
child_samples$drawdt <- as.character(child_samples$drawdt); child_samples$aliquot_type <- as.character(child_samples$aliquot_type)
child_samples$datediff <- as.numeric(child_samples$drawdt2 - child_samples$deldtup)
child_samples <- subset(child_samples, datediff <= timing_window_infant)
child_samples <- ddply(child_samples, .(patid), function(x) {
	i <- which.min(x[, "datediff"])
	x[i,]
})
aliquot.sel <- subset(aliquot.sel, patid %in% rownames(subset(mpatids, cpatid %in% child_samples$patid)))
# greedily assign matches to each preterm-ZDVART sample
match.preterm.zdvart <- data.sel.preterm.zdvart[selected.preterm.zdvart,]
match.preterm.zdvart$match.patid <- NA
match.preterm.zdvart$match.cpatid <- NA
match.preterm.zdvart$match.deldtup <- NA
match.preterm.zdvart$match.delgage <- NA
match.preterm.zdvart$match.GestationalAgeAtCollection <- NA
match.preterm.zdvart$match.drawdt <- NA
match.preterm.zdvart$match.cdrawdt <- NA
match.preterm.zdvart$match.country1 <- NA
for (i in 1:nrow(match.preterm.zdvart)) {
	aliquot.sel2 <- subset(aliquot.sel, country1==match.preterm.zdvart[i, "country1"])
	sel.aliquot <- which.min(abs(aliquot.sel2$GestationalAgeAtCollection - match.preterm.zdvart[i, "GestationalAgeAtCollection"]))
	sel.patid <- as.character(aliquot.sel2[sel.aliquot, "patid"])
	sel.gage <- aliquot.sel2[sel.aliquot, "GestationalAgeAtCollection"]
	match.preterm.zdvart[i, "match.patid"] <- aliquot.sel2[sel.aliquot, "patid"]
	match.preterm.zdvart[i, "match.cpatid"] <- data[aliquot.sel2[sel.aliquot, "patid"], "cpatid"]
	match.preterm.zdvart[i, "match.deldtup"] <- as.character(data[sel.patid, "deldtup"])
	match.preterm.zdvart[i, "match.delgage"] <- data[sel.patid, "delgage"]
	match.preterm.zdvart[i, "match.GestationalAgeAtCollection"] <- aliquot.sel2[sel.aliquot, "GestationalAgeAtCollection"]
	match.preterm.zdvart[i, "match.drawdt"] <- sapply(strsplit(aliquot.sel2[sel.aliquot, "drawdt"], ","), unique)
	aliquot.sel <- subset(aliquot.sel, !(patid %in% sel.patid))
	match.preterm.zdvart[i, "match.cdrawdt"] <- subset(child_samples, patid==match.preterm.zdvart[i, "match.cpatid"])$drawdt
	match.preterm.zdvart[i, "match.country1"] <- as.character(data[sel.patid, "country1"])
}

## add hematocrit values
match.preterm.zdv$matchstrpremom <- sprintf("%s.%s", match.preterm.zdv$patid, match.preterm.zdv$DBS.drawdt)
match.preterm.zdv$matchstrpreinfant <- sprintf("%s.%s", match.preterm.zdv$cpatid, match.preterm.zdv$InfantDBS.drawdt)
match.preterm.zdv$matchstrtermmom <- sprintf("%s.%s", match.preterm.zdv$match.patid, match.preterm.zdv$match.drawdt)
match.preterm.zdv$matchstrterminfant <- sprintf("%s.%s", match.preterm.zdv$match.cpatid, match.preterm.zdv$match.cdrawdt)
match.preterm.zdv$premom_hemaval <- agg[match.preterm.zdv$matchstrpremom, "hemaval"]
match.preterm.zdv$preinfant_hemaval <- agg[match.preterm.zdv$matchstrpreinfant, "hemaval"]
match.preterm.zdv$termmom_hemaval <- agg[match.preterm.zdv$matchstrtermmom, "hemaval"]
match.preterm.zdv$terminfant_hemaval <- agg[match.preterm.zdv$matchstrterminfant, "hemaval"]

match.preterm.zdvart$matchstrpremom <- sprintf("%s.%s", match.preterm.zdvart$patid, match.preterm.zdvart$DBS.drawdt)
match.preterm.zdvart$matchstrpreinfant <- sprintf("%s.%s", match.preterm.zdvart$cpatid, match.preterm.zdvart$InfantDBS.drawdt)
match.preterm.zdvart$matchstrtermmom <- sprintf("%s.%s", match.preterm.zdvart$match.patid, match.preterm.zdvart$match.drawdt)
match.preterm.zdvart$matchstrterminfant <- sprintf("%s.%s", match.preterm.zdvart$match.cpatid, match.preterm.zdvart$match.cdrawdt)
match.preterm.zdvart$premom_hemaval <- agg[match.preterm.zdvart$matchstrpremom, "hemaval"]
match.preterm.zdvart$preinfant_hemaval <- agg[match.preterm.zdvart$matchstrpreinfant, "hemaval"]
match.preterm.zdvart$termmom_hemaval <- agg[match.preterm.zdvart$matchstrtermmom, "hemaval"]
match.preterm.zdvart$terminfant_hemaval <- agg[match.preterm.zdvart$matchstrterminfant, "hemaval"]


## store output
write.table(match.preterm.zdv, file="/Lab_Share/PROMISE/matching_nwcs619.ZDV.031319.txt", quote=F, sep="\t", row.names=F, col.names=T)
write.table(match.preterm.zdvart, file="/Lab_Share/PROMISE/matching_nwcs619.ZDVART.031319.txt", quote=F, sep="\t", row.names=F, col.names=T)

## country of origin data
match.preterm.zdv <- read.table("/Lab_Share/PROMISE/matching_nwcs619.ZDV.030519.txt", header=T, as.is=T, sep="\t")
match.preterm.zdvart <- read.table("/Lab_Share/PROMISE/matching_nwcs619.ZDVART.030519.txt", header=T, as.is=T, sep="\t")
tmp <- melt(match.preterm.zdv[, c("patid", "match.patid")]); tmp$cohort <- "ZDV"
tmp2 <- melt(match.preterm.zdvart[, c("patid", "match.patid")]); tmp2$cohort <- "ZDVART"
df <- rbind(tmp, tmp2)
colnames(df) <- c("Term", "ID", "Cohort"); df$Term <- as.character(df$Term); df$Term[which(df$Term=="patid")] <- "preterm"; df$Term[which(df$Term=="match.patid")] <- "term"; df$ID <- as.character(df$ID)
df$country <- data[df$ID, "country1"]
df$Cohort2 <- paste(df$Term, df$Cohort, sep="-")
table(df[,c("country","Cohort2")])

## EDIT 04/09/2019: Replace ZDVART-term control 601107 with another pair as this sample is not available
ind_to_replace <- which(match.preterm.zdvart$match.patid == "6011070")
i <- ind_to_replace
aliquot.sel2 <- subset(aliquot.sel, country1==match.preterm.zdvart[i, "country1"] & !patid=="6011070")
sel.aliquot <- which.min(abs(aliquot.sel2$GestationalAgeAtCollection - match.preterm.zdvart[i, "GestationalAgeAtCollection"]))
sel.patid <- as.character(aliquot.sel2[sel.aliquot, "patid"])
sel.gage <- aliquot.sel2[sel.aliquot, "GestationalAgeAtCollection"]
match.preterm.zdvart[i, "match.patid"] <- aliquot.sel2[sel.aliquot, "patid"]
match.preterm.zdvart[i, "match.cpatid"] <- data[aliquot.sel2[sel.aliquot, "patid"], "cpatid"]
match.preterm.zdvart[i, "match.deldtup"] <- as.character(data[sel.patid, "deldtup"])
match.preterm.zdvart[i, "match.delgage"] <- data[sel.patid, "delgage"]
match.preterm.zdvart[i, "match.GestationalAgeAtCollection"] <- aliquot.sel2[sel.aliquot, "GestationalAgeAtCollection"]
match.preterm.zdvart[i, "match.drawdt"] <- sapply(strsplit(aliquot.sel2[sel.aliquot, "drawdt"], ","), unique)
aliquot.sel <- subset(aliquot.sel, !(patid %in% sel.patid))
match.preterm.zdvart[i, "match.cdrawdt"] <- subset(child_samples, patid==match.preterm.zdvart[i, "match.cpatid"])$drawdt
match.preterm.zdvart[i, "match.country1"] <- as.character(data[sel.patid, "country1"])

match.preterm.zdvart$matchstrpremom <- sprintf("%s.%s", match.preterm.zdvart$patid, match.preterm.zdvart$DBS.drawdt)
match.preterm.zdvart$matchstrpreinfant <- sprintf("%s.%s", match.preterm.zdvart$cpatid, match.preterm.zdvart$InfantDBS.drawdt)
match.preterm.zdvart$matchstrtermmom <- sprintf("%s.%s", match.preterm.zdvart$match.patid, match.preterm.zdvart$match.drawdt)
match.preterm.zdvart$matchstrterminfant <- sprintf("%s.%s", match.preterm.zdvart$match.cpatid, match.preterm.zdvart$match.cdrawdt)
match.preterm.zdvart$premom_hemaval <- agg[match.preterm.zdvart$matchstrpremom, "hemaval"]
match.preterm.zdvart$preinfant_hemaval <- agg[match.preterm.zdvart$matchstrpreinfant, "hemaval"]
match.preterm.zdvart$termmom_hemaval <- agg[match.preterm.zdvart$matchstrtermmom, "hemaval"]
match.preterm.zdvart$terminfant_hemaval <- agg[match.preterm.zdvart$matchstrterminfant, "hemaval"]


write.table(match.preterm.zdvart, file="/Lab_Share/PROMISE/matching_nwcs619.ZDVART.040919.txt", quote=F, sep="\t", row.names=F, col.names=T)

## EDIT 05/01/2019: add hemoglobin values
match.preterm.zdv <- read.table("/Lab_Share/PROMISE/matching_nwcs619.ZDV.031319.txt", header=T, as.is=T, sep="\t")
match.preterm.zdv$premom_hemocv <- agg.hemo[match.preterm.zdv$matchstrpremom, "hemocv"]
match.preterm.zdv$preinfant_hemocv <- agg.hemo[match.preterm.zdv$matchstrpreinfant, "hemocv"]
match.preterm.zdv$termmom_hemocv <- agg.hemo[match.preterm.zdv$matchstrtermmom, "hemocv"]
match.preterm.zdv$terminfant_hemocv <- agg.hemo[match.preterm.zdv$matchstrterminfant, "hemocv"]
write.table(match.preterm.zdv, file="/Lab_Share/PROMISE/matching_nwcs619.ZDV.050119.txt", quote=F, sep="\t", row.names=F, col.names=T)

match.preterm.zdvart <- read.table("/Lab_Share/PROMISE/matching_nwcs619.ZDVART.040919.txt", header=T, as.is=T, sep="\t")
match.preterm.zdvart$premom_hemocv <- agg.hemo[match.preterm.zdvart$matchstrpremom, "hemocv"]
match.preterm.zdvart$preinfant_hemocv <- agg.hemo[match.preterm.zdvart$matchstrpreinfant, "hemocv"]
match.preterm.zdvart$termmom_hemocv <- agg.hemo[match.preterm.zdvart$matchstrtermmom, "hemocv"]
match.preterm.zdvart$terminfant_hemocv <- agg.hemo[match.preterm.zdvart$matchstrterminfant, "hemocv"]
write.table(match.preterm.zdvart, file="/Lab_Share/PROMISE/matching_nwcs619.ZDVART.050119.txt", quote=F, sep="\t", row.names=F, col.names=T)


