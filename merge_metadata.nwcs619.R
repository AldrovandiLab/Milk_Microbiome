
library(ggplot2)
library(reshape2)
library(plyr)
library(useful)

# merge metadata for NWCS619

mapping <- read.table("/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.081219.txt", header=T, as.is=T, sep="\t")

infant_baseline <- read.csv("/Lab_Share/PROMISE/nwcs619/20190830/INFANT_BASELINE.csv")
rownames(infant_baseline) <- infant_baseline$patid
for (mvar in c("gender", "weight0week", "weight1week")) {
	mapping[, mvar] <- infant_baseline[as.character(mapping$cpatid), mvar]
}

overall <- read.csv("/Lab_Share/PROMISE/nwcs619/20190830/OVERALL_STATUS.csv")
rownames(overall) <- overall$patid
for (mvar in c("ap_onstgage")) {
	mapping[, mvar] <- overall[as.character(mapping$cpatid), mvar]
}

nbw <- read.csv("/Lab_Share/PROMISE/nwcs619/20190830/NBW0012.csv")
rownames(nbw) <- nbw$patid
for (mvar in c("nbclass")) {
	mapping[, mvar] <- nbw[as.character(mapping$cpatid), mvar]
}

arv <- read.csv("/Lab_Share/PROMISE/nwcs619/20181207/ARV_REGIMEN_PROMISE_WITH_HOLDS.csv")
sites <- aggregate(instn ~ patid, arv, function(x) paste(unique(x), collapse=","))
rownames(sites) <- sites$patid
mapping[, "instn.mom"] <- sites[as.character(mapping$patid), "instn"]
mapping[, "instn.infant"] <- sites[as.character(mapping$cpatid), "instn"]

pe6864t <- read.csv("/Lab_Share/PROMISE/nwcs619/20181207/PE6864T.csv")
pe6864t.sel <- subset(pe6864t, primevnt=="Diagnosis" & dx40ap1 %in% c("MALARIA, CONFIRMED", "MALARIA, PROBABLE"))

overall <- read.csv("/Lab_Share/PROMISE/nwcs619/20181207/OVERALL_STATUS.csv")
rownames(overall) <- overall$patid
overall$DaysFromEntryToDelivery <- as.Date(as.character(overall$deldtup), format="%d%B%Y") - as.Date(as.character(overall$ap_randdt1), format="%d%B%Y")
mapping$DaysFromEntryToDelivery <- overall[as.character(mapping$patid), "DaysFromEntryToDelivery"]
mapping$imputed.gage <- mapping$ap_onstgage + as.numeric(mapping$DaysFromEntryToDelivery / 7)
mapping$InfantAgeInDays <- as.numeric(as.Date(mapping$InfantDBS.drawdt, format="%d%B%Y") - as.Date(mapping$deldtup, format="%d%B%Y"))
#mapping$InfantAgeInDaysBinned <- cut(mapping$InfantAgeInDays, breaks=c(-1,7,30,365))
mapping$InfantAgeInDaysBinned <- cut(mapping$InfantAgeInDays, breaks=c(-1,1,3,30,365))

write.table(mapping, file=sprintf("/Lab_Share/PROMISE/nwcs619/nwcs619_Mapping.083019.txt"), quote=F, sep="\t", row.names=F, col.names=T)


