

library(ggplot2)
library(reshape2)
library(plyr)
library(useful)
library(haven)
library(tableone)
library(lubridate)

## match up metadata for NWCS 610
mapping <- read.table("/Lab_Share/PROMISE/nwcs610/combined_nwcs610_ZEBS_Mapping.072722.txt", header=T, as.is=T, sep="\t")
rownames(mapping) <- mapping$Sample.ID

# create study variables
mapping$Study2 <- ifelse(mapping$Study == "ZEBS", "ZEBS", ifelse(mapping$Study %in% c("NWCS610-Aim1", "NWCS610-Aim2", "NWCS610-Aim1,NWCS610-Aim2"), "NWCS610", "Control"))
mapping$Aim1 <- ifelse(mapping$Study %in% c("NWCS610-Aim1", "ZEBS"), "Yes", "No") # some NWCS610-Aim1,NWCS610-Aim2 samples are only Aim2
mapping$Aim2 <- ifelse(mapping$Study %in% c("NWCS610-Aim2", "NWCS610-Aim1,NWCS610-Aim2"), "Yes", "No")

##############################################################################################################
## Pull in original metadata
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
mpatids <- merge(mpatids, tmp[,c(1:4,6)], by="row.names"); rownames(mpatids) <- mpatids[,1]; mpatids <- mpatids[, -1]
data <- mpatids
# recode some regimens
data$first_regimen[which(data$first_regimen %in% c("zdv", "zdv,sd.ftc,sd.tdf,sd.nvp", "zdv,ftc,tdf"))] <- "zdv,sd.nvp,ftc,tdf"
data$first_regimen[which(data$first_regimen %in% c("3tc,zdv,lpv/r,sd.nvp"))] <- "3tc,zdv,lpv/r"
# add country
tmp <- read.csv("/Lab_Share/PROMISE/nwcs610/OVERALL_STATUS.csv", row.names=1)
data$country1 <- tmp[rownames(data), "country1"]
data$ethnicity <- tmp[rownames(data), "ap_race1"]
# add study site
tmp <- read.csv("/Lab_Share/PROMISE/nwcs610/ARV_REGIMEN_PROMISE_WITH_HOLDS.csv")
tmp$instn <- as.character(tmp$instn)
arv_regimen_data <- tmp
tmp2 <- aggregate(instn ~ patid, tmp, function(x) paste(unique(x), collapse=","))
rownames(tmp2) <- as.character(tmp2$patid)
data$instn <- tmp2[rownames(data), "instn"]


##############################################################################################################
## Aim 1: [HIV transmitters, matched non-transmitters from nwcs610 and ZEBS]
##  4 timepoints each [pre, HIV-infection, post1, post2]
##  match by infant age at time of transmission, maternal baseline CD4 and HIV RNA, country of origin, date of breast milk sample

# set Aim 1 timepoint (nwcs610)
metadata <- read.table("/Lab_Share/PROMISE/nwcs610/manually_selected_transmitters.030419.txt", header=T, as.is=T, sep="\t")
metadata$date <- as.Date(metadata$date, format="%Y-%m-%d")
metadata$date2 <- format(metadata$date, format="%Y%m%d")
metadata$SampleID <- sprintf("P%s.%s", metadata$mpatid, metadata$date2)
rownames(metadata) <- metadata$SampleID
mapping$Timepoint <- metadata[rownames(mapping), "timepoint"]
mapping[rownames(metadata), "Aim1"] <- "Yes"

# calculate ARV regimen
arv_regimen_data$regimen_start <- parse_date_time(arv_regimen_data$regimen_start, order="dmy")
arv_regimen_data$regimen_stopdate <- parse_date_time(arv_regimen_data$regimen_stopdate, order="dmy")
for (i in 1:nrow(mapping)) {
	pid <- mapping[i, "Patient.ID"]
	tp <- parse_date_time(mapping[i, "TimePoint.Date"], order="ymd")
	arv.sel <- subset(arv_regimen_data, patid %in% pid)
	j <- which(tp > arv.sel$regimen_start & tp <= arv.sel$regimen_stopdate)
	if (length(j)==0) {
		# if not in the middle of a regimen, then check for prior regimen one day prior
		j <- which(tp-86400 > arv.sel$regimen_start & tp-86400 <= arv.sel$regimen_stopdate)
	}
	mapping[i, "ARVregimen"] <- paste(arv.sel[j, "regimen"], collapse=";")
}
# create drug-specific variables
lut.drug <- c("zdv", "x3tc", "efv", "ftc"); names(lut.drug) <- c("zdv", "3tc", "efv", "ftc")
for (drug in names(lut.drug)) {
	mapping[, lut.drug[drug]] <- unlist(lapply(mapping$ARVregimen, function(x) {
		drug %in% unlist(strsplit(x, ","))
	}))
}

# grab additional metadata (nwcs610)
metadata <- read.table("/Lab_Share/PROMISE/nwcs610/SELECTED.Aim1.030419.txt", header=T, as.is=T, sep="\t")
mapping[, "deldtup2"] <- metadata[as.character(mapping$Patient.ID), "deldtup2"]
mapping[, "cd4bl"] <- metadata[as.character(mapping$Patient.ID), "cd4bl"]
mapping[, "log10vlrna"] <- log10(metadata[as.character(mapping$Patient.ID), "rnacopbl"])
mapping[, "firstregimen"] <- metadata[as.character(mapping$Patient.ID), "first_regimen"]
# metadata on non-transmitters
for (i in 1:nrow(metadata)) {
	id.transmitter <- rownames(metadata)[i]
	id.match1 <- as.character(metadata$Match1[i])
	id.match2 <- as.character(metadata$Match2[i])
	tmp <- subset(mapping, Patient.ID==id.transmitter)
	tmp$date <- as.Date(as.character(tmp$TimePoint.Date), format="%Y%m%d")
	visits <- tmp[order(tmp$date), "Timepoint"]
	# match1
	dates.match1 <- unlist(lapply(strsplit(metadata[i, "Timepoints1"], ","), function(x) as.Date(x, format="%Y-%m-%d")))
	tmp <- subset(mapping, Patient.ID==id.match1)
	tmp$date <- as.Date(as.character(tmp$TimePoint.Date), format="%Y%m%d")
	tmp <- subset(tmp, date %in% dates.match1)
	tmp[order(tmp$date), "Timepoint"] <- c(visits, rep("Post", times=nrow(tmp)-length(visits))) # pad with Post visits as needed
	mapping[rownames(tmp), "Timepoint"] <- tmp$Timepoint
	mapping[rownames(tmp), "Aim1"] <- "Yes"
	mapping[rownames(tmp), "cd4bl"] <- metadata[i, "cd4bl.match1"]
	mapping[rownames(tmp), "log10vlrna"] <- log10(metadata[i, "rnacopbl.match1"])
	mapping[rownames(tmp), "firstregimen"] <- metadata[i, "first_regimen.match1"]
	# match2
	dates.match2 <- unlist(lapply(strsplit(metadata[i, "Timepoints2"], ","), function(x) as.Date(x, format="%Y-%m-%d")))
	tmp <- subset(mapping, Patient.ID==id.match2)
	tmp$date <- as.Date(as.character(tmp$TimePoint.Date), format="%Y%m%d")
	tmp <- subset(tmp, date %in% dates.match2)
	tmp[order(tmp$date), "Timepoint"] <- c(visits, rep("Post", times=nrow(tmp)-length(visits))) # pad with Post visits as needed
	mapping[rownames(tmp), "Timepoint"] <- tmp$Timepoint
	mapping[rownames(tmp), "Aim1"] <- "Yes"
	mapping[rownames(tmp), "cd4bl"] <- metadata[i, "cd4bl.match2"]
	mapping[rownames(tmp), "log10vlrna"] <- log10(metadata[i, "rnacopbl.match2"])
	mapping[rownames(tmp), "firstregimen"] <- metadata[i, "first_regimen.match2"]
}
mapping$deldtup2 <- data[as.character(mapping$Patient.ID), "deldtup2"]
mapping$country <- data[as.character(mapping$Patient.ID), "country1"]
mapping$ethnicity <- data[as.character(mapping$Patient.ID), "ethnicity"]
mapping$delgage <- data[as.character(mapping$Patient.ID), "delgage"]
mapping$parity <- data[as.character(mapping$Patient.ID), "parity1"]
mapping$statdt2 <- data[as.character(mapping$Patient.ID), "statdt2"]

# set Aim 1 timepoint and grab additional metadata (ZEBS)
metadata <- read.table("/Lab_Share/ZEBS/matched_1to2.Aim1.010821.with_metadata.txt", header=T, as.is=T, sep="\t")
rownames(metadata) <- metadata$studyid
zebs <- rownames(subset(mapping, Study=="ZEBS"))
mvars <- c("CD4C", "log10vlrna")
mapping[zebs, "cd4bl"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "CD4C"]
mapping[zebs, "vlrna"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "vlrna"]
mapping[zebs, "log10vlrna"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "log10vlrna"]
mapping[zebs, "deldtup2"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "DOB"]
mapping[zebs, "country"] <- "Zambia"
mapping[zebs, "instn"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "CLINIC"]
mapping[zebs, "lndate"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "lndate"]
mapping[zebs, "lnvis"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "lnvis"]
mapping[zebs, "fpdate"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "fpdate"]
mapping[zebs, "fpvis"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "fpvis"]
mapping[zebs, "lnage"] <- as.numeric(as.Date(mapping[zebs, "lndate"], format="%Y-%m-%d") - as.Date(mapping[zebs, "deldtup2"], format="%Y-%m)-%d"))
mapping[zebs, "fpage"] <- as.numeric(as.Date(mapping[zebs, "fpdate"], format="%Y-%m-%d") - as.Date(mapping[zebs, "deldtup2"], format="%Y-%m)-%d"))
mapping[zebs, "parity"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "PRELIVE"] + metadata[as.character(mapping[zebs, "Patient.ID"]), "PRESTILL"] - metadata[as.character(mapping[zebs, "Patient.ID"]), "PREVMULT"]
mvars <- c("HGB", "CD8C", "CD3C", "PREVPREG", "SEX", "BWT", "CD4_12M", "CD8_12M", "CD3_12M", "rna4m", "rna12m", "cd", "md", "timec", "timem", "bf", "bftime")
for (mvar in mvars) {
	mapping[zebs, mvar] <- metadata[as.character(mapping[zebs, "Patient.ID"]), mvar]
}
mvars <- c("SEX", "cd", "md", "bf")
for (mvar in mvars) {
	mapping[, mvar] <- factor(mapping[, mvar])
}
# populate ARV variables
tmp <- subset(metadata, inarv==1)
metadata[as.character(mapping[zebs, "Patient.ID"]), "CD4C"]
for (i in 1:nrow(tmp)) {
	pid <- rownames(tmp)[i]
	startdate <- parse_date_time(tmp[i, "startdate"], order="ymd")
	tmp.mapping <- subset(mapping, Patient.ID==pid)
	for (j in 1:nrow(tmp.mapping)) {
		tp <- parse_date_time(tmp.mapping[j, "TimePoint.Date"], order="ymd")
		print(sprintf("PID: %s  startdate: %s  Timepoint: %s  Difference: %d days", pid, startdate, tp, startdate-tp))
		if (tp > startdate) {
			# timepoint occurred after start of ARV therapy
			mapping[rownames(tmp)[j], "ARVregimen"] <- "ZEBS ARV"
		}
	}
}


#################


# 5/3/2023: additional metadata from Louise
pplflow_weaning <- as.data.frame(read_sas("/Lab_Share/ZEBS/pplflow_weaning.sas7bdat"))
rownames(pplflow_weaning) <- pplflow_weaning$studyid
mvars <- c("copies0_0", "copies1_0", "copies4_0", "copies4_5")
for (mvar in mvars) {
	mapping[zebs, mvar] <- pplflow_weaning[as.character(mapping[zebs, "Patient.ID"]), mvar]
	mapping[zebs, sprintf("%slog10", mvar)] <- log10(pplflow_weaning[as.character(mapping[zebs, "Patient.ID"]), mvar])
	mapping[zebs, sprintf("%svalid", mvar)] <- !is.na(mapping[zebs, mvar])
	mapping[zebs, sprintf("%sundetectable", mvar)] <- pplflow_weaning[as.character(mapping[zebs, "Patient.ID"]), mvar]==49
}
mvars <- c("ebf0_0", "ebf1_0", "ebf4_0", "ebf4_5", "MODE")
for (mvar in mvars) {
	mapping[zebs, mvar] <- factor(pplflow_weaning[as.character(mapping[zebs, "Patient.ID"]), mvar])
}
mvars <- c("MOMAGE")
for (mvar in mvars) {
	mapping[zebs, mvar] <- pplflow_weaning[as.character(mapping[zebs, "Patient.ID"]), mvar]
}
for (mvar in c("Q2C1", "Q2C2", "Q2C3")) {
	mapping[zebs, mvar] <- pplflow_weaning[as.character(mapping[zebs, "Patient.ID"]), mvar]
} # child death causes
# compute gestational age
pplflow_weaning$DATEE <- as.Date(parse_date_time(pplflow_weaning$DATEE, orders="dmy"))
pplflow_weaning$DOB <- as.Date(pplflow_weaning$DOB, format="%Y-%m-%d")
pplflow_weaning$GWTODAY <- ifelse(pplflow_weaning$GWTODAY==99, NA, pplflow_weaning$GWTODAY)
pplflow_weaning$GWDELIVERY <- pplflow_weaning$GWTODAY + (as.numeric(pplflow_weaning$DOB - pplflow_weaning$DATEE) / 7)
mapping[zebs, "delgage"] <- pplflow_weaning[as.character(mapping[zebs, "Patient.ID"]), "GWDELIVERY"]

#tmp <- melt(metadata[, c("studyid", "Transmitter", "PreDate", "PreAge", "PreVisit", "TransmissionDate", "TransmissionAge", "TransmissionVisit", "Post1Date", "Post1Age", "Post1Visit", "Post2Date", "Post2Age", "Post2Visit")])
tmp <- melt(metadata[, c("studyid",  "PreDate", "TransmissionDate", "Post1Date", "Post2Date")], id.vars=c("studyid"))
lut.tmp <- c("1 Wk", "1 Mo", "4 Mo", "4.5 Mo"); names(lut.tmp) <- c("PreDate", "TransmissionDate", "Post1Date", "Post2Date"); tmp$variable2 <- lut.tmp[tmp$variable]
tmp$variable <- gsub("[12]", "", gsub("Date", "", as.character(tmp$variable)))
tmp$date <- format(as.Date(tmp$value, format="%Y-%m-%d"), "%Y%m%d")
tmp$SampleID <- sprintf("P%s.%s", tmp$studyid, tmp$date)
mapping[tmp$SampleID, "Timepoint"] <- tmp$variable
mapping[tmp$SampleID, "Visit"] <- tmp$variable2

# create Transmitter variable
mapping$Transmitter <- ifelse(mapping$Category == "Transmitter", "Transmitter", NA)
mapping[grep("Control", mapping$Category), "Transmitter"] <- "Control"


##############################################################################################################
## Aim 2: [HIV non-transmitters, n=25 for Maternal TDF+FTC+LVP/r and Infant NVP]
##  4 timepoints each [week 6, week 26, week 50, week 74]
##  match by maternal baseline CD4/HIVRNA, country of origin, and date of randomization +/- 3 months

# set Aim 2 visit (nwcs610)
metadata <- read.table("/Lab_Share/PROMISE/nwcs610/SELECTED.Aim2.030519.txt", header=T, as.is=T, sep="\t")
tmp <- {}
visits <- c("Week06", "Week26", "Week50", "Week74")
for (i in 1:nrow(metadata)) {
	id <- as.character(metadata$id1[i])
	dates <- lapply(strsplit(metadata[i, "selected_aliquots.1"], ","), function(x) as.Date(as.character(x), format="%Y-%m-%d"))[[1]]
	ages <- as.numeric(unlist(strsplit(metadata[i, "ages.1"], ",")))
	df <- data.frame(id=id, date=dates, visit=visits[order(dates)], age=ages, deldtup2=metadata[i, "deldtup2.1"], SampleID=sprintf("P%s.%s", id, format(dates, "%Y%m%d")), country=metadata[i, "country1"], regimen=metadata[i, "pp_trtl1"], cd4bl=metadata[i, "cd4bl.1"], rnacopbl=metadata[i, "rnacopbl.1"])
	tmp <- rbind(tmp, df)
	
	id <- as.character(metadata$id2[i])
	dates <- lapply(strsplit(metadata[i, "selected_aliquots.2"], ","), function(x) as.Date(as.character(x), format="%Y-%m-%d"))[[1]]
	ages <- as.numeric(unlist(strsplit(metadata[i, "ages.2"], ",")))
	df <- data.frame(id=id, date=dates, visit=visits[order(dates)], age=ages, deldtup2=metadata[i, "deldtup2.2"], SampleID=sprintf("P%s.%s", id, format(dates, "%Y%m%d")), country=metadata[i, "country2"], regimen=metadata[i, "pp_trtl2"], cd4bl=metadata[i, "cd4bl.2"], rnacopbl=metadata[i, "rnacopbl.2"])
	tmp <- rbind(tmp, df)
}
rownames(tmp) <- tmp$SampleID
mapping[rownames(tmp), "cd4bl"] <- tmp[, "cd4bl"]
mapping[rownames(tmp), "log10vlrna"] <- log10(tmp[, "rnacopbl"])
mapping[rownames(tmp), "age"] <- tmp[, "age"]
mapping[rownames(tmp), "deldtup2"] <- as.character(tmp[, "deldtup2"])
mapping[rownames(tmp), "country"] <- tmp[, "country"]
mapping[rownames(tmp), "Visit"] <- as.character(tmp[, "visit"])
mapping[rownames(tmp), "Regimen"] <- tmp[, "regimen"]
mapping$Regimen2 <- ifelse(mapping$SampleType != "BMK", NA, ifelse(mapping$Study=="ZEBS", "no ARVs", ifelse(mapping$ARVregimen=="no ARVs", "no ARVs", ifelse(mapping$ARVregimen=="ftc,tdf,lpv/r", "Maternal triple ARV", "Other"))))
mapping$Regimen3 <- ifelse(mapping$SampleType != "BMK", NA, ifelse(mapping$Study=="ZEBS", "no ARVs", ifelse(mapping$ARVregimen=="no ARVs", "no ARVs", "Maternal triple ARV")))

sel <- rownames(subset(mapping, Aim2=="Yes" & is.na(Regimen)))
mapping[sel, "Aim2"] <- "No"



##############################################################################################################
# calculate age as needed
mapping$age <- as.numeric(as.Date(as.character(mapping$TimePoint.Date), format="%Y%m%d") - as.Date(mapping$deldtup2, format="%Y-%m-%d"))

# add institution
mapping$instn <- ifelse(is.na(mapping$instn), data[as.character(mapping$Patient.ID), "instn"], mapping$instn)

# add binned maternal CD4 variables
mapping$cd4bl_binned <- ifelse(mapping$cd4bl > 350, ">350", "<=350")
mapping$cd4bl_binned2 <- ifelse(mapping$cd4bl > 200, ">200", "<=200")

# add age in days relative to transmission
#mapping$DaysRelativeToTransmission <- as.numeric(as.Date(as.character(mapping$TimePoint.Date), format="%Y%m%d") - mapping$statdt2)
a <- parse_date_time(ifelse(mapping$Study2=="NWCS610", as.character(mapping$statdt2), as.character(mapping$fpdate)), order="ymd")
mapping$DaysRelativeToTransmission <- as.numeric(parse_date_time(as.character(mapping$TimePoint.Date), order="ymd") - a, units="days")
mapping$fpage <- as.numeric(a - parse_date_time(as.character(mapping$deldtup2), order="ymd"), units="days")
mapping[which(mapping$Transmitter=="Control"), "fpage"] <- NA

# manually fix PID 6007720 to Pre visits as all occurred prior to infection at statdt2
inds <- which(mapping$Patient.ID=="6007720")
mapping[inds, "Timepoint"] <- "Pre"

# recompute Timepoint2: move Transmission visits to Post except for a few that happen before actual first positive date
mapping$Timepoint2 <- mapping$Timepoint
inds <- which(mapping$Transmitter=="Transmitter" & mapping$Timepoint=="Transmission" & mapping$DaysRelativeToTransmission < 0)
inds2 <- which(mapping$Timepoint=="Transmission")
mapping[inds2, "Timepoint2"] <- "Post"
mapping[inds, "Timepoint2"] <- "Pre"

# 11/8/23: add firstpositive, lastnegative dates, transmission mode for PROMISE transmitters
infectiondate <- read.csv("/Lab_Share/PROMISE/nwcs610/infectiondate.csv")
rownames(infectiondate) <- infectiondate$mpatid
inds <- which(mapping$Study2=="NWCS610" & mapping$Aim1=="Yes" & mapping$Transmitter=="Transmitter")
pids.transmitter <- unique(mapping[inds, "Patient.ID"])
mapping[inds, "fpdate"] <- format(parse_date_time(as.character(infectiondate[as.character(mapping[inds, "Patient.ID"]), "firstpositive"]), order="dmy"), format="%Y-%m-%d")
mapping[inds, "lndate"] <- format(parse_date_time(as.character(infectiondate[as.character(mapping[inds, "Patient.ID"]), "lastnegative"]), order="dmy"), format="%Y-%m-%d")
mapping[inds, "lnage"] <- as.numeric(parse_date_time(mapping[inds, "lndate"], order="ymd") - parse_date_time(mapping[inds, "deldtup2"], order="ymd"))
mapping[inds, "infected_epoch"] <- infectiondate[as.character(mapping[inds, "Patient.ID"]), "infected_epoch"]

# 11/9/23: add ModeOfTransmission based on manual definitions; DaysRelativeToTransmissionMidpoint based on midpoint between lnage/fpage
transmitters.emt <- c("6035760", "6017600", "3037265", "3036777", "3024325")
transmitters.bmt <- c("6034850", "6028820", "6025580", "6021280", "6021000", "6017980", "6012110", "6007720", "6005300", "3035639", "3032974", "3025050", "3024359", "6009290", "6026970") # 6009290 and 6026970 will be filtered out as they are from India/Tanzania
mapping$ModeOfTransmission <- ifelse(mapping$Transmitter=="Control", "Control", NA)
i <- which(mapping$Patient.ID %in% transmitters.emt)
mapping[i, "ModeOfTransmission"] <- "EarlyMucosal"
i <- which(mapping$Patient.ID %in% transmitters.bmt)
mapping[i, "ModeOfTransmission"] <- "BreastMilk"
i <- which(mapping$Patient.ID %in% setdiff(pids.transmitter, c(transmitters.emt, transmitters.bmt)))
mapping[i, "ModeOfTransmission"] <- "Unknown"

lndate <- parse_date_time(as.character(mapping$lndate), order="ymd")
fpdate <- parse_date_time(ifelse(mapping$Study2=="NWCS610", as.character(mapping$statdt2), as.character(mapping$fpdate)), order="ymd")
mapping$DaysRelativeToTransmissionMidpoint <- (as.numeric(parse_date_time(as.character(mapping$TimePoint.Date), order="ymd") - lndate, units="days") + as.numeric(parse_date_time(as.character(mapping$TimePoint.Date), order="ymd") - fpdate, units="days")) / 2


#write.table(mapping, file="/Lab_Share/PROMISE/nwcs610/combined_nwcs610_ZEBS_Mapping.with_metadata.080922.txt", quote=F, sep="\t", row.names=T, col.names=T)
## 04/19/23 - add timing of last negative and first positive visits for ZEBS [lndate, lnvis, lnage, fpdate, fpvis, fpage]
#write.table(mapping, file="/Lab_Share/PROMISE/nwcs610/combined_nwcs610_ZEBS_Mapping.with_metadata.041923.txt", quote=F, sep="\t", row.names=T, col.names=T)
# 09/28/23 - add ARV regimen at time of sample collection (ARVregimen) and first regimen (firstregimen)
#write.table(mapping, file="/Lab_Share/PROMISE/nwcs610/combined_nwcs610_ZEBS_Mapping.with_metadata.092823.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(mapping, file="/Lab_Share/PROMISE/nwcs610/combined_nwcs610_ZEBS_Mapping.with_metadata.101923.txt", quote=F, sep="\t", row.names=T, col.names=T)

output_dir <- "/Lab_Share/PROMISE/nwcs610"
## Table 1 (overall, separately for PROMISE/ZEBS) - Aim 1
# overall
mapping.sel <- subset(mapping, Aim1=="Yes" & Study %in% c("NWCS610-Aim1", "NWCS610-Aim1,NWCS610-Aim2", "ZEBS"))
mapping.sel$country <- droplevels(mapping.sel$country)
demo_vars <- c("Study", "Transmitter","country", "instn", "firstregimen", "cd4bl", "log10vlrna", "delgage", "parity")
demo <- ddply(mapping.sel, .(Patient.ID), function(x) {
	retval <- {}
	for (demo_var in demo_vars) {
		retval <- c(retval, paste(unique(x[!is.na(x[,demo_var]),demo_var]), collapse=";"))
	}
	retval
})
colnames(demo) <- c("Patient.ID", demo_vars)
for (demo_var in c("cd4bl", "log10vlrna", "delgage", "parity")) {
	demo[,demo_var] <- as.numeric(demo[,demo_var])
}
strata_var <- "Transmitter"
write.table(print(CreateTableOne(vars=demo_vars, strata=strata_var, data=demo, smd=T)), file=sprintf("/Lab_Share/PROMISE/nwcs610/Table_1.%s.overall.txt", "Aim1"), quote=F, sep="\t", row.names=T, col.names=T)
# PROMISE
mapping.sel <- subset(mapping, Aim1=="Yes" & Study %in% c("NWCS610-Aim1", "NWCS610-Aim1,NWCS610-Aim2"))
mapping.sel$country <- droplevels(mapping.sel$country)
demo <- ddply(mapping.sel, .(Patient.ID), function(x) {
	retval <- {}
	for (demo_var in demo_vars) {
		retval <- c(retval, paste(unique(x[!is.na(x[,demo_var]),demo_var]), collapse=";"))
	}
	retval
})
colnames(demo) <- c("Patient.ID", demo_vars)
for (demo_var in c("cd4bl", "log10vlrna", "delgage", "parity")) {
	demo[,demo_var] <- as.numeric(demo[,demo_var])
}
strata_var <- "Transmitter"
write.table(print(CreateTableOne(vars=demo_vars, strata=strata_var, data=demo, smd=T)), file=sprintf("/Lab_Share/PROMISE/nwcs610/Table_1.%s.NWCS610.txt", "Aim1"), quote=F, sep="\t", row.names=T, col.names=T)

# PROMISE (by Sample)
aim <- "Aim1"
sample_vars <- c("age", "ARVregimen")
write.table(print(CreateTableOne(vars=sample_vars, strata=c("Timepoint", "Transmitter"), data=mapping.sel, smd=T), noSpaces=T), file=sprintf("%s/Table_1.%s.by_Sample.txt", output_dir, aim), quote=F, sep="\t", row.names=T, col.names=T)
for (tp in unique(mapping.sel$Timepoint)) {
	write.table(print(CreateTableOne(vars=sample_vars, strata=c("Transmitter"), data=subset(mapping.sel, Timepoint==tp), smd=T), noSpaces=T), file=sprintf("%s/Table_1.%s.by_Sample.%s.txt", output_dir, aim, tp), quote=F, sep="\t", row.names=T, col.names=T)
}


# ZEBS
mapping.sel <- subset(mapping, Aim1=="Yes" & Study %in% c("ZEBS"))
demo_vars <- c("country", "instn", "cd4bl", "vlrna", "log10vlrna", "delgage", "MODE", "parity", "lnage", "fpage", "HGB", "CD8C", "CD3C", "PREVPREG", "BWT", "SEX", "CD4_12M", "CD8_12M", "CD3_12M", "rna4m", "rna12m", "cd", "md", "timec", "timem", "Q2C1", "Q2C2", "Q2C3", "bf", "bftime", "copies0_0valid", "copies0_0undetectable", "copies0_0log10", "copies1_0valid", "copies1_0undetectable", "copies1_0log10", "copies4_0valid", "copies4_0undetectable", "copies4_0log10", "copies4_5valid", "copies4_5undetectable", "copies4_5log10", "ebf0_0", "ebf1_0", "ebf4_0", "ebf4_5")
demo <- unique(mapping.sel[,c("Patient.ID", "Study", "Transmitter", demo_vars)])
strata_var <- "Transmitter"
write.table(print(CreateTableOne(vars=demo_vars, strata=strata_var, data=demo, smd=T)), file=sprintf("/Lab_Share/PROMISE/nwcs610/Table_1.%s.ZEBS.txt", "Aim1"), quote=F, sep="\t", row.names=T, col.names=T)

write.table(demo, file="/Lab_Share/PROMISE/nwcs610/ZEBS_subject_data.060623.txt", quote=F, sep="\t", row.names=T, col.names=T)

## Table 1 (overall) - Aim 2
# overall
mapping.sel <- subset(mapping, Aim2=="Yes" & Study %in% c("NWCS610-Aim2", "NWCS610-Aim1,NWCS610-Aim2"))
mapping.sel$country <- droplevels(mapping.sel$country)
demo_vars <- c("Study", "Regimen","country", "instn", "firstregimen", "cd4bl", "log10vlrna", "delgage", "parity")
demo <- ddply(mapping.sel, .(Patient.ID), function(x) {
	retval <- {}
	for (demo_var in demo_vars) {
		retval <- c(retval, paste(unique(x[!is.na(x[,demo_var]),demo_var]), collapse=";"))
	}
	retval
})
colnames(demo) <- c("Patient.ID", demo_vars)
for (demo_var in c("cd4bl", "log10vlrna", "delgage", "parity")) {
	demo[,demo_var] <- as.numeric(demo[,demo_var])
}
strata_var <- "Regimen"
write.table(print(CreateTableOne(vars=demo_vars, strata=strata_var, data=demo, smd=T)), file=sprintf("/Lab_Share/PROMISE/nwcs610/Table_1.%s.overall.txt", "Aim2"), quote=F, sep="\t", row.names=T, col.names=T)


##############################################################################################################
### 11/14/2023: convert SELECTED.Aim1.030419.txt to a more useful (sample-indexed) format

#matching <- read.table("/Lab_Share/PROMISE/nwcs610/SELECTED.Aim1.030419.txt", header=T, as.is=T, sep="\t")
#tmp <- matching[, c("mpatid", "Timepoints", "Match1", "Timepoints1", "Match2", "Timepoints2")]
#out <- {}
#for (i in 1:nrow(tmp)) {
#	tp.transmitter <- parse_date_time(unlist(strsplit(tmp[i, "Timepoints"], ",")), order="ymd")
#	tp.match1 <- parse_date_time(unlist(strsplit(tmp[i, "Timepoints1"], ",")), order="ymd")
#	tp.match2 <- parse_date_time(unlist(strsplit(tmp[i, "Timepoints2"], ",")), order="ymd")
#	
#	sid.transmitter <- sprintf("P%s.%s", tmp[i, "mpatid"], format(tp.transmitter, format="%Y%m%d"))
#	sid.match1 <- sprintf("P%s.%s", tmp[i, "Match1"], format(tp.match1, format="%Y%m%d"))
#	sid.match2 <- sprintf("P%s.%s", tmp[i, "Match2"], format(tp.match2, format="%Y%m%d"))
#	
##	add <- cbind(sid.transmitter, sid.match1, sid.match2)
#	print(sid.transmitter)
#	print(sid.match1)
#	print(sid.match2)
#}










