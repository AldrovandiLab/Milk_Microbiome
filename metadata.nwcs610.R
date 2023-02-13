library(ggplot2)
library(reshape2)
library(plyr)
library(useful)

## match up metadata for NWCS 610
mapping <- read.table("/Lab_Share/PROMISE/nwcs610/combined_nwcs610_ZEBS_Mapping.072722.txt", header=T, as.is=T, sep="\t")
rownames(mapping) <- mapping$Sample.ID

# create study variables
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

# grab additional metadata (nwcs610)
metadata <- read.table("/Lab_Share/PROMISE/nwcs610/SELECTED.Aim1.030419.txt", header=T, as.is=T, sep="\t")
mapping[, "deldtup2"] <- metadata[as.character(mapping$Patient.ID), "deldtup2"]
mapping[, "cd4bl"] <- metadata[as.character(mapping$Patient.ID), "cd4bl"]
mapping[, "log10vlrna"] <- log10(metadata[as.character(mapping$Patient.ID), "rnacopbl"])
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
}
mapping$deldtup2 <- data[as.character(mapping$Patient.ID), "deldtup2"]
mapping$country <- data[as.character(mapping$Patient.ID), "country1"]
mapping$ethnicity <- data[as.character(mapping$Patient.ID), "ethnicity"]
mapping$delgage <- data[as.character(mapping$Patient.ID), "delgage"]
mapping$parity <- data[as.character(mapping$Patient.ID), "parity1"]

# set Aim 1 timepoint and grab additional metadata (ZEBS)
metadata <- read.table("/Lab_Share/ZEBS/matched_1to2.Aim1.010821.with_metadata.txt", header=T, as.is=T, sep="\t")
rownames(metadata) <- metadata$studyid
zebs <- rownames(subset(mapping, Study=="ZEBS"))
mvars <- c("CD4C", "log10vlrna")
mapping[zebs, "cd4bl"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "CD4C"]
mapping[zebs, "log10vlrna"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "log10vlrna"]
mapping[zebs, "deldtup2"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "DOB"]
mapping[zebs, "country"] <- "Zambia"
mapping[zebs, "instn"] <- metadata[as.character(mapping[zebs, "Patient.ID"]), "CLINIC"]

#tmp <- melt(metadata[, c("studyid", "Transmitter", "PreDate", "PreAge", "PreVisit", "TransmissionDate", "TransmissionAge", "TransmissionVisit", "Post1Date", "Post1Age", "Post1Visit", "Post2Date", "Post2Age", "Post2Visit")])
tmp <- melt(metadata[, c("studyid",  "PreDate", "TransmissionDate", "Post1Date", "Post2Date")], id.vars=c("studyid"))
tmp$variable <- gsub("[12]", "", gsub("Date", "", as.character(tmp$variable)))
tmp$date <- format(as.Date(tmp$value, format="%Y-%m-%d"), "%Y%m%d")
tmp$SampleID <- sprintf("P%s.%s", tmp$studyid, tmp$date)
mapping[tmp$SampleID, "Timepoint"] <- tmp$variable

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
mapping[rownames(tmp), "Visit"] <- tmp[, "visit"]
mapping[rownames(tmp), "Regimen"] <- tmp[, "regimen"]

sel <- rownames(subset(mapping, Aim2=="Yes" & is.na(Regimen)))
mapping[sel, "Aim2"] <- "No"



##############################################################################################################
# calculate age as needed
mapping$age <- as.numeric(as.Date(as.character(mapping$TimePoint.Date), format="%Y%m%d") - as.Date(mapping$deldtup2, format="%Y-%m-%d"))

# add institution
mapping$instn <- ifelse(is.na(mapping$instn), data[as.character(mapping$Patient.ID), "instn"], mapping$instn)


write.table(mapping, file="/Lab_Share/PROMISE/nwcs610/combined_nwcs610_ZEBS_Mapping.with_metadata.080922.txt", quote=F, sep="\t", row.names=T, col.names=T)

