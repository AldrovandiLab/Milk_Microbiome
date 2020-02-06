
library(ggplot2)
library(reshape2)
library(plyr)
library(useful)
library(lubridate)
library(dplyr)

# various metadata associated tasks for NWCS 624 (R01 for preterm metabolomics)

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


