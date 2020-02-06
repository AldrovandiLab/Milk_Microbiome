
library(ggplot2)
library(reshape2)
library(plyr)
library(useful)
library(dplyr)
library(MASS)
library(abind)
library(glmmLasso)
library(parallel)
library(gplots)
library(gridExtra)

## simulation for power estimates for NWCS 624 R01
## grouping variables are treatment arm (PI-ART or ZDV) and delivery (preterm/term) or SGA (yes/no)
## 3 timepoints per subject

## basic strategy from Sean:
## 1) sample metabolite values from multivariate normal (mvrnorm from MASS package) for each subject - multivariate here means the 3 timepoints per person; use some arbitrary covariance matrix for sigma, pick mu1 and mu2 based on proposed effect size (draw from N(d_est, sd(d_est) where d_est is estimated d from pilot data); alternatively pick from mu1==mu2==N(0,1) and assume enough metabolites have non-zero effect size?
## 2) randomly select cases and controls (with actual sample size values)
## 3) calculate actual effect size d for each metabolite based on selected cases/controls
## 4) run glmmLasso to select features
## 5) produce result as [d, selected_as_feature]
## 6) repeat 1-5 1000 times
## final output should be a relationship between effect size d (x-axis) and power (% of permutations in which the feature was selected by glmmLasso)

set.seed(624)
out_pdf <- sprintf("/Lab_Share/PROMISE/nwcs624/power_calculation.%s.%s.pdf", "nwcs624", format(Sys.Date(), "%m%d%y"))
pdf(out_pdf, width=12)

num_iter <- 1000 # number of iterations to simulate

## parameters to modify
sample_sizes <- list("PTB.PI-ART"=c(345,345), "PTB.ZDV"=c(162,162), "SGA.PI-ART"=c(449,449), "SGA.ZDV"=c(244,244))
sample_sizes <- lapply(sample_sizes, function(x) floor(x*0.95)) # assume 5% missingness
#n.case <- 25 # number of case subjects
#n.control <- 25 # number of control subjects
num_timepoints <- 3 # number of timepoints per subject
num_metabolites <- 1000 # number of metabolites measured
cormat <- matrix(0.5, ncol=3, nrow=3); diag(cormat) <- 1 # arbitrary covariance matrix
lambda <- 10

# hyperparameters for draws of the desired effect size (IMPORTANT: this controls the distribution of per-metabolite effect sizes)
mu.d <- 0; sd.d <- 0.2; 
#mu.d <- 0
## draw mu.case and mu.control from a hyperparameterized distribution based on the desired effect size
## these parameters do not change across iterations - we keep drawing actual metabolite values based on these parameters
desired_d <- rnorm(n=num_metabolites, mean=mu.d, sd=sd.d); names(desired_d) <- sprintf("metabolite.%04d", 1:num_metabolites)
mu.control <- rep(0, times=num_metabolites); sd.control <- 1
mu.case <- mu.control + desired_d; sd.case <- 1


#for (sd.d in seq(from=0.2, to=0.2, by=0.02)) {
for (subgroupi in names(sample_sizes)) {
	n.case <- sample_sizes[[subgroupi]][1]
	n.control <- sample_sizes[[subgroupi]][2]
	## PUT LOOP HERE
	res <- {}
	
	# for (i in 1:num_iter) {
	res <- mclapply(1:num_iter, function(i) {
		print(i)
		## sample metabolite values for cases/controls based on desired effect size and calculate observed effect size
		data.case <- abind(lapply(1:n.case, function(x) { do.call(rbind, lapply(1:num_metabolites, function(mi) { mvrnorm(n=1, mu=rep(mu.case[mi], times=num_timepoints), Sigma=cormat)})) }), along=0); dimnames(data.case) <- NULL # creates 3d array [subject, metabolite, timepoint]
		data.control <- abind(lapply(1:n.control, function(x) { do.call(rbind, lapply(1:num_metabolites, function(mi) { mvrnorm(n=1, mu=rep(mu.control[mi], times=num_timepoints), Sigma=cormat)})) }), along=0); dimnames(data.control) <- NULL # creates 3d array [subject, metabolite, timepoint]
		observed_mu.case <- apply(data.case, 2, mean); observed_sd.case <- apply(data.case, 2, sd)
		observed_mu.control <- apply(data.control, 2, mean); observed_sd.control <- apply(data.control, 2, sd)
		observed_d <- (abs(observed_mu.case - observed_mu.control)) / sqrt((observed_sd.case^2 + observed_sd.control^2)/2)
		res.iter <- data.frame(iter_i = i, metabolite = sprintf("metabolite.%04d", 1:num_metabolites), proposed_effect_size = abs(desired_d), observed_effect_size = observed_d)

		## run glmmLasso on simulated data (data.case and data.control)
		melt.case <- melt(data.case); colnames(melt.case) <- c("PID", "metabolite", "timepoint", "value")
		melt.case$metabolite <- sprintf("metabolite.%04d", melt.case$metabolite)
		df.case <- dcast(melt.case, PID + timepoint ~ metabolite, value.var="value")
		df.case$Group <- 1
		melt.control <- melt(data.control); colnames(melt.control) <- c("PID", "metabolite", "timepoint", "value")
		melt.control$metabolite <- sprintf("metabolite.%04d", melt.control$metabolite)
		df.control <- dcast(melt.control, PID + timepoint ~ metabolite, value.var="value")
		df.control$Group <- 0; df.control$PID <- df.control$PID + max(df.case$PID) # offset PIDs to not overlap with cases
		df <- rbind(df.case, df.control)
		df$PID <- factor(df$PID)

		metabolites <- setdiff(colnames(df), c("PID", "timepoint", "Group"))
		form <- paste(metabolites, collapse="+")
		mod <- glmmLasso(as.formula(sprintf("Group ~ %s", form)), rnd = list(PID=~1), family=binomial(), lambda=lambda, data=df)
		selected_metabolites <- setdiff(names(which(coefficients(mod)!=0)), "(Intercept)")

		## record metabolites selected in model
		res.iter$selected <- FALSE
		res.iter[match(selected_metabolites, res.iter$metabolite), "selected"] <- TRUE

		## add to results
	#	res <- rbind(res, res.iter)
		res.iter
	}, mc.cores=16)
	res <- do.call(rbind, res); rownames(res) <- {}

	res <- res[order(res$observed_effect_size),]
	write.table(res, file=sprintf("/Lab_Share/PROMISE/nwcs624/power_calculation.%s.%d.%d.txt", subgroupi, n.case, n.control), quote=F, sep="\t", row.names=F, col.names=T)
	#p <- ggplot(res, aes(x=1:nrow(res), y=observed_effect_size, color=selected)) + geom_point() + theme_classic() + scale_color_brewer(palette="Set1")
	#print(p)

	# power curve
	res2 <- ddply(res, .(metabolite), function(x) {
		length(which(x$selected)) / nrow(x)
	})
	colnames(res2) <- c("metabolite", "fraction_detected"); res2$metabolite <- as.character(res2$metabolite)
	rownames(res2) <- res2$metabolite
	res2$proposed_effect_size <- abs(desired_d[rownames(res2)])
	p <- ggplot(res2, aes(x=proposed_effect_size, y=fraction_detected)) + geom_point() + theme_classic() + ggtitle(sprintf("Fraction detected vs proposed effect size (%s, ncase=%d, ncontrol=%d, mu.d=%.2g, sd.d=%.2g)", subgroupi, n.case, n.control, mu.d, sd.d))
	print(p)
	# loess regression and prediction
	span <- 0.4
	mod <- loess(proposed_effect_size ~ fraction_detected, data=res2, span=span)
	predicted_d <- predict(mod, 0.9) # predicted power d at 90% power
	p <- ggplot(res2, aes(x=proposed_effect_size, y=fraction_detected)) + geom_point() + geom_smooth(method="loess", span=span) + theme_classic() + ggtitle(sprintf("Fraction detected vs proposed effect size (90%% power at d=%.4g) (%s, ncase=%d, ncontrol=%d, mu.d=%.2g, sd.d=%.2g)", predicted_d, subgroupi, n.case, n.control, mu.d, sd.d))
	print(p)
	
	# distribution of proposed effect size desired_d
	df <- melt(abs(desired_d))
	p <- ggplot(df, aes(x=value)) + geom_histogram() + theme_classic() + ggtitle(sprintf("Distribution of proposed effect size"))
	print(p)
	df <- as.matrix(summary(abs(desired_d)))
	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Distribution of proposed effect size")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
	print(p)
	df <- melt(quantile(abs(desired_d), probs=seq(from=0,to=1,by=0.05)))
	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Quantiles of proposed effect size")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
	print(p)
	# distribution of actual effect size
	df <- melt(res$observed_effect_size)
	p <- ggplot(df, aes(x=value)) + geom_histogram() + theme_classic() + ggtitle(sprintf("Distribution of observed effect size (%s)", subgroupi))
	print(p)
	
}

dev.off()



