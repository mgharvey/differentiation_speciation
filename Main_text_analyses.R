setwd("~/Differentiation_Speciation/")
getwd()

# Script: Main_text_analyses.r 
# By: Michael G. Harvey 
# Date: 20 February 2017

# A script to replicate tests comparing population divergence to phylogenetic 
# speciation rates presented in main text.

require(BAMMtools)
require(caper)
require(scales)
require(laser)

##################################################

# Get all the data:

##################################################

# Get the taxonomy key relating tips in tree to names used in population genetic datasets:
key <- read.table("taxonomy_key.txt", sep="\t", header=T)

# Get the Jetz et al. tree based on Hackett et al. backbone (w/ lumped species collapsed and species with genetic data only)
birdtree <- read.tree('Hackett_lumped_justgeneticdata.tre')

# Get event data from a BAMM run on the tree above
#ed <- getEventData(birdtree, 'event_data_lumped_hack_vr.txt', burnin=0.1, nsamples=1250) 
#save(ed, file="Hackett_lumped_eventsample.rda")
load("Hackett_lumped_eventsample.rda")

# Subset the event data for only our study species, and get the speciation rates
lumped_subb <- subtreeBAMM(ed, tips = birdtree$tip.label[birdtree$tip.label %in% key$Jetz_Species])
BAMM_spec_rate <- getTipRates(lumped_subb)$lambda.avg

# Get the lumped taxonomy metadata
lumped <- read.table("lumped_data_final.txt", sep="\t", header=T, row.names=2)

# Subset the key for the lumped taxonomy
key.used <- key[key$Jetz_Species %in% names(BAMM_spec_rate),]
key.used <- subset(key.used, !duplicated(key.used$Lumped))

##################################################

# Analysis using BAMM speciation rates and STRAPP:

##################################################

# Calculate population differentiation rate
diff_rate <- log(lumped$Ngroups0.8+lumped$Ngroups.excluded0.8)/lumped$crown.age
names(diff_rate) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]

# Test the correlation
res <- traitDependentBAMM(lumped_subb, traits=diff_rate, method="spearman")

# The observed correlation
res$estimate

# The two-tailed P value
res$p.value

# Relative rates
BAMMrelrate_overall <- sum(diff_rate)/sum(BAMM_spec_rate)
BAMMrelrates <- diff_rate/BAMM_spec_rate[names(diff_rate)]

plot(log(BAMM_spec_rate[names(diff_rate)])~log(diff_rate))
dframe <- data.frame(names(diff_rate), log(BAMM_spec_rate[names(diff_rate)]), log(diff_rate))
colnames(dframe) <- c("Species", "Speciation", "PopDiff")
dframe <- dframe[dframe$PopDiff!="-Inf",]
res <- lm(Speciation~PopDiff, data=dframe)
abline(res)

##################################################

# PGLS using DR statistic:

##################################################

# Get the DR statistics for all species from existing file
DR <- read.table("DR_mean_rates.txt", sep="\t", header=F)
DR_values_all <- DR[,2]
names(DR_values_all) <- DR[,1]
DR_values <- DR_values_all[names(BAMM_spec_rate)]

# Test the correlation
dframe <- data.frame(names(DR_values), log(DR_values), log(diff_rate[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation", "PopDiff")
dframe <- dframe[dframe$PopDiff!="-Inf",]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ PopDiff, data=data)
sum <- summary(test)

# The observed correlations:
sum$coefficients[2,1]

# The P value:
sum$coefficients[2,4]

# Relative rates
DRrelrate_overall <- sum(diff_rate)/sum(DR_values)
DRrelrates <- diff_rate/DR_values[names(diff_rate)]

# DR statistic vs. BAMM speciation rates
res <- lm(log(DR_values)~log(BAMM_spec_rate[names(DR_values)]))
summary(res)

plot(Speciation~PopDiff, data=dframe)
abline(test)

# Correlation between BAMM speciation rate and DR values
dframe <- data.frame(names(DR_values), log(DR_values), log(BAMM_spec_rate[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation_DR", "Speciation_BAMM")
test <- lm(Speciation_BAMM ~ Speciation_DR, data=dframe)
sum <- summary(test)

##################################################

# Tests of false positive rates:

##################################################

SIMS <- 1000
vv <- vcv.phylo(as.phylo(lumped_subb))
pvalues_bammcorr <- numeric(SIMS) # pvalues from BAMM permutation test
pvalues_uncorrected <- numeric(SIMS) # pvalues from simple Spearman test
pvalues_uncorrected2 <- numeric(SIMS) # pvalues from simple Spearman test with DR stat
pvalues_pgls <- numeric(SIMS) # pvalues from PGLS
spec_BAMM <- BAMM_spec_rate[rownames(vv)]
spec_DR <- DR_values[rownames(vv)]

for (i in 1:SIMS){
	cat(i, '\n')
	traits <- rmvnorm(1, sigma=vv)[1,]
	names(traits) <- rownames(vv)
	tmp_bamm <- traitDependentBAMM(lumped_subb, traits)
	tmp_spear <- cor.test(spec_BAMM, traits)
	tmp_spear2 <- cor.test(spec_DR, traits)
	dframe <- data.frame(names(spec_DR), traits, spec_DR)
	data <- comparative.data(data=dframe, phy=birdtree, names.col="names.spec_DR.")
	tmp_pgls <- pgls(spec_DR ~ traits, data=data)
	pvalues_bammcorr[i] <- tmp_bamm$p.value
	pvalues_pgls[i] <- summary(tmp_pgls)$coefficients[2,4]
	pvalues_uncorrected[i] <- tmp_spear$p.value
	pvalues_uncorrected2[i] <- tmp_spear2$p.value;
}

# The effective Type I error rate from the BAMM permutation test:
sum(pvalues_bammcorr <= 0.05)/SIMS

# The effective Type I error rate from the PGLS:
sum(pvalues_pgls <= 0.05)/SIMS

# The effective Type I error rate from ignoring the BAMM covariance structure in rates:
sum(pvalues_uncorrected <= 0.05)/SIMS

# The effective Type I error rate from ignoring the BAMM covariance structure in rates:
sum(pvalues_uncorrected2 <= 0.05)/SIMS

##################################################

# Tests with raw number of bGMYC clusters:

##################################################

# Pull out number of bGMYC clusters
xvals <- lumped$Ngroups0.8+lumped$Ngroups.excluded0.8
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]

# BAMM/STRAPP
res <- traitDependentBAMM(lumped_subb, traits=xvals, method="spearman")
res$estimate
res$p.value

# DR/PGLS
dframe <- data.frame(names(DR_values), log(DR_values), log(xvals[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation", "XVals")
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ XVals, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

##################################################

# Sensitivity tests:

##################################################

# Split taxonomy
splittree <- read.tree('Hackett_split_justgeneticdata.tre')
#ed <- getEventData(splittree, 'event_data_split_hack_vr.txt', burnin=0.1, nsamples=1250) 
#save(ed, file="Hackett_split_eventsample.rda")
load("Hackett_split_eventsample.rda")
split_subb <- subtreeBAMM(ed, tips = splittree$tip.label[splittree$tip.label %in% key$Jetz_Species])
BAMM_split_spec_rate <- getTipRates(split_subb)$lambda.avg
split <- read.table("split_data_final.txt", sep="\t", header=T, row.names=1)
key.used.split <- key[key$Jetz_Species %in% names(BAMM_split_spec_rate),]
key.used.split <- subset(key.used.split, !duplicated(key.used$Split))
key.used.split <- subset(key.used.split, !duplicated(key.used$Jetz_Species))
# Calculate population differentiation rate
split_diff_rate <- log(split$Ngroups0.8total)/split$crown.age
names(split_diff_rate) <- key.used.split$Jetz_Species[match(rownames(split), key.used.split$Split)]
res <- traitDependentBAMM(split_subb, traits=split_diff_rate, method="spearman")
res$estimate
res$p.value

# Time threshold taxonomy
timetree <- read.tree('Hackett_time_justgeneticdata.tre')
#ed <- getEventData(timetree, 'event_data_time_hack_vr.txt', burnin=0.1, nsamples=1250) 
#save(ed, file="Hackett_time_eventsample.rda")
load("Hackett_time_eventsample.rda")
time_subb <- subtreeBAMM(ed, tips = timetree$tip.label[timetree$tip.label %in% key$Jetz_Species])
BAMM_time_spec_rate <- getTipRates(time_subb)$lambda.avg
key.used.time <- key[key$Jetz_Species %in% names(BAMM_time_spec_rate),]
key.used.time <- subset(key.used.time, !duplicated(key.used$Lumped))
time <- read.table("lumped_data_final.txt", sep="\t", header=T, row.names=2)
time <- time[match(rownames(time), key.used.time$Lumped),]
# Calculate population differentiation rate
time_diff_rate <- log(time$Ngroups0.8+time$Ngroups.excluded0.8)/time$crown.age
names(time_diff_rate) <- key.used.time$Jetz_Species[match(rownames(time), key.used.time$Lumped)]
res <- traitDependentBAMM(time_subb, traits=time_diff_rate, method="spearman")
res$estimate
res$p.value

# Lower (0.7) PP threshold for clustering
xvals <- log(lumped$Ngroups0.7+lumped$Ngroups.excluded0.7)/lumped$crown.age
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
res <- traitDependentBAMM(lumped_subb, traits=xvals, method="spearman")
res$estimate
res$p.value

# Higher (0.9) PP threshold for clustering
xvals <- log(lumped$Ngroups0.9+lumped$Ngroups.excluded0.9)/lumped$crown.age
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
res <- traitDependentBAMM(lumped_subb, traits=xvals, method="spearman")
res$estimate
res$p.value

# Stem age rate
xvals <- log(lumped$Ngroups0.8+lumped$Ngroups.excluded0.8)/lumped$stem.age
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
res <- traitDependentBAMM(lumped_subb, traits=xvals, method="spearman")
res$estimate
res$p.value

# Removal of 20% of samples
xvals <- log(lumped$Ngroups0.8_w0.2removed)/lumped$crown.age_0.2removed
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
res <- traitDependentBAMM(lumped_subb, traits=xvals, method="spearman")
res$estimate
res$p.value

# Removal of 40% of samples (not significant)
xvals <- log(lumped$Ngroups0.8_w0.4removed)/lumped$crown.age_0.4removed
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
res <- traitDependentBAMM(lumped_subb, traits=xvals, method="spearman")
res$estimate
res$p.value

# Removal of single-individual populations
xvals <- log(lumped$Ngroups0.8)/lumped$crown.age
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
res <- traitDependentBAMM(lumped_subb, traits=xvals, method="spearman")
res$estimate
res$p.value

# Population differentiation rate with moderate (eps=0.45) extinction
xvals <- lambda.stem.ms01((lumped$Ngroups0.8+lumped$Ngroups.excluded0.8), lumped$crown.age, eps = 0.45)$lambda
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
res <- traitDependentBAMM(lumped_subb, traits=xvals, method="spearman")
res$estimate
res$p.value

# Population differentiation rate with high (eps=0.90) extinction
xvals <- lambda.stem.ms01((lumped$Ngroups0.8+lumped$Ngroups.excluded0.8), lumped$crown.age, eps = 0.9)$lambda
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
res <- traitDependentBAMM(lumped_subb, traits=xvals, method="spearman")
res$estimate
res$p.value

# Crown age vs. speciation rate
xvals <- lumped$crown.age
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
res <- traitDependentBAMM(lumped_subb, traits=xvals, method="spearman")
res$estimate
res$p.value

##################################################

# Multivariate analysis:

##################################################

order <- as.vector(key$Lumped[match(names(DR_values), key$Jetz_Species)])
lumped2 <- lumped[match(order, rownames(lumped)),]
dframe <- data.frame(names(DR_values), log(DR_values), log(diff_rate[names(DR_values)]), abs(lumped2$lat.mid), lumped2$PC1, lumped2$PC2, log(lumped2$PC1_var), lumped2$PC2_var, log(lumped2$range.size), log(lumped2$MeanTarsusLength), log(lumped2$Tarsus_var), log(lumped2$MeanHandWingIndex), lumped2$Kipps_var, lumped2$Dichromatism_presence, lumped2$migratory.distance)
colnames(dframe) <- c("Species", "Speciation", "PopDiff", "lat", "pc1", "pc2", "pc1v", "pc2v", "rangesize", "tarsus", "tarsusv", "kipps", "kippsv", "dichrom_pa", "migration")
dframe <- dframe[dframe$PopDiff!="-Inf",]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
nrow(dframe)
test1 <- pgls(Speciation ~ PopDiff + rangesize + lat + migration + pc1 + pc2 + tarsus + kipps + dichrom_pa, data=data)
test2 <- pgls(Speciation ~ PopDiff + lat + migration + pc1 + pc2 + tarsus + kipps + dichrom_pa, data=data)
test3 <- pgls(Speciation ~ PopDiff + rangesize + migration + pc1 + pc2 + tarsus + kipps + dichrom_pa, data=data)
test4 <- pgls(Speciation ~ PopDiff + rangesize + lat + pc1 + pc2 + tarsus + kipps + dichrom_pa, data=data)
test5 <- pgls(Speciation ~ PopDiff + rangesize + lat + migration + pc2 + tarsus + kipps + dichrom_pa, data=data)
test6 <- pgls(Speciation ~ PopDiff + rangesize + lat + migration + pc1 + tarsus + kipps + dichrom_pa, data=data)
test7 <- pgls(Speciation ~ PopDiff + rangesize + lat + migration + pc1 + pc2 + kipps + dichrom_pa, data=data)
test8 <- pgls(Speciation ~ PopDiff + rangesize + lat + migration + pc1 + pc2 + tarsus + dichrom_pa, data=data)
test9 <- pgls(Speciation ~ PopDiff + rangesize + lat + migration + pc1 + pc2 + tarsus + kipps, data=data)
test10 <- pgls(Speciation ~ rangesize + lat + migration + pc1 + pc2 + tarsus + kipps, data=data)
summary(test1)
test <- pgls(Speciation ~ PopDiff, data=data)
summary(test)
AIC(test1, test2, test3, test4, test5, test6, test7, test8, test9, test10)

test10$aicc-test1$aicc
test9$aicc-test1$aicc
test8$aicc-test1$aicc
test7$aicc-test1$aicc
test6$aicc-test1$aicc
test5$aicc-test1$aicc
test4$aicc-test1$aicc
test3$aicc-test1$aicc
test2$aicc-test1$aicc

##################################################

# Tropical vs. temperate analyses:

##################################################

# Temperate correlation
lumpedTemp <- lumped[abs(lumped$lat.mid) >= 23.43706,] # For examining Temp subset
Temp_diff_rate <- log(lumpedTemp$Ngroups0.8+lumpedTemp$Ngroups.excluded0.8)/lumpedTemp$crown.age
names(Temp_diff_rate) <- key.used$Jetz_Species[match(rownames(lumpedTemp), key.used$Lumped)]
Tempres <- traitDependentBAMM(lumped_subb, traits=Temp_diff_rate, method="spearman")
Tempres$estimate
Tempres$p.value

# Tropical correlation
lumpedTrop <- lumped[abs(lumped$lat.mid) < 23.43706,] # For examining Trop subset
Trop_diff_rate <- log(lumpedTrop$Ngroups0.8+lumpedTrop$Ngroups.excluded0.8)/lumpedTrop$crown.age
names(Trop_diff_rate) <- key.used$Jetz_Species[match(rownames(lumpedTrop), key.used$Lumped)]
Tropres <- traitDependentBAMM(lumped_subb, traits=Trop_diff_rate, method="spearman")
Tropres$estimate
Tropres$p.value	

# Using 1000 subsamples of same size as Temperate data set
lumpedTrop <- lumped[abs(lumped$lat.mid) < 23.43706,] # For examining Temp subset
BAMMslopes <- vector()
BAMMpvals <- vector()
for(i in 1:1000) {
	ilumpedTrop <- lumpedTrop[sample(1:nrow(lumpedTrop), nrow(lumpedTemp)),]
	Trop_diff_rate <- log(ilumpedTrop$Ngroups0.8+ilumpedTrop$Ngroups.excluded0.8)/ilumpedTrop$crown.age
	names(Trop_diff_rate) <- key.used$Jetz_Species[match(rownames(ilumpedTrop), key.used$Lumped)]
	res <- traitDependentBAMM(lumped_subb, traits=Trop_diff_rate, method="spearman")
	BAMMslopes <- c(BAMMslopes, res$estimate)
	BAMMpvals <- c(BAMMpvals, res$p.value)	
}
length(BAMMslopes[BAMMslopes < Tempres$estimate])
length(BAMMpvals[BAMMpvals > Tempres$p.value])

# 1000 permutations among species of groups the same size as the initial analysis (110/63)
lumpedTrop <- lumped[abs(lumped$lat.mid) < 23.43706,] # For examining Temp subset
Tempslopes <- vector()
Tropslopes <- vector()
Temppvals <- vector()
Troppvals <- vector()
for(i in 1:1000) {
	ilumpedTemp <- lumped[sample(1:nrow(lumped), nrow(lumpedTemp)),]
	ilumpedTrop <- lumped[sample(1:nrow(lumped), nrow(lumpedTrop)),]
	iTemp_diff_rate <- log(ilumpedTemp$Ngroups0.8+ilumpedTemp$Ngroups.excluded0.8)/ilumpedTemp$crown.age
	iTrop_diff_rate <- log(ilumpedTrop$Ngroups0.8+ilumpedTrop$Ngroups.excluded0.8)/ilumpedTrop$crown.age
	names(iTemp_diff_rate) <- key.used$Jetz_Species[match(rownames(ilumpedTemp), key.used$Lumped)]
	names(iTrop_diff_rate) <- key.used$Jetz_Species[match(rownames(ilumpedTrop), key.used$Lumped)]
	iTempres <- traitDependentBAMM(lumped_subb, traits=iTemp_diff_rate, method="spearman")
	iTropres <- traitDependentBAMM(lumped_subb, traits=iTrop_diff_rate, method="spearman")
	Tempslopes <- c(Tempslopes, iTempres$estimate)
	Tropslopes <- c(Tropslopes, iTropres$estimate)
	Temppvals <- c(Temppvals, iTempres$p.value)	
	Troppvals <- c(Troppvals, iTropres$p.value)	
	print(i)
}
mean(Tempslopes)
mean(Tropslopes)
mean(Temppvals)
mean(Troppvals)
sloperatios <- (Tropslopes/Tempslopes)
length(sloperatios[sloperatios > (Tropres$estimate/Tempres$estimate)])

# PGLS of DR statistic
lat <- abs(lumped$lat.mid)
names(lat) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
dframe <- data.frame(names(DR_values), log(DR_values), log(diff_rate[names(DR_values)]), lat[names(DR_values)])
colnames(dframe) <- c("Species", "Speciation", "PopDiff", "Lat")
preTempdframe <- dframe[abs(dframe$Lat) >= 23.43706,]
preTropdframe <- dframe[abs(dframe$Lat) < 23.43706,]
Tempdframe <- preTempdframe[preTempdframe$PopDiff!="-Inf",]
Tropdframe <- preTropdframe[preTropdframe$PopDiff!="-Inf",]
Tempdata <- comparative.data(data=Tempdframe, phy=birdtree, names.col="Species")
Tropdata <- comparative.data(data=Tropdframe, phy=birdtree, names.col="Species")
Temptest <- pgls(Speciation ~ PopDiff, data=Tempdata)
Troptest <- pgls(Speciation ~ PopDiff, data=Tropdata)
Tempsum <- summary(Temptest)
Tropsum <- summary(Troptest)
Tempsum$coefficients[2,1]
Tropsum$coefficients[2,1]
Tempsum$coefficients[2,4]
Tropsum$coefficients[2,4]

# PGLS of DR statistic resampling the tropical data
DRslopes <- vector()
DRpvals <- vector()
for(i in 1:1000) {
	iTropdframe <- preTropdframe[sample(1:nrow(preTropdframe), nrow(preTempdframe)),]
	iTropdframe <- iTropdframe[iTropdframe$PopDiff!="-Inf",]
	iTropdata <- comparative.data(data=iTropdframe, phy=birdtree, names.col="Species")
	iTroptest <- pgls(Speciation ~ PopDiff, data=iTropdata)
	iTropsum <- summary(iTroptest)
	DRslopes <- c(DRslopes, iTropsum$coefficients[2,1])
	DRpvals <- c(DRpvals, iTropsum$coefficients[2,4])
}
length(DRslopes[DRslopes < Tempsum$coefficients[2,1]])
length(DRpvals[DRpvals > Tempsum$coefficients[2,4]])

# 1000 permutations among species of groups the same size as the initial analysis (110/63) for DR statistic/PGLS
Tempslopes <- vector()
Tropslopes <- vector()
Temppvals <- vector()
Troppvals <- vector()
for(i in 1:1000) {
	iTempdframe <- dframe[sample(1:nrow(dframe), nrow(preTempdframe)),]
	iTropdframe <- dframe[sample(1:nrow(dframe), nrow(preTropdframe)),]
	iTempdata <- comparative.data(data=iTempdframe, phy=birdtree, names.col="Species")
	iTropdata <- comparative.data(data=iTropdframe, phy=birdtree, names.col="Species")
	iTemptest <- pgls(Speciation ~ PopDiff, data=iTempdata)
	iTroptest <- pgls(Speciation ~ PopDiff, data=iTropdata)
	iTempsum <- summary(iTemptest)
	iTropsum <- summary(iTroptest)
	Tempslopes <- c(Tempslopes, iTempsum$coefficients[2,1])
	Tropslopes <- c(Tropslopes, iTropsum$coefficients[2,1])
	Temppvals <- c(Temppvals, iTempsum$coefficients[2,4])
	Troppvals <- c(Troppvals, iTropsum$coefficients[2,4])
	print(i)
}
mean(Tempslopes)
mean(Tropslopes)
mean(Temppvals)
mean(Troppvals)
sloperatios <- (Tropslopes/Tempslopes)
length(sloperatios[sloperatios > (Tropsum$coefficients[2,1]/Tempsum$coefficients[2,1])])
pvalratios <- (Troppvals/Temppvals)
length(pvalratios[pvalratios < (Tropsum$coefficients[2,4]/Tempsum$coefficients[2,4])])

# Differentiation and speciation rates vs. latitude
lat <- abs(lumped$lat.mid)
names(lat) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
dframe <- data.frame(names(DR_values), log(BAMM_spec_rate[names(DR_values)]), log(DR_values), log(diff_rate[names(DR_values)]), log(lat[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation_BAMM", "Speciation_DR", "PopDiff", "Lat")
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation_BAMM ~ Lat, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]
test <- pgls(Speciation_DR ~ Lat, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]
dframe <- dframe[dframe$PopDiff!="-Inf",]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(PopDiff ~ Lat, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

# Relative rate versus latitude
dframe <- data.frame(names(DR_values), log(diff_rate[names(DR_values)]/BAMM_spec_rate[names(DR_values)]), log(diff_rate[names(DR_values)]/DR_values), lat[names(DR_values)])
colnames(dframe) <- c("Species", "Relrates_BAMM", "Relrates_DR", "Lat")
dframe1 <- dframe[dframe$Relrates_BAMM!="-Inf",]
data <- comparative.data(data=dframe1, phy=birdtree, names.col="Species")
test <- pgls(Relrates_BAMM ~ Lat, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]
plot(Relrates_BAMM ~ Lat, data=dframe1)
dframe2 <- dframe[dframe$Relrates_DR!="-Inf",]
data <- comparative.data(data=dframe2, phy=birdtree, names.col="Species")
test <- pgls(Relrates_DR ~ Lat, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]
plot(Relrates_DR ~ Lat, data=dframe2)

# Compare relative rate distribution between temperate and tropical species
lumpedTemp <- lumped[abs(lumped$lat.mid) >= 23.43706,] # For examining Temp subset
lumpedTrop <- lumped[abs(lumped$lat.mid) < 23.43706,] # For examining Temp subset
Temp_diff_rate <- log(lumpedTemp$Ngroups0.8+lumpedTemp$Ngroups.excluded0.8)/lumpedTemp$crown.age
Trop_diff_rate <- log(lumpedTrop$Ngroups0.8+lumpedTrop$Ngroups.excluded0.8)/lumpedTrop$crown.age
names(Temp_diff_rate) <- key.used$Jetz_Species[match(rownames(lumpedTemp), key.used$Lumped)]
names(Trop_diff_rate) <- key.used$Jetz_Species[match(rownames(lumpedTrop), key.used$Lumped)]
temp_spec <- BAMM_spec_rate[names(Temp_diff_rate)]
trop_spec <- BAMM_spec_rate[names(Trop_diff_rate)]
Temprelrates <- as.numeric(Temp_diff_rate)/as.numeric(temp_spec)
Troprelrates <- as.numeric(Trop_diff_rate)/as.numeric(trop_spec)
par(mfrow=c(2,1))
hist(log(Temprelrates), xlim=c(-2,5), breaks=14)
hist(log(Troprelrates), xlim=c(-2,5), breaks=7)
var.test(log(Temprelrates), log(Troprelrates))

# Compare relative rate distribution between temperate and tropical species - DR statistic
dframe2 <- data.frame(names(DR_values), DR_values, diff_rate[names(DR_values)], lat[names(DR_values)])
colnames(dframe2) <- c("Species", "Speciation", "PopDiff", "Lat")
Tempdframe2 <- dframe2[abs(dframe2$Lat) >= 23.43706,]
Tropdframe2 <- dframe2[abs(dframe2$Lat) < 23.43706,]
TemprelratesDR <- Tempdframe2$PopDiff/Tempdframe2$Speciation
TroprelratesDR <- Tropdframe2$PopDiff/Tropdframe2$Speciation
par(mfrow=c(2,1))
hist(log(TemprelratesDR), xlim=c(-1,5), breaks=12)
hist(log(TroprelratesDR), xlim=c(-1,5), breaks=12)
var.test(log(TemprelratesDR), log(TroprelratesDR))








