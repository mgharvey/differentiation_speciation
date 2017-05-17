setwd("~/Desktop/Feb_2017")
getwd()

# Script: Supplementary_analyses.r 
# By: Michael G. Harvey 
# Date: 20 February 2017

# A script to replicate tests comparing population divergence to phylogenetic 
# speciation rates presented in supplementary materials.

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

# Calculate population differentiation rate
diff_rate <- log(lumped$Ngroups0.8+lumped$Ngroups.excluded0.8)/lumped$crown.age
names(diff_rate) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]

# Get the DR statistics for all species from existing file
DR <- read.table("DR_mean_rates.txt", sep="\t", header=F)
DR_values_all <- DR[,2]
names(DR_values_all) <- DR[,1]
DR_values <- DR_values_all[names(BAMM_spec_rate)]

##################################################

# Compare number of differentiated populations to range size and age:

##################################################

res <- lm(log(lumped$Ngroups0.8+lumped$Ngroups.excluded0.8)~log(lumped$range.size))
summary(res)
res <- lm(log(lumped$Ngroups0.8+lumped$Ngroups.excluded0.8)~log(lumped$crown.age))
summary(res)
res <- lm(log(lumped$Ngroups0.8+lumped$Ngroups.excluded0.8)~log(lumped$nSamples))
summary(res)


##################################################

# Compare number of differentiated populations to range size and age:

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



##################################################

# Evaluate a one-tailed version of the STRAPP test:

##################################################

diff_rate <- log(lumped$Ngroups0.8+lumped$Ngroups.excluded0.8)/lumped$crown.age
names(diff_rate) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
# The two-tailed P value
res <- traitDependentBAMM(lumped_subb, traits=diff_rate, method="spearman")
res$p.value
# The one-tailed P value
res <- traitDependentBAMM(lumped_subb, traits=diff_rate, method="spearman", two.tailed=FALSE, traitorder="positive")
res$p.value

##################################################

# Sensitivity analyses with PGLS of DR statistic:

##################################################

# Split taxonomy
splittree <- read.tree('Hackett_split_justgeneticdata.tre')
DR_split <- read.table("DR_split_mean_rates.txt", sep="\t", header=F)
DR_split_values_all <- DR_split[,2]
names(DR_split_values_all) <- DR_split[,1]
load("Hackett_split_eventsample.rda")
split_subb <- subtreeBAMM(ed, tips = splittree$tip.label[splittree$tip.label %in% key$Jetz_Species])
BAMM_split_spec_rate <- getTipRates(split_subb)$lambda.avg
DR_split_values <- DR_split_values_all[names(BAMM_split_spec_rate)]
split <- read.table("split_data_final.txt", sep="\t", header=T, row.names=1)
key.used.split <- key[key$Jetz_Species %in% names(BAMM_split_spec_rate),]
key.used.split <- subset(key.used.split, !duplicated(key.used$Split))
key.used.split <- subset(key.used.split, !duplicated(key.used$Jetz_Species))
split_diff_rate <- log(split$Ngroups0.8total)/split$crown.age
names(split_diff_rate) <- key.used.split$Jetz_Species[match(rownames(split), key.used.split$Split)]
dframe <- data.frame(names(DR_split_values), log(DR_split_values), log(split_diff_rate[names(DR_split_values)]))
colnames(dframe) <- c("Species", "Speciation", "PopDiff")
dframe <- dframe[dframe$PopDiff!="-Inf",]
dframe <- dframe[!is.na(dframe$PopDiff),]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ PopDiff, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

# Time threshold taxonomy
timetree <- read.tree('Hackett_time_justgeneticdata.tre')
DR_time <- read.table("DR_time_mean_rates.txt", sep="\t", header=F)
DR_time_values_all <- DR_time[,2]
names(DR_time_values_all) <- DR_time[,1]
load("Hackett_time_eventsample.rda")
time_subb <- subtreeBAMM(ed, tips = timetree$tip.label[timetree$tip.label %in% key$Jetz_Species])
BAMM_time_spec_rate <- getTipRates(time_subb)$lambda.avg
DR_time_values <- DR_time_values_all[names(BAMM_time_spec_rate)]
key.used.time <- key[key$Jetz_Species %in% names(BAMM_time_spec_rate),]
key.used.time <- subset(key.used.time, !duplicated(key.used$Lumped))
time <- read.table("lumped_data_final.txt", sep="\t", header=T, row.names=2)
time2 <- time[rownames(time) %in% key.used.time$Lumped,]
time_diff_rate <- log(time$Ngroups0.8+time$Ngroups.excluded0.8)/time$crown.age
names(time_diff_rate) <- key.used.time$Jetz_Species[match(rownames(time), key.used.time$Lumped)]
dframe <- data.frame(names(DR_time_values), log(DR_time_values), log(time_diff_rate[names(DR_time_values)]))
colnames(dframe) <- c("Species", "Speciation", "PopDiff")
dframe <- dframe[dframe$PopDiff!="-Inf",]
dframe <- dframe[!is.na(dframe$PopDiff),]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ PopDiff, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

# Lower (0.7) PP threshold for clustering
xvals <- log(lumped$Ngroups0.7+lumped$Ngroups.excluded0.7)/lumped$crown.age
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
dframe <- data.frame(names(DR_values), log(DR_values), log(xvals[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation", "XVals")
dframe <- dframe[dframe$XVals!="-Inf",]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ XVals, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

# Higher (0.9) PP threshold for clustering
xvals <- log(lumped$Ngroups0.9+lumped$Ngroups.excluded0.9)/lumped$crown.age
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
dframe <- data.frame(names(DR_values), log(DR_values), log(xvals[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation", "XVals")
dframe <- dframe[dframe$XVals!="-Inf",]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ XVals, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

# Stem age rate
xvals <- log(lumped$Ngroups0.8+lumped$Ngroups.excluded0.8)/lumped$stem.age
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
dframe <- data.frame(names(DR_values), log(DR_values), log(xvals[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation", "XVals")
dframe <- dframe[dframe$XVals!="-Inf",]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ XVals, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

# Removal of 20% of samples
xvals <- log(lumped$Ngroups0.8_w0.2removed)/lumped$crown.age_0.2removed
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
dframe <- data.frame(names(DR_values), log(DR_values), log(xvals[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation", "XVals")
dframe <- dframe[dframe$XVals!="-Inf",]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ XVals, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

# Removal of 40% of samples 
xvals <- log(lumped$Ngroups0.8_w0.4removed)/lumped$crown.age_0.4removed
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
dframe <- data.frame(names(DR_values), log(DR_values), log(xvals[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation", "XVals")
dframe <- dframe[dframe$XVals!="-Inf",]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ XVals, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

# Removal of single-individual populations
xvals <- log(lumped$Ngroups0.8)/lumped$crown.age
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
dframe <- data.frame(names(DR_values), log(DR_values), log(xvals[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation", "XVals")
dframe <- dframe[dframe$XVals!="-Inf",]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ XVals, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

# Population differentiation rate with moderate (eps=0.45) extinction
xvals <- lambda.stem.ms01((lumped$Ngroups0.8+lumped$Ngroups.excluded0.8), lumped$crown.age, eps = 0.45)$lambda
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
dframe <- data.frame(names(DR_values), log(DR_values), log(xvals[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation", "XVals")
dframe <- dframe[dframe$XVals!="-Inf",]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ XVals, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

# Population differentiation rate with high (eps=0.90) extinction
xvals <- lambda.stem.ms01((lumped$Ngroups0.8+lumped$Ngroups.excluded0.8), lumped$crown.age, eps = 0.9)$lambda
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
dframe <- data.frame(names(DR_values), log(DR_values), log(xvals[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation", "XVals")
dframe <- dframe[dframe$XVals!="-Inf",]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ XVals, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

# Crown age vs. speciation rate
xvals <- lumped$crown.age
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
dframe <- data.frame(names(DR_values), log(DR_values), log(xvals[names(DR_values)]))
colnames(dframe) <- c("Species", "Speciation", "XVals")
dframe <- dframe[dframe$XVals!="-Inf",]
data <- comparative.data(data=dframe, phy=birdtree, names.col="Species")
test <- pgls(Speciation ~ XVals, data=data)
sum <- summary(test)
sum$coefficients[2,1]
sum$coefficients[2,4]

##################################################

# Subspecies analysis:

##################################################

# Subspecies vs. population differentiation
res <- lm(log(lumped$Ngroups0.8+lumped$Ngroups.excluded0.8)~log(lumped$nSubspecies))
summary(res)
plot(log(lumped$Ngroups0.8+lumped$Ngroups.excluded0.8)~log(lumped$nSubspecies))

# BAMM speciation rate versus number of subspecies
xvals <- lumped$nSubspecies
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
res <- traitDependentBAMM(lumped_subb, traits=xvals, method="spearman")
res$estimate
res$p.value

# BAMM speciation rate versus rate of subspecies formation
xvals <- log(lumped$nSubspecies)/lumped$crown.age
names(xvals) <- key.used$Jetz_Species[match(rownames(lumped), key.used$Lumped)]
res <- traitDependentBAMM(lumped_subb, traits=xvals, method="spearman")
res$estimate
res$p.value
