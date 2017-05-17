setwd("~/Desktop/Feb_2017/")
getwd()

require(ape)

DR_statistic <- function(x, return.mean = FALSE){
	rootnode <- length(x$tip.label) + 1
	sprates <- numeric(length(x$tip.label))
	for (i in 1:length(sprates)){
		node <- i
		index <- 1
		qx <- 0
		while (node != rootnode){
			el <- x$edge.length[x$edge[,2] == node]
			node <- x$edge[,1][x$edge[,2] == node]			
			qx <- qx + el* (1 / 2^(index-1))			
			index <- index + 1
		}
		sprates[i] <- 1/qx
	}
	if (return.mean){
		return(mean(sprates))		
	}else{
		names(sprates) <- x$tip.label
		return(sprates)
	}
}

path = "~/Desktop/Feb_2017/100_complete_trees/"
file.names <- dir(path, pattern =".tre")

outfile.names <- NULL
for (i in 1:100) {
	outfile.names[i] <- paste("~/Desktop/Feb_2017/DR_results_lumped/", i, ".txt", sep="")
}

# Tree used for BAMM analyses
birdtree <- read.tree('Hackett_lumped_justgeneticdata.tre')

# Get a list of tips to drop
key <- read.table("~/Desktop/Feb_2017/taxonomy_key.txt", sep="\t", header=T)
lumped <- read.table("~/Desktop/Feb_2017/lumped_data_final.txt", header=T, sep="\t", row.names=2)
tips_to_drop <- NULL
tips_to_keep <- NULL
for (i in 1:nrow(lumped)) { 
	jetz_nom <- as.character(key$Jetz_Species[!is.na(match(key$Lumped, rownames(lumped)[i]))]) # Get the Jetz species names
	if(length(jetz_nom) > 1) {
		keep <- jetz_nom[jetz_nom %in% birdtree$tip.label] # Keep one tip (the one in the tree used for BAMM analyses)
		tips_to_keep <- c(tips_to_keep, keep) 
		tips_to_drop <- c(tips_to_drop, jetz_nom[jetz_nom != keep]) # Drop the rest
	} else {
		tips_to_keep <- c(tips_to_keep, jetz_nom) 
	}	
}

for (i in 1:length(file.names)) {
	pretree <- read.tree(paste0(path, file.names[i]))	
	unique(tips_to_drop)
	tree <- drop.tip(pretree, tips_to_drop) # Drop tips collapsed in lumped taxonomy
	DR <- DR_statistic(tree)
	write(paste(names(DR), DR, sep="\t"), file=outfile.names[i], sep="\t")
	print(i)
}
