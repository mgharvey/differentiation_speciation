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
	outfile.names[i] <- paste("~/Desktop/Feb_2017/DR_results_split/", i, ".txt", sep="")
}

for (i in 1:length(file.names)) {
	tree <- read.tree(paste0(path, file.names[i]))	
	DR <- DR_statistic(tree)
	write(paste(names(DR), DR, sep="\t"), file=outfile.names[i], sep="\t")
	print(i)
}
