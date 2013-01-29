# Including TNRS in Datelife
install_github("taxize_", "ropensci", ref="release")
library(taxize)

# testing with one tree
library(ape); library(stringr)
load("/Users/scottmac2/phyloorchard/pkg/branches/data/AlfaroEtAl2009.rda")
# tree <- read.tree(file="~/github/sac/datelife/data/alfaro2009.rda")
alfarotree <- AlfaroEtAl2009_tree
alfarotree_tips <- str_replace_all(alfarotree$tip.label, "_", " ")
out <- tnrs(alfarotree_tips, getpost="POST")[,1:5] # get data

system.time(out <- tnrs(alfarotree_tips, getpost="POST")[,1:5])
system.time(out <- tnrs(alfarotree_tips, getpost="POST", splitby=50, sleep=1)[,1:5])

registerDoMC(cores=4)
system.time( out <- llply(slice(alfarotree_tips, 75), function(x) tnrs(x, getpost="POST")[,1:5], .parallel=TRUE)  )

# check for exact matches
library(doMC)
registerDoMC(cores=4)

#' 
#' @param phylo Phylogeny object, class phylo
#' @param source_ Source to match names to. One of ncbi, iPlant
#' @param splitby Length of chunks by which to split species list
#' @examples \dontrun{
#' replacenames()
#' }
replacenames <- function(phylo, source_ = "NCBI", splitby = 100){
	tree_tips <- str_replace_all(phylo$tip.label, "_", " ") # replace underscores
	registerDoMC(cores=4)
	out <- ldply(slice(tree_tips, splitby), function(x) tnrs(x, getpost="POST")[,1:5], .parallel=TRUE) # get TNRS data
	
	temp <- out[out$sourceId %in% source_, ]
	out$submittedName[!out$submittedName %in% temp$submittedName]
	
	foo2 <- function(x){ 
		if(x$submittedName %in% x$acceptedName){
			return(as.character(x$submittedName))
		} else
		{
			return(as.character(x$acceptedName))
		}
	} 
	registerDoMC(4)
	temp2 <- ddply(temp, .(submittedName), foo2, .parallel=TRUE)
	order_ <- sapply(temp$sub, function(x) grep(x, temp2$submittedName))
	temp3 <- temp2[order_,]
	phylo$tip.label <- temp3
	return(phylo)
}
tomatch <- replacenames(out) # data.frame to match against input data.frame


alfarotree_tips_orig <- alfarotree_tips

alfarotree_tips_orig %in% alfarotree_tips


# parallelize = machine has 4 cores - may be able to add more cores
# How much will we be able to parallelize taxosaurus calls?
library(doMC)
registerDoMC(cores=4)


# 