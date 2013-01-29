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
#' replacenames(alfarotree, source_="NCBI", splitby=100)
#' }
replacenames <- function(phylo, source_ = "NCBI", splitby = 100){
	tree_tips <- str_replace_all(phylo$tip.label, "_", " ") # replace underscores
	registerDoMC(cores=4)
	out <- ldply(slice(tree_tips, splitby), function(x) tnrs(x, getpost="POST")[,1:4], .parallel=TRUE) # get TNRS data
	
	temp <- out[out$sourceId %in% source_, ]
	
	foo2 <- function(x){  # function to grab submitted name if matched or accepted if not
		if(x$submittedName %in% x$acceptedName){
			return(as.character(x$submittedName))
		} else
		{
			return(as.character(x$acceptedName))
		}
	} 
	registerDoMC(4)
	temp2 <- ddply(temp, .(submittedName), foo2, .parallel=TRUE)
	
	nosourcematch <- unique(out$submittedName[!out$submittedName %in% temp$submittedName]) # no match to source
	temp22 <- rbind(temp2, data.frame(submittedName=nosourcematch, V1=nosourcematch)) # add in species for which there was no match
	
	# replace spaces with underscores in both columns
	temp22$submittedName <- str_replace_all(temp22$submittedName, " ", "_")
	temp22$V1 <- str_replace_all(temp22$V1, " ", "_")
	
	notnrsmatch <- phylo$tip.label[!phylo$tip.label %in% as.character(temp22$submittedName)] # no match at all
	temp33 <- rbind(temp22, data.frame(submittedName=notnrsmatch, V1=notnrsmatch)) # add notnrsmatches to data.frame
	
# 	order_ <- sapply(temp$sub, function(x) grep(x, temp22$submittedName))
	
	order_ <- sapply(phylo$tip.label, function(x) grep(x, as.character(temp33$submittedName)))
	
	temp3 <- temp33[order_,]
	phylo$tip.label <- temp3$V1
	return(phylo)
}

# parallelize = machine has 4 cores - may be able to add more cores
# How much will we be able to parallelize taxosaurus calls?
library(doMC)
registerDoMC(cores=4)


# 