# Including TNRS in Datelife

#' Function to break up a vector by certain N sized chunks.
#' 
#' @import plyr
#' @param input Input character vector
#' @pararm by Number by which to split the input character vector
slice <- function(input, by = 2) {
	starts <- seq(1, length(input), by)
	tt <- lapply(starts, function(y) input[y:(y + (by - 1))])
	llply(tt, function(x) x[!is.na(x)])
}
	
#' Function to check names on a species list or phylogeny.
#' 
#' @import taxize
#' @param phylo Phylogeny object, class phylo
#' @param source_ Source to match names to. One of ncbi, iPlant
#' @param splitby Length of chunks by which to split species list
#' @examples \dontrun{
#' # A phylogeny as input
#' library(doMC)
#' out <- checknames(obj=BergmannEtAl2012_tree, source_="NCBI", splitby=200)
#' out <- checknames(obj=AlfaroEtAl2009_tree, source_="NCBI", splitby=100)
#' out <- checknames(obj=EastmanEtAlUnpublished_tree, source_="NCBI", splitby=100)
#' out <- checknames(obj=HeathEtAl2012_tree[[1]], source_="NCBI")
#' out <- checknames(obj=JaffeEtAl2011_tree, source_="NCBI", splitby=100)
#' out <- checknames(obj=Chaetodontidae2011[[1]][[1]], source_="NCBI", splitby=100)
#' 
#' # A character vector as input (of species names)
#' out <- checknames(obj=AlfaroEtAl2009_tree$tip.label, source_="NCBI", splitby=50)
#' out <- checknames(obj=AlfaroEtAl2009_tree$tip.label[1:20], source_="NCBI")
#' }
checknames <- function(obj, source_ = "NCBI", splitby = NULL)
{	
	if(class(obj)=="phylo"){
		tree_tips <- str_replace_all(obj$tip.label, "_", " ") # replace underscores
		orig_tips <- obj$tip.label
	} else
	{
		tree_tips <- str_replace_all(obj, "_", " ") # replace underscores
		orig_tips <- tree_tips
	}
	
	if(is.null(splitby)){
		out <- tnrs(tree_tips, getpost="POST")[,1:4] # get TNRS data
	} else
	{
		registerDoMC(cores=4)
		out <- ldply(slice(tree_tips, splitby), function(x) tnrs(x, getpost="POST")[,1:4], .parallel=TRUE) # get TNRS data
	}
	
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
	
	notnrsmatch <- orig_tips[!orig_tips %in% as.character(temp22$submittedName)] # no match at all
	temp33 <- rbind(temp22, data.frame(submittedName=notnrsmatch, V1=notnrsmatch)) # add notnrsmatches to data.frame
	order_ <- sapply(orig_tips, function(x) match(x, temp33$submittedName), USE.NAMES=F)
	temp3 <- temp33[order_,] # reorder data.frame to order of tip.labels
	
	# Keep names that match genera, but not those that don't match genera
	foooo <- function(x){ 
		if(as.character(x$old) %in% as.character(x$new)) { 0 } else { 1 }
	}
	df_ <- data.frame(old=orig_tips, new=temp3$V1)
	df_2 <- ddply(df_, .(old), foooo)
	df_2_merge <- merge(df_, df_2, by="old")
	df_3 <- droplevels(df_2_merge[df_2_merge$V1 == 1, ])
	
	genusmatch <- function(x){
		one <- str_split(x$old, "_")[[1]]
		two <- str_split(x$new, "_")[[1]]
		ifelse(one[[1]] %in% two[[1]], as.character(x$new), as.character(x$old))
	}
	asdf <- adply(df_3, 1, genusmatch) 
	
	if(!nrow(asdf)==0){
		order_2 <- sapply(temp3[temp3$submittedName %in% asdf$old, "submittedName"], function(x) match(x, as.character(asdf$old)), USE.NAMES=F)
		asdf2 <- asdf[order_2,]
		temp3[temp3$submittedName %in% asdf$old, "V1"] <- asdf2$V1
	} else
	{
		temp3 <- temp3
	}
	
	if(class(obj)=="phylo"){
		obj$tip.label <- temp3$V1 # assign new names to tip.labels on phylo object
		return(obj)
	} else
	{
		obj <- temp3$V1
		return(obj)
	}
}

# Check species names against Taxosaurus
# Run and write trees to directory with new file name, just appending "_new" to the end
# Run the function replacenames across all trees
# install_github("taxize_", "ropensci", ref="release")
library(taxize); library(ape); library(stringr)
# treefilenames <- c("AlfaroEtAl2009.rda","BergmannEtAl2012.rda","BinindaEmondsEtAl2007.rda",
# 		"EastmanEtAlUnpublished.rda","HeathEtAl2012.rda","JaffeEtAl2011.rda","megaTree.rda",
# 		"Oaks2011.rda","PyronWiens2011.rda","RaboskyEtAlUnpublished.rda","SantiniEtAl2009.rda",
# 		"SmithEtAl2011.rda","TreeBase.rda","ZhangandWake2009.rda")
treefilenames <- dir("/Users/scottmac2/phyloorchard3")

# Load all trees into the workspace
# CHANGE PATH!!!!
llply(treefilenames, function(x) load(paste("/Users/scottmac2/phyloorchard3/",x,sep=""), .GlobalEnv))

# Get tree names in the workspace
trees <- sapply(treefilenames, function(x) str_replace(x, ".rda", "_tree"), USE.NAMES=F)

# Run em all
checknames_safe <- plyr::failwith(NULL, checknames)
l_ply(trees, checknames_safe)