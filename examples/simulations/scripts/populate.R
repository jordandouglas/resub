if (!require(rjson , quietly = T)) install.packages("rjson")
library(rjson)


BURNIN = 0
NSIMS = 100


getTrees = function(fileName){

	file.in = readLines(fileName, warn=F)
	trees = file.in[grep("tree STATE_", file.in)]
	trees = gsub(".+[=] ", "", trees)
	
	# Translate
	trans = file.in[(grep("Translate", file.in)+1):(grep("^;$", file.in)-1)[1]]
	trans = gsub("\t", "", trans)
	trans = gsub("^ +", "", trans)
	trans = gsub(",", "", trans)
	indexes = sapply(strsplit(trans, " "), function(ele) ele[1])
	labels  = sapply(strsplit(trans, " "), function(ele) ele[2])
	
	
	for (i in 1:length(indexes)){
	
		trees = gsub(paste0("[(]", indexes[i], "[:]"), paste0("(", labels[i], ":"), trees)
		trees = gsub(paste0(",", indexes[i], "[:]"), paste0(",", labels[i], ":"), trees)
		trees = gsub(paste0("[(]", indexes[i], "[[]"), paste0("(", labels[i], "["), trees)
		trees = gsub(paste0(",", indexes[i], "[[]"), paste0(",", labels[i], "["), trees)
	
	}
	
	trees

}




# Read parameters
truth.df = read.table("truth/truth.log", header=T, sep = "\t")

# Read species tree
species.trees = getTrees("truth/truth.trees")



# Postprocessing burnin
burnin.start = floor(BURNIN*nrow(truth.df))
if (burnin.start <= 0) burnin.start = 1
include = seq(from = burnin.start, to = nrow(truth.df), length = NSIMS)
include = floor(include)



# Ensure equal counts of each model (100 counts)
useResub = rep(c(TRUE, FALSE), NSIMS/2)
aaRefine = rep(c(TRUE, TRUE, FALSE, FALSE), NSIMS/4)
aaExpand = rep(c(TRUE, FALSE), each=NSIMS/2)


# Prepare folders
for (simnum in 1:NSIMS){

	rownum = include[simnum]
	
	
	#tmp
	#simnum = simnum - 1

	f = paste0("templates/rep", simnum)
	cat(paste(f, "\n"))
	dir.create(f, showWarnings=F) 
	
	
	
	# Create json
	JSON = list()
	
	
	# Prior variables
	sub.df = truth.df[rownum,]
	for(x in colnames(sub.df)){
		val = sub.df[,x]
		JSON[[x]] = val
	}	
	
	
	# Tree
	JSON[["tree"]] = gsub("&", "&amp;", species.trees[rownum])
	write(JSON[["tree"]], paste0(f, "/trueTree.newick"))


	
	# Branch rates
	branchRates = as.numeric(truth.df[rownum,grep("branchRates", colnames(truth.df))])
	JSON[["branchRates.all"]] = paste0(branchRates, collapse = " ")

	aa.frequencies = as.numeric(truth.df[rownum,grep("aa.frequencies", colnames(truth.df))])
	JSON[["aa.frequencies"]] = paste0(aa.frequencies, collapse = " ")



	ntaxa = truth.df[rownum,"tree.ntaxa"]
	JSON[["ntaxa"]] = ntaxa
	JSON[["taxonRange"]] = paste0(1:ntaxa, collapse=",")



	# Resub state 1 and 2
	cherries = c("WY", "EQ", "IV", "DN", "DK", "FH", "AG", "SG", "PT")
	cherry = sample(cherries, 1)
	JSON[["state1"]] = strsplit(cherry, "")[[1]][1]
	JSON[["state2"]] = strsplit(cherry, "")[[1]][2]


	# Epoch boundary
	#treeHeight = JSON[["tree.height"]]
	#JSON[["transitionHeight"]] = rbeta(1, 2, 2) * treeHeight # Beta(4,4) distribution on tree height

	JSON[["useResub"]] = useResub[simnum]
	JSON[["aaExpand"]] = aaExpand[simnum]
	JSON[["aaRefine"]] = aaRefine[simnum]


	JSON_str = as.character(rjson::toJSON(JSON, indent=1))
	write(JSON_str, paste0(f, "/var.json"))

}









