## Input: phenotype list --> Output: lists of children/parent for each non-leaf/root node
#aya43@sfu.ca 20170131

## root directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)



## input directories
meta_dir = paste0(result_dir,"/meta")
meta_cell_dir = paste(meta_dir, "/cell", sep="")

## output directories
meta_cell_child_dir = paste(meta_dir, "/cell_child",sep="") #specifies a phenotypes children
meta_cell_child_ind_dir = paste(meta_dir, "/cell_child_ind",sep="")
meta_cell_child_names_dir = paste(meta_dir, "/cell_child_names",sep="")
meta_cell_childpn_dir = paste(meta_dir, "/cell_childpn",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
meta_cell_childpn_ind_dir = paste(meta_dir, "/cell_childpn_ind",sep="")
meta_cell_childpn_names_dir = paste(meta_dir, "/cell_childpn_names",sep="")
meta_cell_parent_dir = paste(meta_dir, "/cell_parent",sep="") #specifies a phenotypes parents
meta_cell_parent_ind_dir = paste(meta_dir, "/cell_parent_ind",sep="")
meta_cell_parentpn_dir = paste(meta_dir, "/cell_parentpn",sep="") #specifies a phenotypes parents and splits them into +/- (only for when both -/+ exists)
meta_cell_parentpn_ind_dir = paste(meta_dir, "/cell_parentpn_ind",sep="")

## libraries
library(foreach)
library(doMC)
source("~/projects/IMPC/code/_funcdist.R")
source("~/projects/IMPC/code/_funcAlice.R")

## cores
no_cores = detectCores()-1
registerDoMC(no_cores)









## options
options(stringsAsFactors=F)
options(device="cairo")
options(na.rm=T)

## prepare data
meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
# colnames(meta_cell) = c("phenotype","phenocode","phenolevel")
# save(meta_cell,file=paste0(meta_cell_dir, pcp,".Rdata"))










start = Sys.time()

start1 = Sys.time()

cat("\ncreating children indices ")
pc = getphenoChild(meta_cell, no_cores=no_cores)

meta_cell_child = pc$phenoChild
meta_cell_child_ind = pc$phenoChild_ind
meta_cell_child_names = pc$phenoChild_names
meta_cell_childpn = pc$phenoChildpn
meta_cell_childpn_ind = pc$phenoChildpn_ind
meta_cell_childpn_names = pc$phenoChildpn_names

#save
save(meta_cell_child, file=paste0(meta_cell_child_dir, ".Rdata"))
save(meta_cell_child_ind, file=paste0(meta_cell_child_ind_dir, ".Rdata"))
save(meta_cell_child_names, file=paste0(meta_cell_child_names_dir, ".Rdata"))
save(meta_cell_childpn, file=paste0(meta_cell_childpn_dir, ".Rdata"))
save(meta_cell_childpn_ind, file=paste0(meta_cell_childpn_ind_dir, ".Rdata"))
save(meta_cell_childpn_names, file=paste0(meta_cell_childpn_names_dir, ".Rdata"))

TimeOutput(start1)








start1 = Sys.time()

cat("\ncreating parents indices ")
pp = getphenoParent(meta_cell, meta_cell_childpn, meta_cell_childpn_ind, no_cores)

meta_cell_parent = pp$phenoParent
meta_cell_parent_ind = pp$phenoParent_ind
meta_cell_parentpn = pp$phenoParentpn
meta_cell_parentpn_ind = pp$phenoParentpn_ind

#save
save(meta_cell_parent, file=paste0(meta_cell_parent_dir, ".Rdata"))
save(meta_cell_parent_ind, file=paste0(meta_cell_parent_ind_dir, ".Rdata"))
save(meta_cell_parentpn, file=paste0(meta_cell_parentpn_dir, ".Rdata")) 
save(meta_cell_parentpn_ind, file=paste0(meta_cell_parentpn_ind_dir, ".Rdata"))

TimeOutput(start1)

TimeOutput(start)





