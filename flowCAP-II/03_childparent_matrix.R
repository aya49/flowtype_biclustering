# aya43@sfu.ca 20170131
# creates list of children for each non-leaf node and creates children/parent list[[matrices]] (-/+ are only for phenotypes where both -,+ data exists)

## root directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)


## input directories
meta_dir = paste0(result_dir,"/meta")
meta_cell_dir = paste(meta_dir, "/cell", sep="")
meta_cell_child_dir = paste(meta_dir, "/cell_child",sep="") #specifies a phenotypes children
meta_cell_child_ind_dir = paste(meta_dir, "/cell_child_ind",sep="")
meta_cell_childpn_dir = paste(meta_dir, "/cell_childpn",sep="") #specifies a phenotypes children and splits them into +/- (only for when both -/+ exists)
meta_cell_childpn_ind_dir = paste(meta_dir, "/cell_childpn_ind",sep="")
meta_cell_parent_dir = paste(meta_dir, "/cell_parent",sep="") #specifies a phenotypes parents
meta_cell_parent_ind_dir = paste(meta_dir, "/cell_parent_ind",sep="")
meta_cell_parentpn_dir = paste(meta_dir, "/cell_parentpn",sep="") #specifies a phenotypes parents and splits them into +/- (only for when both -/+ exists)
meta_cell_parentpn_ind_dir = paste(meta_dir, "/cell_parentpn_ind",sep="")

feat_dir = paste(result_dir, "/feat", sep=""); dir.create(feat_dir, showWarnings=F)
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")

## output directories
feat_file_edge_pnratio_dir = paste(feat_dir, "/file-edge-pnratio",sep="")
feat_file_edge_prop_dir = paste(feat_dir, "/file-edge-prop",sep="")
feat_file_cell_entropychild_dir = paste(feat_dir, "/file-cell-entropychild",sep="")
feat_file_cell_entropyparent_dir = paste(feat_dir, "/file-cell-entropyparent",sep="")

## libraries
library(stringr)
library(entropy)
library(foreach)
library(doMC)
source("code/_funcAlice.R")





## cores
no_cores = 3#detectCores() - 1
registerDoMC(no_cores)




## options
options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

create_child_entropy = T
create_parent_entropy = T



writecsv = F




start = Sys.time()




start1 = Sys.time()

#get list of children for each non-leaf node & save
cat("\ncreating child matrix")

m = get(load(paste0(feat_file_cell_count_dir,".Rdata")))
meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
meta_cell_child = get(load(paste0(meta_cell_child_dir,".Rdata")))
meta_cell_child_ind = get(load(paste0(meta_cell_child_ind_dir,".Rdata")))




## child proportion --------------------------------------------

start1 = Sys.time()

#mlist=list()
cat("; childprop")
childprop = foreach(i = 1:length(meta_cell_child)) %dopar% { #for each phenotype
  #for (i in 1:length(meta_cell_child)) {
  parent = m[,meta_cell_child_ind[i]]
  
  # Child proportion matrix
  children = m[,unlist(meta_cell_child[[i]])]
  children[which(children<1)] = 1
  parent0 = parent
  parent0[which(parent<1)] = 1
  childprop = exp(children/parent0)
  if (is.null(dim(childprop))) {
    childprop = matrix(childprop,ncol=1)
    colnames(childprop) = colnames(m)[unlist(meta_cell_child[[i]])]
  }
  
  return(childprop)
}
names(childprop) = meta_cell$phenotype[meta_cell_child_ind]
for (i in 1:length(childprop)) {
  colnames(childprop[[i]]) = paste0(names(childprop)[i], "_", colnames(childprop[[i]]))
}
childprop = Reduce("cbind",childprop)
save(childprop, file=paste0(feat_file_edge_prop_dir,".Rdata"))
if (writecsv) write.csv(childprop, file=paste0(feat_file_edge_prop_dir,".csv"))


TimeOutput(start1)




## child pn ratio + child entropy --------------------------------------------

start1 = Sys.time()


cat(", childratio + childentropy")
meta_cell_childpn = get(load(paste0(meta_cell_childpn_dir, ".Rdata")))
meta_cell_childpn_ind = get(load(paste0(meta_cell_childpn_ind_dir, ".Rdata")))

mlist = foreach(i = 1:length(meta_cell_childpn)) %dopar% { #for each phenotype
  parent = m[,meta_cell_childpn_ind[i]]
  
  # P/N ratio matrix
  pnratio = NULL
  
  pos = m[,meta_cell_childpn_ind[i]]-m[,meta_cell_childpn[[i]][[1]]]
  neg = m[,meta_cell_childpn[[i]][[1]]]
  pnratio = pos/neg ##get rid of 0, Inf
  pnratio[pos==0 & neg==0] = 1
  pnratio[pos==0 & neg!=0] = 1/neg[pos==0 & neg!=0]
  pnratio[pos!=0 & neg==0] = pos[pos!=0 & neg==0]
  pnratio = log(pnratio)
  
  if (is.null(dim(pnratio))) {
    pnratio = matrix(pnratio,ncol=1)
    colnames(pnratio) = colnames(m)[meta_cell_childpn[[i]][[1]]]
  } 
  if (sum(parent==0)>0) pnratio[which(parent==0),] = rep(0,length(meta_cell_childpn[[i]][[1]]))
  
  # Entropy matrix
  en = rep(0,nrow(m))
  
  parent = m[,meta_cell_childpn_ind[i]]
  children = m[,unlist(meta_cell_childpn[[i]])]
  if (length(meta_cell_childpn[[i]])==1) children = matrix(children,ncol=1)
  
  
  no_child = length(meta_cell_childpn[[i]][[1]])
  non0parents = which(parent>0)
  en[non0parents] = sapply(non0parents, function(x) return(entropy(children[x,]/parent[x])/no_child)) #average entropy over # of markers added
  
  #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
  return(list(entropy=en, pnratio=pnratio)) #ratio = +child_prop / -child_prop
}


pnratio = lapply(mlist, function(x) x$pnratio)
en = lapply(mlist, function(x) x$en)

rm(mlist)

feat_file_cell_entropychild = Reduce('cbind',en)
names(pnratio) = colnames(feat_file_cell_entropychild) = meta_cell$phenotype[meta_cell_childpn_ind]
for (i in 1:length(pnratio)) {
  colnames(pnratio[[i]]) = paste0(names(pnratio)[i], "_", colnames(pnratio[[i]]))
}
pnratio = Reduce("cbind",pnratio)

save(pnratio, file=paste0(feat_file_edge_pnratio_dir, ".Rdata"))
save(feat_file_cell_entropychild, file=paste0(feat_file_cell_entropychild_dir, ".Rdata"))
if (writecsv) write.csv(Reduce('cbind',pnratio), file=paste0(feat_file_edge_pnratio_dir, ".csv"))
if (writecsv) write.csv(feat_file_cell_entropychild, file=paste0(feat_file_cell_entropychild_dir, ".csv"))




TimeOutput(start1)





## parent entropy --------------------------------------------

start1 = Sys.time()


cat(", parententropy")

meta_cell_parentpn = get(load(paste0(meta_cell_parentpn_dir, ".Rdata")))
meta_cell_parentpn_ind = get(load(paste0(meta_cell_parentpn_ind_dir, ".Rdata")))


feat_file_cell_entropyparent = foreach(i = 1:length(cell_parentpn), .combine='cbind') %dopar% { #for each phenotype
  # Entropy matrix
  en = rep(0,nrow(m))
  
  childr = m[,meta_cell_parentpn_ind[i]]
  parent = m[,meta_cell_parentpn[[i]]]
  if (length(meta_cell_parentpn[[i]])==1) parent = matrix(parent,ncol=1)
  
  no_parent = length(meta_cell_parentpn[[i]])
  non0childrs = which(childr>0)
  en[non0childrs] = sapply(non0childrs, function(x) return(entropy(childr[x]/parent[x,])/no_parent)) #average entropy over # of markers added
  
  #mlist[[i]] = list(pnratio=pnratio, childprop=childprop)
  return(en) #ratio = +child_prop / -child_prop
}

colnames(feat_file_cell_entropyparent) = meta_cell$phenotype[meta_cell_parentpn_ind]
save(feat_file_cell_entropyparent, file=paste0(feat_file_cell_entropyparent_dir, ".Rdata"))
if (writecsv) write.csv(feat_file_cell_entropyparent, file=paste0(feat_file_cell_entropyparent_dir, ".csv"))







TimeOutput(start1)



TimeOutput(start)





