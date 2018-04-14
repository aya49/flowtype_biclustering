## Input: original count matrix --> Output: normalized count matrix
#aya43@sfu.ca 20151228

## root directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)



## input directories
meta_dir = paste0(result_dir,"/meta")
meta_cell_dir = paste(meta_dir, "/cell", sep="")
meta_file_dir = paste(meta_dir, "/file", sep="")
feat_dir = paste(result_dir, "/feat", sep="")
feat_file_cell_count_dir = paste(feat_dir, "/file-cell-count", sep="")

## output directories
feat_file_cell_countAdj_dir = paste(feat_dir, "/file-cell-countAdj", sep="")
feat_file_cell_countAdjLog_dir = paste(feat_dir, "/file-cell-countAdjLog", sep="")
norm_dir = paste(result_dir, "/cell_count_norm",sep=""); dir.create(norm_dir,showWarnings=F)
norm_factor_dir = paste(norm_dir, "/norm_factor", sep=""); dir.create(norm_factor_dir,showWarnings=F) #plot of norm factor for each file
norm_factor_diff_log_dir = paste(norm_dir, "/norm_factor_diff_dens_logged", sep="")

## libraries
library(stringr)
library(pracma)
library(foreach)
library(doMC)
# library(flowDensity)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")

## cores
no_cores = detectCores()-1
registerDoMC(no_cores)











## options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

writecsv = T

id_col = "fileName"
target_col = "aml" #column with control/experiment
control = "normal" #control value in target_col column

split_col = "tube" #AML ONLY #column to split by -- FlowCAP-II normalizes only between files of same tube/panel; set to NULL if don't split

cutoff = c(Inf) #c(.6) #if TMM-peak>cutoff, then apply peak instead of TMM; run this script and look at norm_fdiffplot plot to determine this number
layer_norm = 4 #0 #calculate TMM using only phenotypes in this layer; set to 0 if do for all layers
cellCountThres = 1200 #don't use phenotypes with cell count lower than cellCountThres




#Prepare data
meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
meta_cell = get(load(paste0(meta_cell_dir,".Rdata")))
feat_file_cell_count0 = get(load(paste0(feat_file_cell_count_dir,".Rdata")))










start = Sys.time()

#calculate TMM
fdiff0 = f0 = rep(0,nrow(meta_file0))
split_unique = NULL; if (!is.null(split_col)) split_unique = unique(meta_file0[,split_col])

for (ti in split_unique) {
  start1 = Sys.time()
  
  #split by ti
  tube_ind = rep(T,nrow(meta_file0)); if(!is.null(ti)) tube_ind = meta_file0[,split_col]==ti
  feat_file_cell_count = feat_file_cell_count0[tube_ind,]
  meta_file = meta_file0[tube_ind,]
  
  #prepare feat_file_cell_counts
  x = x0 = as.matrix(feat_file_cell_count)[,-1] #take out total cell count
  if (layer_norm>0) x = as.matrix(x0[,colnames(x0)%in%meta_cell$phenotype[meta_cell$phenolevel==layer_norm] & sapply(1:ncol(x0), function(y) any(x0[,y]>cellCountThres))])
  lib.size = feat_file_cell_count[,1]
  refColumn = which.min(abs( lib.size - median(lib.size[grepl(control,meta_file[,target_col])]) )) #reference column: median total count out of all control files
  
  #prepare plot paths/titles
  loop.ind = 1:nrow(x)
  pngnames = sapply(loop.ind, function(i) {
    pngname = paste0(norm_factor_dir,"/")
    if (!is.null(ti)) pngname = paste0(pngname,split_col,"-",ti,"_")
    pngname = paste0(pngname, meta_file[i,target_col], "_", meta_file[i,id_col],".png")
  })
  mains = sapply(loop.ind, function(i) paste0("mean count vs. ln fold change:\n", meta_file[i,target_col]," over refColumn ", meta_file[refColumn,target_col], "___layer-",layer_norm))
  
  ## calculate absolute count TMM, mostly taken from TMM
  fresult = tmm(x,x0,lib.size,refColumn,cutoff=Inf,plotimg=T,pngnames=pngnames,mains=mains,no_cores=no_cores,samplesOnCol=F)
  f0[tube_ind] = fresult$f
  fdiff0[tube_ind] = fresult$fdiff
  
  #plot difference between TMM and peak for all files
  pngname = ifelse(is.null(ti), paste0(norm_factor_dir,"/all.png"), paste0(norm_factor_dir,"/all_",split_col,"-",ti,".png"))
  png(file=pngname , width=700, height=700)
  plot(sort(abs(fdiff0[tube_ind])), cex=.4, ylim=c(0,3), main="cell-count-norm-factor_f_diff_from_peak_abs")
  lines(sort(abs(fdiff0[tube_ind])), col="blue")
  graphics.off()
  
  TimeOutput(start1)
}

feat_file_cell_countAdj = sapply(c(1:nrow(feat_file_cell_count0)), function(x) {feat_file_cell_count0[x,]*f0[x]})
feat_file_cell_countAdj = t(feat_file_cell_countAdj)
colnames(feat_file_cell_countAdj) = colnames(feat_file_cell_count0)
rownames(feat_file_cell_countAdj) = rownames(feat_file_cell_count0)

#phenotype on cols
feat_file_cell_countAdjlog = log(feat_file_cell_countAdj)
feat_file_cell_countAdjlog[which(feat_file_cell_countAdjlog<0)] = 0

#save
save(f0, file=paste0(norm_factor_dir,".Rdata"))
save(fdiff0, file=paste0(norm_factor_diff_log_dir,".Rdata"))
save(feat_file_cell_countAdj, file=paste0(feat_file_cell_countAdj_dir,".Rdata"))
if (writecsv) write.csv(feat_file_cell_countAdj, file=paste0(feat_file_cell_countAdj_dir,".csv"), row.names=F)
save(feat_file_cell_countAdjlog, file=paste0(norm_factor_diff_log_dir,".Rdata"))
if (writecsv) write.csv(feat_file_cell_countAdjlog, file=paste0(norm_factor_diff_log_dir,".csv"), row.names=F)

TimeOutput(start)






