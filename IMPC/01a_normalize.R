## Input: original features --> Output: distance matrices
# aya43@sfu.ca 20151228

## root directory
root = "~/projects/IMPC"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")
controlL = c("+_+|+_Y","+_+|+_Y","WildType","WildType","WildType") #control value in target_col column
ci = 1; panel = panelL[ci]; centre = centreL[ci]

result_dir = paste0("result/", panelL, "/", centreL); suppressWarnings(dir.create (result_dir))


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
norm_factor_diff_dir = paste(norm_dir, "/norm_factor_diff", sep="")

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
target_col = "gene" #column with control/experiment #IMPC ONLY
control = controlL[ci] #IMPC ONLY

split_col = NULL

cutoff = Inf#c(.56)#c(.6) #run this script and look at plot to determine this number; determines when normalization factor is too far off from peak
layer_norm = 4 #layer at which to do normalization, should have larger cell populations; set to 0 if do all
cellCountThres = c(200)#,200,500,500,500) #insignificant if count under


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

#save
save(feat_file_cell_countAdj, file=paste0(feat_file_cell_countAdj_dir,".Rdata"))
if (writecsv) write.csv(feat_file_cell_countAdj, file=paste0(feat_file_cell_countAdj_dir,".csv"), row.names=T)
save(f0, file=paste0(norm_factor_dir,".Rdata"))
if (writecsv) write.csv(f0, file=paste0(norm_factor_dir,".csv"), row.names=T)
save(fdiff0, file=paste0(norm_factor_diff_dir,".Rdata"))
if (writecsv) write.csv(fdiff0, file=paste0(norm_factor_diff_dir,".csv"), row.names=T)

TimeOutput(start)

# #max change in slope
# cutoff = sort(abs(fdiff))[which.max(diff(sort(abs(fdiff))))]

# #Find inflection point of fdiff curve to define cut-off for when to switch to highest density method
# x = seq(1,length(fdiff))
# y = sort(abs(fdiff))
# #The exact inflection point is ip=5.0
# #Because of the total symmetry we expect inflection point near the middle of x-range:
# A=findiplist(x,y,0);A;#Our expectation came true
# #Let's make some ESE and EDE iterations and plot them:
# a=findipiterplot(x,y,0,TRUE,TRUE,TRUE);
# a$first;#Show first solution
# a$BESE;#Show ESE iterations
# a$CRESE;#Show cross ESE iterations from EDE iterations
# a$esm;#Show all ESE iterations
# a$BEDE;#Show EDE iterations
# a$CREDE;#Show cross EDE iterations from ESE iterations
# a$edm;#Show all EDE iterations
# a$aesmout;#Statistics and 95%c c.i. for ESE
# a$besmout;#Statistics and 95%c c.i. for cross EsE
# a$aedmout;#Statistics and 95%c c.i. for EDE
# a$bedmout;#Statistics and 95%c c.i. for cross EDE
# a$esmout;#Statistics and 95%c c.i. for all ESE results
# a$edmout;#Statistics and 95%c c.i. for all EDE results
# a$ipall;#Statistics and 95%c c.i. for results
# #Close the 2 previously opened devises





TimeOutput(start) #21min parallel, 5 centres






