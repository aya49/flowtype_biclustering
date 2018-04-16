## Input: bicluster --> Output: bicluster scores
# aya43@sfu.ca 20180328

## root directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

## input directories
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste(meta_dir, "/file", sep="")
feat_dir = paste(result_dir, "/feat", sep="")

## output directories
biclust_dir = paste(result_dir,  "/biclust", sep=""); dir.create (biclust_dir,showWarnings=F)
biclust_source_dir = paste(biclust_dir,  "/source", sep=""); dir.create (biclust_source_dir,showWarnings=F)
biclust_score_dir = paste(biclust_dir,  "/score", sep="")


## libraries
library(biclust)
library(NMF)
library(pheatmap)
library(foreach)
library(doMC)
library(stringr)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
source("~/projects/IMPC/code/_bayesianbiclustering.R")

#Setup Cores
no_cores = 3#detectCores()-3
registerDoMC(no_cores)








## options for script
attributes = c("aml")

control = "normal" #control value in target_col column for each centre
id_col = "fileName"
target_col = "aml"
order_cols = c("specimen","aml")
split_col = "tube"

# bcmethods = c("plaid","CC","bimax","BB-binary","nmf-nsNMF","nmf-lee","nmf-brunet","CC")
# onlysigBB = T #only evaluate significant BB-binary clusters
# Kr = 20; Kc = 20 #number of row and column biclusters for binary bayesian biclustering
# pval_thres = .01 #pval_threshold for choosing biclusters out of all bayesian biclusters
# min_iter = 100 #min number of iterations for BB-binary (B2PS)
# nmf_thres = .05 # * max contribution: threshold at which a row/col can be considered a significant contribution to a factor in nmf

#data paths
clust_paths = list.files(path=biclust_source_dir,pattern=".Rdata", full.names=T)
                     # list.files(path=clust_source_dir,pattern=".Rdata", full.names=T))
clust_paths = gsub(".Rdata","",clust_paths)
feat_count = "file-cell-countAdj"
















start = Sys.time()

centre = paste0(panel," ",centre)
cat("\n",centre,sep="")


mc = get(load(paste0(feat_dir,"/", feat_count,".Rdata")))
meta_file = get(load(paste0(meta_file_dir,".Rdata")))

score_list = foreach(clust_path=clust_paths) %dopar% {
  cat("\n", clust_path, " ",sep="")
  start2 = Sys.time()
  
  ## prep clusters / label
  rowclust = rowclust0 = read.csv(paste0(clust_path,"_rowclust.csv"), row.names=1)
  rowclust = rowclust[,1]; 
  rowlabel = rowlabel0 = read.csv(paste0(clust_path,"_rowlabel.csv"), row.names=1)
  rowlabel = rowlabel[,1]; 
  # rowclust1 = rep(NA,length(rowlabel)) #s.t. each cluster is labeled by what majority of its real contents
  # tubes0 = unique(rowclust)[order(table(rowclust))]
  # for (tubei in tubes0) {
  #   tci = which(rowclust==tubei) #index of cluster tubei
  #   tubej = Mode(rowlabel[tci]) #tubei is label taking up majority of cluster tubei
  #   rowclust1[tci] = tubej ## MAJORITY IN 2+ classes?
  # }
  # names(rowclust) = names(rowclust1) = rownames(rowclust0)
  names(rowclust) = rownames(rowclust0)
  names(rowlabel) = rownames(rowlabel0)
  
  rowclust_df = model.matrix(~ factor(rowclust) - 1); colnames(rowclust_df) = sort(unique(rowclust)); rownames(rowclust_df) = names(rowclust)
  # rowclust1_df = model.matrix(~ factor(rowclust1) - 1); colnames(rowclust1_df) = sort(unique(rowclust1)); rownames(rowclust1_df) = names(rowclust1)
  rowlabel_df = model.matrix(~ factor(rowlabel) - 1); colnames(rowlabel_df) = sort(unique(rowlabel)); rownames(rowlabel_df) = names(rowlabel)
  rownames(rowlabel_df) = names(rowlabel)
  colnames(rowlabel_df) = laname = sort(unique(rowlabel))
  
  
  ## score
  score = NULL
  
  
  ## external validation F1 (classification & clustering)
  # if (ncol(rowlabel_df)==ncol(rowclust1_df)) F11 = F.measure.single.over.classes(rowlabel_df,rowclust1_df)$average[-6]; names(F1) = paste0(names(F1),"_1")
  # f11c = f.measure.comembership(rowlabel,rowclust1); names(f11c) = paste0(names(f11c), "_co_1")
  # r1 = adjustedRand(rowlabel,rowclust1); names(r1) = paste0(names(r1), "_1")
  # score = c(score, f11c, r1)
  
  
  ## external validation f1 (clustering)
  f1c = f.measure.comembership(rowlabel,rowclust); names(f1c) = paste0(names(f1c), "_co")
  r = adjustedRand(rowlabel,rowclust)
  score = c(score,f1c,r)
  
  
  # ## internal validation NCA (distance)
  # cl = la0
  # if (!"NCA"%in%names(fm[[colnam]][[dindname]][[cltype]][[par]]) | overwritef) {
  #   if (length(unique(cl))==1) { nca = NA
  #   } else { nca = NCA_score(as.matrix(d[[dindname]]), rowclust_df)$p }
  #   fm[[colnam]][[dindname]][[cltype]][[par]]["NCA"] = nca
  # }
  # 
  # 
  # ## internal validation silmed (distance & clustering)
  # if (!"silmed"%in%names(fm[[colnam]][[dindname]][[cltype]][[par]]) | overwritef) {
  #   if (length(unique(cl))==1) { sil = NA
  #   } else { sil = median(silhouette(cl,d[[dindname]])[,3]) }
  #   fm[[colnam]][[dindname]][[cltype]][[par]]["silmed"] = sil
  # }
  # 
  # if (length(unique(cl))==1) {
  #   score = rep(NA,9)
  #   names(score) = c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")
  # } else { 
  #   score = unlist(cluster.stats(as.dist(d[[dindname]]),cl))[c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")] 
  # }
  # fm[[colnam]][[dindname]][[cltype]][[par]] = append(fm[[colnam]][[dindname]][[cltype]][[par]], score)
  # 
  # 
  # ## internal validation (adjusted clustering)
  # if (length(unique(rowclust1))==1) { sil = NA
  # } else { sil = median(silhouette(rowclust1,d[[dindname]])[,3]) }
  # fm[[colnam]][[dindname]][[cltype]][[par]]["silmed_1"] = sil
  # 
  # if (length(unique(rowclust1))==1) {
  #   score = rep(NA,9)
  #   names(score) = c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")
  # } else {
  #   score = unlist(cluster.stats(as.dist(d[[dindname]]),rowclust1))[c("average.between","average.within","avg.silwidth","pearsongamma","dunn","dunn2","entropy","wb.ratio","ch")] 
  # }
  # names(score) = paste0(names(score),"_1")
  # fm[[colnam]][[dindname]][[cltype]][[par]] = append(fm[[colnam]][[dindname]][[cltype]][[par]], score)
  
  write.csv(score,file=paste0(clust_path,"_score.csv"))
  return(score)
  TimeOutput(start2)
}


















score_table = Reduce("rbind",score_list)

clust_files = fileNames(clust_paths)
clust_files_attr = str_split(clust_files,"_")
clust_files_table = t(sapply(clust_files_attr, function(x) {
  clust_file = c(x[1], x[2], str_split(x[3],"-")[[1]][2], gsub("layer","",x[4]), str_split(x[5],"-")[[1]][2])
}))
colnames(clust_files_table) = c("method","feature","splitby","layer","count-thresh")

clust_files_table_final = cbind(clust_files_table, score_table)
write.csv(clust_files_table_final, file = paste0(biclust_score_dir,".csv"))


TimeOutput(start)





