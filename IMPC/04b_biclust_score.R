## Input: bicluster --> Output: bicluster scores
# aya43@sfu.ca 20180328

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
meta_file_dir = paste(meta_dir, "/file", sep="")
feat_dir = paste(result_dir, "/feat", sep="")
network_dir = paste(result_dir, "/network", sep=""); dir.create(network_dir, showWarnings=F)
# network_dist_dir = paste(network_dir, "/dist", sep=""); dir.create(network_dir, showWarnings=F)

## output directories
biclust_dir = paste(result_dir,  "/biclust", sep=""); dir.create (biclust_dir,showWarnings=F)
biclust_source_dir = paste(biclust_dir,  "/source", sep=""); dir.create (biclust_source_dir,showWarnings=F)
biclust_plot_dir = paste(biclust_dir,  "/plot", sep=""); dir.create (biclust_source_dir,showWarnings=F)
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
no_cores = detectCores()-1
registerDoMC(no_cores)








## options for script
network_socre = T #for IMPC / gene clusters
gene_col = "gene"

readcsv = F

control = str_split(controlL,"[|]")[[1]]
control2 = "[+]_[+]|[+]_Y"
time_col = "date"
id_col = "fileName"
target_col = "gene"
order_cols = c("barcode","sample","specimen","date")

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

netdist_paths = list.files(path=network_dir, full.names=T, recursive=T, pattern="dist_")















start = Sys.time()

#get meta file
if (readcsv) {
  # mc = read.csv(paste0(feat_dir,"/", feat_count,".csv"),row.names=1, check.names=F)
  meta_file = read.csv(paste0(meta_file_dir,".csv"),check.names=F)
} else {
  # mc = get(load(paste0(feat_dir,"/", feat_count,".Rdata")))
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
}

#get network file
netdists = list()
for (netdist_path in netdist_paths) {
  netdist = netdist0 = read.table(netdist_path,  check.names=F, sep="\t")
  netdist = netdist0[-1,-1]
  netgenes = netdist0[-1,1]
  if (all(str_sub(netgenes,1,1)==str_sub(netgenes[1],1,1))) netgenes = str_sub(netgenes,2,-1)
  rownames(netdist) = colnames(netdist) = netgenes
  netdists[[fileNames(netdist_path)]] = apply(as.matrix(netdist), 2, as.numeric)
}

scorelength = 5 + 1 + length(netdist_paths)*3 #total score length = f measure + network x 3 nca, silhouette, number of genes tested


score_list = foreach(clust_path=clust_paths) %dopar% {
  tryCatch ({
    cat("\n", clust_path, " ",sep="")
    start2 = Sys.time()
    
    ## prep clusters / label
    rowclust = rowclust0 = read.csv(paste0(clust_path,"_rowclust.csv"), row.names=1)
    rowclust = rowclust[,1]; 
    names(rowclust) = rownames(rowclust0)
    mm = get_feat_matrix2(clust_fileName=NULL, feat_dir=NULL, meta_file=meta_file, id_col=id_col, row_names=names(rowclust), col_names=NULL, getcsv=readcsv)
    sm = mm$sm
    rowlabel = sm[,target_col]
    names(rowlabel) = names(rowclust)
    # rowlabel = rowlabel0 = read.csv(paste0(clust_path,"_rowlabel.csv"), row.names=1)
    # rowlabel = rowlabel[,1]; 
    # names(rowlabel) = rownames(rowlabel0)
    
    # rowclust1 = rep(NA,length(rowlabel)) #s.t. each cluster is labeled by what majority of its real contents
    # tubes0 = unique(rowclust)[order(table(rowclust))]
    # for (tubei in tubes0) {
    #   tci = which(rowclust==tubei) #index of cluster tubei
    #   tubej = Mode(rowlabel[tci]) #tubei is label taking up majority of cluster tubei
    #   rowclust1[tci] = tubej ## MAJORITY IN 2+ classes?
    # }
    # names(rowclust) = names(rowclust1) = rownames(rowclust0)
    
    if (length(unique(rowclust))==1) {
      return(NA)
      # rowclust_df = matrix(rep(1,length(rowclust)),ncol=1)
      # rownames(rowclust_df) = names(rowclust)
      # colnames(rowclust_df) = rowclust[1]
    } else {
      rowclust_df = model.matrix(~ factor(rowclust) - 1); colnames(rowclust_df) = sort(unique(rowclust)); rownames(rowclust_df) = names(rowclust)
    }
    # rowclust1_df = model.matrix(~ factor(rowclust1) - 1); colnames(rowclust1_df) = sort(unique(rowclust1)); rownames(rowclust1_df) = names(rowclust1)
    rowlabel_df = model.matrix(~ factor(rowlabel) - 1); colnames(rowlabel_df) = sort(unique(rowlabel)); rownames(rowlabel_df) = names(rowlabel)
    rownames(rowlabel_df) = names(rowlabel)
    colnames(rowlabel_df) = laname = sort(unique(rowlabel))
    
    
    ## score
    score = length(unique(rowclust))
    
    
    ## external validation F1 (classification & clustering)
    # if (ncol(rowlabel_df)==ncol(rowclust1_df)) F11 = F.measure.single.over.classes(rowlabel_df,rowclust1_df)$average[-6]; names(F1) = paste0(names(F1),"_1")
    # f11c = f.measure.comembership(rowlabel,rowclust1); names(f11c) = paste0(names(f11c), "_co_1")
    # r1 = adjustedRand(rowlabel,rowclust1); names(r1) = paste0(names(r1), "_1")
    # score = c(score, f11c, r1)
    
    
    # ## internal validation NCA (distance)
    # nca = NCA_score(netdist, rowclust_df)$p
    


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
    
    
    
    
    ## external validation f1 (clustering)
    f1c = f.measure.comembership(rowlabel,rowclust); names(f1c) = paste0(names(f1c), "_co")
    # r = adjustedRand(rowlabel,rowclust)
    # score = c(score,f1c,r)
    score = c(score,unlist(f1c))
    
    
    
    ## internal validation of genes using gene to gene interaction network
    
    #format IMPC label genes
    gene0 = sm[,gene_col]
    names(gene0) = sm[,id_col]
    gene0 = gene0[!grepl(" KO$",gene0)]
    gene0 = gene0[!grepl("failed", gene0)]
    gene0 = gene0[!grepl(control2, gene0)]
    
    #continue if there is non contorl gene in biclusters
    if (length(gene0)==0) {
      length(score) = scorelength
      return(score)
    }
    
    gene = str_split(gene0,"_")
    gene1 = sapply(gene, function(x) x[1])
    gene2 = sapply(gene, function(x) x[2])
    names(gene) = names(gene1) = names(gene2) = names(gene0)
    cat("\ntwo copies of genes are the same for all mice: ")
    identical(gene1,gene2)
    cat("\nonly one copy of gene KO-ed for genes: ")
    # sort(unique(gene1[gene2=="+"]))
    cat("\n(one copy of gene KO-ed / all) for genes: ")
    for (g in unique(gene1[gene2=="+"])) cat("\n(",sum(gene1==g & gene2=="+"),"/",sum(gene1==g),") ", g, sep="")
    # sort(unique(gene1[gene2=="Y"]))
    cat("\n(one copy of gene=Y / all) for genes: ")
    for (g in unique(gene1[gene2=="Y"])) cat("\n(",sum(gene1==g & gene2=="Y"),"/",sum(gene1==g),") ", g, sep="")
    cat("\ndelete ")
    # gene = gene1[gene2!="+" & gene2!="Y"] #delete 118 genes because mostly 1/2 or 1/1 -/+
    gene = gene1
    gene = str_split(gene,"[(]") #take out "(b)" at the end of gene name, number of genes dont decrease
    gene = sapply(gene, function(x) x[1])
    rowgene = toupper(gene)
    names(rowgene) = names(gene1)
    
    
    
    #get dist matrix between genes
    for (netdisti in names(netdists)) {
      netdist = netdists[[netdisti]]
      netgenes = colnames(netdist)
      evalgenes = intersect(rowgene, netgenes)
      
      rowgene1 = rowgene[rowgene%in%evalgenes]
      rowclust1 = rowclust[match(names(rowgene1),names(rowclust))]
      rowclust1 = rowclust1[rowclust1>0]
      if (length(rowclust1)==0) next
      if (length(unique(rowclust1))==1) next
      rowgene1 = rowgene1[match(names(rowclust1),names(rowgene1))]
      
      rowclust1_df = model.matrix(~ factor(rowclust1) - 1); colnames(rowclust1_df) = sort(unique(rowclust1)); rownames(rowclust1_df) = names(rowclust1)
      netdist_rowcol = match((rowgene1),netgenes)
      netdist1 = netdist[netdist_rowcol,netdist_rowcol]
      netdist1 = apply(as.matrix(netdist1), 2, as.numeric)
      rownames(netdist1) = colnames(netdist1) = names(rowgene1)
      
      #metrics
      nca = NCA_score(netdist1, rowclust1)$p
      names(nca) = paste0("NCA_net-",netdisti,".g-",length(unique(rowgene1)))
      score = append(score,nca)
      
      silhouette = silhouette(rowclust1,netdist1)
      pngame = paste0(biclust_plot_dir,"/",fileNames(clust_path),"_silhouette.net-",netdisti,".g-",length(unique(rowgene1)),".png")
      png(pngame)
      plot(silhouette)
      graphics.off()
      sil = median(silhouette[,3])
      names(sil) = paste0("median-silhouette_net-",netdisti,".g-",length(unique(rowgene1)))
      score = append(score,sil)
      
      score = append(score, length(unique(rowgene1)))
    }
    
    
    length(score) = scorelength
    
    
    write.csv(score,file=paste0(biclust_score_dir, "/", clust_path, "_score.csv"))
    TimeOutput(start2)
  }, error = function(err) { cat(paste("error:  ",err)); return(NA) })
  return((score))
}

















## put scores into a table

score_list0 = score_list
scorelength = sapply(score_list, function(x) length(x))

keepind = !is.na(score_list) & scorelength>1
score_list = score_list[keepind]

scorelength = sapply(score_list, function(x) length(x))

if (!all(scorelength==max(scorelength))) {
  for (si in which(scorelength<max(scorelength))) {
    a = score_list[[si]]
    length(a) = max(scorelength)
    score_list[[si]] = a
  }
}
score_table = Reduce("rbind",score_list)

clust_files = fileNames(clust_paths[keepind])
clust_files_attr = str_split(clust_files,"_")
clust_files_table = t(sapply(clust_files_attr, function(x) 
  c(x[1], x[2], str_split(x[3],"-")[[1]][2], gsub("layer","",x[4]), str_split(x[5],"-")[[1]][2])
))
colnames(clust_files_table) = c("method","feature","splitby","layer","count-thresh")

clust_files_table_final = cbind(clust_files_table, score_table)
write.csv(clust_files_table_final, file = paste0(biclust_score_dir,"_number-of-cluster_f-measures_NCA_Silhouette_number-of-genes.csv"))

TimeOutput(start)



