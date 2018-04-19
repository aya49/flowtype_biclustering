## Input: bicluster --> Output: bicluster plots
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

## output directories
biclust_dir = paste(result_dir,  "/biclust", sep=""); dir.create (biclust_dir,showWarnings=F)
biclust_source_dir = paste(biclust_dir,  "/source", sep=""); dir.create (biclust_source_dir,showWarnings=F)
biclust_plot_dir = paste(biclust_dir,  "/plot", sep=""); dir.create (biclust_plot_dir,showWarnings=F) #path to store plots


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
no_cores = detectCores()-3
registerDoMC(no_cores)








## options for script
readcsv = T

overwrite = T #overwrite biclust?
writecsv = F

good_count = 3
good_sample = 3


plot_size = c(500,500)
plot_size_bar = c(1300,2000)
attributes = c("gene","gender","date") #interested attribute to plot

control = str_split(controlL,"[|]")[[1]]
time_col = "date"
id_col = "fileName"
target_col = "gene"
order_cols = c("barcode","sample","specimen","date")
split_col = NULL

# bcmethods = c("plaid","CC","bimax","BB-binary","nmf-nsNMF","nmf-lee","nmf-brunet","CC")
# onlysigBB = T #only evaluate significant BB-binary clusters
# Kr = 20; Kc = 20 #number of row and column biclusters for binary bayesian biclustering
# pval_thres = .01 #pval_threshold for choosing biclusters out of all bayesian biclusters
# min_iter = 100 #min number of iterations for BB-binary (B2PS)
# nmf_thres = .05 # * max contribution: threshold at which a row/col can be considered a significant contribution to a factor in nmf

#biclustering result paths
clust_paths = list.files(path=biclust_source_dir,pattern=".Rdata", full.names=T)
# , list.files(path=clust_source_dir,pattern=".Rdata", full.names=T))
clust_paths = gsub(".Rdata","",clust_paths)
feat_count = "file-cell-countAdj"
















start = Sys.time()



if (readcsv) {
  mc = read.csv(paste0(feat_dir,"/", feat_count,".csv"),row.names=1, check.names=F)
  meta_file = read.csv(paste0(meta_file_dir,".csv"),check.names=F)
} else {
  mc = get(load(paste0(feat_dir,"/", feat_count,".Rdata")))
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
}

a = foreach(clust_path=clust_paths) %dopar% {
  cat("\n", clust_path, " ",sep="")
  start2 = Sys.time()
  
  ## prep clusters / label
  bc0 = get(load(paste0(clust_path,".Rdata")))
  bc = bc0$source
  if (is.null(bc0)) next
  if (bc@Number == 0) next
  rowclust = bc0$rowclust
  colclust = bc0$colclust
  
  # score = score0 = read.csv(paste0(clust_path,"_score.csv"), row.names=1)
  # score = score[,1]
  # names(score) = rownames(score0)
  
  
  
  
  
  
  ## plot ---------------------------------
  clust_path_plot = gsub(biclust_source_dir,biclust_plot_dir,clust_path)
  
  ## nmf special plot to see factors if needed
  if (grepl("nmf|GrNMF|fabia",clust_path)) {
    rowxfactor = bc@info$rowxfactor
    factorxcol = bc@info$factoxcol
    png(paste0(clust_path_plot, "_rowxfactor.png", sep=""), height=500*2, width=600*2)
    par(mfrow=c(2,2))
    plot(sort(rowxfactor[,1]),type="l", main="contribution of rows to each factor; each line = factor" )
    for (fi in 2:ncol(rowxfactor)) {
      lines(sort(rowxfactor[,fi]))
    }
    plot(sort(factorxcol[1,]),type="l", main="contribution of cols to each factor; each line = factor" )
    for (fi in 2:nrow(factorxcol)) {
      lines(sort(factorxcol[fi,]))
    }
    
    aheatmap(rowxfactor, main="row x factor")
    aheatmap(factorxcol, main="factor x col")
    graphics.off()
  }
  
  
  
  ## pretty heatmap
  
  # get original feature matrix and meta file
  mm = get_feat_matrix(fileNames(clust_path), feat_dir, mc, meta_file, id_col, target_col, control, order_cols, good_count, good_sample)
  m = mm$m
  sm = mm$sm
  
  # prepare col/row annotation
  if (grepl("BB-binary",clust_path)) {
    colcluster_temp = bc@info$patient.clusters
    rowcluster_temp = bc@info$transcript.clusters
    colcluster_temp[!colcluster_temp%in%bc0$colclust] = 0
    rowcluster_temp[!rowcluster_temp%in%bc0$rowclust] = 0
  } else {
    colcluster_temp = unlist(apply(bc@NumberxCol, 2, function(x) {
      clust = which(x)
      if (length(clust) == 0) return(0)
      return(max(clust))
    }))
    rowcluster_temp = unlist(apply(bc@RowxNumber, 1, function(x) {
      clust = which(x)
      if (length(clust) == 0) return(0)
      return(max(clust))
    }))
  }
  
  
  # plot
  tryCatch ({
    for (attribute in attributes) {
      
      col_annot = data.frame(colcluster=factor(colcluster_temp))
      row_annot = data.frame(rowcluster=factor(rowcluster_temp), class=factor(sm[,attribute]))
      rownames(col_annot)=names(colclust)
      
      row_order = order(row_annot$rowcluster, decreasing=T)
      col_order = order(col_annot$colcluster, decreasing=T)
      annotation_row = row_annot[row_order,]
      annotation_col = data.frame(col_annot[col_order,])
      
      mh = as.matrix(m)[row_order,col_order]
      rownames(mh) = 1:nrow(mh)
      rownames(annotation_col) = colnames(mh)
      rownames(annotation_row) = 1:nrow(mh)
      
      pheatmap(mh, main=paste0(attribute),#,"; precision, recall, f measure comembership:\n",paste0(c(f1$p,f1$r,f1$f_comember),collapse=", ")),
               annotation_row = annotation_row, annotation_col = annotation_col,
               show_rownames = F, 
               show_colnames=F,
               cluster_cols = T, cluster_rows = F,
               # cellwidth = 3, cellheight = 3, 
               # fontsize = 3, 
               filename = paste0(clust_path_plot, "_pheatmap-",attribute,"-sorted.pdf"))
      pheatmap(mh, main=paste0(attribute),#,"; precision, recall, f measure comembership:\n",paste0(c(f1$p,f1$r,f1$f_comember),collapse=", ")),
               annotation_row = annotation_row, annotation_col = annotation_col,
               show_rownames = F, 
               show_colnames=F,
               cluster_cols = F, cluster_rows = F,
               # cellwidth = 3, cellheight = 3, 
               # fontsize = 3, 
               filename = paste0(clust_path_plot, "_pheatmap-",attribute,"-sorted.png"))
      
      pheatmap(mh, main=paste0(attribute),#,"; precision, recall, f measure comembership:\n",paste0(c(f1$p,f1$r,f1$f_comember),collapse=", ")),
               annotation_row = annotation_row, annotation_col = annotation_col,
               show_rownames = F, 
               # show_colnames=F,
               cluster_cols = F, cluster_rows = F,
               # cellwidth = 3, cellheight = 3, 
               fontsize = 3, 
               filename = paste0(clust_path_plot, "_pheatmap-",attribute,"-sorted-smallfont.pdf"))
    }
  }, error = function(err) { cat(paste("pheatmap error:  ",err)); return(T) })
  # # Specifying clustering from distance matrix
  # drows = dist(test, method = "minkowski")
  # dcols = dist(t(test), method = "minkowski")
  # pheatmap(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)
  
  
  
  # png(paste0(clust_path_plot, "_bar.png", sep=""), height=plot_size_bar[1], width=plot_size_bar[2])
  # par(mar=c(2,6,3,2))
  # try({
  #   plotclust(bc,as.matrix(m))
  # })
  # graphics.off()
  
  # png(paste0(clust_path_plot, "_bubble.png", sep=""), height=plot_size_bar[1], width=plot_size_bar[2])
  # par(mar=c(20,20,20,20))
  # try({
  #   bubbleplot(as.matrix(m),bc)
  # })
  # graphics.off()
  
  
  ## plot heatmaps
  tryCatch ({
    png(paste0(clust_path_plot, "_heatmap0.png", sep=""), height=plot_size_bar[1], width=plot_size_bar[2])
    par(mar=c(5,3,6,5))
    try({
      heatmapBC(as.matrix(m),bc, order=T, local=T, outside=T)
    })
    graphics.off()
    
    rowcolno = ceiling(sqrt(bc@Number))
    png(paste0(clust_path_plot, "_heatmap.png", sep=""), height=plot_size[1]*rowcolno, width=plot_size[2]*rowcolno)
    par(mar=c(25,50,20,50))
    par(mfrow=rep(rowcolno,2))
    for (BCi in 1:bc@Number) {
      drawHeatmap(as.matrix(m),bc,BCi,plotAll=T)
    }
    graphics.off()
  }, error = function(err) { cat(paste("heatmap error:  ",err)); return(T) })
  
  
  
  ## plot row clusters against different sample attributes
  tryCatch({
    
    for (attri in attributes) {
      # png(paste0(clust_path_plot, "_stats_",attri,".png", sep=""), height=plot_size[1]*rowcolno, width=plot_size[2]*rowcolno)
      # par(mar=c(10,5,3,2), mfrow=rep(rowcolno,2))
      attri_valuesL = list()
      attri_topics = c()
      for (BCi in 1:bc@Number) {
        attri_values = sm[bc@RowxNumber[,BCi],attri]
        attri_valuesL[[BCi]] = attri_valuesT = table(attri_values)
        attri_topics = union(attri_topics,names(attri_valuesT))
        # try({
        #   barplot(attri_valuesT, las=2, xlab=attri, ylab="# of samples in bicluster with attribute on x axis", main=paste0(length(attri_values)," samples in bicluster ", BCi))
        # })
      }
      # graphics.off()
      
      attri_topics = sort(attri_topics)
      attri_valuesM = sapply(attri_valuesL, function(x) {
        sapply(attri_topics, function(y) {
          if (y %in% names(x)) return (x[y])
          return (0)
        })
      })
      attri_valuesM = t(attri_valuesM)
      if (all(colnames(attri_valuesM)==colnames(attri_valuesM)[1])) attri_valuesM = t(attri_valuesM)
      colnames(attri_valuesM) = attri_topics
      rownames(attri_valuesM) = c(1:nrow(attri_valuesM))
      avm = attri_valuesM/20
      
      png(paste0(clust_path_plot, "_stats_",attri,"_vs.png", sep=""), height=plot_size[1], width=plot_size[2])
      par(mar=c(10,5,5,8), xpd=TRUE)
      
      colour = rainbow(ncol(attri_valuesM))
      barplot(t(attri_valuesM), xlab="bicluster",ylab="# of samples", col=colour, main=paste0("# of samples in biclusters with attribute of ",attri,"\n",attribute))#,"; precision, recall, f measure comembership:\n",paste0(c(f1$p,f1$r,f1$f_comember),collapse=", ")))
      legend("topright",legend=colnames(attri_valuesM),fill=colour,inset=c(-.2,0))
      # plot(row(avm), col(avm),
      #   cex=avm,
      #   xlim=c(0.5,nrow(avm)+0.5), ylim=c(0.5,ncol(avm)+0.5),
      #   axes=FALSE, ann=FALSE
      # )
      # text(row(avm), col(avm), paste0(attri_valuesM, "\nsamples"), col="brown", pos=3)
      # axis(1,at=1:nrow(avm),labels=rownames(avm),cex.axis=0.8)
      # axis(2,at=1:ncol(avm),labels=colnames(avm),cex.axis=0.8)
      # title(xlab="bicluster",ylab=attri)
      # box()
      graphics.off()
    }
  }, error = function(err) { cat(paste("attribute plot error:  ",err)); return(T) })
  
  TimeOutput(start2)
}

TimeOutput(start)