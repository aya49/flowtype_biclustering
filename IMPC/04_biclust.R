## Input: original features --> Output: bicluster & plots
# aya43@sfu.ca 20161220

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
overwrite = T #overwrite biclust?
writecsv = F

good_count = 3
good_sample = 3


plot_size = c(500,500)
plot_size_bar = c(1300,2000)
attributes = c("gene","gender","date")

cellCountThres = c(200) #(needed if sample x cell pop matrices are used) insignificant if count under
good_samples = c(0,3)
control = str_split(controlL,"[|]")[[1]]
time_col = "date"
id_col = "fileName"
target_col = "gene"
order_cols = c("barcode","sample","specimen","date")
split_col = NULL

bcmethods = c("plaid","CC","bimax","BB-binary","nmf-nsNMF","nmf-lee","nmf-brunet","CC")
onlysigBB = T #only evaluate significant BB-binary clusters
Kr = 20; Kc = 20 #number of row and column biclusters for binary bayesian biclustering
pval_thres = .01 #pval_threshold for choosing biclusters out of all bayesian biclusters
min_iter = 100 #min number of iterations for BB-binary (B2PS)
nmf_thres = .05 # * max contribution: threshold at which a row/col can be considered a significant contribution to a factor in nmf


#data paths
feat_types = list.files(path=feat_dir,pattern=".Rdata")
feat_count = "file-cell-countAdj"
#feat_types = list.files(path=result_dir,pattern=glob2rx("matrix*.Rdata"))
# feat_types = c("PvalTRIM_CountAdjBH"                       ,"PvalTRIM_CountAdjbonferroni"               ,"PvalTRIM_CountAdjBY"                       ,"PvalTRIM_CountAdjPEERBH"                  
#                 ,"PvalTRIM_CountAdjPEERbonferroni"           ,"PvalTRIM_CountAdjPEERBY"                   ,"PvalTRIM_CountAdjPEER"                     ,"PvalTRIM_CountAdj"       
#                 ,"Pval_CountAdjBH"                           ,"Pval_CountAdjbonferroni"                   ,"Pval_CountAdjBY"                           ,"Pval_CountAdjPEERBH"                      
#                 ,"Pval_CountAdjPEERbonferroni"               ,"Pval_CountAdjPEERBY"                       ,"Pval_CountAdjPEER"                         ,"Pval_CountAdj"                            
#                 ,"Child_entropyTRIM_CountAdjBH"              ,"Child_entropyTRIM_CountAdjbonferroni"      ,"Child_entropyTRIM_CountAdjBY"              ,"Child_entropyTRIM_CountAdjPEERBH"         
#                 ,"CountAdjPEER"                              ,"CountAdj"                                  ,"LogFold_CountAdjPEER"                      ,"LogFold_CountAdj"                         
#                 ,"LogFoldTRIM_CountAdjBH"                    ,"LogFoldTRIM_CountAdjbonferroni"            ,"LogFoldTRIM_CountAdjBY"                    ,"LogFoldTRIM_CountAdjPEERBH"               
#                 ,"LogFoldTRIM_CountAdjPEERbonferroni"        ,"LogFoldTRIM_CountAdjPEERBY"                ,"LogFoldTRIM_CountAdjPEER"                  ,"LogFoldTRIM_CountAdj"                     
#                 ,"Child_entropyTRIM_CountAdjPEERbonferroni"  ,"Child_entropyTRIM_CountAdjPEERBY"          ,"Child_entropyTRIM_CountAdjPEER"            ,"Child_entropyTRIM_CountAdj"               
#                 ,"Parent_entropyTRIM_CountAdjBH"             ,"Parent_entropyTRIM_CountAdjbonferroni"     ,"Parent_entropyTRIM_CountAdjBY"             ,"Parent_entropyTRIM_CountAdjPEERBH"        
#                 ,"Parent_entropyTRIM_CountAdjPEERbonferroni" ,"Parent_entropyTRIM_CountAdjPEERBY"         ,"Parent_entropyTRIM_CountAdjPEER"           ,"Parent_entropyTRIM_CountAdj"              
#                 ,"Child_pnratioTRIM_CountAdjBH"              ,"Child_pnratioTRIM_CountAdjbonferroni"      ,"Child_pnratioTRIM_CountAdjBY"              ,"Child_pnratioTRIM_CountAdjPEERBH"         
#                 ,"Child_pnratioTRIM_CountAdjPEERbonferroni"  ,"Child_pnratioTRIM_CountAdjPEERBY"          ,"Child_pnratioTRIM_CountAdjPEER"            ,"Child_pnratioTRIM_CountAdj"               
#                 ,"Child_propTRIM_CountAdjBH"                 ,"Child_propTRIM_CountAdjbonferroni"         ,"Child_propTRIM_CountAdjBY"                 ,"Child_propTRIM_CountAdjPEERBH"            
#                 ,"Child_propTRIM_CountAdjPEERbonferroni"     ,"Child_propTRIM_CountAdjPEERBY"             ,"Child_propTRIM_CountAdjPEER"               ,"Child_propTRIM_CountAdj"                  
#                 ,"Parent_contrib_CountAdjPEER"               ,"Parent_contrib_CountAdj"                   ,"Parent_contribTRIM_CountAdjBH"             ,"Parent_contribTRIM_CountAdjbonferroni"    
#                 ,"Parent_contribTRIM_CountAdjBY"             ,"Parent_contribTRIM_CountAdjPEERBH"         ,"Parent_contribTRIM_CountAdjPEERbonferroni" ,"Parent_contribTRIM_CountAdjPEERBY"        
#                 ,"Parent_contribTRIM_CountAdjPEER"           ,"Parent_contribTRIM_CountAdj"               ,"Parent_effort_CountAdjPEER"                ,"Parent_effort_CountAdj"                   
#                 ,"Parent_effortTRIM_CountAdjBH"              ,"Parent_effortTRIM_CountAdjbonferroni"      ,"Parent_effortTRIM_CountAdjBY"              ,"Parent_effortTRIM_CountAdjPEERBH"         
#                 ,"Parent_effortTRIM_CountAdjPEERbonferroni"  ,"Parent_effortTRIM_CountAdjPEERBY"          ,"Parent_effortTRIM_CountAdjPEER"            ,"Parent_effortTRIM_CountAdj"               
# )
# matrix_count = c("CountAdj") #(needed if sample x cell pop matrices are used) count matrix, to get rid of low cell count cells
















start = Sys.time()

centre = paste0(panel," ",centre)
cat("\n",centre,sep="")


mc = get(load(paste0(feat_dir,"/", feat_count,".Rdata")))
meta_file = get(load(paste0(meta_file_dir,".Rdata")))

f1scores0 = foreach(feat_type=feat_types) %dopar% {
  f1scores = NULL
  tryCatch({
    cat("\n", feat_type, " ",sep="")
    start2 = Sys.time()
    
    ## upload and prep feature matrix
    m0 = as.matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata"))))
    layers = 0
    countThres = 0
    if (str_split(feat_type,"-")[[1]][2]=="cell") {
      layers = c(1,2,4,max(unique(sapply(unlist(str_split(colnames(m0),"_")), function(x) str_count(x,"[+-]")))))
      countThres = cellCountThres
    }
    
    for (layer in layers) {
      #trim matrix
      mm = trimMatrix(m0,TRIM=T, mc=mc, sampleMeta=meta_file, sampleMeta_to_m1_col=id_col, target_col=target_col, control=control, order_cols=order_cols, colsplitlen=NULL, k=layer, countThres=countThres, goodcount=good_count, good_sample=good_sample)
      if (is.null(mm)) next
      m_ordered = mm$m
      meta_file_ordered = mm$sm
      
      #split up analysis by tube etc.
      if (is.null(split_col)) {
        split_ind = list(all = 1:nrow(meta_file_ordered))
      } else {
        split_ids = unique(meta_file_ordered[,split_col])
        split_ids = split_ids[!is.na(split_ids)]
        split_ind = lapply(split_ids, function(split_id) which(meta_file_ordered[,split_col]==split_id) )
        names(split_ind) = split_ids
      }
      
      for (tube in names(split_ind)) {
        m = m_ordered[split_ind[[tube]],]
        if (!sum(m_ordered!=0)>0) next
        sm = meta_file_ordered_split = meta_file_ordered[split_ind[[tube]],]
        
        #check matrix; if matrix is a list, merge; if matrix is empty, skip
        if (is.null(m)) next
        
        #for each biclustering method
        for (bcmethod in bcmethods) { #requires binary matrix
          if (bcmethod=="BB-binary" & !grepl("TRIM",feat_type)) next
          
          dname0 = paste("/",bcmethod, "_", feat_type, "_splitby-", tube, "_layer", str_pad(k, 2, pad = "0"), "_countThres-", countThres, sep = "")
          dname = paste0(biclust_source_dir, ,dname0)
          
          ## bicluster ---------------------------------------
          if (overwrite | !file.exists(paste0(dname,".Rdata"))) {
            if (bcmethod == "plaid") d = biclust(as.matrix(m), method=BCPlaid(), row.release=.3,col.release=.7, back.fit=10, verbose = F)
            if (bcmethod == "CC") d = biclust(as.matrix(m), method=BCCC(), number=Kr)
            if (bcmethod == "Xmotifs") d = biclust(as.matrix(m), method=BCXmotifs(), number=Kr, ns=50, nd=500, alpha=10)
            if (bcmethod == "spectral") d = biclust(as.matrix(m), method=BCSpectral(), numberOfEigenvalues=10)
            if (bcmethod == "bimax") d = biclust(as.matrix(m), method=BCBimax(),number=Kr)
            if (bcmethod == "quest") d = biclust(as.matrix(m), method=BCQuest(), number=Kr, ns=50)
            
            if (bcmethod == "fabia") { next #don't do, do nmf instead
              db = fabia(as.matrix(abs(m)), p=Kr,alpha=0.01,cyc=max(ncol(m)/10,500),spl=0,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0)
            }
            
            db = NULL
            try({
              if (bcmethod == "nmf-brunet") db = nmf(as.matrix(abs(m)), rank=Kr, method="brunet")#, nrun=10, method=list("lee", "brunet", "nsNMF"))
              if (bcmethod == "nmf-lee") db = nmf(as.matrix(abs(m)), rank=Kr, method="lee")
              if (bcmethod == "nmf-nsNMF") db = nmf(as.matrix(abs(m)), rank=Kr, method="nsNMF")
            })
            
            
            
            
            ## adjust format / special plot if needed
            
            if (grepl("nmf",bcmethod) & !is.null(db)) {
              db$rowxfactor = basis(db)
              db$factorxcol = coef(db)
              
              png(paste0(dname, "_rowxfactor.png", sep=""), height=500*2, width=600*2)
              par(mfrow=c(2,2))
              plot(sort(db$rowxfactor[,1]),type="l", main="contribution of rows to each factor; each line = factor" )
              for (fi in 2:ncol(db$rowxfactor)) {
                lines(sort(db$rowxfactor[,fi]))
              }
              plot(sort(db$factorxcol[1,]),type="l", main="contribution of cols to each factor; each line = factor" )
              for (fi in 2:nrow(db$factorxcol)) {
                lines(sort(db$factorxcol[fi,]))
              }
              
              aheatmap(db$rowxfactor, main="row x factor")
              aheatmap(db$factorxcol, main="factor x col")
              graphics.off()
              
              rthres = max(db$rowxfactor) * nmf_thres
              cthres = max(db$factorxcol) * nmf_thres
              
              d = biclust(array(0,dim=c(2,2)), method=BCPlaid()) #get the framwork
              d@RowxNumber = array(F, dim=dim(db$rowxfactor)) 
              for(ri in 1:nrow(db$rowxfactor)) {
                mi = which.max(db$rowxfactor[ri,])
                if (max(db$rowxfactor[ri,])>rthres) d@RowxNumber[ri,mi] = T
              }
              d@NumberxCol = array(F, dim=dim(db$factorxcol)) 
              for(ri in 1:ncol(db$factorxcol)) {
                mi = which.max(db$factorxcol[,ri])
                if (max(db$factorxcol[,ri])>cthres) d@NumberxCol[mi,ri] = T
              }
              d@Number = nrow(db$factorxcol)
              
              d@info = list(rowxfactor=db$rowxfactor,factoxcol=db$factorxcol)
            } else if (is.null(db)) {
              next
            }
            
            if (bcmethod == "BB-binary") {
              mbinary = m
              mbinary[mbinary != 0] = 1 #make matrix binary (for p values TRIM only)
              db = B2PS(as.matrix(mbinary), sideData=NULL, Kt=Kr, Kp=Kc, iterations=max(ncol(mbinary)/20,min_iter), alpha_p = 1, alpha_t = 1, alpha_e = 1, alpha_sd = 1)
              # THETA is trasnposed!
              db$theta = t(db$Theta[,,2])
              
              #get significant biclusters
              theta = db$theta
              theta = theta[sort(unique(db$transcript.clusters)),sort(unique(db$patient.clusters))]
              db$p = array(1,dim = c(nrow(db$theta), ncol(db$theta)))
              for (i in unique(db$transcript.clusters)) {
                for (j in unique(db$patient.clusters)) {
                  db$p[i,j] = t.test.single(as.vector(theta),db$theta[i,j])
                  if (is.na(db$p[i,j])) db$p[i,j] = 1
                }
              }
              
              # convert format
              d = biclust(array(0,dim=c(2,2)), method=BCPlaid())
              d@info = db
              if (onlysigBB) {
                d@info$cid = cid = which(db$p<pval_thres, arr.ind=T)
                if (nrow(cid) == 0) {
                  if (length(theta)>=6 & !all(db$p==1)) {
                    d@info$cid = cid = which(db$p<=sort(db$p)[4], arr.ind=T)
                  } else {
                    save(NULL, file=paste0(dname,".Rdata"))
                    next
                  }
                }
              } else {
                d@info$cid = cid = which(db$p<=1, arr.ind=T)
              }
              
              # cidmatrix = data.frame(row=rep(0,nrow(cid)), col=rep(0,nrow(cid)))
              # cidmatrix[cid] = 1
              # d@info$sigcluster = cidmatrix
              d@Number = nrow(cid)
              d@RowxNumber = array(F, dim=c(nrow(m),nrow(cid))); colnames(d@RowxNumber) = rep(1,ncol(d@RowxNumber))
              d@NumberxCol = array(F, dim=c(nrow(cid), ncol(m))); rownames(d@NumberxCol) = rep(1,nrow(d@NumberxCol))
              for (i in 1:nrow(cid)) { #column = cluster
                rowclust = cid[i,1]
                colclust = cid[i,2]
                colnames(d@RowxNumber)[i] = rownames(d@NumberxCol)[i] = db$p[rowclust,colclust]
                rows = db$transcript.clusters==rowclust
                cols = db$patient.clusters==colclust
                d@RowxNumber[rows, i] = T
                d@NumberxCol[i, cols] = T
                # cat("\n", sum(db$transcript.clusters==rowclust), ", ", sum(db$patient.clusters==colclust), sep="")
              }
              d@info$rowclust = d@info$transcript.clusters
            }
            
            if (nrow(d@RowxNumber) != nrow(m) | ncol(d@RowxNumber) != d@Number) d@RowxNumber = t(d@RowxNumber)
            if (ncol(d@NumberxCol) != ncol(m) | nrow(d@NumberxCol) != d@Number) d@NumberxCol = t(d@NumberxCol)
            
            ## get clustering & labels!
            
            if (bcmethod != "BB-binary") {
              clrow = rowxcluster_to_cluster(d@RowxNumber)
              clcol = clusterxcol_to_cluster(d@NumberxCol)
            } 
            
            rowlabel = as.numeric(factor(sm[,target_col]))
            d@info$f1 = f.measure.comembership(la,cl)
            
            d0 = d
            
            names(rowclust) = names(rowlabel) = rownames(m)
            names(colclust) = colnames(m)
            d = list(source=d0,rowclust=rowclust,colclust=colclust,rowlabel=rowlabel)
            
            save(d, file=paste0(dname,".Rdata")); if (writecsv) write.csv(as.matrix(d), file=paste0(checkm(d,dname),".Rdata"))
          } else {
            d0 = get(load(paste0(dname,".Rdata")))
            if (is.null(d0)) next
            d = d0$source
          }
          
          # #temporary
          # if (nrow(d@RowxNumber) != nrow(m) | ncol(d@RowxNumber) != d@Number) d@RowxNumber = t(d@RowxNumber)
          # if (ncol(d@NumberxCol) != ncol(m) | nrow(d@NumberxCol) != d@Number) d@NumberxCol = t(d@NumberxCol)
          # if (is.null(d0$rowclust) & bcmethod!="BB-binary") d0$rowclust = rowxcluster_to_cluster(d@RowxNumber)
          # try({
          #   # get hard clustering; put samples into larger cluster
          #   la = as.numeric(factor(sm[,target_col]))
          #   if (bcmethod=="BB-binary") {
          #     cl = d0$rowclust = d@info$transcript.clusters
          #   } else {
          #     cl = d0$rowclust
          #   }
          #   d@info$f1 = f.measure.comembership(la,cl)
          # })
          # save(d, file=paste0(dname,".Rdata")); if (writecsv) write.csv(as.matrix(d), file=paste0(checkm(d,dname),".Rdata"))
          
          
          
          ## plot ---------------------------------
          if (d@Number > 0) {
            try ({
              f1 = d@info$f1
              f1scores = rbind(f1scores,c(feature=feat_type,countThres=countThres,layer=k,method=bcmethod,unlist(f1)))
            })
            
            #pretty heatmap
            try ({#bcmethod=="BB-binary") {
              if (bcmethod=="BB-binary") {
                colcluster_temp = d0@info$patient.clusters
                rowcluster_temp = d0@info$transcript.clusters
                colcluster_temp[!colcluster_temp%in%d0$colclust] = 0
                rowcluster_temp[!rowcluster_temp%in%d0$rowclust] = 0
              } else {
                colcluster_temp = unlist(apply(d@NumberxCol, 2, function(x) {
                  clust = which(x)
                  if (length(clust) == 0) return(0)
                  return(max(clust))
                }))
                rowcluster_temp = unlist(apply(d@RowxNumber, 1, function(x) {
                  clust = which(x)
                  if (length(clust) == 0) return(0)
                  return(max(clust))
                }))
              }
              
              for (attribute in attributes) {
                
                col_annot = data.frame(colcluster=factor(colcluster_temp))
                row_annot = data.frame(rowcluster=factor(rowcluster_temp), class=factor(sm[,attribute]))
                rownames(col_annot)=colnames(m)
                
                row_order = order(row_annot$rowcluster, decreasing=T)
                col_order = order(col_annot$colcluster, decreasing=T)
                annotation_row = row_annot[row_order,]
                annotation_col = data.frame(col_annot[col_order,])
                
                mh = as.matrix(m)[row_order,col_order]
                rownames(mh) = 1:nrow(mh)
                rownames(annotation_col) = colnames(mh)
                rownames(annotation_row) = 1:nrow(mh)
                
                pheatmap(mh, main=paste0(attribute,"; precision, recall, f measure comembership:\n",paste0(c(f1$p,f1$r,f1$f_comember),collapse=", ")),
                         annotation_row = annotation_row, annotation_col = annotation_col,
                         show_rownames = F, 
                         show_colnames=F,
                         cluster_cols = T, cluster_rows = F,
                         # cellwidth = 3, cellheight = 3, 
                         # fontsize = 3, 
                         filename = paste0(dname, "_pheatmap-",attribute,"-sorted.pdf"))
                pheatmap(mh, main=paste0(attribute,"; precision, recall, f measure comembership:\n",paste0(c(f1$p,f1$r,f1$f_comember),collapse=", ")),
                         annotation_row = annotation_row, annotation_col = annotation_col,
                         show_rownames = F, 
                         show_colnames=F,
                         cluster_cols = F, cluster_rows = F,
                         # cellwidth = 3, cellheight = 3, 
                         # fontsize = 3, 
                         filename = paste0(dname, "_pheatmap-",attribute,"-sorted.png"))
                
                pheatmap(mh, main=paste0(attribute,"; precision, recall, f measure comembership:\n",paste0(c(f1$p,f1$r,f1$f_comember),collapse=", ")),
                         annotation_row = annotation_row, annotation_col = annotation_col,
                         show_rownames = F, 
                         # show_colnames=F,
                         cluster_cols = F, cluster_rows = F,
                         # cellwidth = 3, cellheight = 3, 
                         fontsize = 3, 
                         filename = paste0(dname, "_pheatmap-",attribute,"-sorted-smallfont.pdf"))
                
              }
              # # Specifying clustering from distance matrix
              # drows = dist(test, method = "minkowski")
              # dcols = dist(t(test), method = "minkowski")
              # pheatmap(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)
              
            })
            
            
            #save row/col as csv
            try({
              dgene = d@RowxNumber; rownames(dgene) = sm[,attributes[1]]
              dgene = dgene[apply(dgene, 1, function(x) all(!x)),]
              write.csv(dgene,file=paste0(dname,"_row.csv"))
            })
            try ({
              dcol = d@NumberxCol
              if (ncol(dcol)==d@Number) dcol = t(dcol)
              colnames(dcol) = colnames(m)
              dcol = dcol[,apply(dcol, 2, function(x) all(!x))]
              write.csv(dcol,file=paste0(dname,"_col.csv"))
            })
            
            # png(paste0(dname, "_bar.png", sep=""), height=plot_size_bar[1], width=plot_size_bar[2])
            # par(mar=c(2,6,3,2))
            # try({
            #   plotclust(d,as.matrix(m))
            # })
            # graphics.off()
            
            # png(paste0(dname, "_bubble.png", sep=""), height=plot_size_bar[1], width=plot_size_bar[2])
            # par(mar=c(20,20,20,20))
            # try({
            #   bubbleplot(as.matrix(m),d)
            # })
            # graphics.off()
            
            png(paste0(dname, "_heatmap0.png", sep=""), height=plot_size_bar[1], width=plot_size_bar[2])
            par(mar=c(5,3,6,5))
            try({
              heatmapBC(as.matrix(m),d, order=T, local=T, outside=T)
            })
            graphics.off()
            
            rowcolno = ceiling(sqrt(d@Number))
            png(paste0(dname, "_heatmap.png", sep=""), height=plot_size[1]*rowcolno, width=plot_size[2]*rowcolno)
            par(mar=c(25,50,20,50))
            par(mfrow=rep(rowcolno,2))
            for (BCi in 1:d@Number) {
              try({
                drawHeatmap(as.matrix(m),d,BCi,plotAll=T)
              })
            }
            graphics.off()
            
            for (attri in attributes) {
              png(paste0(dname, "_stats_",attri,".png", sep=""), height=plot_size[1]*rowcolno, width=plot_size[2]*rowcolno)
              par(mar=c(10,5,3,2), mfrow=rep(rowcolno,2))
              attri_valuesL = list()
              attri_topics = c()
              for (BCi in 1:d@Number) {
                attri_values = sm[d@RowxNumber[,BCi],attri]
                attri_valuesL[[BCi]] = attri_valuesT = table(attri_values)
                attri_topics = union(attri_topics,names(attri_valuesT))
                try({
                  barplot(attri_valuesT, las=2, xlab=attri, ylab="# of samples in bicluster with attribute on x axis", main=paste0(length(attri_values)," samples in bicluster ", BCi))
                })
              }
              graphics.off()
              
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
              
              png(paste0(dname, "_stats_",attri,"_vs.png", sep=""), height=plot_size[1], width=plot_size[2])
              par(mar=c(10,5,5,8), xpd=TRUE)
              
              colour = rainbow(ncol(attri_valuesM))
              barplot(t(attri_valuesM), xlab="bicluster",ylab="# of samples", col=colour, main=paste0("# of samples in biclusters with attribute of ",attri,"\n",attribute,"; precision, recall, f measure comembership:\n",paste0(c(f1$p,f1$r,f1$f_comember),collapse=", ")))
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
          }
          # }
          
        }
        
        
      } #layer
    } #countThres
    TimeOutput(start2)
  }, error = function(err) { cat(paste("ERROR:  ",err)); return(NULL) })
  return(f1scores)
}
f1scores0 = Reduce("rbind",f1scores0)
write.csv(f1scores0,file=paste0(biclust_dir,"/scores.csv"))



TimeOutput(start)




