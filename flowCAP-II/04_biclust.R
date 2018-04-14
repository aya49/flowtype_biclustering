## Input: original features --> Output: distance matrices
# aya43@sfu.ca 20161220

#Directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir,  "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir,  "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir,  "/matrix", sep="")

#Output
biclust_dir = paste(result_dir,  "/biclust_p6_F4-CD34Pos", sep=""); suppressWarnings(dir.create (biclust_dir))

#Libraries/Functions
library(biclust)
library(fabia)
library(foreach)
library(doMC)
library(stringr)
library(pheatmap)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
source("~/projects/IMPC/code/_bayesianbiclustering.R")

#Setup Cores
no_cores = 1#detectCores()-1
registerDoMC(no_cores)








#Options for script
overwrite = T #overwrite biclust?
writecsv = F

plot_size = c(500,500)
plot_size_bar = c(1300,2000)
attributes = c("aml")

cellCountThres = c(2000) #(needed if sample x cell pop matrices are used) insignificant if count under
good_samples = c(3)
control = "normal" #control value in target_col column for each centre
target_col = "aml"
# order_cols = c("barcode","sample","specimen","date")

bcmethods = c("plaid","CC","bimax","BB-binary","nmf-nsNMF","nmf-lee","nmf-brunet","CC")
#,"quest", "CC", "spectral", "Xmotifs", have to change this manually in function...
onlysigBB = T
Kr = 6; Kc = 20 #number of row and column biclusters for binary bayesian biclustering
pthres = .01 #pthreshold for choosing biclusters out of all bayesian biclusters
miniter = 100
nmf_thres = .05 # * max contribution: threshold at which a row/col can be considered a significant contribution to a factor

tube = 6 #which panel to use (flowCAP-II only); any of tubes 1-7

#data paths
#matrix_type = list.files(path=result_dir,pattern=glob2rx("matrix*.Rdata"))
matrix_type = c("Pval_CountAdj","PvalTRIM_CountAdj","CountAdj","CountAdjTRIM_CountAdj","Prop","PropTRIM_Prop")
matrix_count = c("CountAdj") #(needed if sample x cell pop matrices are used) count matrix, to get rid of low cell count cells










start = Sys.time()

#load count matrix and sample meta
mc = get(load(paste0(matrix_dir, matrix_count,".Rdata")))
sampleMeta = get(load(sampleMeta_dir))

f1scores = NULL

#for each feature matrix
# fe = foreach(mcp=matrix_type) %dopar% {
for (mcp in matrix_type) {
  tryCatch({
    cat("\n", mcp, " ",sep="")
    start2 = Sys.time()
    
    ## upload and prep feature matrix
    if (file.exists(paste0(matrix_dir, mcp,".Rdata"))) { m1 = get(load(paste0(matrix_dir, mcp,".Rdata")))
    } else if (file.exists(paste0(matrix_dir, mcp,".csv"))) { m1 = read.csv(paste0(matrix_dir, mcp,".csv", header=T)) } 
    
    m0 = m1
    
    #get matrix colnames (/celltype features)
    m0cn = colnames(m0); if (is.null(m0cn)) m0cn = names(m0)
    
    #get feature layers if feature names represent cell types
    colsplitlen = cell_type_layers(m0cn); k0 = 0; if (!is.null(colsplitlen)) k0 = c(1,2,3,4,max(colsplitlen)) #1,4, # how many layers to consider i.e. k=max(phenolevel) only
    
    #get to-delete low count phenotype indices; CountAdj should be first one
    for (countThres in cellCountThres) { cat("\ncountThres: ",countThres," >",sep="")
      
      #get to-delete high no of marker phenotypes
      for (k in k0) { cat(" level",k," ",sep="")
        
        mm = trimMatrix(m0, TRIM=grepl("TRIM",mcp), mc=mc,
                        sampleMeta=sampleMeta, sampleMeta_to_m1_col="fileName", target_col=target_col, control=str_split(control,"[|]")[[1]], order_cols=NULL, 
                        colsplitlen=NULL, k=k, countThres=countThres, goodcount=1, good_sample=1)
        m = mm$m
        sm = mm$sm
        
        #check matrix; if matrix is a list, merge; if matrix is empty, skip
        if (is.null(m)) next
        if (is.null(dim(m))) { m = Reduce('cbind',m); if (all(m==0)) next
        } else { if (all(m==0)) next }
        
        
        #trim matrix rows by tube
        m = m[sm$tube==tube,]
        sm = sm[sm$tube==tube,]
        
        
        rownames(m) = sm[,target_col]
        
        
        #for each biclustering method
        for (bcmethod in bcmethods) {
          if (bcmethod=="BB-binary" & !grepl("TRIM",mcp)) next
          #biclustering file name
          # dname = paste(biclust_dir, "/", mcp, "_", bcmethod, "_layer", str_pad(k, 2, pad = "0"), "_countThres-", countThres,"_KOsamplecountover-",good_sample, sep = "")
          dname = paste(biclust_dir, "/", mcp, "_", bcmethod, "_layer", str_pad(k, 2, pad = "0"), "_countThres-", countThres, sep = "")
          
          ## calculate biclusters
          if (overwrite | !file.exists(paste0(dname,".Rdata"))) {
            if (bcmethod == "plaid") d = biclust(as.matrix(m), method=BCPlaid(), row.release=.3,col.release=.7, back.fit=10, verbose = F)
            if (bcmethod == "CC") d = biclust(as.matrix(m), method=BCCC(), number=Kr)
            if (bcmethod == "Xmotifs") d = biclust(as.matrix(m), method=BCXmotifs(), number=Kr, ns=50, nd=500, alpha=10)
            if (bcmethod == "spectral") d = biclust(as.matrix(m), method=BCSpectral(), numberOfEigenvalues=10)
            if (bcmethod == "bimax") d = biclust(as.matrix(m), method=BCBimax(),number=Kr)
            if (bcmethod == "quest") d = biclust(as.matrix(m), method=BCQuest(), number=Kr, ns=50)
            
            if (bcmethod == "fabia") {
              next
              db = fabia(as.matrix(abs(m)), p=Kr,alpha=0.01,cyc=max(ncol(m)/10,500),spl=0,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0)
            }
            
            db = NULL
            try({
              if (bcmethod == "nmf-brunet") db = nmf(as.matrix(abs(m)), rank=Kr, method="brunet")#, nrun=10, method=list("lee", "brunet", "nsNMF"))
              if (bcmethod == "nmf-lee") db = nmf(as.matrix(abs(m)), rank=Kr, method="lee")
              if (bcmethod == "nmf-nsNMF") db = nmf(as.matrix(abs(m)), rank=Kr, method="nsNMF")
            })
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
              
              d = biclust(array(0,dim=c(2,2)), method=BCPlaid())
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
              db = B2PS(as.matrix(mbinary), sideData=NULL, Kt=Kr, Kp=Kc, iterations=max(ncol(mbinary)/20,miniter), alpha_p = 1, alpha_t = 1, alpha_e = 1, alpha_sd = 1)
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
              
              #convert format
              d = biclust(array(0,dim=c(2,2)), method=BCPlaid())
              d@info = db
              if (onlysigBB) {
                d@info$cid = cid = which(db$p<pthres, arr.ind=T)
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
            
            if (bcmethod != "BB-binary") d@info$rowclust = rowxcluster_to_cluster(d@RowxNumber)
            
            try({
              # get hard clustering; put samples into larger cluster
              la = as.numeric(factor(sm[,target_col]))
              cl = d@info$rowclust
              d@info$f1 = f.measure.comembership(la,cl)
            })
            
            
            save(d, file=paste0(dname,".Rdata")); if (writecsv) write.csv(as.matrix(d), file=paste0(checkm(d,dname),".Rdata"))
          } else {
            d = get(load(paste0(dname,".Rdata")))
            if (is.null(d)) next
          }
          
          #temporary
          if (nrow(d@RowxNumber) != nrow(m) | ncol(d@RowxNumber) != d@Number) d@RowxNumber = t(d@RowxNumber)
          if (ncol(d@NumberxCol) != ncol(m) | nrow(d@NumberxCol) != d@Number) d@NumberxCol = t(d@NumberxCol)
          if (is.null(d@info$rowclust) & bcmethod!="BB-binary") d@info$rowclust = rowxcluster_to_cluster(d@RowxNumber)
          try({
            # get hard clustering; put samples into larger cluster
            la = as.numeric(factor(sm[,target_col]))
            if (bcmethod=="BB-binary") {
              cl = d@info$rowclust = d@info$transcript.clusters
            } else {
              cl = d@info$rowclust
            }
            d@info$f1 = f.measure.comembership(la,cl)
          })
          save(d, file=paste0(dname,".Rdata")); if (writecsv) write.csv(as.matrix(d), file=paste0(checkm(d,dname),".Rdata"))
          
          
          
          ## plot
          if (d@Number > 0) {
            
            try ({
              f1 = d@info$f1
              f1scores = rbind(f1scores,c(feature=mcp,countThres=countThres,layer=k,method=bcmethod,unlist(f1)))
            })
            
            #pretty heatmap
            try ({#bcmethod=="BB-binary") {
              if (bcmethod=="BB-binary") {
                colcluster_temp = d@info$patient.clusters
                rowcluster_temp = d@info$transcript.clusters
                colcluster_temp[!colcluster_temp%in%d@info$cid[,2]] = 0
                rowcluster_temp[!rowcluster_temp%in%d@info$cid[,1]] = 0
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
  }, error = function(err) { cat(paste("ERROR:  ",err)); #return(T)
  })
  #return(F)
}

write.csv(f1scores,file=paste0(biclust_dir,"/scores.csv"))
TimeOutput(start)




