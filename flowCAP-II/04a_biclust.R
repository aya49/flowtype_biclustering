## Input: original features --> Output: biclusters
# aya43@sfu.ca 20161220

#Directory
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
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

overwrite = T #overwrite biclust?
writecsv = F

good_count = 3
good_sample = 3


plot_size = c(500,500)
plot_size_bar = c(1300,2000)

attributes = c("aml")

plot_size = c(500,500)
plot_size_bar = c(1300,2000)
attributes = c("aml")

cellCountThres = c(2000) #(needed if sample x cell pop matrices are used) insignificant if count under
control = "normal" #control value in target_col column for each centre
id_col = "fileName"
target_col = "aml"
order_cols = c("specimen","aml")
split_col = "tube"

bcmethods = c("plaid","CC","bimax","BB-binary","nmf-nsNMF","nmf-lee","nmf-brunet","CC")
#,"quest", "CC", "spectral", "Xmotifs", have to change this manually in function...
onlysigBB = T
Kr = 6; Kc = 20 #number of row and column biclusters for binary bayesian biclustering
pthres = .01 #pthreshold for choosing biclusters out of all bayesian biclusters
miniter = 100
nmf_thres = .05 # * max contribution: threshold at which a row/col can be considered a significant contribution to a factor

# tube = 6 #which panel to use (flowCAP-II only); any of tubes 1-7

#data paths
feat_types = list.files(path=feat_dir,pattern=".Rdata")
feat_types = gsub(".Rdata","",feat_types)
feat_count = "file-cell-countAdj"
#feat_types = list.files(path=result_dir,pattern=glob2rx("matrix*.Rdata"))
















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
    
    for (k in layers) {
      #trim matrix
      mm = trimMatrix(m0,TRIM=T, mc=mc, sampleMeta=meta_file, sampleMeta_to_m1_col=id_col, target_col=target_col, control=control, order_cols=order_cols, colsplitlen=NULL, k=k, countThres=countThres, goodcount=good_count, good_sample=good_sample)
      if (is.null(mm)) next
      m_ordered = mm$m
      meta_file_ordered = mm$sm
      if (is.null(mm)) next
      
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
        
        
        #for each biclustering method
        for (bcmethod in bcmethods) { #requires binary matrix
          if (bcmethod=="BB-binary" & !grepl("pval[A-z]*TRIM",feat_type)) next
          
          bcname0 = paste("/",bcmethod, "_", feat_type, "_splitby-",split_col,".", tube, "_layer", str_pad(k, 2, pad = "0"), "_countThres-", countThres, sep = "")
          bcname = paste0(biclust_source_dir, bcname0)
          
          ## bicluster ---------------------------------------
          if (overwrite | !file.exists(paste0(bcname,".Rdata"))) {
            if (bcmethod == "plaid") bc = biclust(as.matrix(m), method=BCPlaid(), row.release=.3,col.release=.7, back.fit=10, verbose = F)
            if (bcmethod == "CC") bc = biclust(as.matrix(m), method=BCCC(), number=Kr)
            if (bcmethod == "Xmotifs") bc = biclust(as.matrix(m), method=BCXmotifs(), number=Kr, ns=50, nd=500, alpha=10)
            if (bcmethod == "spectral") bc = biclust(as.matrix(m), method=BCSpectral(), numberOfEigenvalues=10)
            if (bcmethod == "bimax") bc = biclust(as.matrix(m), method=BCBimax(),number=Kr)
            if (bcmethod == "quest") bc = biclust(as.matrix(m), method=BCQuest(), number=Kr, ns=50)
            
            if (bcmethod == "fabia") { next #don't do, do nmf instead
              bcb = fabia(as.matrix(abs(m)), p=Kr,alpha=0.01,cyc=max(ncol(m)/10,500),spl=0,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0)
            }
            
            bcb = NULL
            try({
              if (bcmethod == "nmf-brunet") bcb = nmf(as.matrix(abs(m)), rank=Kr, method="brunet")#, nrun=10, method=list("lee", "brunet", "nsNMF"))
              if (bcmethod == "nmf-lee") bcb = nmf(as.matrix(abs(m)), rank=Kr, method="lee")
              if (bcmethod == "nmf-nsNMF") bcb = nmf(as.matrix(abs(m)), rank=Kr, method="nsNMF")
            })
            
            ## adjust format / special plot if needed
            if (grepl("nmf",bcmethod)) {
              if (!is.null(bcb)) {
                bcb$rowxfactor = basis(bcb)
                bcb$factorxcol = coef(bcb)
                
                rthres = max(bcb$rowxfactor) * nmf_thres
                cthres = max(bcb$factorxcol) * nmf_thres
                
                bc = biclust(array(0,dim=c(2,2)), method=BCPlaid()) #get the framwork
                bc@RowxNumber = array(F, dim=dim(bcb$rowxfactor)) 
                for(ri in 1:nrow(bcb$rowxfactor)) {
                  mi = which.max(bcb$rowxfactor[ri,])
                  if (max(bcb$rowxfactor[ri,])>rthres) bc@RowxNumber[ri,mi] = T
                }
                bc@NumberxCol = array(F, dim=dim(bcb$factorxcol)) 
                for(ri in 1:ncol(bcb$factorxcol)) {
                  mi = which.max(bcb$factorxcol[,ri])
                  if (max(bcb$factorxcol[,ri])>cthres) bc@NumberxCol[mi,ri] = T
                }
                bc@Number = nrow(bcb$factorxcol)
                
                bc@info = list(rowxfactor=bcb$rowxfactor,factoxcol=bcb$factorxcol)
              } else { next }
            }
            
            if (bcmethod == "BB-binary") {
              mbinary = m
              mbinary[mbinary != 0] = 1 #make matrix binary (for p values TRIM only)
              bcb = B2PS(as.matrix(mbinary), sideData=NULL, Kt=Kr, Kp=Kc, iterations=max(ncol(mbinary)/20,min_iter), alpha_p = 1, alpha_t = 1, alpha_e = 1, alpha_sd = 1)
              # THETA is trasnposed!
              bcb$theta = t(bcb$Theta[,,2])
              
              #get significant biclusters
              theta = bcb$theta
              theta = theta[sort(unique(bcb$transcript.clusters)),sort(unique(bcb$patient.clusters))]
              bcb$p = array(1,dim = c(nrow(bcb$theta), ncol(bcb$theta)))
              for (i in unique(bcb$transcript.clusters)) {
                for (j in unique(bcb$patient.clusters)) {
                  bcb$p[i,j] = t.test.single(as.vector(theta),bcb$theta[i,j])
                  if (is.na(bcb$p[i,j])) bcb$p[i,j] = 1
                }
              }
              
              # convert format
              bc = biclust(array(0,dim=c(2,2)), method=BCPlaid())
              bc@info = bcb
              if (onlysigBB) {
                bc@info$cid = cid = which(bcb$p<pval_thres, arr.ind=T)
                if (nrow(cid) == 0) {
                  if (length(theta)>=6 & !all(bcb$p==1)) {
                    bc@info$cid = cid = which(bcb$p<=sort(bcb$p)[4], arr.ind=T)
                  } else {
                    save(NULL, file=paste0(bcname,".Rdata"))
                    next
                  }
                }
              } else {
                bc@info$cid = cid = which(bcb$p<=1, arr.ind=T)
              }
              
              # cidmatrix = data.frame(row=rep(0,nrow(cid)), col=rep(0,nrow(cid)))
              # cidmatrix[cid] = 1
              # bc@info$sigcluster = cidmatrix
              bc@Number = nrow(cid)
              bc@RowxNumber = array(F, dim=c(nrow(m),nrow(cid))); colnames(bc@RowxNumber) = rep(1,ncol(bc@RowxNumber))
              bc@NumberxCol = array(F, dim=c(nrow(cid), ncol(m))); rownames(bc@NumberxCol) = rep(1,nrow(bc@NumberxCol))
              for (i in 1:nrow(cid)) { #column = cluster
                rowclust = cid[i,1]
                colclust = cid[i,2]
                colnames(bc@RowxNumber)[i] = rownames(bc@NumberxCol)[i] = bcb$p[rowclust,colclust]
                rows = bcb$transcript.clusters==rowclust
                cols = bcb$patient.clusters==colclust
                bc@RowxNumber[rows, i] = T
                bc@NumberxCol[i, cols] = T
                # cat("\n", sum(bcb$transcript.clusters==rowclust), ", ", sum(bcb$patient.clusters==colclust), sep="")
              }
              rowclust = bc@info$transcript.clusters
              colclust = bc@info$patient.clusters
            }
            
            if (nrow(bc@RowxNumber) != nrow(m) | ncol(bc@RowxNumber) != bc@Number) bc@RowxNumber = t(bc@RowxNumber)
            if (ncol(bc@NumberxCol) != ncol(m) | nrow(bc@NumberxCol) != bc@Number) bc@NumberxCol = t(bc@NumberxCol)
            
            
            
            ## get & save clustering & labels! ----------------------------------
            
            if (bcmethod != "BB-binary") {
              rowclust = rowxcluster_to_cluster(bc@RowxNumber)
              colclust = clusterxcol_to_cluster(bc@NumberxCol)
            } 
            rowlabel = sm[,target_col]
            
            names(rowclust) = names(rowlabel) = rownames(m)
            names(colclust) = colnames(m)
            f1 = f.measure.comembership(rowlabel,rowclust)
            
            bc0 = list(source=bc,rowclust=rowclust,colclust=colclust,rowlabel=rowlabel,scores=f1)
            save(bc0, file=paste0(bcname,".Rdata"))
            write.csv(rowclust, file=paste0(bcname,"_rowclust.csv"))
            write.csv(colclust, file=paste0(bcname,"_colclust.csv"))
            write.csv(rowlabel, file=paste0(bcname,"_rowlabel.csv"))
            
            # #save row/col as csv
            # try({
            #   bcgene = bc@RowxNumber; rownames(bcgene) = sm[,attributes[1]]
            #   bcgene = bcgene[apply(bcgene, 1, function(x) all(!x)),]
            #   write.csv(bcgene,file=paste0(bcname,"_row.csv"))
            # })
            # try ({
            #   bccol = bc@NumberxCol
            #   if (ncol(bccol)==bc@Number) bccol = t(bccol)
            #   colnames(bccol) = colnames(m)
            #   bccol = bccol[,apply(bccol, 2, function(x) all(!x))]
            #   write.csv(bccol,file=paste0(bcname,"_col.csv"))
            # })
            
          }
          
          # #temporary
          # if (nrow(bc@RowxNumber) != nrow(m) | ncol(bc@RowxNumber) != bc@Number) bc@RowxNumber = t(bc@RowxNumber)
          # if (ncol(bc@NumberxCol) != ncol(m) | nrow(bc@NumberxCol) != bc@Number) bc@NumberxCol = t(bc@NumberxCol)
          # if (is.null(bc0$rowclust) & bcmethod!="BB-binary") bc0$rowclust = rowxcluster_to_cluster(bc@RowxNumber)
          # try({
          #   # get hard clustering; put samples into larger cluster
          #   la = as.numeric(factor(sm[,target_col]))
          #   if (bcmethod=="BB-binary") {
          #     cl = bc0$rowclust = bc@info$transcript.clusters
          #   } else {
          #     cl = bc0$rowclust
          #   }
          #   bc@info$f1 = f.measure.comembership(la,cl)
          # })
          # save(bc, file=paste0(bcname,".Rdata")); if (writecsv) write.csv(as.matrix(bc), file=paste0(checkm(bc,bcname),".Rdata"))
          
          rowlabel = as.numeric(factor(sm[,target_col]))
          
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




