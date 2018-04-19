## Input: original features --> Output: biclusters
# aya43@sfu.ca 20161220

## root directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

## input directories
meta_dir = paste0(result_dir,"/meta") # meta files directory
meta_file_dir = paste(meta_dir, "/file", sep="") #meta for rows (sample)
meta_cell_child_names_dir = paste(meta_dir, "/cell_child_names", sep="") #lists children cell population of each cell population; only used for graph regularized non matrix factorization
feat_dir = paste(result_dir, "/feat", sep="") #feature files directory

## output directories
biclust_dir = paste(result_dir,  "/biclust", sep=""); dir.create (biclust_dir,showWarnings=F)
biclust_source_dir = paste(biclust_dir,  "/source", sep=""); dir.create (biclust_source_dir,showWarnings=F) # path to save biclusterings: row clusterings, row labels, column colusterings

## libraries
library(biclust)
library(NMF)
library(GrNMF) #library(devtools); install_github("jstjohn/GrNMF")
library(pheatmap)
library(foreach)
library(doMC)
library(stringr)
library(Matrix)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")
source("~/projects/IMPC/code/_bayesianbiclustering.R")

## setup Cores for parallel processing (parallelized for each feature)
no_cores = 3#detectCores()-3
registerDoMC(no_cores)








## options for script
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

readcsv = F #read features as csv or Rdata
overwrite = T #overwrite biclust?
# writecsv = F

good_count = 1 #trim matrix; only keep col/rows that meet criteria for more than 3 elements
good_sample = 3 #trim matrix; only keep rows that are a part of a class with more than 3 samples


cellCountThres = c(2000) #a cell is insignificant if count under cell CountThres so delete -- only for matrices that have cell populations as column names
target_col = "aml" #the interested column in meta_file
control = "normal" #control value in target_col column
id_col = "fileName" #the column in meta_file matching rownames in feature matrices
order_cols = NULL #if matrix rows should be ordered by a certain column
split_col = "tube" # if certain rows in matrices should be analyzed in isolation, split matrix by this column in meta_file

bcmethods = c("plaid","CC","bimax","BB-binary","nmf-nsNMF","nmf-lee","nmf-brunet","CC","GrNMF-0","GrNMF-1","GrNMF-5","GrNMF-10") #biclustering methods; GrNMF-<weight of graph regularization>
#,"quest", "CC", "spectral", "Xmotifs", have to change this manually in function...
onlysigBB = T #only extract significant or all biclusters from BB-binary B2PS biclustering?
pval_thres = .05
Kr = 6; Kc = 20 #number of row and column biclusters for binary bayesian biclustering
pthres = .01 #pthreshold for choosing biclusters out of all bayesian biclusters
min_iter = 100 #BB-binary
sig_biclust_thres = .025 # * max contribution: threshold at whifeature matrices
qthres = .15 # quantile of nmf type methods factors; how large does factor effect need to be for bicluster to be significant


#feature matrix paths
if (readcsv) {
  feat_types = list.files(path=feat_dir,pattern=".csv")
  feat_types = gsub(".csv","",feat_types)
} else {
  feat_types = list.files(path=feat_dir,pattern=".Rdata")
  feat_types = gsub(".Rdata","",feat_types)
}
feat_count = "file-cell-countAdj" # cell count features used to trim matrix
#feat_types = list.files(path=result_dir,pattern=glob2rx("matrix*.Rdata"))
















start = Sys.time()

# read meta_file and count matrix (count matrix only to trim feature matrix of small cell populations -- only for matrices that have cell populations as column names)
if (readcsv) {
  mc = read.csv(paste0(feat_dir,"/", feat_count,".csv"),row.names=1, check.names=F)
  meta_file = read.csv(paste0(meta_file_dir,".csv"),check.names=F)
} else {
  mc = get(load(paste0(feat_dir,"/", feat_count,".Rdata")))
  meta_file = get(load(paste0(meta_file_dir,".Rdata")))
}
meta_cell_child_names = get(load(paste0(meta_cell_child_names_dir,".Rdata")))

## for each feature
a = foreach(feat_type=feat_types) %dopar% {
  tryCatch({
    cat("\n", feat_type, " ",sep="")
    start2 = Sys.time()
    
    ## upload and prep feature matrix
    if (readcsv) {
      m0 = as.matrix(read.csv(paste0(feat_dir,"/", feat_type,".csv"),row.names=1, check.names=F))
    } else {
      m0 = as.matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata"))))
    }
    if (!rownames(m0)[1]%in%meta_file[,id_col]) {
      cat("\nskipped: ",feat_type,", matrix rownames must match fileName column in meta_file","\n", sep="")
    }
    
    ## does feature matrix have cell populations on column names?
    layers = 0
    countThres = 0
    colhascell = ifelse(str_split(feat_type,"-")[[1]][2]=="cell",T,F)
    if (colhascell) {
      layers = c(1,2,4,max(unique(sapply(unlist(str_split(colnames(m0),"_")[[1]]), function(x) str_count(x,"[+-]")))))
      countThres = cellCountThres
    }
    
    ## for each layer, trim feature matrix accordingly
    for (k in layers) {
      #trim matrix
      mm = trimMatrix(m0,TRIM=T, mc=mc, sampleMeta=meta_file, sampleMeta_to_m1_col=id_col, target_col=target_col, control=control, order_cols=order_cols, colsplitlen=NULL, k=k, countThres=countThres, goodcount=good_count, good_sample=good_sample)
      if (is.null(mm)) next
      m_ordered = mm$m
      meta_file_ordered = mm$sm
      if (all(meta_file_ordered[,target_col]==meta_file_ordered[1,target_col])) next
      
      ## split up analysis of feature matrix rows by split_col
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
        if (length(unique(sm[,target_col]))<=1) next
        
        ## for each biclustering method
        for (bcmethod in bcmethods) { #requires binary matrix
          if (bcmethod=="BB-binary" & !grepl("pval[A-z]*TRIM",feat_type)) next
          
          # where to save biclustering
          bcname0 = paste("/",bcmethod, "_", feat_type, "_splitby-",split_col,".", tube, "_layer", str_pad(k, 2, pad = "0"), "_countThres-", countThres, sep = "")
          bcname = paste0(biclust_source_dir, bcname0)
          
          # bicluster
          if (overwrite | !file.exists(paste0(bcname,".Rdata"))) {
            if (bcmethod == "plaid") bc = biclust(as.matrix(m), method=BCPlaid(), row.release=.3,col.release=.7, back.fit=10, verbose = F)
            if (bcmethod == "CC") bc = biclust(as.matrix(m), method=BCCC(), number=Kr)
            if (bcmethod == "Xmotifs") bc = biclust(as.matrix(m), method=BCXmotifs(), number=Kr, ns=50, nd=500, alpha=10)
            if (bcmethod == "spectral") bc = biclust(as.matrix(m), method=BCSpectral(), numberOfEigenvalues=10)
            if (bcmethod == "bimax") bc = biclust(as.matrix(m), method=BCBimax(),number=Kr)
            if (bcmethod == "quest") bc = biclust(as.matrix(m), method=BCQuest(), number=Kr, ns=50)
            if (bc@Number==0) next
            
            bcb = NULL
            if (grepl("GrNMF",bcmethod) & colhascell & !grepl("_",colnames(m)[1])) {
              #build binary relation graph between features
              cellpops = colnames(m)
              medge = matrix(0,nrow=length(cellpops),ncol=length(cellpops))
              for (cellpop_ind in 1:ncol(m)) {
                cellpop = colnames(m)[cellpop_ind]
                cind = which(names(meta_cell_child_names)==cellpop)
                if (length(cind)==0) next
                children = meta_cell_child_names[[cind]]
                children_ind = colnames(m) %in% unlist(children)
                medge[cellpop_ind,children_ind] = medge[children_ind,cellpop_ind] = 1
              }
              # medge = Matrix(medge,sparse=T)
              # bicluster
              bcb = grnmf(t(as.matrix(m)), medge, k=Kr, lambda_multiple=as.numeric(str_split(bcmethod,"-")[[1]][2]), n_iter=max(ncol(m),1000), converge=1e-06, dynamic_lambda=T)
              # bcb = grnmf(t(m), medge, k=Kr, lambda_multiple=1, n_iter=1000, converge=1e-06, dynamic_lambda=T)
              if (is.null(bcb)) next
              
              bcb$rowxfactor = bcb$V
              bcb$factorxcol = t(bcb$U)
              if (all(is.na(bcb$rowxfactor))) next
            }
            
            # unused
            if (bcmethod == "fabia") {
              bcb = fabia(as.matrix(abs(m)), p=Kr,alpha=0.01,cyc=max(ncol(m),1000),spl=0,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0)
              if (is.null(bcb)) next
              
              bcb$rowxfactor = bcb@L; if(all(bcb$rowxfactor==0)) next
              bcb$factorxcol = bcb@Z
            }
            
            tryCatch({
              if (grepl("nmf",bcmethod)) {
                bcb = nmf(as.matrix(abs(m)), rank=Kr, method=str_split(bcmethod,"-")[[1]][2])#, nrun=10, method=list("lee", "brunet", "nsNMF"))
                if (is.null(bcb)) next
                
                bcb$rowxfactor = basis(bcb)
                bcb$factorxcol = coef(bcb)
              } 
            }, error = function(err) { cat(paste("nmf error:  ",err)); bcb = NULL })
            
            
            
            
            
            
            # adjust format of biclustering to match output of biclust()
            if (grepl("nmf|GrNMF|fabia",bcmethod)) {
              if (is.null(bcb)) next
              
              # threshold to determine significant bicluster
              rthres = quantile(bcb$rowxfactor,qthres)
              cthres = quantile(bcb$factorxcol,qthres)
              
              # get the framework of biclust output to put nmf output into
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
                
                bc@Number = nrow(bcb$factorxcol)
                
                bc@info = list(rowxfactor=bcb$rowxfactor,factoxcol=bcb$factorxcol)
              } 
            }
            
            
            if (bcmethod == "BB-binary") {
              # create binary matrix as input into BB-binary
              mbinary = m
              mbinary[mbinary != 0] = 1 #make matrix binary (for p values TRIM only)
              bcb = B2PS(as.matrix(mbinary), sideData=NULL, Kt=Kr, Kp=Kc, iterations=max(ncol(mbinary)/20,min_iter), alpha_p = 1, alpha_t = 1, alpha_e = 1, alpha_sd = 1)
              # THETA is trasnposed!
              bcb$theta = t(bcb$Theta[,,2])
              
              # get significant biclusters only
              theta = bcb$theta
              theta = theta[sort(unique(bcb$transcript.clusters)),sort(unique(bcb$patient.clusters))]
              bcb$p = array(1,dim = c(nrow(bcb$theta), ncol(bcb$theta)))
              for (i in unique(bcb$transcript.clusters)) {
                for (j in unique(bcb$patient.clusters)) {
                  bcb$p[i,j] = t.test.single(as.vector(theta),bcb$theta[i,j])
                  if (is.na(bcb$p[i,j])) bcb$p[i,j] = 1
                }
              }
              
              # adjust format of biclustering to match output of biclust()
              bc = biclust(array(0,dim=c(2,2)), method=BCPlaid())
              bc@info = bcb
              if (onlysigBB) {
                bc@info$cid = cid = which(bcb$p<pval_thres, arr.ind=T)
                if (nrow(cid) == 0) {
                  if (length(theta)>=6 & !all(bcb$p==1)) {
                    bc@info$cid = cid = which(bcb$p<=sort(bcb$p)[4], arr.ind=T)
                  } else {
                    # save(NULL, file=paste0(bcname,".Rdata"))
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
            
            
            
            
            
            ## get & save clustering & labels!
            
            if (bcmethod != "BB-binary") {
              rowclust = rowxcluster_to_cluster(bc@RowxNumber)
              colclust = clusterxcol_to_cluster(bc@NumberxCol)
            } 
            # rowlabel = sm[,target_col]
            
            # names(rowclust) = names(rowlabel) = rownames(m)
            names(rowclust) = rownames(m)
            names(colclust) = colnames(m)
            # f1 = f.measure.comembership(rowlabel,rowclust)
            
            # save biclustering
            bc0 = list(source=bc,rowclust=rowclust,colclust=colclust)
            # bc0 = list(source=bc,rowclust=rowclust,colclust=colclust,rowlabel=rowlabel)
            save(bc0, file=paste0(bcname,".Rdata"))
            write.csv(rowclust, file=paste0(bcname,"_rowclust.csv"))
            write.csv(colclust, file=paste0(bcname,"_colclust.csv"))
            # write.csv(rowlabel, file=paste0(bcname,"_rowlabel.csv"))
            
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
          
          # rowlabel = as.numeric(factor(sm[,target_col]))
          
          # }
          
        }
        
        
      } 
    } #layer
    TimeOutput(start2)
  }, error = function(err) { cat(paste("ERROR:  ",err)); return(T) })
  return(F)
}




TimeOutput(start)




