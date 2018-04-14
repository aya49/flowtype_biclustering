## Input: original normalized count matrix --> Output: phenodeviance features
# aya43@sfu.ca 20151228

#Directory
root = "~/projects/flowCAP-II"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/sampleMeta.Rdata", sep="")
phenoParent_dir = paste(result_dir, "/phenoParent.Rdata",sep="")
phenoParent_ind_dir = paste(result_dir, "/phenoParent_ind.Rdata",sep="")
matrix_dir = paste(result_dir, "/matrix", sep="")

#Output
phenoMeta_trim_dir = paste(result_dir, "/phenoMeta", sep="")
sampleMeta_trim_dir = paste(result_dir, "/sampleMeta", sep="")
matrixPval_dir = paste(result_dir, "/matrixPval",sep="")
matrixPvalTRIM_dir = paste(result_dir, "/matrixPvalTRIM",sep="")
matrixLogFold_dir = paste(result_dir, "/matrixLogFold",sep="")
matrixLogFoldTRIM_dir = paste(result_dir, "/matrixLogFoldTRIM",sep="")
matrixMaxCountAdj_dir = paste(result_dir, "/matrixMaxCountAdj",sep="")
matrixMaxCountAdjTRIM_dir = paste(result_dir, "/matrixMaxCountAdjTRIM",sep="")
matrixKOCountAdj_dir = paste(result_dir, "/matrixKOCountAdj",sep="")

#Libraries/Functions
library(Matrix)
library(stringr)
library(foreach)
library(doMC)
library(pracma)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")

#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)









#Options for script
#adjust = "BY" #"BH" #pvalue adjustment
#test = "wilcox" #pvalue test
cellCountThres = 200 #count threshold -- insignificant if count too low

pvalThres = .025 #p value significant if larger
sampleCountThres = 3 #only compare if >=3 samples per experiment is available; not really necessary for single file p value...
wtcount = 70 #min 70 wt used to compare with ko; not used here

splitby = "tube"
target_col = "aml" #column with control/experiment
control = "normal"

qdefault = .9 #log pvalue, anything above qp (percentile) is set to qp
alpha = .5 #threshold for what is largest change in slope; used to set threshold to prevent extreme values
logplim = 30; #largest pvalue allowable; change according to data

matrix_type = c("CountAdj")#, "Prop")
matrix_count = "Count"



#Prepare data
sampleMeta = get(load(sampleMeta_dir))
phenoMeta = get(load(phenoMeta_dir))
m0 = get(load(paste0(matrix_dir, matrix_count,".Rdata")))








start = Sys.time()

for (mcp in matrix_type) {
  cat("\n",mcp,"\t")
  start1 = Sys.time()
  
  #load matrices
  m = get(load(paste0(matrix_dir, mcp,".Rdata")))
  m[m0<=cellCountThres] = 0
  
  #define wildtypes
  g = getGTindex(sampleMeta$aml, control, sampleCountThres, sampleMeta$fileName, all=T)
  ftGT = g$attlist; ftWTIndex = g$controlIndex; ftKOGT = g$exp; ftKOIndex = g$expIndex; ftKOgoodGTIndex = g$goodexpIndex; rm(g) #Index of unique KO genotypes that has 3+ samples available ftKOIndex(ftKOGTIndexIndex)
  
  #split by tubes
  rowcombos = NULL
  for (i in length(ftKOIndex):1) {
    kotube = sampleMeta[ftKOIndex[[i]],splitby]
    wttube = sampleMeta[ftWTIndex,splitby]
    rowcombos[[i]][[1]] = ftWTIndex[wttube==kotube]
    rowcombos[[i]][[2]] = ftKOIndex[[i]] 
  }
  
  #calculate p value for every phenotype and control experiment combo
  cat(paste("getting pValues of ", ncol(m), " possible cell populations for ", length(rowcombos), " genotypes ", sep="")) #3iTCell specific
  
  loop.ind = 1:ncol(m)
  #for (k in loop.ind){
  result = foreach(k = loop.ind) %dopar% { #for each phenotype
    pvalcol = rep(0,length(rowcombos))
    logfold = rep(0,length(rowcombos))
    maxcount = rep(0,length(rowcombos))
    kocount = rep(0,length(rowcombos))
    for (j in 1:length(rowcombos)) { #for each experiment file
      #adjust values; only calculate p value if cell counts are over threshold
      c1 = as.numeric(m0[ rowcombos[[j]][[1]],k ]); c1 = c1[c1>0]
      c2 = as.numeric(m0[ rowcombos[[j]][[2]],k ]); c2 = c2[c2>0]
      c10 = 0; if (length(c1)>0) c10 = exp(mean(log(c1)))
      c20 = 0; if (length(c2)>0) c20 = exp(mean(log(c2)))
      if (c10>=cellCountThres | c20>=cellCountThres) { #set as 1 if both WT and KO medians (to avoid outler influence) < cellCountThreshold then change isn't significant
        # if (test=="wilcox") { pvalcol[j] = wilcox.test(compare1, compare2)$p.value 
        # } else if (test=="ttest") { try({ pvalcol[j] = t.test(compare1, compare2)$p.value }) } 
        compare1 = as.numeric(m[ rowcombos[[j]][[1]],k ]); compare1. = compare1[compare1>0]
        compare2 = as.numeric(m[ rowcombos[[j]][[2]],k ]); compare2. = compare2[compare2>0]
        compare10 = 0; if (length(compare1.)>0) compare10 = exp(mean(log(compare1.)))
        compare20 = 0; if (length(compare2.)>0) compare20 = exp(mean(log(compare2.)))
        
        ## pval
        t = (compare20-mean(compare1)) / (sd(compare1) )# / log(max(mean(compare1),mean(compare2))))
        p = 2*pt(-abs(t), df=length(compare1)-1)
        #pvalcol[j] = -log(t.test.single(compare1, compare2))
        pvalcol[j] = -log(p)
        
        ## logfold
        logfold[j] = log(compare20/compare10)
        if (compare10==0 & compare20>0) logfold[j] = log(compare20)
        if (compare10>0 & compare20==0) logfold[j] = log(1/compare10)
        if (compare10==0 & compare20==0) logfold[j] = 0
        
        ## counts
        maxcount[j] = max(compare20,compare10)
        kocount[j] = compare20
      }
    }
    return(list(pvalcol=pvalcol,logfold=logfold,maxcount=maxcount, kocount=kocount))
  }
  
  matrixPval = lapply(1:length(result), function(i) return(result[[i]]$pvalcol))
  matrixLogFold = lapply(1:length(result), function(i) return(result[[i]]$logfold))
  matrixMaxCountAdj = lapply(1:length(result), function(i) return(result[[i]]$maxcount))
  matrixKOCountAdj = lapply(1:length(result), function(i) return(result[[i]]$kocount))
  
  TimeOutput(start1)
  
  matrixPval = foreach(k=1:length(matrixPval),.combine="cbind") %dopar% { return(matrixPval[[k]]) }
  matrixLogFold = foreach(k=1:length(matrixLogFold),.combine="cbind") %dopar% { return(matrixLogFold[[k]]) }
  matrixMaxCountAdj = foreach(k=1:length(matrixMaxCountAdj),.combine="cbind") %dopar% { return(matrixMaxCountAdj[[k]]) }
  matrixKOCountAdj = foreach(k=1:length(matrixKOCountAdj),.combine="cbind") %dopar% { return(matrixKOCountAdj[[k]]) }
  
  colnames(matrixPval) = colnames(matrixLogFold) = colnames(matrixMaxCountAdj) = colnames(matrixKOCountAdj) = colnames(m)
  rownames(matrixPval) = rownames(matrixLogFold) = rownames(matrixMaxCountAdj) = rownames(matrixKOCountAdj) = ftKOGT
  
  #get rid of extreme values
  matrixPval[is.nan(matrixPval)] = 0
  matrixPval[matrixPval==Inf] = 10^(ceiling(log(max(matrixPval[matrixPval!=Inf]),10)))
  
  matrixPvalTRIM = matrixPval 
  matrixLogFoldTRIM = matrixLogFold 
  matrixMaxCountAdjTRIM = matrixMaxCountAdj 
  
  TimeOutput(start1)
  
  #adjust if doing Wilcox
  # matrixPvalAdj = foreach(i=1:ncol(matrixPval), .combine='cbind') %dopar% { return(p.adjust(matrixPval[,i], method=adjust)) }
  # matrixPvalFULL = matrixPvalAdj
  
  #trim extreme and insignificant pval
  trimind = abs(matrixPval) <= -log(pvalThres)
  
  matrixPvalTRIM[trimind] = 0
  d = density(matrixPvalTRIM[(!trimind)&matrixPvalTRIM<logplim])
  p = as.vector(findpeaks(d$y,npeaks=1))[2]
  qp0 = seq(from=.5,to=.95,by=.05)
  qp1 = sapply(qp0,function(x) return(quantile(matrixPvalTRIM[matrixPvalTRIM>d$x[p]],x)))
  qp2 = diff(qp1) #slope
  qp3 = diff(qp2) #acceleration
  q = qdefault
  if (sum(qp3>=alpha)) q = qp1[names(qp3)[which(qp3>=alpha)[1]+1]]
  shrinkind = matrixPvalTRIM>q
  matrixPval[shrinkind] = matrixPvalTRIM[shrinkind] = q
  matrixPval[matrixPval>0] = matrixPval[matrixPval>0]/q
  matrixPvalTRIM[!trimind] = matrixPvalTRIM[!trimind]/q
  
  matrixLogFoldTRIM[trimind] = 0
  matrixMaxCountAdjTRIM[trimind] = 0
  
  trimRowIndex <- apply(matrixPvalTRIM[,-1], 1, function(x) all(x==0))==T
  trimColIndex <- apply(matrixPvalTRIM[-1,], 2, function(x) all(x==0))==T
  
  matrixPvalTRIM = Matrix(matrixPvalTRIM[!trimRowIndex,!trimColIndex],sparse=T)
  matrixLogFoldTRIM = Matrix(matrixLogFoldTRIM[!trimRowIndex,!trimColIndex],sparse=T)
  matrixMaxCountAdjTRIM = Matrix(matrixMaxCountAdjTRIM[!trimRowIndex,!trimColIndex],sparse=T)
  
  sampleMeta_trim = sampleMeta[!trimRowIndex,]
  phenoMeta_trim = phenoMeta[!trimColIndex,]
  
  #make pval negative if decrease in cells
  neg = matrixLogFold<0
  matrixPval[neg] = -matrixPval[neg]
  matrixPvalTRIM[neg] = -matrixPvalTRIM[neg]
  
  TimeOutput(start1)
  
  #save
  save(sampleMeta_trim, file=paste0(sampleMeta_trim_dir, "_", mcp,".Rdata"))
  save(phenoMeta_trim, file=paste0(phenoMeta_trim_dir, "_", mcp,".Rdata"))
  
  save(matrixPval, file=paste0(matrixPval_dir, "_", mcp,".Rdata")); write.csv(matrixPval, file=paste0(matrixPval_dir, "_", mcp,".csv"))
  save(matrixPvalTRIM, file=paste0(matrixPvalTRIM_dir, "_", mcp,".Rdata"))
  a = matrixPvalTRIM; rownames(a) = sampleMeta[match(rownames(a),sampleMeta$fileName),target_col]
  write.table(as.matrix(a), file=paste0(matrixPvalTRIM_dir, "_", mcp,".csv"), sep=",", col.names=F)
  
  save(matrixLogFold, file=paste0(matrixLogFold_dir, "_", mcp,".Rdata")); write.csv(matrixLogFold, file=paste0(matrixLogFold_dir, "_", mcp,".csv"))
  save(matrixLogFoldTRIM, file=paste0(matrixLogFoldTRIM_dir, "_", mcp,".Rdata"))
  b = matrixLogFoldTRIM; rownames(b) = sampleMeta[match(rownames(b),sampleMeta$fileName),target_col]
  write.table(as.matrix(b), file=paste0(matrixLogFoldTRIM_dir, "_", mcp,".csv"), sep=",", col.names=F)
  
  save(matrixMaxCountAdj, file=paste0(matrixMaxCountAdj_dir, "_", mcp,".Rdata"))
  save(matrixMaxCountAdjTRIM, file=paste0(matrixMaxCountAdjTRIM_dir, "_", mcp,".Rdata"))
  save(matrixKOCountAdj, file=paste0(matrixKOCountAdj_dir,"_",mcp,".Rdata"))
}

cat("\nTime taken to calculate p values & barcode matrices is: ",TimeOutput(start), "\n", sep="") #3iTcell ~40min






