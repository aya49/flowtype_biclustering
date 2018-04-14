# Plot: plot samples & divide dates up for p value calculation
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
stat_dir = paste(result_dir, "/stats", sep=""); dir.create(stat_dir, showWarnings=F)
cp_dir = paste(stat_dir, "/changepoint", sep=""); dir.create(cp_dir, showWarnings=F)
single_dir = paste(stat_dir, "/singlephen", sep=""); dir.create(single_dir, showWarnings=F)


## libraries
library(stringr)
library(colorspace)
library(changepoint) # library(proxy)
library(FKF)
library(fastcluster)
library(dendextend)
library(circlize)
library(Rtsne)
library(MASS)
library(RDRToolbox)
library(scater, quietly = TRUE)
library(foreach)
library(doMC)
library(lubridate) #if there are date variables
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")



## cores
no_cores = 8#detectCores() - 1
registerDoMC(no_cores)





## options
options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

countThres = 2000 #columns/rows must have at least goodCount number of elements that is greater than countThres, else delete
good_count = 3
good_sample = 3 #need more than good_sample samples per class for class to be included

interested_cols = c("gene", "gender", "date", "colony", "strain", "fur", "birth_date", "specimen", "sample") #sampleMeta columns to plot
interested_cont_cols = c("date","birth_date") #continuous variables
target_col = "gene" #column with control/experiment
id_col = "fileName"
split_col = NULL #NULL if split by nothing
time_col = "date"
order_cols = c("barcode","sample","specimen","date")

control = str_split(controlL[ci],"[|]")[[1]]

dofullPCA = F #do PCA for all cell popoulations not just first level
plotpc = 5 #number of pca pc to plot

doISO = F #do ISO feature reduction
iso_k = 3 # number of neighbours
mds_type = c("iso", "mds")

doTsne = T #do Tsne feature reduction
theta=.5 #parameter for Tsne

doHC = F #do hierarchical clustering
link = c("ward.D", "ward.D2", "mcquitty") #, "single", "complete") # , "median", "average", "centroid"

methods = c("BinSeg","AMOC","PELT") #changepoint analysis; AMOC at most one change, PELT is poisson (quick exact), BinSeg (quick, approx), SegNeigh (slow, exact)
usemethod = "AMOC"

feat_types = c("file-cell-countAdj","file-cell-prop")
# feat_types = c("CountAdj", "Child_entropyTRIM_CountAdjBH"              ,"Child_entropyTRIM_CountAdjbonferroni"      ,"Child_entropyTRIM_CountAdjBY"              ,"Child_entropyTRIM_CountAdjPEERBH"         
#                ,"Child_entropyTRIM_CountAdjPEERbonferroni"  ,"Child_entropyTRIM_CountAdjPEERBY"          ,"Child_entropyTRIM_CountAdjPEER"            ,"Child_entropyTRIM_CountAdj"               
#                ,"CountAdjPEER"                                                                ,"LogFold_CountAdjPEER"                      ,"LogFold_CountAdj"                         
#                ,"LogFoldTRIM_CountAdjBH"                    ,"LogFoldTRIM_CountAdjbonferroni"            ,"LogFoldTRIM_CountAdjBY"                    ,"LogFoldTRIM_CountAdjPEERBH"               
#                ,"LogFoldTRIM_CountAdjPEERbonferroni"        ,"LogFoldTRIM_CountAdjPEERBY"                ,"LogFoldTRIM_CountAdjPEER"                  ,"LogFoldTRIM_CountAdj"                     
#                ,"Parent_entropyTRIM_CountAdjBH"             ,"Parent_entropyTRIM_CountAdjbonferroni"     ,"Parent_entropyTRIM_CountAdjBY"             ,"Parent_entropyTRIM_CountAdjPEERBH"        
#                ,"Parent_entropyTRIM_CountAdjPEERbonferroni" ,"Parent_entropyTRIM_CountAdjPEERBY"         ,"Parent_entropyTRIM_CountAdjPEER"           ,"Parent_entropyTRIM_CountAdj"              
#                ,"Pval_CountAdjBH"                           ,"Pval_CountAdjbonferroni"                   ,"Pval_CountAdjBY"                           ,"Pval_CountAdjPEERBH"                      
#                ,"Pval_CountAdjPEERbonferroni"               ,"Pval_CountAdjPEERBY"                       ,"Pval_CountAdjPEER"                         ,"Pval_CountAdj"                            
#                ,"PvalTRIM_CountAdjBH"                       ,"PvalTRIM_CountAdjbonferroni"               ,"PvalTRIM_CountAdjBY"                       ,"PvalTRIM_CountAdjPEERBH"                  
#                ,"PvalTRIM_CountAdjPEERbonferroni"           ,"PvalTRIM_CountAdjPEERBY"                   ,"PvalTRIM_CountAdjPEER"                     ,"PvalTRIM_CountAdj"       )
feat_count = c("file-cell-countAdj")












start = Sys.time()

meta_file = get(load(paste0(meta_file_dir,".Rdata")))
mc = get(load(paste0(feat_dir,"/",feat_count,".Rdata")))

#order samples by date, exclude genotypes with less than 3 samples

a = foreach (feat_type=feat_types) %dopar% { cat("\n  ", feat_type, ": ")
  start2 = Sys.time()
  
  m0 = as.matrix(get(load(paste0(feat_dir,"/", feat_type,".Rdata"))))
  layers = c(1,2,4,max(unique(sapply(unlist(str_split(colnames(m0),"_")), function(x) str_count(x,"[+-]")))))
  
  for (layer in layers) {
    #trim matrix
    mm = trimMatrix(m0,TRIM=T, mc=mc, sampleMeta=meta_file, sampleMeta_to_m1_col=id_col, target_col=target_col, control=control, order_cols=order_cols, colsplitlen=NULL, k=layer, countThres=countThres, goodcount=good_count, good_sample=good_sample)
    if (is.null(mm)) next
    m_ordered = mm$m
    meta_file_ordered = mm$sm
    
    ## get interested columns
    interested_col_ind = which(colnames(meta_file_ordered)%in%interested_cols)
    
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
      meta_file_ordered_split = meta_file_ordered[split_ind[[tube]],]
      meta_file_ordered_split_factor = as.data.frame(sapply(meta_file_ordered_split, function(x) as.numeric(factor(x, ordered=T))))
      
      
      ## get interested columns
      uniquecols = apply(meta_file_ordered_split_factor, 2, function(x) nrow(meta_file_ordered_split_factor)==length(unique(x)))
      meta_file_ordered_split_factor_attonly = meta_file_ordered_split_factor[,!uniquecols]
      interested_cols_ = match(interested_cols,colnames(meta_file_ordered_split_factor_attonly))
      interested_cols_ = interested_cols_[!is.na(interested_cols_)]
      
      
      
      g = getGTindex(meta_file_ordered_split[,target_col], control, good_sample)
      ko = g$expIndex
      wt = g$controlIndex #wildtype
      
      #prepare plot colours by date
      ts = meta_file_ordered_split_factor[,time_col]
      tswt = ts[wt]
      tsc = c(1,1+which(diff(ts)!=0))
      tscwt = c(1,1+which(diff(tswt)!=0))
      tcolour = heat.colors(length(unique(ts))+25)[1:length(unique(ts))]
      
      
      ## time plots -- do only for cell populations in the first layer -----------------------------------------
      if (layer==1) {
        sm = meta_file_ordered_split_factor_attonly
        cols = seq(from=1,to=ncol(m),by=2)
        
        
        ## create kalman filtering line plots for all files, for use in plots

        cat("Kalman filtering for all; ")
        kmff = kmf(m, cols)
        fkffitall = kmff$fkffitall
        statsfitall = kmff$statsfitall
        
        
        ## plot all single phenotypes onto an image ------------------------------------------------
        cat("single phen;")
        pngname <- paste0(stat_dir, "/time-1layer-all_", feat_type, "_split-",tube,".png")
        png(filename=pngname, width=3*1000, height=length(cols)*400)
        par(mfrow=c(length(cols),3), mar=c(5,5,5,5), cex.axis=1.5)
        for (i in cols) {
          y <- as.numeric(m[,i])
          ylim <- c(min(m),max(m))
          
          mvalueall <- cpt.mean(fkffitall[[i]][wt],method=usemethod, Q=20, penalty="MBIC", minseglen=5)
          plot(mvalueall, main=paste0(colnames(m)[i], "; kalman filtered on WT only; heat colours = days since ", min(sm[,time_col]), "; vertical lines = first sample on day"), cex.axis=2)
          abline(v=tsc, col="#DCDCDC")
          points(y, col=tcolour[ts], cex=.4)
          points(wt ,y[wt], cex=.4, col="black")#pch=19, 
          lines(wt, statsfitall[[i]][wt], col = "green")
          lines(wt, fkffitall[[i]][wt], col = "blue")
          legend("top", c("Actual datapoints (blue=WT)", "Local level (StructTS)", "Local level (fkf)"), col = c("red", "green", "blue"), lty = 1)
          legend.col(col=tcolour, lev=ts)
          
          for (j in 1:2) {
            if (j==1) { #everything to same scale
              plot(y, main=paste0(colnames(m)[i], "; heat colours = days since ", min(sm[,time_col]), "; vertical lines = first sample on day"), xlab="date", ylab=feat_type, col=tcolour[ts], pch=19,cex=1.5, ylim=ylim)
            } else {
              plot(y, main=paste0(colnames(m)[i], "; heat colours = days since ", min(sm[,time_col]), "; vertical lines = first sample on day"), xlab="date", ylab=feat_type, col=tcolour[ts], pch=19,cex=1.5)
            }
            abline(v=tsc, col="#DCDCDC")
            points(wt ,y[wt], col="black")
            lines(wt, statsfitall[[i]][wt], col = "green")
            lines(wt, fkffitall[[i]][wt], col = "blue")
            legend("top", c("Actual datapoints (clack=WT)", "Local level (StructTS)", "Local level (fkf)"), col = c("red", "green", "blue"), lty = 1)
            legend.col(col=tcolour, lev=ts)
          }
        }
        graphics.off()
        
        
        
        ## plot one single phenotypes and its changepoints ------------------------------------------------
        cat(" changepoint; ")
        for (i in cols) {
          y <- as.numeric(m[,i])
          ywt <- y[wt]
          
          
          pngname <- paste0(cp_dir[ci], "/time-1layer-",colnames(m)[i],"_", feat_type, "_split-",tube, ".png")
          png(filename=pngname, width=length(methods)*800, height=(1+3)*400)
          layout(matrix(c(1,1,1, 2:(3*length(methods)+1)),ncol=length(methods),byrow=T))
          par(mar=rep(5,4))
          
          mvalueall <- cpt.mean(fkffitall[[i]][wt],method=usemethod, Q=20, penalty="MBIC", minseglen=5)
          plot(mvalueall, main=paste0(colnames(m)[i], "; kalman filtered on WT only; heat colours = days since ", min(sm[,time_col]), "; vertical lines = first sample on day"), cex.axis=2)
          abline(v=tsc, col="#DCDCDC")
          points(y, col=tcolour[ts], cex=.4)
          points(wt ,ywt, cex=.4, col="black")#pch=19, 
          lines(wt, statsfitall[[i]][wt], col = "green")
          lines(wt, fkffitall[[i]][wt], col = "blue")
          legend("top", c("Actual datapoints (blue=WT)", "Local level (StructTS)", "Local level (fkf)"), col = c("red", "green", "blue"), lty = 1)
          legend.col(col=tcolour, lev=ts)
          
          ylim <- c(min(ywt), max(ywt))
          for (j in 1:length(methods)) {
            mvalue <- cpt.mean(fkffitall[[i]][wt],method=methods[j], Q=20, penalty="MBIC", minseglen=5)
            if (methods[j]==usemethod & i==1 & feat_type==feat_count) { mvaluewt = mvalue }
            plot(mvalue, main=paste("WT only; mean change: ",methods[j], "; Penalty MBIC; ",colnames(m)[i], sep=""), ylim=ylim, cex.axis=2)
            points(ywt, col=tcolour[tswt], pch=19,cex=.4)
            lines(statsfitall[[i]][wt], col = "green")
            lines(fkffitall[[i]][wt], col = "blue")
            legend.col(col=tcolour, lev=ts)
          }
          for (j in 1:length(methods)) {
            vvalue <- cpt.var(diff(fkffitall[[i]][wt]), method=methods[j], penalty="MBIC")
            plot(vvalue, main=paste("variance change: ",methods[j],sep=""), cex.axis=2)
            #points(wt ,ywt, col=tcolour[tswt], pch=19,cex=.4)
            legend.col(col=tcolour, lev=ts)
          }
          for (j in 1:length(methods)) {
            mvvalue <- cpt.meanvar(diff(fkffitall[[i]][wt]), method=methods[j], penalty="MBIC")
            plot(mvvalue, main=paste("variance/mean change: ",methods[j],sep=""), cex.axis=2)
            #points(wt ,ywt, col=tcolour[tswt], pch=19,cex=.4)
            legend.col(col=tcolour, lev=ts)
          }
          graphics.off()
        }
        
        

      }
      
      
      
      
      
      
      ## pca analysis ------------------------------------------------
      if (doISO) { cat("iso; ")
        fit <- Isomap(m,k=iso_k)
        save(fit)
      }
      cat("pca; ")
      pc <- prcomp(m)
      
      #pca scatterplot
      pngname = paste0(stat_dir, "/pca-iso_", feat_type, "_split-",tube,"_layer-",str_pad(layer, 2, pad = "0"),".png")
      if (length(split_ind)>1) pngname = gsub(".png",paste0("_splitby-",split_col,"-",tube, ".png"),pngname)
      png(filename=pngname, width=length(interested_col_ind)*400, height=(1+doISO+plotpc)*400)
      layout(matrix(c(rep(1,length(interested_col_ind)), 2:(((2*plotpc)+doISO)*length(interested_col_ind)+1)),ncol=length(interested_col_ind),byrow=T))
      par(cex=1)
      
      plot(pc$sdev^2/sum(pc$sdev^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", ylab = "variance explained")
      for (i in 1:plotpc) {
        for (col in interested_col_ind) {
          colname = colnames(meta_file_ordered_split)[col]
          coloursm <- meta_file_ordered_split_factor[,col]
          plot(pc$x[,i+1], pc$x[,i], col = coloursm, main = paste0("PCA ", colname), xlab = paste0("PC_",i+1), ylab = paste0("PC_",i))
          points(0, 0, pch = 3, cex = 4, lwd = 4)
        }
        for (col in interested_col_ind) {
          colname = colnames(meta_file_ordered_split)[col]
          coloursm <- meta_file_ordered_split_factor[,col]
          
          attribute = meta_file_ordered_split_factor[,col]
          attributen = meta_file_ordered_split[,col]
          attributenames = sort(unique(attributen))
          
          cor = cor(attribute, pc$x[,i])
          
          if (colname%in%interested_cont_cols) {
            plot(attribute, pc$x[,i], col=coloursm, main=paste0("PCA ", colname," Pearson Cor = ", cor), xlab = colname, ylab = paste0("PC_",i))
          } else {
            xy_boxplot = lapply(attributenames, function(x) pc$x[attribute==x,i])
            boxplot(xy_boxplot, lwd = 1, outline=F, ylim=c(min(pc$x[,i]),max(pc$x[,i])),
                    main = paste0("PCA ", colname," Pearson Cor (ok if 2 var values) = ", cor),
                    xaxt = 'n', xlab=colname, ylab = paste0("PC_",i)) #,xaxt ='n', xlab=testcol
            axis(1, at=1:length(attributenames), labels=attributenames)
            points(jitter(attribute, factor=1), pc$x[,i], col = coloursm)
          }
        }
      }
      if (doISO) {
        for (col in interested_col_ind) {
          coloursm <- meta_file[,col]
          plot(fit$dim2, col = coloursm, main = paste0("ISO ", colname))
          points(0, 0, pch = 3, cex = 4, lwd = 4)
        }
      }
      
      graphics.off()
    }
    
    
    
    #same as SVD on centred data
    # cx <- sweep(cbind(sm,m[,]), 2, colMeans(x), "-")
    # sv <- svd(cx)
    
    
  }
  
  cat("\n centre ", centre, " ",TimeOutput(start2)," \n",sep="")
}

TimeOutput(start)





## See dates without WT
#Harwell
# control <- "FR-FCM-ZYCB-WildType_01_All"
# KOdays <- table(sampleMeta[!sampleMeta$gene%in%control,]$date)
# colnames(sampleCountThres) <- table(sampleMeta[sampleMeta$gene%in%control,]$date)
# noWTdays <- KOdays [!(names(KOdays) %in% names(WTdays))]
# 
# 
# 
# noWTdaysInd <- which(as.character(sampleMeta2$date)%in%names(noWTdays))
# plot()








## AMI evaluation sci.kit.learn ==FAST on python, not good on R NMI

## library (fpc)

#set.seed(4634)
#face <- rFace(300,dMoNo=2,dNoEy=0)
#grface <- as.integer(attr(face,"grouping"))
#plotcluster(face,grface==1)
#plotcluster(face,grface, clnum=1, method="vbc")








