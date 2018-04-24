# PEER coundounding factor analysis
# aya43@sfu.ca 20170924

## root directory
root = "~/projects/IMPC"
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")
controlL = c("+_+|+_Y","+_+|+_Y","WildType","WildType","WildType") #control value in target_col column
ci = 1; panel = panelL; centre = centreL

result_dir = paste0("result/", panelL, "/", centreL); suppressWarnings(dir.create (result_dir))


## input directories
meta_dir = paste0(result_dir,"/meta")
meta_file_dir = paste(meta_dir, "/file", sep="")
feat_dir = paste(result_dir, "/feat", sep="")

## output directories
scale = T
peer_dir = paste(result_dir, "/PEER", sep="")
if (scale) peer_dir = paste0(peer_dir, "_scale", sep="")
dir.create(peer_dir, showWarnings=F)
resid_dir = paste(peer_dir, "/peer_residual_analysis.txt", sep="")


## libraries
library(stringr)
library(colorspace)
library(lubridate) #if there are date variables
library(peer)
library(qtl)
library(DMwR)
library(arules)
library(foreach)
library(doMC)
source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")


no_cores = detectCores()-1
registerDoMC(no_cores)







## options
options(stringsAsFactors=FALSE)
options(device="cairo")
options(na.rm=T)

writecsv = T

countThres = 2000 #columns/rows must have at least goodCount number of elements that is greater than countThres, else delete
good_count = 3
good_sample = 3 #need more than good_sample samples per class for class to be included

layers = c(5)
factors = c(15,5,10)

wtonly = "" #"_WTonly" if only analyzing wildtypes, else ""; code modified, just leave this! do both!
# layerequal = c(T,F) # only include cell populations in the layer (T) or it and all above (F)
le = T #only cell pops in one layer; F means cell pops in that layer and above
lei = "-all"; if (le) lei = "-layerbylayer"

interested_cols = c("gene", "date", "birth_date", "gender", "colony", "strain", "fur") #meta_file columns to plot against generated covriates
interested_covs_cols = list(NULL,c("gene"),c("gender"),c("gene","gender"),c("gene", "date", "gender"),c("gene", "date", "birth_date", "gender", "colony", "strain", "fur")) #meta_file columns to input into model as known covariates
interested_cont_cols = c("date","birth_date") #continuous variables
target_col = "gene" #column with control/experiment
id_col = "fileName"
split_col = NULL #NULL if split by nothing
order_cols = c("barcode","sample","specimen","date")
time_col = "date"

control = str_split(controlL,"[|]")[[1]]

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






replaceModel = F #if true, run the whole thing, else use existing predictions
# dateRestriction = c("-dateRestrict","2015-08-25","2016-03-03") #if "", don't consider, if "_dateRestrict", cut samples into date interval
dateRestriction = c("")
niter = 10000

scale_center = T
scale_scale = T

plotCountThres = 500 #how high must mean cell population count be ot get plotted for before and after plot
plotCountNo = 4 #how many of those cell population counts to plot














start = Sys.time()

#Prepare data
meta_file0 = get(load(paste0(meta_file_dir,".Rdata")))
meta_file0$gene[grepl("[+]_Y",meta_file0$gene)] = "+_+"
mc = get(load(paste0(feat_dir,"/",feat_count,".Rdata")))

#for (ci in 4:1) {
cat("\n",paste0(panel," ",centre))


for (feat_type in feat_types) { cat("\n  ", feat_type, ": ") #for each type of feature
  start2 = Sys.time()
  
  # ## get interested_cols columns
  # interested_colsCols = which(colnames(sm)%in%interested_cols)
  
  m0 = get(load(paste0(feat_dir, "/", feat_type,".Rdata")))
  meta_file1 = meta_file0[match(rownames(m0),meta_file0[,id_col]),]
  
  if (dateRestriction[1]!="") {
    dateind = ymd(meta_file1[,time_col]) > ymd(dateRestriction[2]) & ymd(meta_file1[,time_col]) < ymd(dateRestriction[3]) & m0[,5] < 250000 & m0[,5] > 75000 #set according to graphs
    m0 = m0[dateind,]
    meta_file1 = meta_file1[dateind,]
  }
  
  #if wildtype only, clip meta_file and feature matrix
  if (wtonly!="") {
    wtind = grepl(paste(control,collapse="|"),meta_file1[,target_col])
    m0 = m0[wtind,]
    meta_file1 = meta_file1[wtind,]#!colnames(meta_file)%in%target_col]
  }
  
  
  
  ## get only 1st layer phenotypes
  phenoLevel = sapply(unlist(str_split(colnames(m0),"_")), function(x) str_count(x,"[+-]"))
  
  ## peer analysis ------------------------------------------------
  
  
  
  
  ## no covariates included (meta_file not used)
  for (interested_cov_col in interested_covs_cols) { cat(", cov")
    compare_only = F; if (identical(interested_cov_col,target_col) & !wtonly=="") compare_only = T
    
    # keep only needed columns in meta_file_ordered_colonly (interested_cols)
    meta_file_colonly = meta_file1[,colnames(meta_file1)%in%interested_cols]
    interested_colnames = colnames(meta_file_colonly)
    
    # create known covariate matrix
    noCov = T
    covs = NULL
    if (!is.null(interested_cov_col) | compare_only) {
      col_ind = match(interested_cov_col,interested_colnames)
      col_ind = col_ind[!is.na(col_ind)]
      if (length(col_ind)>0) {
        covs1 = meta_file_colonly[,col_ind]
        if (is.null(dim(covs1))) { covs1 = matrix(covs1,ncol=1) }
        covs = sapply(1:ncol(covs1), function(x) as.numeric(factor(covs1[,x], ordered=T)))
        if (is.null(dim(covs))) { covs = matrix(covs,ncol=1) }
        colnames(covs) = colnames(meta_file_colonly)[col_ind] 
        noCov = F
      }
    }
    
    
    fnamepartt = ""; if (!noCov) fnamepartt = paste0("_knowncov-",paste(interested_cov_col,collapse="."))
    fnamepure = paste0(feat_type,".PEER-",wtonly,dateRestriction[1],lei)
    fname0 = paste0(peer_dir, "/", fnamepure)
    fname1 = paste0(fname0,fnamepartt)
    
    
    modelfactors0 = foreach (layer = unique(phenoLevel)) %dopar% {
      #for (layer in layers) { cat("layer ",layer,", ")
      # if (le==F & layer==1) next()
      if (layer == 0) return(NA)
      
      
      
      # trim feature matrix
      if (layer == 1) le = F
      mm = trimMatrix(m0,TRIM=T, mc=mc, sampleMeta=meta_file1, sampleMeta_to_m1_col=id_col, target_col=target_col, control=control, order_cols=order_cols, colsplitlen=NULL, k=layer, konly=le, countThres=countThres, goodcount=good_count, good_sample=good_sample)
      if (is.null(mm)) return(NA)
      if (layer == 1 & lei == "_only") le = T
      Y = mm$m
      meta_file_ordered_colonly = mm$sm[,match(colnames(meta_file_colonly),colnames(mm$sm))]
      
      
      
      fname = paste0(fname1,"_layer", str_pad(layer, 2, pad = "0"),"_countThres-",countThres,"_n.",ncol(Y),fnamepartt)
      
      
      
      # compare predictions given different number of factors to find
      colour = rainbow(length(factors)) 
      
      
      
      
      # prepare list of models ((list of the predicted factors from PEER))
      modelexists = T
      if (compare_only) {
        modelsALL = get(load(paste0(gsub("_WTonly","",fname),".Rdata")))
        modelsALL = Filter(Negate(is.null), modelsALL)
        models = get(load(paste0(fname0,".Rdata")))
        models = Filter(Negate(is.null), models)
        modelsfactors = as.numeric(intersect(names(models),names(modelsALL)))
        modelsfactors = modelsfactors[!is.na(modelsfactors)]
        if (length(modelsfactors)==0) next
      } else if (!replaceModel & file.exists(paste0(fname,".Rdata"))) {
        models = get(load(paste0(fname,".Rdata")))
        if (length(models)==0) next()
        modelsfactors = as.numeric(names(models))
        #if models was already made, just replot it, likely plotting's changed that's why i'm coming back; and it doesn't take long
      } else {
        models = list()
        modelsfactors = factors
        modelexists = F
      }
      
      # for different factors 
      for (factori in 1:length(modelsfactors)) {
        kfactor = modelsfactors[factori]
        
        
        if (compare_only) {
          AlphaALL = modelsALL[[factori]]$Alpha
          #factors:
          XALL = modelsALL[[factori]]$X
          #weights:
          WALL = modelsALL[[factori]]$W
        }
        if (modelexists) {
          Alpha = models[[factori]]$Alpha
          AlphaPlot = which(Alpha!=0)
          #factors:
          X = models[[factori]]$X
          #weights:
          W = models[[factori]]$W
        } else {
          # build model
          # set data and parameters
          model = PEER()
          if (scale) {
            PEER_setPhenoMean(model, scale(as.matrix(Y),center=scale_center,scale=scale_scale)) # data for inference - note the as.matrix() !
          } else {
            PEER_setPhenoMean(model, as.matrix(Y)) # data for inference - note the as.matrix() !
          }
          # set priors (these are the default settings of PEER)
          PEER_setPriorAlpha(model,0.1,0.1)
          PEER_setPriorEps(model,0.1,10)
          PEER_setNmax_iterations(model,niter)
          PEER_setNk(model, kfactor) #number of factor for learning
          PEER_setAdd_mean(model, TRUE) #doesn't work otherwise; adds one factor that is all 1 with weights as mean count?
          if (!noCov) PEER_setCovariates(model, as.matrix(covs)) # covariates (e.g batch, RNA quality, ...) - not the as.matrix()!
          
          # perform inference
          PEER_update(model)
          
          #investigate results
          #ARD parameters
          Alpha = PEER_getAlpha(model)
          Alpha[!is.finite(Alpha)] = 0
          #factors:
          X0 = X = PEER_getX(model)
          #delete factors and fix
          AlphaPlot = which(Alpha!=0 & sapply(1:ncol(X), function(x) length(unique(X[,x]))>1))
          if (length(AlphaPlot)==0) next() #skip if no valid factors found
          AlphaPlot = AlphaPlot[order(Alpha[AlphaPlot,],decreasing=F)]
          Alpha = Alpha[AlphaPlot,]
          if (is.null(dim(Alpha))) Alpha = matrix(Alpha,ncol=1)
          X = X[,AlphaPlot]
          if (is.null(dim(X))) X = matrix(X,ncol=1)
          #weights:
          W = PEER_getW(model)[,AlphaPlot]
          if (is.null(dim(W))) W = matrix(W,ncol=1)
          #get corrected dataset:
          Yc = PEER_getResiduals(model)
          dimnames(Yc) = dimnames(Y)
          YcUnscaled = NULL
          if (scale) {
            YcUnscaled = unscale(Yc,scale(as.matrix(Y),center=scale_center,scale=scale_scale))
            dimnames(YcUnscaled) = dimnames(as.matrix(Y))
          }
          
          models[[as.character(kfactor)]] = list(X=X, Alpha=Alpha, W=W, Yc=Yc, YcUnscaled=YcUnscaled, scale_scale=scale_scale, scale_center=scale_center)
        }
        
        ## plot inverse variance of factors - in this case, we expect a natural elbow where there are 5 active factors, as 5 were simulated
        # Alpha[!is.finite(Alpha)] = 0
        # if (any(Alpha!=0)) {
        #   tryCatch({
        #     lines(1.0 / Alpha, type="l", col=colour[factori])
        #   }, error = function(err) {
        #     plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops)", main="PEER inference", type="l",col=colour[factori])#,xlim=nrow(Alpha),ylim=c(0,max(Alpha[is.finite(Alpha)])))       
        #     legend("topright",legend=factors,text.col=colour,title="no. factors to predict")
        #   })
        # }
        
        
        
        
        
        # #compare PEER on all data vs WTonly only if all data is given the gene known covariate; don't change order of WTonlyL, this way, we'd have already done WT only for comparison
        # if (compare_only) {
        #   meta_file_ordered_colonlytemp = meta_file_ordered_colonly
        #   wtind = grepl(paste(control,collapse="|"),meta_file_ordered_colonly[,meta_file_ordered_colonly])
        #   meta_file_ordered_colonly = XALL[wtind,]
        #   colnames(meta_file_ordered_colonly) = AlphaALL
        # }
        # 
        # ## does learned factors matter
        # if (wtonly=="" & !compare_only) {
        #   png(filename=paste0(fname, "_factor-",kfactor, "_andWTOnlyFactorPlots2.png"), width=length(AlphaPlot)*500*2, height=ncol(meta_file_ordered_colonly)*300)
        #   layout(matrix(c(1,rep(2,(length(AlphaPlot)*2)-1), 3:(ncol(meta_file_ordered_colonly)*length(AlphaPlot)*2+2)),ncol=length(AlphaPlot)*2,byrow=T))
        # } else {
        #   if (compare_only) {
        #     png(filename=paste0(fname, "_factor-",kfactor, "_compareWTandALL.png"), width=length(AlphaPlot)*500, height=ncol(meta_file_ordered_colonly)*300)
        #   } else {
        #     png(filename=paste0(fname, "_factor-",kfactor, ".png"), width=length(AlphaPlot)*500, height=ncol(meta_file_ordered_colonly)*300)
        #   }
        #   layout(matrix(c(1,rep(2,length(AlphaPlot)-1), 3:(ncol(meta_file_ordered_colonly)*length(AlphaPlot)+2)),ncol=length(AlphaPlot),byrow=T))
        # }
        # 
        # par(mar=c(3,3,3,3)) #mfrow=c(length(AlphaPlot),ncol(meta_file_ordered_colonly))
        # 
        # meta_file_ordered_colonly.1 = sapply(meta_file_ordered_colonly, function(x) as.numeric(factor(x, ordered=T))) #numeric
        # 
        # 
        # ## do plotting again!
        # if (compare_only) {
        #   if (max(1/Alpha)<max(1/as.numeric(colnames(meta_file_ordered_colonly)))) {
        #     ## plot inverse variance of factors - in this case, we expect a natural elbow where there are 5 active factors, as 5 were simulated
        #     plot(1.0 /  as.numeric(colnames(meta_file_ordered_colonly)),xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops)", main="PEER inference", type="l",col=colour[factori])#,xlim=nrow(Alpha),ylim=c(0,max(Alpha[is.finite(Alpha)])))       
        #     lines(1.0 / Alpha,col="black",lty=2)
        #     legend("topright",legend=factors,text.col=colour,title="no. factors to predict")
        #     
        #     plot(1.0 /  as.numeric(colnames(meta_file_ordered_colonly)),xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops) (dashed line = aLL files)", main="PEER inference", type="l",col=colour[factori])#,xlim=nrow(Alpha),ylim=c(0,max(Alpha[is.finite(Alpha)]))) 
        #     lines(1.0 /Alpha,col="black",lty=2)
        #     legend("topright",legend=factors,text.col=colour,title="no. factors to predict")
        #   } else {
        #     ## plot inverse variance of factors - in this case, we expect a natural elbow where there are 5 active factors, as 5 were simulated
        #     plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops)", main="PEER inference", type="l",col=colour[factori])#,xlim=nrow(Alpha),ylim=c(0,max(Alpha[is.finite(Alpha)])))       
        #     lines(1.0 / as.numeric(colnames(meta_file_ordered_colonly)),col="black",lty=2)
        #     legend("topright",legend=factors,text.col=colour,title="no. factors to predict")
        #     
        #     plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops) (dashed line = aLL files)", main="PEER inference", type="l",col=colour[factori])#,xlim=nrow(Alpha),ylim=c(0,max(Alpha[is.finite(Alpha)]))) 
        #     lines(1.0 / as.numeric(colnames(meta_file_ordered_colonly)),col="black",lty=2)
        #     legend("topright",legend=factors,text.col=colour,title="no. factors to predict")
        #   }
        # } else {
        #   ## plot inverse variance of factors - in this case, we expect a natural elbow where there are 5 active factors, as 5 were simulated
        #   plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops)", main="PEER inference", type="l",col=colour[factori])#,xlim=nrow(Alpha),ylim=c(0,max(Alpha[is.finite(Alpha)])))       
        #   legend("topright",legend=factors,text.col=colour,title="no. factors to predict")
        #   
        #   plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops)", main="PEER inference", type="l",col=colour[factori])#,xlim=nrow(Alpha),ylim=c(0,max(Alpha[is.finite(Alpha)])))       
        #   legend("topright",legend=factors,text.col=colour,title="no. factors to predict")
        # }
        # 
        # 
        # ## plot factors against variables if they do
        # 
        # #heat colours for continuous values
        # colourH = heat.colors(max(sapply(1:ncol(meta_file_ordered_colonly), function(x)length(unique(meta_file_ordered_colonly[,x])))))
        # 
        # for (col in 1:ncol(meta_file_ordered_colonly)) { #for each interested_cols known factor
        #   for (ai in 1:length(AlphaPlot)) { #for each known/learned factor
        #     for (ploti in 1:2) {
        #       if (ploti==2 & !wtonly=="") break()
        #       if (ploti==2 & wtonly=="" & !compare_only) {
        #         meta_file_ordered_colonlytemp0 = meta_file_ordered_colonly
        #         Xtemp0 = X
        #         wtind = grepl(paste(control,collapse="|"),meta_file_ordered_colonlytemp0[,target_col])
        #         meta_file_ordered_colonly = meta_file_ordered_colonlytemp0[wtind,]
        #         X = Xtemp0[wtind,]
        #       }
        #       
        #       cont = colnames(meta_file_ordered_colonly)[col]%in%interested_cont_cols #is interested_cols factor continuous
        #       
        #       # interested_cols factor values to plot
        #       x_value = meta_file_ordered_colonly[,col]
        #       coloursm = x_split = as.integer(factor(x_value))
        #       coloursm[is.na(coloursm)] = "gray" #grey for non-finite values
        #       x_splitnames = sort(unique(meta_file_ordered_colonly[,col]))
        #       
        #       if (length(unique(X[,ai]))==1) { #should only happen in ploti=2 when we compare gene factor
        #         cort = "NA"
        #       } else {
        #         cort = specify_decimal(cor.test(x_split, X[,ai])$p.value,4)
        #       }
        #       cor = specify_decimal(cor(x_split, X[,ai]),4)
        #       
        #       # try({
        #       if (cont | compare_only) {
        #         
        #         coloursm = colourH[as.integer(coloursm)]
        #         if (grepl("date",colnames(meta_file_ordered_colonly)[col])) {
        #           x_value = as_date(meta_file_ordered_colonly[,col])
        #           x_splitnames = sort(unique(x_value))
        #         }
        #         if (!compare_only) mainn = paste0("Factor ",ai," and ", colnames(meta_file_ordered_colonly)[col]," Pearson Cor = ", cor,"; p value=",cort)
        #         if (compare_only) mainn = paste0("Factor ",ai," and ALL.", col," Pearson Cor = ", cor,"; p value=",cort)
        #         plot(x_value, X[,ai], col = coloursm,
        #              main = mainn,
        #              #ylim=c(min(X[is.finite(X[,ai]),ai]),max(X[X[,ai],ai]))+1,
        #              xlab = colnames(meta_file_ordered_colonly)[col], ylab = paste0("factor ",ai))
        #       } else {
        #         if (length(unique(x_split))>3 & length(unique(X[,ai]))>1) {
        #           if (length(unique(X[,ai]))>1.5*length(unique(x_split))) {
        #             Xai = as.numeric(factor(discretize(X[,ai],categories=2*length(unique(x_split)))))
        #           } else {
        #             Xai = X[,ai]
        #           }
        #           cort = specify_decimal(chisq.test(x_split, Xai)$p.value,4)
        #           cort = paste0(cort," (Chi2)")
        #         } 
        #         if (!compare_only) mainn = paste0("Factor ",ai," and ", colnames(meta_file_ordered_colonly)[col]," Pearson Cor (ok if 2 var values) = ", cor,"; p value=",cort)
        #         if (compare_only) mainn = paste0("Factor ",ai," and ALL.", col," Pearson Cor (ok if 2 var values) = ", cor,"; p value=",cort)
        #         
        #         xy_boxplot = lapply(sort(unique(x_split)), function(x) X[x_split==x,ai])
        #         boxplot(xy_boxplot, lwd = 1, outline=F, #ylim=c(min(X[is.finite(X[,ai]),ai]),max(X[X[,ai],ai])),
        #                 main = mainn,
        #                 xaxt = 'n', xlab = colnames(meta_file_ordered_colonly)[col], ylab = paste0("factor",ai), ylab = "Covariate value") #,xaxt ='n', xlab=testcol
        #         axis(1, at=1:length(x_splitnames), labels=x_splitnames)
        #         points(jitter(x_split, factor=1), X[,ai], col = coloursm, cex=.5, pch=16)
        #       }
        #       # })
        #       
        #       if (ploti==2 & wtonly=="" & !compare_only) {
        #         meta_file_ordered_colonly = meta_file_ordered_colonlytemp0
        #         X = Xtemp0
        #       }
        #       
        #     }
        #     
        #   }
        #   
        # }
        # if (compare_only) {
        #   meta_file_ordered_colonly = meta_file_ordered_colonlytemp
        # }
        # graphics.off()
        # 
        
        
        
        
        
        
      }
      # if (!modelexists & length(models)>0 & !compare_only) {
      #   # save(models,file=paste0(fname,".Rdata"))
      # } else {
      #   next()
      # }
      
      
      
      # ## compare old data with new fixed data
      # if (wtonly=="") {
      #   png(paste0(fname,"_beforeAfter.png"), width=length(models)*700*2, height=ncol(meta_file_ordered_colonly)*300*2+1)
      #   par(mfcol=c(ncol(meta_file_ordered_colonly)*2+1,length(models)*2))
      # } else {
      #   png(paste0(fname,"_beforeAfter.png"), width=length(models)*700, height=ncol(meta_file_ordered_colonly)*300*2+1)
      #   par(mfcol=c(ncol(meta_file_ordered_colonly)*2+1,length(models)))
      # }
      # par(mar=c(10,3,3,3))
      # for (factori in 1:length(models)) {
      #   kfactor = names(models)[factori]
      #   Alpha = models[[factori]]$Alpha
      #   W = models[[factori]]$W
      #   YcUnscaled = models[[factori]]$YcUnscaled
      #   cellpops0 = order(colMeans(as.matrix(Y))) #ranking of smallest count to largest 1,2,3...
      #   cellpops0[which(colMeans(as.matrix(Y))>plotCountThres)] = -Inf
      #   if (sum(cellpops0!=-Inf)==0) next()
      #   cellpops1 = order(rowSums(abs(W)))
      #   cellpops = order(cellpops0+cellpops1, decreasing=T)[1:min(sum(cellpops0!=-Inf),plotCountNo)]
      #   colours = rainbow(length(cellpops))
      #   
      #   for (ploti in 1:2) {
      #     if (ploti==2) {
      #       if (wtonly!="") next()
      #       wtind = grepl(paste(control,collapse="|"),meta_file_ordered[,target_col])
      #       
      #       meta_file_ordered_colonlytemp = meta_file_ordered_colonly
      #       meta_file_ordered_colonly = meta_file_ordered_colonlytemp[wtind,]
      #       
      #       Ytemp = Y
      #       Y = as.matrix(Ytemp)[wtind,]
      #       colnames(Y) = colnames(as.matrix(Ytemp))
      #       
      #       YcUnscaledtemp = YcUnscaled
      #       YcUnscaled = YcUnscaledtemp[wtind,]
      #     }
      #     
      #     ## plot inverse variance of factors - in this case, we expect a natural elbow where there are 5 active factors, as 5 were simulated
      #     if (ploti==2) {
      #       plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops)", main=paste0("WT only\nPEER inference; factors=", kfactor), type="l")
      #     } else {
      #       plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance (factor inverse variance over cell pops)", main=paste0("PEER inference; factors=", kfactor), type="l")
      #     }
      #     for (coli in 1:ncol(meta_file_ordered_colonly)) {
      #       #BEFORe
      #       plot(jitter(as.numeric(factor(meta_file_ordered_colonly[,coli]))), as.matrix(Y)[,cellpops[1]], cex=.5, col=colours[1], xlab='n', ylab="Normalized Cell Count", main=paste0("[",colnames(meta_file_ordered_colonly)[coli],"] BEFORE \n(most effected cell pops with countAdj>",plotCountThres,")"), xaxt = 'n', ylim=c(min(as.matrix(Y)[,cellpops]),max(as.matrix(Y)[,cellpops])+(.2*max(as.matrix(Y)[,cellpops])-(min(as.matrix(Y)[,cellpops])))))
      #       par(las=2)
      #       axis(1, at=1:length(unique(meta_file_ordered_colonly[,coli])), labels=sort(unique(meta_file_ordered_colonly[,coli])))
      #       pvar = c(var(as.matrix(Y)[,cellpops[1]]))
      #       if (length(cellpops)>1) {
      #         for (cellpopi in 1:length(cellpops)) {
      #           points(jitter(as.numeric(factor(meta_file_ordered_colonly[,coli]))), as.matrix(Y)[,cellpops[cellpopi]], cex=.5, col=colours[cellpopi])
      #           pvar = append(pvar,var(as.matrix(Y)[,cellpops[cellpopi]]))
      #         }
      #       }
      #       legend("topright",legend=paste0(colnames(as.matrix(Y))[cellpops]," var=",specify_decimal(pvar,3)), fill=colours)
      #       
      #       #AFTER
      #       plot(jitter(as.numeric(factor(meta_file_ordered_colonly[,coli]))), YcUnscaled[,cellpops[1]], cex=.5, col=colours[1], xlab='n', ylab="Normalized Cell Count", main=paste0("[",colnames(meta_file_ordered_colonly)[coli],"] AFTER \n(most effected cell pops with countAdj>",plotCountThres,")"), xaxt = 'n', ylim=c(min(YcUnscaled[,cellpops]),max(YcUnscaled[,cellpops])+(.2*max(YcUnscaled[,cellpops])-(min(YcUnscaled[,cellpops])))))
      #       par(las=2)
      #       axis(1, at=1:length(unique(meta_file_ordered_colonly[,coli])), labels=sort(unique(meta_file_ordered_colonly[,coli])))
      #       pvar = c(var(YcUnscaled[,cellpops[1]]))
      #       if (length(cellpops)>1) {
      #         for (cellpopi in 1:length(cellpops)) {
      #           points(jitter(as.numeric(factor(meta_file_ordered_colonly[,coli]))), YcUnscaled[,cellpops[cellpopi]], cex=.5, col=colours[cellpopi])
      #           pvar = append(pvar,var(YcUnscaled[,cellpops[cellpopi]]))
      #         }
      #       }
      #       legend("topright",legend=paste0(colnames(as.matrix(Y))[cellpops]," var=",specify_decimal(pvar,3)), fill=colours)
      #       
      #     }
      #     
      #     if (ploti==2) {
      #       meta_file_ordered_colonly = meta_file_ordered_colonlytemp
      #       Y = Ytemp
      #       YcUnscaled = YcUnscaled
      #     }
      #     
      #   }
      #   
      # }
      # graphics.off()
      
      return(models)
      
    }
    modelfactors = modelfactors0
    names(modelfactors) = unique(phenoLevel)
    save(modelfactors, file=paste0(fname1,".Rdata"))
    
    factors1 = factors
    keepind = rep(T,length(modelfactors))
    for (modelsn in 1:length(modelfactors)) {
      if (is.null(modelfactors[[modelsn]]) | is.na(modelfactors[[modelsn]])) {
        keepind[modelsn] = F
        next()
      }
      factors1 = intersect(factors1,names(modelfactors[[modelsn]]))
    }
    modelfactors = modelfactors[keepind]
    
    Yc0 = list()
    for (factor in factors1) {
      Yc0i = lapply(names(modelfactors), function(xi) {
        x = modelfactors[[xi]]
        if (is.null(x[[factor]]$YcUnscaled)) y = x[[factor]]$Yc
        y = x[[factor]]$YcUnscaled
        if (!le & as.numeric(xi)>1) {
          phenolevel = sapply(unlist(str_split(colnames(y),"_")), function(z) str_count(z,"[+-]"))
          return(y[,phenolevel==as.numeric(xi)])
        }
        return(y)
      })
      Yc0im = Reduce('cbind',Yc0i)
      if(!le) {
        Yc0imYc = modelfactors[[length(modelfactors)]][[factor]]$YcYcUnscaled
        if (is.null(Yc0imYc)) Yc0imYc = modelfactors[[length(modelfactors)]][[factor]]$Yc
        save(Yc0imYc, file=paste0(feat_dir,"/",fnamepure,".Rdata"))
        if (writecsv) write.csv(Yc0imYc, file=paste0(feat_dir,"/",fnamepure,"-all.csv"))
      }
      save(Yc0im, file=paste0(feat_dir,"/",fnamepure,".Rdata"))
      if (writecsv) write.csv(Yc0im, file=paste0(feat_dir,"/",fnamepure,"-layerbylayer.csv"))
    }
    
  }
  cat("\n centre ", centre, " ",TimeOutput(start2)," \n",sep="")
}

TimeOutput(start)




# #try a whole bunch of residual parameters
# sink(file=resid_dir,append=T)
# cat("\n\n\n======================================================\nLayer",layer,": plot factor weight variances for a large set of maximum number K of factors, and see if K has an effect\n")
# # plot factor weight variances for a large set of maximum number K of factors, and see if K has an effect
# for(pa in c(0.0001, 0.1, 1000)){
#   for(pb in c(0.1,10,1000)){
#     model = get_simple_PEER_object() # see simple_unsupervised_demo for how it is constructed
#     PEER_setPriorEps(model,0.1, pb);
#     PEER_update(model)
#     cat(paste("\nEps pa=", pa, "pb=", pb, "mean(residuals^2)=",mean((PEER_getResiduals(model))**2)))
#   }
# }
# sink()











# ## See dates without WT
# #Harwell
# control <- "FR-FCM-ZYCB-WildType_01_All"
# KOdays <- table(meta_file[!meta_file$gene%in%control,]$date)
# colnames(sampleCountThres) <- table(meta_file[meta_file$gene%in%control,]$date)
# noWTdays <- KOdays [!(names(KOdays) %in% names(WTdays))]
# 
# 
# 
# noWTdaysInd <- which(as.character(meta_file_ordered$date)%in%names(noWTdays))
# plot()








## AMI evaluation sci.kit.learn ==FAST on python, not good on R NMI

## library (fpc)

#set.seed(4634)
#face <- rFace(300,dMoNo=2,dNoEy=0)
#grface <- as.integer(attr(face,"grouping"))
#plotcluster(face,grface==1)
#plotcluster(face,grface, clnum=1, method="vbc")








