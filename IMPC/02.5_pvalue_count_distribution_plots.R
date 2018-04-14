## Input: matrixPval matrices --> Output: plots of distribution of p values & significant cell population counts
# aya43@sfu.ca 20171201

#Directory
root = "~/projects/IMPC"
result_dir = "result"; suppressWarnings(dir.create (result_dir))
setwd(root)

panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")

#Options
options(stringsAsFactors=FALSE)
# options(device="cairo")
options(na.rm=T)

#Input
phenoMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/phenoMeta.Rdata", sep="")
sampleMeta_dir = paste(result_dir, "/", panelL, "/", centreL, "/sampleMeta.Rdata", sep="")
matrix_dir = paste(result_dir, "/", panelL, "/", centreL, "/matrix", sep="")

#Output
plot_dir = paste(result_dir, "/", panelL, "/", centreL, "/plots", sep=""); for(i in 1:length(plot_dir)) { suppressWarnings(dir.create(plot_dir[i])) }
plot_html_dir = paste0(plot_dir,"/pval_count_distributions"); for(i in 1:length(plot_html_dir)) { suppressWarnings(dir.create(plot_html_dir[i])) }


#Libraries/Functions
library(Matrix)
library(stringr)
# library(foreach)
# library(doMC)
# library(pracma)
source("~/projects/IMPC/code/_funcAlice.R")
# source("~/projects/IMPC/code/_funcdist.R")

#Setup Cores
# no_cores = detectCores()-1
# registerDoMC(no_cores)









#Options for script

# matrix_type = c("CountAdjPEER","CountAdj")
adjustL = c("none","BH","BY","bonferroni") #"BH" #pvalue adjustment; "none" if don't adjust
target_col = "gene"
controlL = c("[+]_[+]|[+]_Y","[+]_[+]|[+]_Y","WildType","WildType","WildType") #control value in target_col column

matrix_type = c("CountAdj","CountAdjPEER")#"CountAdj", "Prop") #must be same as that of 02_plot.R

sampleno = 600 # sample sampleno of cells in matrix to plot



start = Sys.time()

for (ci in 1:length(paste0(panelL,centreL))) {
  centre = paste0(panelL," ",centreL)[ci]
  sampleMeta = get(load(sampleMeta_dir[ci]))
  phenoMeta = get(load(phenoMeta_dir[ci]))
  
  mm = mpt = list()
  for (mcp in matrix_type){
    cat("\n", mcp, sep="")
    start1 = Sys.time()
    
    mm[[mcp]] = mpt[[mcp]] = list()
    mm[[mcp]][["original"]] = m = get(load(paste0(matrix_dir[ci],mcp,".Rdata")))
    for (adjust in adjustL) {
      cat("\n- ",adjust,": ",sep="")
      
      mm[[mcp]][[adjust]] = mpt[[mcp]][[adjust]] = list()
      adjustname = adjust; if (adjust=="none") adjustname = ""
      
      # get trimmed p value matricies (all, WT, KO)
      mpt[[mcp]][[adjust]][["all"]] = mpvalt = get(load(paste0(matrix_dir[ci],"PvalTRIM_",mcp,adjustname,".Rdata")))
      wtind = grepl(controlL[ci], sampleMeta[match(rownames(mpvalt),sampleMeta$fileName), target_col])
      mpvalt1 = mpvalt[wtind,]
      mpvalt2 = mpvalt[!wtind,]
      mpt[[mcp]][[adjust]][["WT"]] = mpvalt1[,apply(mpvalt1,2,function(x) !all(x==0))]
      mpt[[mcp]][[adjust]][["KO"]] = mpvalt2[,apply(mpvalt2,2,function(x) !all(x==0))]
      
      cat("all (",nrow(mpt[[mcp]][[adjust]][["all"]])," x ",ncol(mpt[[mcp]][[adjust]][["all"]]),"; ",sum(mpt[[mcp]][[adjust]][["all"]]!=0),
          "), WT (",nrow(mpt[[mcp]][[adjust]][["WT"]])," x ",ncol(mpt[[mcp]][[adjust]][["WT"]]),"; ",sum(mpt[[mcp]][[adjust]][["WT"]]!=0),
          "), KO (",nrow(mpt[[mcp]][[adjust]][["KO"]])," x ",ncol(mpt[[mcp]][[adjust]][["KO"]]),"; ",sum(mpt[[mcp]][[adjust]][["KO"]]!=0),")",sep="")
      
      # get trimmed count matrices (all, WT, KO)
      for (mname in names(mpt[[mcp]][[adjust]]))
        mm[[mcp]][[adjust]][[mname]] = mm[[mcp]][["original"]][,match(colnames(mpt[[mcp]][[adjust]][[mname]]), colnames(mm[[mcp]][["original"]]))]
      
      # get html plots of count distribution of significant cell populations
      for (mname in append("original",names(mpt[[mcp]][[adjust]])) ) {
        if (mname=="original") {
          mi = mi0 = mm[[mcp]][[mname]]
        }  else {
          mi = mi0 = mm[[mcp]][[adjust]][[mname]]
        }
        pm = phenoMeta[match(colnames(mi0),phenoMeta$phenotype),]
        mi[mi<1] = 1
        mi = log(mi)
        val = unlist(lapply(unique(pm$phenolevel), function(x) sample(mi[,pm$phenolevel==x],min(sampleno,sum(pm$phenolevel==x)))))
        group = unlist(lapply(unique(pm$phenolevel), function(x) rep(paste0("layer ",x,": ",sum(mi0[,pm$phenolevel==x]<1),"/",sum(pm$phenolevel==x)*nrow(mi0)," m cell pops x n files = 0"),min(sampleno,sum(pm$phenolevel==x)))))
        dens_plot(val=val,group=as.character(group),filename=paste0(root, "/", plot_html_dir,"/counts_sample-",sampleno,"_",mcp,"_",mname,adjustname,".html"),title=paste0("Distribution of ln(normalized cell count) \nin different layers (<",sampleno," sampled) \n(if not original, matrix is trimmed)"))
      } # all, WT, KO
    } # adjust
  } # matrix_type
  TimeOutput(start1)
}
TimeOutput(start)


# mpvalt1 = mpvalt[grepl(controlL[ci], sampleMeta[match(rownames(mpvalt),sampleMeta$fileName), target_col]),]
# mppvalt1= mppvalt[grepl(controlL[ci], sampleMeta[match(rownames(mppvalt),sampleMeta$fileName), target_col]),]
# mpvalt2 = mpvalt[!grepl(controlL[ci], sampleMeta[match(rownames(mpvalt),sampleMeta$fileName), target_col]),]
# mppvalt2 = mppvalt[!grepl(controlL[ci], sampleMeta[match(rownames(mppvalt),sampleMeta$fileName), target_col]),]

# mvalt$mm1 = mvalt$mm[grepl(controlL[ci], sampleMeta[match(rownames(mvalt$mm),sampleMeta$fileName), target_col]),]
# mvalt$mpp1 = mvalt$mpp[grepl(controlL[ci], sampleMeta[match(rownames(mvalt$mpp),sampleMeta$fileName), target_col]),]
# mvalt$mm2 = mvalt$mm[!grepl(controlL[ci], sampleMeta[match(rownames(mvalt$mm),sampleMeta$fileName), target_col]),]
# mvalt$mpp2 = mvalt$mpp[!grepl(controlL[ci], sampleMeta[match(rownames(mvalt$mpp),sampleMeta$fileName), target_col]),]

# mvalt$mpvalt1 = mpvalt1 = mpvalt1[,apply(mpvalt1,2,function(x) !all(x==0))] 
# mvalt$mpvalt2 = mpvalt2 = mpvalt2[,apply(mpvalt2,2,function(x) !all(x==0))] 
# mvalt$mppvalt1 = mppvalt1 = mppvalt1[,apply(mppvalt1,2,function(x) !all(x==0))] 
# mvalt$mppvalt2 = mppvalt2 = mppvalt2[,apply(mppvalt2,2,function(x) !all(x==0))] 

# how many cell pop/ files are significat
# dim(mpvalt)
# sum(mpvalt!=0)
# dim(mppvalt)
# sum(mppvalt!=0)
# dim(mvalt$mpvalt1)
# sum(mvalt$mpvalt1!=0)
# dim(mvalt$mppvalt1)
# sum(mvalt$mppvalt1!=0)
# dim(mvalt$mpvalt2)
# sum(mvalt$mpvalt2!=0)
# dim(mvalt$mppvalt2)
# sum(mvalt$mppvalt2!=0)


# for (pe in c("","PEER")) {
#   if (pe=="") m = mvalt$mm
#   if (pe=="PEER") m = mvalt$mpp
#   for (nmi in names(mvalt)) {
#     # pe is mcp, as is nmi
#     if (pe=="" & grepl("mpp",nmi)) next()
#     if (pe=="PEER" & !grepl("mpp",nmi)) next()
#     pm = phenoMeta[match(colnames(mvalt[[nmi]]),phenoMeta$phenotype),]
#     mi = mi0 = m[,match(pm$phenotype,colnames(m))]
#     mi[mi<1] = 1
#     mi = log(mi)
#     val = unlist(lapply(unique(pm$phenolevel), function(x) sample(mi[,pm$phenolevel==x],min(sampleno,sum(pm$phenolevel==x)))))
#     group = unlist(lapply(unique(pm$phenolevel), function(x) rep(paste0("layer ",x,": ",sum(mi0[,pm$phenolevel==x]<1),"/",sum(pm$phenolevel==x)*nrow(mi0)," m cell pops x n files = 0"),min(sampleno,sum(pm$phenolevel==x)))))
#     dens_plot(val=val,group=as.character(group),binwidth=20,filename=paste0("~/projects/IMPC/result/P1/Sanger_SPLEEN/plots/counts_sample-",sampleno,"_",nmi,pe,adjust,".html"),title=paste0("Distribution of ln(normalized cell count) in different layers (<","sampleno"," sampled)"))
#   }
# }


