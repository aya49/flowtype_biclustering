# use countadj to add a group column to meta_file; to group samples by date
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
# options(device="cairo")
options(na.rm=T)

countThres = 2000 #columns/rows must have at least goodCount number of elements that is greater than countThres, else delete
good_count = 3
good_sample = 3 #need more than good_sample samples per class for class to be included

target_col = "gene" #column with control/experiment
id_col = "fileName"
split_col = NULL #NULL if split by nothing
time_col = "date"
order_cols = c("barcode","sample","specimen","date")

control = str_split(controlL[ci],"[|]")[[1]]

methods = c("AMOC") #changepoint analysis; AMOC at most one change, PELT is poisson (quick exact), BinSeg (quick, approx), SegNeigh (slow, exact)
usemethod="AMOC" #changepoint analysis; use this method to determine changepoints for pvalue calculation

feat_type = c("file-cell-countAdj")
feat_count = c("file-cell-countAdj")












start = Sys.time()

meta_file = get(load(paste0(meta_file_dir,".Rdata")))
m0 = get(load(paste0(feat_dir,"/",feat_type,".Rdata")))
mc = get(load(paste0(feat_dir,"/",feat_count,".Rdata")))

mm = trimMatrix(m0,TRIM=T, mc=mc, sampleMeta=meta_file, sampleMeta_to_m1_col=id_col, target_col=target_col, control=control, order_cols=order_cols, colsplitlen=NULL, countThres=countThres)
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


g = getGTindex(meta_file_ordered[,target_col], control, good_sample)
ko = g$expIndex
wt = g$controlIndex #wildtype

#prepare plot colours by date
meta_file_ordered_factor = as.data.frame(sapply(meta_file_ordered, function(x) as.numeric(factor(x, ordered=T))))
ts = meta_file_ordered_factor[,time_col]
tswt = ts[wt]
tsc = c(1,1+which(diff(ts)!=0))
tscwt = c(1,1+which(diff(tswt)!=0))
tcolour = heat.colors(length(unique(ts))+25)[1:length(unique(ts))]


sm = meta_file_ordered
m = m_ordered
# cols = seq(from=1,to=ncol(m),by=2)


## create kalman filtering line plots for all files, for use in plots

cat("Kalman filtering for all; ")
kmff = kmf(m, 1)
fkffitall = kmff$fkffitall
statsfitall = kmff$statsfitall


## plot one single phenotypes and its changepoints ------------------------------------------------
cat(" changepoint; ")
i = 1
y <- as.numeric(m[,i])
ywt <- y[wt]

ylim <- c(min(ywt), max(ywt))
mvaluewt <- cpt.mean(fkffitall[[i]][wt],method=usemethod, Q=20, penalty="MBIC", minseglen=5)


## group days by change in mean; incomplete, this only works for AMOC, large amounts of dates grouped into one.
if (feat_type==feat_count & nrow(sm)==nrow(meta_file)) {
  group = rep(0,nrow(meta_file))
  enddate=0
  for (i in 1:length(mvaluewt@cpts)) {
    startdate = min(meta_file[meta_file[,time_col]>enddate,time_col])
    enddate = sm[wt[mvaluewt@cpts[i]],time_col]
    group[meta_file[,time_col]>=startdate & meta_file[,time_col]<=enddate] <- i
  }
  meta_file$group = group
  save(meta_file, file=paste0(meta_file_dir,".Rdata"))
}


cat("\n centre ", centre, " ",TimeOutput(start)," \n",sep="")





