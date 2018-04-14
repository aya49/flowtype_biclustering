### Matrix (sample x phenotype) Preprocessing ########################################

##Input: matrix / list of matrices with columnname = cell population;;; layer to prune at, each cell population must have atleast goodcount number of samples with more than countThres cells, each class must have more thatn good_sample amount of samples or else prune
##Output: trimmed matrix / list of matrices
trimMatrix <- function(m0,TRIM=T,mc=NULL,sampleMeta,sampleMeta_to_m1_col="fileName", target_col=NULL ,control="+_+", order_cols=NULL, colsplitlen=NULL, k=NULL, konly=F,countThres=0, goodcount=1, good_sample=1) {
  
  require(stringr)  
  #get matrix colnames/rownames (sample/celltype features)
  m0cn = colnames(m0)
  m0cn = sapply(str_split(m0cn,"_"), function(x) ifelse(length(x)>1, x[2], x[1]) )
  m0rn = rownames(m0)
  
  #get feature layers if feature names represent cell types
  colsplitlen = cell_type_layers(m0cn)
  if (is.null(colsplitlen)) print("missing cell type features or column names; columns untrimmed.")
  
  #get to-delete low count phenotype indices; CountAdj should be first one
  colIndexC = rep(T,length(m0cn))
  rowIndexC = rep(T,length(m0rn))
  if(countThres>0 & !is.null(mc)) {
    mcordercol = match(m0cn,colnames(mc))
    mcordercol = mcordercol[!is.na(mcordercol)]
    mcorderrow = match(m0rn,rownames(mc))
    mcorderrow = mcorderrow[!is.na(mcorderrow)]
    mc0 = mc[mcorderrow,mcordercol]
    colIndexC = apply(as.matrix(mc0), 2, function(x) sum(x>=countThres)>=goodcount)
    rowIndexC = apply(as.matrix(mc0), 1, function(x) sum(x>=countThres)>=goodcount)
  }
  
  #get to-delete high no of marker phenotypes
  colIndexL = rep(T,length(m0cn))
  if (!is.null(k) | k == 0) {
    if (k<max(colsplitlen)) {
      if (konly) {
        colIndexL = colsplitlen == k
      } else {
        colIndexL = colsplitlen <= k
      }
    } 
  }
    
      
  
  #trim matrix based on cell count and layer
  # if (m0list) {
  #   m = m0[colIndexC & colIndexL]
  #   if (TRIM & F) {
  #     m = lapply(m, function(y0) {
  #       colIndex0 = apply(y0, 2, function(x) sum(x!=0)>0 )
  #       # rowIndex0 = apply(y, 1, function(x) sum(x!=0)>0 )
  #       # y = y[rowIndex0,colIndex0]
  #       y = y0[,colIndex0]
  #       if (length(colIndex0)==1) {
  #         y = matrix(y,ncol=1)
  #         rownames(y) = rownames(y0)
  #         colnames(y) = colnames(y0)[colIndex0]
  #       }
  #       if (sum(dim(y))>0) return(y)
  #       return(NULL)
  #     })
  #     m = Filter(Negate(is.null), m)
  #   }
  #   
  # } else {
  m = m0[rowIndexC, colIndexC & colIndexL]
  if (TRIM) {
    colIndex0 = apply(m, 2, function(x) sum(x!=0)>0 )
    rowIndex0 = apply(m, 1, function(x) sum(x!=0)>0 )
    m = m[rowIndex0,colIndex0]
  }
  # }
  
  #check matrix; if matrix is a list, merge; if matrix is empty, skip
  if (is.null(m)) return (NULL)
  if (all(m==0)) return (NULL)
  
  
  sm = NULL
  if (!is.null(sampleMeta)) {
    # if (m0list) {
    #   smorder = match(rownames(m[[1]]),sampleMeta[[sampleMeta_to_m1_col]])
    # } else {
    smorder = match(rownames(m),sampleMeta[[sampleMeta_to_m1_col]])
    # }
    smorder = smorder[!is.na(smorder)]
    sm = sampleMeta[smorder,]
    
    #order samples by date etc.
    if (!is.null(order_cols)) {
      for (order_col in order_cols) { # different centres have different filename colnames
        if (length(grep(order_col, colnames(sm)))>0) sm = sm[order(sm[,order_col]),]
      }
      # if (m0list) {
      #   mrowind = match(sm[,sampleMeta_to_m1_col],rownames(m[[1]]))
      #   morder = mrowind[!is.na(mrowind)]
      #   m = lapply(m, function(x) {
      #     xo = x[morder,]
      #     if (is.null(dim(xo))) xo = matrix(xo, ncol=1, dimnames=list(names(xo),colnames(x))) 
      #     return(xo)
      #   })
      # } else {
      mrowind = match(sm[,sampleMeta_to_m1_col],rownames(m))
      m = m[mrowind[!is.na(mrowind)],]
      # }
    }
    
    #exclude genotypes with less than 3 samples
    if (good_sample>1 & !is.null(target_col)) {
      g = getGTindex(sm[,target_col], control, good_sample)
      goodind = as.vector(c(unlist(g$expIndex[g$goodexpIndex]), g$controlIndex))
      goodind = sort(goodind)
      sm = sm[goodind,]
      # if(m0list) {
      #   m = lapply(m, function(x) {
      #     xo =  x[goodind,]
      #     if (is.null(dim(xo))) xo = matrix(xo, ncol=1, dimnames=list(names(xo),colnames(x))) 
      #     return(xo)
      #   })
      # } else {
      m = m[goodind,]
      # }
    }
  }
  
  return(list(m=m,sm=sm))
}







##Input: matrix
##Output: if features (rownames) have +/- symbols, returns corresponding feature layers; else returns NULL
cell_type_layers <- function(cell_types) {
  colsplit = str_split(cell_types,"[+-]")
  colsplitlen = sapply(colsplit, function(x) length(x)) - 1
  return(colsplitlen)
}






##Input: matrix with rows in order of date/time etc.
##Output: kalman filter for each column in matrix
kmf <- function(m,cols) {
  fkffitall = NULL
  statsfitall = NULL
  for (i in cols) {
    yall = as.numeric(m[,i])
    if (length(unique(yall))==1) {
      fkffitall[[i]] <- statsfitall[[i]] <- yall
    } else {
      ## Set constant parameters:
      dt <- ct <- matrix(0) 
      Zt <- Tt <- matrix(1)
      a0 <- yall[1]           # Estimation of the first sample count
      P0 <- matrix(100)     # Variance of 'a0'
      ## Estimate parameters 23min if TS
      fit.fkf <- optim(c(HHt = var(yall, na.rm = TRUE) * .5,
                         GGt = var(yall, na.rm = TRUE) * .5),
                       fn = function(par, ...) -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
                       yt = rbind(yall), a0 = a0, P0 = P0, dt = dt, ct = ct,
                       Zt = Zt, Tt = Tt, check.input = FALSE)
      ## Filter Nile data with estimated parameters:
      fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fit.fkf$par[1]), GGt = matrix(fit.fkf$par[2]), yt = rbind(yall))
      fkffitall[[i]] <- fkf.obj$att[1,]
      ## Compare with the stats' structural time series implementation: 5min if TS
      statsfitall[[i]] <- fitted(StructTS(yall, type = "level"))
    }
  }
  return(list(fkffitall=fkffitall,statsfitall=statsfitall))
}














### Dist Score ########################################################################

## Input: Distance matrix, class list (numeric)
## Output: NCA score
NCA_score <- function(x,y,delta=1,fast=F,doUnderflow=T) {
  require(Brobdingnag)
  require(foreach)
  require(doMC)
  
  registerDoMC(no_cores)
  
  #preprocess matrix
  x = as.matrix(x)
  if (nrow(x)!=length(y)) {
    cat("x and y dimensions don't match")
    return(NULL) #check if x pairs with y
  }
  power = 2
  if (fast) power = 1
  if (!fast) delta = 1
  xe = -(x^power)/delta ## PREVENT UNDERFLOW
  underflow = F
  if (length(which(xe<(-700)))>0) {
    if (!doUnderflow) {
      cat(length(which(xe<700)), "values too big")
      return(NULL)
    }
    underflow = T
    diag(xe) = -Inf
  } else {
    xe = exp(xe)
    diag(xe) = 0
  }
  
  #pij = prob of xi picking xj as neighbour (fill top half)
  if (underflow) {
    pij = list()
    for (i in 1:length(y)) {
      pij[[i]] = brob(xe[i,])/sum(brob(xe[i,-i])) ## PREVENT UNDERFLOW
    }
  } else {
    pij = matrix(0, nrow=nrow(xe), ncol=ncol(xe))
    for (i in 1:length(y)) {
      pij[i,] = xe[i,]/sum(xe[i,-i])
    }
  }
  
  
  #pi = prob of classifying xi correctly; piy = prob of classifying points of class y correctly
  yt = table(y)
  yf = as.numeric(as.factor(y)) # as factor
  pyt = rep(0,length(yt))
  pi = rep(0,length(y))
  yi = lapply(1:length(yt), function(yn) return(which(y==names(yt)[yn])) )
  
  
  if (underflow) {
    # pi = sum(pij[[1]][yi[[yf[1]]]])
    # if (length(yf)==2) pi = cbrob(pi,sum(pij[[2]][yi[[yf[2]]]]))
    # if (length(yf)>2) {
    #   for (i in 2:length(yf)) {
    #     pi = cbrob(pi,sum(pij[[i]][yi[[yf[i]]]]))
    #   }
    # }
    # for (yn in 1:length(yt)) {
    #   pyt[yn] = as.numeric(sum(pi[yi[[yn]]]))
    # }
    for (yn in 1:length(yt)) { #for each class yn
      pi[yi[[yn]]] = sapply(yi[[yn]], function(i) as.numeric(sum(pij[[i]][yi[[yn]]])) )
      pyt[yn] = sum(pi[yi[[yn]]])
    }
  } else {
    for (yn in 1:length(yt)) { #for each class yn
      pi[yi[[yn]]] = sapply(yi[[yn]], function(i) sum(pij[i,yi[[yn]]]) )
      pyt[yn] = sum(pi[yi[[yn]]])
    }
  }
  names(pyt) = names(yt)
  return(list(pi=pi,yt=yt,pyt=pyt,p=sum(pyt))) #p is total score
}


## Input: 2 vector of labels
## Output: Precision, Recall, Specificity, F, Accuracy
f.measure.comembership = function(la,cl) {
  require(clusteval)
  ftpn = comembership_table(la,cl)
  tn = ftpn$n_00
  tp = ftpn$n_11
  fn = ftpn$n_10
  fp = ftpn$n_01
  
  p=tp/(tp+fp)
  r=tp/(tp+fn)
  sp=tn/(tn+fp)
  return(list(p=p, r=r, sp=sp, f_comember=2*p*r/(p+r), a=(tp+tn)/(sum(unlist(ftpn)))))
}










## Input: Clt, a matrix
## Output: Outputs column inds of duplicate columns
duplicateindM = function(clt) {
  delcol_dup = duplicated(t(clt))
  col_dup = 1:ncol(clt); col_dup[delcol_dup] = sapply(which(delcol_dup),function(x) {
    for (xi in which(!delcol_dup)) {
      if (identical(clt[,xi],clt[,x]) ) { a = xi; break }
    }
    return(a)
  })
  return(col_dup)
}









### Feature Generation ########################################################################



## Input: phenoMeta ($phenolevel, $phenotype, $phenocode)
## Output: Index of list of children
getphenoChild <- function (phenoMeta, no_cores=1) {
  require(foreach)
  require(doMC)
  
  registerDoMC(no_cores)
  
  notFinalLevel = which(phenoMeta$phenolevel!=max(phenoMeta$phenolevel))
  ppcc = foreach (i=1:nrow(phenoMeta)) %dopar% {
    if (!sum(notFinalLevel%in%i)>0) return(list(pc=NULL,pcpn=NULL))
    pheni = unlist(strsplit(phenoMeta$phenocode[i],""))
    phenind = which(pheni!="0")
    zeroind = setdiff(1:length(pheni),phenind)
    childphenocode = as.vector(sapply(zeroind, function(x) { pi1=pi2=pheni; pi1[x]="1"; pi2[x] = "2"; return(c(paste(pi1,collapse=""),paste(pi2,collapse=""))) } ))
    childrenind = match(childphenocode, phenoMeta$phenocode)
    childrenind = childrenind[!is.na(childrenind)]
    if (length(childrenind)>0) { #if non-leaf node; this condition is only used when full hierarchy not available; otherwise, can take this out
      childsum = sapply(phenoMeta$phenocode[childrenind], function(x) sum(as.numeric(unlist(strsplit(x,"")))) ) #sum children phenocodes up, larger ones have +, smaller ones have - in new marker, split them up
      childplus = which(childsum==max(childsum))
      #split between positive and negative
      # phenoChild[[i]] = list()
      # phenoChild[[i]][[1]] = childrenind[-childplus]
      # phenoChild[[i]][[2]] = childrenind[childplus]
      pc = list()
      pc[[1]] = childrenind[-childplus]
      pc[[2]] = childrenind[childplus]
      
      ## if there are missing counterparts
      nopm1 = gsub("[+-]","",phenoMeta$phenotype[pc[[1]]])
      nopm2 = gsub("[+-]","",phenoMeta$phenotype[pc[[2]]])
      # unique1 = which(nopm1%in%setdiff(nopm1,nopm2))
      # unique2 = which(nopm2%in%setdiff(nopm2,nopm1))
      # if (length(unique1)>0) {
      # }
      # if (length(unique2)>0) {
      # }
      negposintersect1 = which(nopm1%in%intersect(nopm1,nopm2))
      negposintersect2 = which(nopm2%in%intersect(nopm2,nopm1))
      if (length(negposintersect1)>0) {
        # phenoChildpn[[i]] = list()
        # phenoChildpn[[i]][[1]] = phenoChild[[i]][[1]][negposintersect1]
        # phenoChildpn[[i]][[2]] = phenoChild[[i]][[2]][negposintersect2]
        pcnp = list()
        pcnp[[1]] = pc[[1]][negposintersect1]
        pcnp[[2]] = pc[[2]][negposintersect2]
        
        return(list(pc=pc, pcnp=pcnp))
      } else {
        return(list(pc=pc, pcnp=NULL))
      }
    } else {
      return(list(pc=NULL, pcnp=NULL))
    }
  }
  
  phenoChild = list()
  phenoChildpn = list()
  for (i in 1:length(ppcc)) {
    if (!is.null(ppcc[[i]]$pc)) phenoChild[[i]] = ppcc[[i]]$pc
    if (!is.null(ppcc[[i]]$pc)) phenoChildpn[[i]] = ppcc[[i]]$pcnp
  }
  
  phenoChild_ind <- which(vapply(phenoChild, Negate(is.null), NA))
  phenoChildpn_ind <- which(vapply(phenoChildpn, Negate(is.null), NA))
  
  phenoChild = phenoChild[phenoChild_ind]
  phenoChildpn = phenoChildpn[phenoChildpn_ind]
  
  names(phenoChild) = phenoMeta$phenotype[phenoChild_ind]
  names(phenoChildpn) = phenoMeta$phenotype[phenoChildpn_ind]
  
  return(list(phenoChild=phenoChild, phenoChild_ind=phenoChild_ind, phenoChildpn=phenoChildpn, phenoChildpn_ind=phenoChildpn_ind))
}

getphenoParent_phenoChild <- function(phenoChild,phenoChild_ind, phenoChildpn=NULL, phenoChildpn_ind=NULL, no_cores=1) {
  
  require(foreach)
  require(doMC)
  
  registerDoMC(no_cores)
  
  phenoParent = NULL
  phenoParentpn = NULL
  phenoParent_ind = NULL
  phenoParentpn_ind = NULL
  for(i in 1:length(phenoChild)) {
    parind = phenoChild_ind[i]
    children = unlist(phenoChild[[i]])
    childrennp = NULL
    childrenboth = NULL
    childrenonly = NULL
    if (sum(phenoChildpn_ind==phenoChild_ind[i])>0) childrennp = unlist(phenoChildpn[[which(phenoChildpn_ind==phenoChild_ind[i])]])
    if (length(childrennp)>0) {
      childrenboth = intersect(children, childrennp)
      childrenonly = setdiff(children, childrennp)
      if (length(childrenboth)>0) {
        for (ind in childrenboth) {
          if (!(sum(phenoParent_ind==ind)>0)) {
            phenoParent_ind = append(phenoParent_ind, ind)
            phenoParent[[ind]] = c(parind)
          } else {
            a = append(phenoParent[[ind]], parind)
            phenoParent[[ind]] = a
          }
          if (!(sum(phenoParentpn_ind==ind)>0)) {
            phenoParentpn_ind = append(phenoParent_ind, ind)
            phenoParentpn[[ind]] = c(parind)
          } else {
            phenoParentpn[[ind]] = append(phenoParentpn[[ind]], parind)
          }
        }
      }
      if (length(childrenonly)>0) {
        for (ind in childrenonly) {
          if (!(sum(phenoParent_ind==ind)>0)) {
            phenoParent_ind = append(phenoParent_ind, ind)
            phenoParent[[ind]] = c(parind)
          } else {
            phenoParent[[ind]] = append(phenoParent[[ind]], parind)
          }
        }
      }
    } else {
      for (ind in children) {
        if (!(sum(phenoParent_ind==ind)>0)) {
          phenoParent_ind = append(phenoParent_ind, ind)
          phenoParent[[ind]] = c(parind)
        } else {
          phenoParent[[ind]] = append(phenoParent[[ind]], parind)
        }
      }
    }
  }
  return(list(phenoParent=phenoParent[which(vapply(phenoChild, Negate(is.null), NA))], phenoParent_ind=sort(phenoParent_ind), phenoParentpn=phenoParentpn[which(vapply(phenoChild, Negate(is.null), NA))], phenoParentpn_ind=sort(phenoParentpn_ind)))  
}

getphenoParent <- function(phenoMeta, phenoChildpn, phenoChildpn_ind, no_cores) {
  require(foreach)
  require(doMC)
  
  registerDoMC(no_cores)
  
  phenoParent = foreach (i=(1:nrow(phenoMeta))) %dopar% {
    pheni = unlist(strsplit(phenoMeta$phenocode[i],""))
    phenind = which(pheni!="0")
    #enumerate all possible parents
    parentphenocode = sapply(phenind, function(x) { pi=pheni; pi[x]="0"; return(paste(pi,collapse="")) } )
    parentind = match(parentphenocode, phenoMeta$phenocode)
    parentind = parentind[!is.na(parentind)]
    if (length(parentind)>0) return(parentind)
    return(NULL)
  }
  
  phenoParent_ind = which(vapply(phenoParent, Negate(is.null), NA))
  phenoParent = phenoParent[phenoParent_ind]
  names(phenoParent) = phenoMeta$phenotype[phenoParent_ind]
  
  #delete parents from child whom doesn't have a twin under that parent.
  phenoParentpn = phenoParent
  phenoParentpn_ind = phenoParent_ind
  delind = c()
  for (i in 1:length(phenoParent)) {
    child = phenoParent_ind[i]
    parents = phenoParent[[i]]
    delparents = c()
    for (j in 1:length(phenoParent[[i]])) {
      parent = phenoParent[[i]][j]
      ind = (phenoChildpn_ind==parent)
      if (sum(ind)>0) {
        if (!sum(unlist(phenoChildpn[[which(ind)]])==child)>0) delparents = append(delparents, j)
      }
    }
    if (length(delparents)>0) parents = parents[-delparents]
    if (!length(parents)>0) {
      delind = append(delind, i)
    } else {
      phenoParentpn[[i]] = parents
    }
  }
  
  if (length(delind)>0) {
    phenoParentpn[delind] = NULL
    phenoParentpn_ind = phenoParentpn_ind[-delind]
  }
  
  return(list(phenoParent=phenoParent, phenoParent_ind=phenoParent_ind, phenoParentpn=phenoParentpn, phenoParentpn_ind=phenoParentpn_ind))
  
}




## Input: feature matrix (row = sample, col = phenotypes), phenocode or phenotype (prefer phenocode)
## Output: feature matrix with proportions of parent
prop_parent <- function(m, phen, list=F, par=T, no_cores=detectCores()-3) {
  if (ncol(m)!=length(phen)) {
    cat("ncol(matrix) doesn't match length(phen)")
    return(NULL)
  }
  
  #list of child phenotypes for each non-leaf phenotype
  children = phen_children(phen)
  parentind = which(!is.null(children))
  
  if (par) {
    mplist = foreach(i = 1:length(children), .combine='list') %dopar% { #for each phenotype
      if (is.null(children[[i]])) return(NULL)
      return(m[,children[[i]]]/m[,i])
    }
  } else {
    mplist = list()
    for (i in parentind) {
      mplist[[i]] = m[,children[[i]]]/m[,i]
    }
  }
  if (list) return(list(mp=mplist,parentind=parentind,children=children))
  mp = do.call(cbind, mplist)
  return(list(mp=mp,parentind=parentind,children=children))
}




## IN PROGRESS
## Input: normalized count matrix m, phenoChild list of children, phenoChild_ind indices of parents in phenoChild, whether this list is negative only (with positive counterpart), weighting function
child_entrop <- function(m, phenoChild, phenoChild_ind, pc_overlapneg=T) {
  #get entropy for each cell population
  require(entropy)
  m[m<1] = 1 #prevent Inf, NaN's
  for (i in 1:(nrow(m)-1)) {
    for (j in (i+1):nrow(m)) { #all combination of rows
      ratio = m[,i]/m[,j] #10/5=2, 5/10=.5;; log(2)=.69, log(.5)=-.69
      growth = m[,i]-m[,j]
      weight = matrix(0,nrow=2,ncol=length(ratio)) #row 1 = negative, row 2 = positive (more than parent)
      for (s in 1:length(phenoChild)) { # all features
        l = log(ratio[phenoChild[[s]]]/ratio[phenoChild_ind[s]])
      }
    }
  }
}



## Input: normalized count matrix m, phenoChild list of children, phenoChild_ind indices of parents in phenoChild, whether this list is negative only (with positive counterpart), weighting function
child_entropeasy <- function(m, phenoChild, phenoChild_ind, pc_overlapneg=T) { # use only overlapping indices, easier to normalize
  #get entropy for each cell population
  require(entropy)
  m[m<1] = 1 #prevent Inf, NaN's
  for (i in 1:(nrow(m)-1)) {
    for (j in (i+1):nrow(m)) { #all combination of rows
      ratio = m[,i]/m[,j] #10/5=2, 5/10=.5;; log(2)=.69, log(.5)=-.69
      growth = m[,i]-m[,j]
      weight = matrix(0,nrow=2,ncol=length(ratio)) #row 1 = negative, row 2 = positive (more than parent)
      for (s in 1:length(phenoChild)) { # all features
        l = log(ratio[phenoChild[[s]]]/ratio[phenoChild_ind[s]])
      }
    }
  }
}


# JS Divergence; 0<JSD<log(2)

# input: normalized count matrix m, phenoChild list of children, phenoChild_ind indices of parents in phenoChild, whether this list is negative only (with positive counterpart), weighting function
child_JSD <- function(m, phenoChildpn, phenoChildpn_ind) { # use only overlapping indices, easier to normalize
  #get entropy for each cell population
  m[m<1] = 1 #prevent Inf, NaN's
  for (i in 1:(length(phenoChildpn_ind)-1)) {
    for (j in (i+1):length(phenoChildpn_ind)) { #all combination of rows
      coli = phenoChildpn_ind[i]
      colj = phenoChildpn_ind[j]
      dist = rep(0,nrow(m)) #jsd for all features
      dist[m[,coli]==0 & m[,colj]==0] = 0
      dist[(m[,coli]!=0 & m[,colj]==0) | m[,coli]==0 & m[,colj]!=0] = log(2)
      non0parents = (m[,coli]!=0 & m[,colj]!=0) #parents that are not 000
      
      parenti = m[,coli]
      childreni = cbind(m[,unlist(phenoChildpn[[i]][[1]])])
      parentj = m[,colj]
      childrenj = cbind(m[,unlist(phenoChildpn[[j]][[1]])])
      
      dist[non0parents] = sapply(non0parents, function(x) return(jsd(
        mapply(c,as.list(childreni[x,]),as.list(parenti[x]-childreni[x,]),SIMPLIFY=F),
        mapply(c,as.list(childrenj[x,]),as.list(parentj[x]-childrenj[x,]),SIMPLIFY=F) )))
      
    }
  }
}







## Input: matrix (row = objects)
## Output: distance matrix










### Pvalue #########################################################################

## Significance; take note to put in same size control
## tstat = (a-mean(control))/se
## se = sd(control)/sqrt(length(control))
t.test.single <- function(control,a) { #input a vector and a single number
  se = sd(control) #/ sqrt(length(control))
  tstat = (a-mean(control)) / se
  p = 2*pt(-abs(tstat), df=length(control)-1) # 2 tailed (less stringent)
  return(p)
}


## find lv such that b= exp(log(a,lv))
getlv = function(a,b=100) return(exp(log(a)/(log(b))))









### Normalize ######################################################################

## Input: x matrix (phenotypes on columns, will convert in function) -- x0 is optional, plots according to x0 so if you want to plot more cell populations than those in x
## Output: f = normalize factors per sample; fidff = difference between f and peak of count ratio per sample
tmm <- function(x,x0=NULL,lib.size,refColumn,cutoff=Inf,plotimg=T,pngnames=NULL,mains=NULL,no_cores=detectCores()-1,samplesOnCol=F,Acutoff=-10) {
  require(foreach)
  require(doMC)
  registerDoMC(no_cores)
  
  if (is.null(x0)) x0 = x
  if (!samplesOnCol) { x = t(x); x0 = t(x0) }
  
  ## Taken from TMM 
  
  #f <- rep(NA,ncol(x))
  #fdiff <- rep(NA,ncol(x)) #diff between density peak and value (note: logged)
  ref <- x[,refColumn]
  refn <- lib.size[refColumn]
  doWeighting <- F
  logratioTrim <- .3
  sumTrim <- 0.05
  minlogR <- 1e-6 #min value of log2((obs/obsn)/(ref/refn))
  
  
  ff <- foreach(i=1:ncol(x), .combine = list, .maxcombine = ncol(x), .multicombine = T) %dopar% {
    #for(i in ncol(x):1) { cat(i," ",sep="")
    obs <- x[,i]
    obsn <- lib.size[i]
    #logR <- log2((obs/obsn)/(ref/refn))
    logR <- log2(obs/ref)			#log ratio of expression, accounting for library size
    absE <- (log2(obs/obsn) + log2(ref/refn))/2	#absolute expression
    v <- (obsn-obs)/obsn/obs + (refn-ref)/refn/ref	 #estimated asymptotic variance
    
    #remove infinite values, cutoff based on A
    fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
    logR <- logR[fin]
    absE <- absE[fin]
    v <- v[fin]
    
    if(max(abs(logR)) < minlogR) { return(list(f=1, fdiff=0)) # f[i] <- 1
    } else {
      
      #taken from the original mean() function
      n <- length(logR)
      loL <- floor(n * logratioTrim) + 1
      hiL <- n + 1 - loL
      loS <- floor(n * sumTrim) + 1
      hiS <- n + 1 - loS
      
      #keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
      #a fix from leonardo ivan almonacid cardenas, since rank() can return
      #non-integer values when there are a lot of ties
      keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)
      
      if(doWeighting) {
        fi <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
      } else { fi <- mean(logR[keep], na.rm=TRUE) } #f[i] <- mean(logR[keep], na.rm=TRUE) }
      
      #Results will be missing if the two libraries share no features with positive counts
      #In this case, return unity
      #if(is.na(f[i])) f[i] <- 0
      if(is.na(fi)) fi <- 0
      
      #check if close to peak; if not, switch to peak
      d <- density(log2((obs)/ref), na.rm=T)
      p <- as.matrix(findpeaks(d$y)); if(ncol(p)==1) p <- t(p)
      p1 <- d$x[p[which.max(p[,1]),2]]
      #fdiff[i] <- p1-f[i]
      fdiffi <- p1-fi
      
      if (plotimg) {
        pngname <- pngnames[i]
        png (file=pngname , width=700, height=1800)
        par(mfrow=c(3,1), mar=(c(5, 5, 4, 2) + 0.1))
        
        #plot(d); abline(v=f[i], col="red"); abline(v=p1, col="blue"); 
        plot(d); abline(v=fi, col="red"); abline(v=p1, col="blue"); 
        
        plot((x0[,i]+x0[,refColumn])/2, log(x0[,i]/x0[,refColumn]), cex=.5, main=paste(mains[i],": f=",fi, sep=""))
        #abline(h=f[i], col="red")
        abline(h=fi, col="red")
      }
      
      #if f[i] too far from peak
      #if (abs(f[i]-p1)>cutoff) {
      if (abs(fi-p1)>cutoff) {
        abline(h=p1, col="blue")
        #f[i] <- p1
        fi = p1
      }
      
      #f[i] <- 1/2^f[i]
      fi <- 1/2^fi
      
      #plot((matrixCount[,i]+matrixCount[,refColumn])/2, log2((matrixCount[,i]*f[i])/matrixCount[,refColumn]), cex=.5, main=paste("AFTER CHANGE: mean count vs. log2 fold change: ", sampleMeta$gene[i]," over refColumn ", sampleMeta$gene[refColumn],": f=",f[i], sep=""))
      if (plotimg) {
        plot((x0[,i]+x0[,refColumn])/2, log((x0[,i]*fi)/x0[,refColumn]), cex=.5, main=paste(mains[i],": f=",fi, sep=""))
        abline(h=0, col="red")
        dev.off()
      }
      
      return(list(f=fi, fdiff=fdiffi))
      
    }
  }
  #multiple of 1
  rm(x)
  
  f <- rep(NA,length(ff))
  fdiff <- rep(NA,length(ff)) #diff between density peak and value (note: logged)
  for (i in 1:length(ff)) {
    f[i] <- ff[[i]]$f
    try({ fdiff[i] <- ff[[i]]$fdiff })
  }
  
  return(list(f=f,fdiff=fdiff))
}








## Classify/Cluster ############################################

## Input: distance matrix, parameter values
## Output: table of labels/cluster assigned to each item in distance matrix

knntable = function(dd,knn,label) {
  testind = which(is.na(label))
  clt = t(matrix(sapply(testind, function(x) {
    topkind = order(dd[x,])
    topkind = topkind[-(topkind==x)]
    topkfn = colnames(dd)[topkind]
    topkla = label[match(topkfn,rownames(dd))]
    topkla = topkla[!is.na(topkla)]
    return(sapply(knn, function(y) return(Mode(topkla[1:y])) ))
  }),nrow=length(knn)))
  if (!length(knn)>1) clt = as.matrix(clt,ncol=1)
  colnames(clt) = knn
  rownames(clt) = rownames(dd)[testind]
  return(clt)
}


pamtable = function(dd,nclass,pamtries) {
  require(cluster)
  pp0t = sapply(1:length(pamtries), function(j) pam(dd, k=(nclass),diss=T)$clustering)
  if (is.null(dim(pp0t))) pp0t = matrix(pp0t,ncol=1)
  colnames(pp0t) = rep("none",ncol(pp0t))
  rownames(pp0t) = rownames(dd)
  return(pp0t)
}


spectable = function(mm,nclass,label,methods=c("rbf"),tries=1,savedist=NULL,savesim=NULL,replace='kerns') {
  #methods=c("rbf","poly","vanilla","tanh","laplace","bessel","anova","spline")
  require(kernlab)
  require(stringr)
  tryCatch({
    
    parlist = pp0t = NULL
    for (method in methods) {
      parlist0 = sil = NCA = cl = sim = dist = NULL
      for (i in 1:tries) {
        sp = specc(as.matrix(mm),kernel=paste0(method,"dot"),centers=nclass)
        cl[[i]] = sp@.Data
        parlist0[[i]] = unlist(sp@kernelf@kpar)
        sim[[i]] = kernelMatrix(kernel=match.fun(paste0(method,"dot"))(as.numeric(parlist0[[i]])),x=as.matrix(mm))
        dist[[i]] = get_graphd(sim[[i]])
        sil = append(sil,median(silhouette(label,dist[[i]])[,3]))
        NCA = append(NCA,NCA_score(as.matrix(dist[[i]]), label, doUnderflow=doUnderflow)$p)
      }
      maxsilind = which.max(sil)
      pp0t = cbind(pp0t, cl[[maxsilind]])
      parlist = c(parlist, paste0(paste0(names(sp@kernelf@kpar),"-",signif(parlist0[[maxsilind]],4),collapse="_"),"_silmed-",signif(max(sil),5),"_NCA-",signif(NCA[maxsilind],5),collapse=""))
      s = sim[[maxsilind]]
      d = dist[[maxsilind]]
      if (!is.null(savesim)) save(s,file=paste0(gsub(replace,paste0("kern_",method,"_",paste0(names(sp@kernelf@kpar),"-",signif(parlist0[[maxsilind]],4),collapse="_")),savesim),"_simmatrix.Rdata"))
      if (!is.null(savedist)) save(d,file=paste0(gsub(replace,paste0("kern_",method,"_",paste0(names(sp@kernelf@kpar),"-",signif(parlist0[[maxsilind]],4),collapse="_")),savesim),"_dist.Rdata"))
    }
    colnames(pp0t) = parlist
    rownames(pp0t) = rownames(mm)
    return(pp0t)
  }, error = function(e) {
    return(NULL)
  })
  
}

spec1table = function(sim,nclass) {
  require(kernlab)
  sim = sim/max(sim)
  diag(sim) = 1
  pp0t = matrix(specc(as.kernelMatrix(sim),centers=nclass)@.Data, ncol=1)
  colnames(pp0t) = "none"
  rownames(pp0t) = rownames(sim)
  return(pp0t)
}

hctable = function(dd,nclass,links=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")) {
  pp0t = sapply(links,function(x) return(cutree(hclust(as.dist(dd),method=x),k=nclass)))
  colnames(pp0t) = links
  rownames(pp0t) = rownames(dd)
  return(pp0t)
}


lvtable = function(sim,rwThres) {
  require(igraph)
  pp0t = sapply(rwThres, function(rwt) {
    d2 = sim
    tops = quantile(as.vector(sim),rwt)
    d2[d2<tops] = 0
    gr = graph_from_adjacency_matrix(d2, weighted=T,mode='undirected', diag=F)
    return(cluster_louvain(gr)$membership)
  })
  colnames(pp0t) = rwThres
  rownames(pp0t) = rownames(sim)
  return(pp0t)
}

dctable = function(mm,alpha=.85,nu=seq(0.0, 1.0, by=0.05)) {
  require(densitycut)
  clt = DensityCut(X=mm, alpha=alpha, nu=nu, show.plot=F)$cluster
  clt = matrix(clt,ncol=1)
  rownames(clt) = rownames(mm)
  colnames(clt) = "none"
  return(clt)
}

dc1table = function(dd,k=3,alpha=.85,nu=seq(0.0, 1.0, by=0.05)) {
  require(densitycut)
  require(FastKNN)
  
  #create knn.ind, knn.dist
  dd = as.matrix(dd)
  ki = t(sapply(1:nrow(dd), function(x) k.nearest.neighbors(x,dd,k=k)))
  kd = t(sapply(1:nrow(dd), function(x) dd[x,ki[x,]]))
  rownames(ki) = rownames(kd) = rownames(dd)
  colnames(ki) = colnames(kd) = 1:k
  
  clt = DensityCut(knn.dist=kd,knn.index=ki, alpha=alpha, nu=nu, show.plot=F)$cluster
  clt = matrix(clt,ncol=1)
  rownames(clt) = rownames(dd)
  colnames(clt) = "none"
  return(clt)
}


rw1table = function(sim,rwThres) {
  require(igraph)
  gr = graph_from_adjacency_matrix(sim, weighted=T,mode='undirected', diag=F)
  
  pp0t = sapply(rwThres, function(rwt) {
    tops = quantile(as.vector(sim),rwt)
    gr1 = delete.edges(gr, which(E(gr)$weight<tops))
    return(components(gr1)$membership)
  })
  colnames(pp0t) = rwThres
  rownames(pp0t) = rownames(sim)
  return(pp0t)
  
}








##Input: Row x Bicluster (TRUE/FALSE)
##Output: cluster array for each row; put row into largest cluster its in
rowxcluster_to_cluster <- function(rowxcluster) {
  clusters_size = apply(rowxcluster, 2, function(x) sum(x) )
  rclusters = apply(rowxcluster, 1, function(x) {
    cl = which(x)
    if (length(cl) == 0) cl = 0
    if (length(cl) > 1) cl = cl[which.max(clusters_size[cl])]
    return(cl)
  })
}

##Input: Bicluster x Col (TRUE/FALSE)
##Output: cluster array for each row; put row into largest cluster its in
clusterxcol_to_cluster <- function(clusterxcol) {
  clusters_size = apply(clusterxcol, 1, function(x) sum(x) )
  rclusters = apply(clusterxcol, 2, function(x) {
    cl = which(x)
    if (length(cl) == 0) cl = 0
    if (length(cl) > 1) cl = cl[which.max(clusters_size[cl])]
    return(cl)
  })
}



