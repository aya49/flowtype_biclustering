B2PS <- function(data, sideData=NULL, Kp, Kt, iterations, alpha_p = 1, alpha_t = 1, alpha_e = 1, alpha_sd = 1, conv_amount = 15, conv_upbound = 0){
  if(min(min(data))==0) data = data + 1
  withSideData=!is.null(sideData)
  if(withSideData && min(min(sideData))==0) sideData = sideData + 1
  Np = ncol(data)
  Nt = nrow(data)
  logKp=log(Kp)
  logKt=log(Kt)
  logDBLMax=log(.Machine$double.xmax)
  # Clustering
  kp = c()
  kt = c()
  npINkp = c()
  ntINkt = c()
  # Bicluster Counters & Parameters
  nesINkpkt = array(0,dim = c(Kp, Kt, 2))
  neINpFORkt = array(0,dim = c(Np, Kt, 2))
  neINtFORkp = array(0,dim = c(Nt, Kp, 2))
  # Side Data Counters
  nvsINktkt = array(0,dim = c(Kt, Kt, 2))
  nvINtFORkt = array(0,dim = c(Nt, Kt, 2))
  
  # Random Initialization
  kp = sample(c(1:Kp),Np,replace = T)
  kt = sample(c(1:Kt),Nt,replace = T)
  for(i in 1:Kp){
    npINkp[i]=length(which(kp==i))
    if(npINkp[i]!=0){
      for(j in 1:Kt){
        if(length(which(kt==j))==0) next
        nesINkpkt[i,j,1]=length(which(data[which(kt==j),which(kp==i)]==1))
        nesINkpkt[i,j,2]=length(which(data[which(kt==j),which(kp==i)]==2))
      }
      for(t in 1:Nt){
        neINtFORkp[t,i,1]= length(which(data[t,which(kp==i)]==1))
        neINtFORkp[t,i,2]= length(which(data[t,which(kp==i)]==2))
      }
    }
  }
  for(j in 1:Kt){
    ntINkt[j]=length(which(kt==j))
    if(ntINkt[j]!=0){
      if(withSideData){
        for(jj in 1:Kt){
          nvsINktkt[j,jj,1]=length(which(sideData[which(kt==j),which(kt==jj)]==1))
          nvsINktkt[j,jj,2]=length(which(sideData[which(kt==j),which(kt==jj)]==2))
        }
        for(t in 1:Nt){
          nvINtFORkt[t,j,1]=length(which(sideData[t,which(kt==j)]==1))
          nvINtFORkt[t,j,2]=length(which(sideData[t,which(kt==j)]==2))
        }
      }
      for(p in 1:Np){
        neINpFORkt[p,j,1]=length(which(data[which(kt==j),p]==1))
        neINpFORkt[p,j,2]=length(which(data[which(kt==j),p]==2))
      }
    }
  }
  
  
  # Sampling
  p.diffL = rep(0,conv_amount)
  t.diffL = rep(0,conv_amount)
  for(ite in 1:iterations){
    print(paste0("Iteration: ",ite))
    model.prob = 0
    ## Transcript Clusters
    t.diff = 0
    for(t in 1:Nt){
      ### take out + variables
      tempkt=kt[t]
      ntINkt[tempkt]=ntINkt[tempkt]-1
      nesINkpkt[,tempkt,] = nesINkpkt[,tempkt,] - neINtFORkp[t,,]
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] - (2-data[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] - (data[t,]-1)
      if (withSideData){
        nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]-nvINtFORkt[tempkt,]
        nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]-nvINtFORkt[tempkt,]
        nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
        for(tt in 1:Nt)
        {
          nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]-ifelse(sideData[tt,t]==1,1,0)
          nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]-ifelse(sideData[tt,t]==2,1,0)
        }
      }
      
      ### compute probability
      Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
      t.prob=log(ntINkt + alpha_t / Kt)+
        colSums(log(Theta[,,1])*matrix(rep(neINtFORkp[t,,1],Kt),nrow = Kp,ncol = Kt,byrow = F))+
        colSums(log(Theta[,,2])*matrix(rep(neINtFORkp[t,,2],Kt),nrow = Kp,ncol = Kt,byrow = F))
      
      ### put back
      t.prob.max=max(t.prob)
      t.prob.max = t.prob.max - (logDBLMax - logKt)
      t.prob = exp(t.prob - t.prob.max)
      total.prob = sum(t.prob)
      tempkt=sample(1:Kt,1,prob = t.prob/total.prob)
      if (kt[t] != tempkt)
        t.diff=t.diff+1
      model.prob = model.prob + log(t.prob[tempkt] / total.prob)
      
      kt[t]=tempkt
      ntINkt[tempkt]=ntINkt[tempkt]+1
      nesINkpkt[,tempkt,] = nesINkpkt[,tempkt,] + neINtFORkp[t,,]
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] + (2-data[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] + (data[t,]-1)
      if (withSideData){
        nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]+nvINtFORkt[tempkt,]
        nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]+nvINtFORkt[tempkt,]
        nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
        for(tt in 1:Nt)
        {
          nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]+ifelse(sideData[tt,t]==1,1,0)
          nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]+ifelse(sideData[tt,t]==2,1,0)
        }
      }
    }#t
    
    ## Patient Clusters
    p.diff = 0
    for(p in 1:Np){
      ### take out
      tempkp=kp[p]
      npINkp[tempkp]=npINkp[tempkp]-1
      nesINkpkt[tempkp,,] = nesINkpkt[tempkp,,] - neINpFORkt[p,,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] - (2-data[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] - (data[,p]-1)
      
      ### compute probability
      reverseExp=matrix(rep(2-data[,p],Kp),nrow = Kp,ncol = Nt,byrow = T)
      Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
      p.prob=rowSums(log(abs(reverseExp-Theta[,kt,2])))+ #add other data
        log(npINkp + alpha_p / Kp)
      p.prob.max=max(p.prob)
      p.prob.max = p.prob.max - (logDBLMax - logKp)
      p.prob = exp(p.prob - p.prob.max)
      total.prob = sum(p.prob)
      tempkp=sample(1:Kp,1,prob = p.prob/total.prob)
      if (kp[p] != tempkp)
        p.diff=p.diff+1
      model.prob = model.prob + log(p.prob[tempkp] / total.prob)
      
      ### put back
      kp[p]=tempkp
      npINkp[tempkp]=npINkp[tempkp]+1
      nesINkpkt[tempkp,,] = nesINkpkt[tempkp,,] + neINpFORkt[p,,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] + (2-data[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] + (data[,p]-1)
    }#p
    
    print("Difference:")
    print(paste0("p: ",p.diff))
    print(paste0("t: ",t.diff))
    
    ## test convergence
    p.diffL = append(p.diffL[-1], p.diff)
    t.diffL = append(t.diffL[-1], t.diff)
    # if (iterations == conv_amount) {
    #   p.diffLd = diff(p.diffL)
    #   t.diffLd = diff(t.diffL)
    #   p.diffLdsum = sum(abs(p.diffLd))
    #   t.diffLdsum = sum(abs(t.diffLd))
    # } else if (iterations > conv_amount) {
    #   p.diffLdsum = p.diffLdsum - abs(p.diffLd[1]) + 
    #   p.diffLd = append(p.diffLd[-1], p.diff - p.diffL[conv_amount])
    #   t.diffLd = append(t.diffLd[-1], t.diff - t.diffL[conv_amount])
    # }
    if (all(p.diffL < conv_upbound) & all(t.diffL < conv_upbound)) break
    
  }#ite
  
  Pi=(npINkp+alpha_p/Kp)/(Np+alpha_p)
  Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
  model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
             bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt,
             Pi=Pi,Theta=Theta,
             sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
  return(model)
}
