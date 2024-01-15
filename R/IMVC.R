#'@title Integrated Mean Variance Correlation
#'@description This function is used to calculate the integrated mean variance correlation between two vectors
#'@param y is a numeric vector
#'@param x is a numeric vector
#'@param K is the number of quantile levels
#'@param NN is the number of B spline basis, default is 3
#'@param type is an indicator for measuring linear or nonlinear correlation, "linear" represents linear correlation and "nonlinear" represents linear or nonlinear correlation using B splines
#'@importFrom splines bs
#'@importFrom quantreg rq
#'@export
IMVC<-function(y,x,K,NN=3,type){
  new_env1 = new.env()
  new_env1$x = x
  new_env1$y = y
  n=length(new_env1$y)
  int_tau=sapply(1:K,function(i) (i/(K+1)))
  wei=rep(1/K,K)
  xx=x
  if(type=="nonlinear")
  {
    B_x<-bs(xx,df=NN)
    B_x=as.matrix(B_x)
    rec_int=c()
    for (l in 1:length(int_tau)) {
      beta_con=rq(y~B_x,tau=int_tau[l])$coeff
      xx_new=cbind(rep(1,n),B_x)
      fitted_con=xx_new%*%as.matrix(beta_con)
      ind1=c()
      for (t in 1:n) {
        ind1[t]=((length(which(y<=fitted_con[t]))/n)-int_tau[l])^2
      }
      rec_int[l]=sum(ind1)/n
    }
    IMVcor=6*sum(wei*rec_int)
  }
  if(type=="linear")
  {
    B_x<-xx
    B_x=as.matrix(B_x)
    rec_int=c()
    for (l in 1:length(int_tau)) {
      beta_con=rq(y~B_x,tau=int_tau[l])$coeff
      xx_new=cbind(rep(1,n),B_x)
      fitted_con=xx_new%*%as.matrix(beta_con)
      ind1=c()
      for (t in 1:n) {
        ind1[t]=((length(which(y<=fitted_con[t]))/n)-int_tau[l])^2
      }
      rec_int[l]=sum(ind1)/n
    }
    IMVcor=6*sum(wei*rec_int)
  }
  IMVcor
}

#'@title Integrated Mean Variance Correlation Based Screening
#'@description This function is used to select important features using integrated mean variance correlation
#'@param y is the response vector
#'@param x is the covariate matrix
#'@param K is the number of quantile levels
#'@param d is the size of selected variables
#'@param NN is the number of B spline basis, default is 3
#'@param type is an indicator for measuring linear or nonlinear correlation, "linear" represents linear correlation and "nonlinear" represents linear or nonlinear correlation using B splines
#'@importFrom splines bs
#'@importFrom quantreg rq
#'@export
IMVCS<-function(y,x,K,d,NN=3,type){
  new_env2 = new.env()
  new_env2$x = x
  new_env2$y = y
  n=length(y)
  p=ncol(x)
  int_tau=sapply(1:K,function(i) (i/(K+1)))
  wei=rep(1/K,K)
  scr_con=c()
  if(type=="nonlinear"){
    for (j in 1:p) {
      xx=x[,j]
      B_x<-bs(xx,df=NN)
      B_x=as.matrix(B_x)
      rec_int=c()
      for (l in 1:length(int_tau)) {
        beta_con=rq(y~B_x,tau=int_tau[l])$coeff
        xx_new=cbind(rep(1,n),B_x)
        fitted_con=xx_new%*%as.matrix(beta_con)
        ind1=c()
        for (t in 1:n) {
          ind1[t]=((length(which(y<=fitted_con[t]))/n)-int_tau[l])^2
        }
        rec_int[l]=sum(ind1)/n
      }
      scr_con[j]=6*sum(wei*rec_int)
    }
  }
  if(type=="linear"){
    for (j in 1:p) {
      xx=x[,j]
      B_x<-xx
      B_x=as.matrix(B_x)
      rec_int=c()
      for (l in 1:length(int_tau)) {
        beta_con=rq(y~B_x,tau=int_tau[l])$coeff
        xx_new=cbind(rep(1,n),B_x)
        fitted_con=xx_new%*%as.matrix(beta_con)
        ind1=c()
        for (t in 1:n) {
          ind1[t]=((length(which(y<=fitted_con[t]))/n)-int_tau[l])^2
        }
        rec_int[l]=sum(ind1)/n
      }
      scr_con[j]=6*sum(wei*rec_int)
    }
  }
  sel=order(scr_con,decreasing = T)[1:d]
  sel
}

#'@title Integrated Mean Variance Correlation Based Hypothesis Test
#'@description This function is used to test significance of linear or nonlinear correlation using integrated mean variance correlation
#'@param y is the response vector
#'@param x is the univariate covariate vector
#'@param K is the number of quantile levels
#'@param num_per is the number of permutation times
#'@param NN is the number of B spline basis, default is 3
#'@param type is an indicator for measuring linear or nonlinear correlation, "linear" represents linear correlation and "nonlinear" represents linear or nonlinear correlation using B splines
#'@importFrom splines bs
#'@importFrom quantreg rq
#'@importFrom expm sqrtm
#'@importFrom CompQuadForm liu
#'@export
IMVCT<-function(x,y,K,num_per,NN=3,type){
  n=length(y)
  new_env4 = new.env()
  new_env4$x = x
  new_env4$y = y
  int_tau=sapply(1:K,function(i) (i/(K+1)))
  wei=rep(1/K,K)
  xx=x
  if(type=="nonlinear"){
    B_x<-bs(xx,df=NN)
    B_x=as.matrix(B_x)
    rec_int=c()
    for (l in 1:length(int_tau)) {
      beta_con=rq(y~B_x,tau=int_tau[l])$coeff
      xx_new=cbind(rep(1,n),B_x)
      fitted_con=xx_new%*%as.matrix(beta_con)
      ind1=c()
      for (t in 1:n) {
        ind1[t]=((length(which(y<=fitted_con[t]))/n)-int_tau[l])^2
      }
      rec_int[l]=sum(ind1)/n
    }
    IMVT=6*n*sum(wei*rec_int)
    
    pr_imv=c()
    new_env4$pr_imv = pr_imv
    for (o in 1:num_per) {
      pr_sim=sample(1:n,n,replace =FALSE)
      y_pe=y[pr_sim]
      rec_int=c()
      for (l in 1:length(int_tau)) {
        beta_con=rq(y_pe~B_x,tau=int_tau[l])$coeff
        xx_new=cbind(rep(1,n),B_x)
        fitted_con=xx_new%*%as.matrix(beta_con)
        ind1=c()
        for (t in 1:n) {
          ind1[t]=((length(which(y_pe<=fitted_con[t]))/n)-int_tau[l])^2
        }
        rec_int[l]=sum(ind1)/n
      }
      pr_imv[o]=6*n*sum(wei*rec_int)
    }
    p.value=length(which(pr_imv>=IMVT))/num_per 
  }
  if(type=="linear"){
    B_x<-xx
    B_x=as.matrix(B_x)
    rec_int=c()
    for (l in 1:length(int_tau)) {
      beta_con=rq(y~B_x,tau=int_tau[l])$coeff
      xx_new=cbind(rep(1,n),B_x)
      fitted_con=xx_new%*%as.matrix(beta_con)
      ind1=c()
      for (t in 1:n) {
        ind1[t]=((length(which(y<=fitted_con[t]))/n)-int_tau[l])^2
      }
      rec_int[l]=sum(ind1)/n
    }
    IMVT=6*n*sum(wei*rec_int)
    
    lam_k1=(int_tau*(1-int_tau))/K
    lam_k2=6*lam_k1
    sig_Z=diag(1,K)
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        sig_Z[i,j]=sig_Z[j,i]=(min(int_tau[i],int_tau[j])*
                                 (1-max(int_tau[i],int_tau[j])))/sqrt(int_tau[i]*
                                                                        int_tau[j]*(1-int_tau[i])*(1-int_tau[j]))
      }
    }
    sig_sqr=sqrtm(sig_Z)
    AA=diag(lam_k2)
    BB=sig_sqr%*%AA%*%sig_sqr
    value_BB=eigen(BB)$values
    lambda=value_BB
    p.value=liu(IMVT, lambda, h = rep(1, length(lambda)),
        delta = rep(0, length(lambda)))
  }
  p.value
}



#'@title Integrated Mean Variance Correlation Based FDR Control
#'@description This function is used for FDR control with integrated mean variance correlation
#'@param y is the response vector
#'@param x is the covariate matrix
#'@param K is the number of quantile levels
#'@param NN is the number of B spline basis, default is 3
#'@param numboot is the size of bootstrap samples
#'@param timeboot is the number of bootstrap times for computing standard deviation of the IMVC
#'@param true_signal is the true active set
#'@param null_method is the estimation method for proportion of true null hypotheses. Choices are "lfdr", "mean", "hist" or "convest"
#'@param alpha is the nominal FDR level
#'@importFrom splines bs
#'@importFrom quantreg rq
#'@importFrom GGMridge getEfronp
#'@importFrom limma propTrueNull
#'@return A list of FDP, power and selected variables
#'@export
IMVCFDR<-function(y,x,K,NN=3,numboot,timeboot,true_signal,null_method,alpha){
  new_env5 = new.env()
  new_env5$x = x
  new_env5$y = y
  n=length(y)
  p=ncol(x)
  int_tau=sapply(1:K,function(i) (i/(K+1)))
  wei=rep(1/K,K)
  scr_con=c()
  for (j in 1:p) {
    xx=x[,j]
    B_x<-bs(xx,df=NN)
    B_x=as.matrix(B_x)
    rec_int=c()
    for (l in 1:length(int_tau)) {
      beta_con=rq(y~B_x,tau=int_tau[l])$coeff
      xx_new=cbind(rep(1,n),B_x)
      fitted_con=xx_new%*%as.matrix(beta_con)
      ind1=c()
      for (t in 1:n) {
        ind1[t]=((length(which(y<=fitted_con[t]))/n)-int_tau[l])^2
      }
      rec_int[l]=sum(ind1)/n
    }
    scr_con[j]=6*sum(wei*rec_int)
  }
  
  
  scr_sig=c()
  for (j in 1:p) {
    boot_rec=c()
    for (boot in 1:timeboot) {
      boot_s=sample(n,numboot,replace = FALSE)
      xx=x[boot_s,j]
      yy=y[boot_s]
      B_x<-bs(xx,df=NN)
      B_x=as.matrix(B_x)
      rec_int=c()
      for (l in 1:length(int_tau)) {
        beta_con=rq(yy~B_x,tau=int_tau[l])$coeff
        xx_new=cbind(rep(1,numboot),B_x)
        fitted_con=xx_new%*%as.matrix(beta_con)
        ind1=c()
        for (t in 1:numboot) {
          ind1[t]=((length(which(yy<=fitted_con[t]))/numboot)-int_tau[l])^2
        }
        rec_int[l]=sum(ind1)/numboot
      }
      boot_rec[boot]=6*sum(wei*rec_int)
    }
    scr_sig[j]=sd(boot_rec)
  }
  
  neww_z=scr_con/scr_sig
  index_bor=rbinom(p,1,0.5)
  index_bor[which(index_bor==0)]=-1
  new_z=index_bor*neww_z
  
  Efron=getEfronp(new_z)
  F0muu=Efron$mu0hat
  F0sig=Efron$sigma0hat
  
  
  conver_p=2*pnorm(-abs(new_z),mean =F0muu,sd=F0sig)
  pi_est0=1-propTrueNull(conver_p, method=null_method)
  if(pi_est0==0||pi_est0>(length(true_signal)/p)||pi_est0<((length(true_signal)/p)/2)){pi_est=length(true_signal)/p}
  else{pi_est=pi_est0}
  
  bandw1=bw.nrd0(scr_sig)
  bandw2=bw.nrd0(new_z)
  
  GauK<-function(u){
    gau=(1/sqrt(2*pi))*exp(-u^2/2)
    gau
  }
  
  f0star=c()
  for (k in 1:p) {
    f0star[k]=sum(sapply(1:p,function(i) (1/bandw1)*GauK((scr_sig[k]-scr_sig[i])/bandw1)*((1/bandw2)*
                                                                                            GauK((new_z[k]-new_z[i])/bandw2))))/sum(sapply(1:p,function(i) (1/bandw1)*GauK((scr_sig[k]-scr_sig[i])/bandw1)))
  }
  
  
  weight0=c()
  for (k in 1:p) {
    weight0[k]=1-min((((1-pi_est)*dnorm(new_z[k],mean =F0muu,sd=F0sig ))/f0star[k]),1)
  }
  
  f1_ini=c()
  for (k in 1:p) {
    f1_ini[k]=sum(sapply(1:p,function(i) (1/bandw1)*weight0[i]*GauK((scr_sig[k]-scr_sig[i])/bandw1)*((1/bandw2)*
                                                                                                       GauK((new_z[k]-new_z[i])/bandw2))))/sum(sapply(1:p,function(i) (1/bandw1)*weight0[i]*GauK((scr_sig[k]-scr_sig[i])/bandw1)))
  }
  
  weight1=c()
  for (k in 1:p) {
    TT1=((1-pi_est)*dnorm(new_z[k],mean =F0muu,sd=F0sig))/((1-pi_est)*
                                                             dnorm(new_z[k],mean =F0muu,sd=F0sig)+pi_est*f1_ini[k])
    weight1[k]=1-TT1
  }
  
  f1_update=c()
  for (k in 1:p) {
    f1_update[k]=sum(sapply(1:p,function(i) (1/bandw1)*weight1[i]*GauK((scr_sig[k]-scr_sig[i])/bandw1)*((1/bandw2)*
                                                                                                          GauK((new_z[k]-new_z[i])/bandw2))))/sum(sapply(1:p,function(i) (1/bandw1)*weight1[i]*GauK((scr_sig[k]-scr_sig[i])/bandw1)))
  }
  
  LFT=c()
  for (k in 1:p) {
    LFT[k]=((1-pi_est)*dnorm(new_z[k],mean =F0muu,sd=F0sig))/((1-pi_est)*
                                                                dnorm(new_z[k],mean =F0muu,sd=F0sig)+pi_est*f1_update[k])
  }
  
  order_LFT=order(LFT,decreasing = FALSE)
  
  LFTrec=c()
  for (k in 1:p) {
    LFTrec[k]=sum(LFT[order_LFT][1:k])/k
  }
  
  indlft=max(which(LFTrec<alpha))
  name_result=order_LFT[1:indlft]
  sett=true_signal
  FDR=length(setdiff(name_result,sett))/length(name_result)
  Power=length(intersect(sett,name_result))/length(sett)
  return(list(selected=name_result,FDR=FDR,Power=Power))
}