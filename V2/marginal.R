diffExp<-function(X,y,m_X,m_W,m_G,m_I){
  
  X_base<-X[,c(1:(m_X+m_W))]
  ma<-rep(1:m_G,c(rep(1,m_G)))
  ma<-c(ma,rep(1:m_G,rep(1,m_I)))
  X_G<-X[,-c(1:(m_X+m_W))]
  X_G<-split(X_G,ma)
  X_G<-lapply(X_G,function(x) x=matrix(x,nrow=n))
  
  coeff<-matrix(0,ncol=m_X+m_W+2,nrow=m_G)
  pval<-matrix(1,ncol=m_X+m_W+2,nrow=m_G)
  
  for(i in seq_along(X_G)){
    X<-cbind(X_base,X_G[[i]])
    colnames(X)<-c(paste("base",c(1:m_X),sep=""),"T","G_main","G_inter")
    L<-lm(y~-1+X)
    coeff[i,]<-L$coefficients
    pval[i,which(is.na(coeff[i,])==F)]<-summary(L)$coefficients[,4]
  }
  colnames(coeff)<-names(L$coefficients)
  colnames(pval)<-names(L$coefficients)
  
  return(list("coeff"=coeff,"pval"=pval))
}



pval_G<-NULL
pval_I<-NULL

for(i in 1:t){
  print(i)
  
  true_beta<-sim_beta2(m_X,m_W,m_G,main_zero,inter_zero,bit=T)
  
  
  X<-sim_X2(m_X,m_W,m_G)
  y<-X%*%true_beta+rnorm(n,sd=sdErr)                          
  
  true_beta<-split_beta2(true_beta,m_X,m_W,m_G,m_I)
  
  rst<-diffExp(X,y,m_X,m_W,m_G,m_I)
  
  pval_G<-c(pval_G,rst$pval[which(true_beta$G!=0),4])
  pval_I<-c(pval_I,rst$pval[which(true_beta$I!=0),5])
}

hist(pval_G,n=20,main = "Distribution of p-values for true main effect")
hist(pval_I,n=20,main = "Distribution of p-values for true interaction effect")
