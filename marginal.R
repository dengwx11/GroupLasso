diffExp<-function(X,y,m_X,m_W,m_G,m_I){
  
  X_base<-X[,c(1:(m_X+m_W))]
  ma<-rep(1:m_G,c(rep(1,m_G)))
  ma<-c(ma,rep(1:m_G,rep(1,m_I)))
  X_G<-X[,-c(1:(m_X+m_W))]
  X_G<-split(X_G,ma)
  X_G<-lapply(X_G,function(x) x=matrix(x,nrow=n))
  
  coeff<-matrix(0,ncol=1+m_X+m_W+2,nrow=m_G)
  pval<-matrix(0,ncol=1+m_X+m_W+2,nrow=m_G)
  
  for(i in seq_along(X_G)){
    X<-cbind(X_base,X_G[[i]])
    colnames(X)<-c(paste("base",c(1:m_X),sep=""),"T","G_main","G_inter")
    L<-lm(y~X)
    coeff[i,]<-L$coefficients
    pval[i,]<-summary(L)$coefficients[,4]
  }
  colnames(coeff)<-names(L$coefficients)
  colnames(pval)<-names(L$coefficients)
  
  return(list("coeff"=coeff,"pval"=pval))
}

diffExp(X,y,m_X,m_W,m_G,m_I)

