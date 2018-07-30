library(partitions)
library(MASS)

cov_block<-function(p,rho,block_num){
  block<-rmultinom(n = 1, size = p, prob = rep(1/block_num, block_num))
  blockstart<-cumsum(c(1,block))[-(block_num+1)]
  blockend<-cumsum(block)

  fun<-function(i,j){ (max(which(i>=blockstart))==max(which(j>=blockstart)))*rho/max(abs(i-j),0) } 
  corrmat<-outer(1:p, 1:p , Vectorize(fun) )
  diag(corrmat)<-1
  
  return(corrmat)
}

sim_X<-function(m_W,m_G,sigma,n){
  X<-mvrnorm(n,rep(0,m_G),sigma)
  W<-(sample(c(0,1),n,replace = T)*2-1)
  X<-cbind(W,X,W*X)
  return(X)
  
}


sim_beta<-function(m_X,m_W,m_G,main_zero,inter_zero,bit=TRUE,hier=TRUE){
  if(m_X!=0){
    beta_X<-matrix(rnorm(m_X),m_X,1)
  } else{ beta_X<-NULL}
  if(m_W!=0){
    beta_W<-matrix(rnorm(m_W),m_W,1)
  } else{ beta_W<-NULL }
  beta_G<-matrix(rnorm(m_G),m_G,1)
  zero_main<-sample(1:m_G,main_zero,replace = F)
  zero<-sample(1:m_G,inter_zero,replace = F)
  beta_I<-matrix(rnorm(m_G),m_G,1)
  
  if(bit){
    if(m_X!=0){
      beta_X<-beta_X+sign(beta_X)*0.4
    }
    if(m_W!=0){
      beta_W<-beta_W+sign(beta_W)*0.4
    }
    beta_G<-beta_G+sign(beta_G)*0.4
    beta_I<-beta_I+sign(beta_I)*0.4
  }
  
  if(hier){
    beta_G[intersect(zero_main,zero)]<-0   
  } else{
    beta_G[zero_main]<-0
  }
  beta_I[zero]<-0   # Onl three are nonzero
  
  beta<-rbind(beta_X,beta_W,beta_G,beta_I)
  return(beta)
}

