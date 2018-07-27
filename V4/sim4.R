library(partitions)


cov_block<-function(p,rho,block_num){
  block<-rmultinom(n = 1, size = p, prob = rep(1/block_num, block_num))
  blockstart<-cumsum(c(1,block))[-(block_num+1)]
  blockend<-cumsum(block)

  fun<-function(i,j){ (max(which(i>=blockstart))==max(which(j>=blockstart)))*rho/max(abs(i-j),0) } 
  corrmat<-outer(1:p, 1:p , Vectorize(fun) )
  diag(corrmat)<-1
  
  return(corrmat)
}

