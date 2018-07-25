
# Simulation Setup
t=100 # iteration times

n=1000
m_X=5
m_W=1
m_G=400
m_I=m_G
p<-m_X+m_W+m_G+m_I

main_zero=floor(m_G*0.25)
inter_zero=floor(m_G*0.15)

main_nonzero=floor(m_G*0.1)
inter_nonzero=floor(m_G*0.05)
inter_nonzero=floor(m_G*0.1)


# Design matrix
sim_X<-function(m_X,m_W,m_G){
  X0<-matrix(rnorm(n*m_X),n,m_X)
  W<-(sample(c(0,1),nrow(X0),replace = T)*2-1)
  G<-matrix(rnorm(n*m_G),n,m_G)
  X<-cbind(X0,W,G,W*G)
  return(X)
}

# Coefficients
sim_beta<-function(m_X,m_W,m_G,main_zero,inter_zero,bit=TRUE,hier=TRUE){
  beta_X<-matrix(rnorm(m_X),m_X,1)
  beta_W<-matrix(rnorm(m_W),m_W,1)
  beta_G<-matrix(rnorm(m_G),m_G,1)
  zero_main<-sample(1:m_G,main_zero,replace = F)
  zero<-sample(1:m_G,inter_zero,replace = F)
  beta_I<-matrix(rnorm(m_G),m_G,1)
  
  if(bit){
    beta_X<-beta_X+sign(beta_X)*0.4
    beta_W<-beta_W+sign(beta_W)*0.4
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

sim_beta_const<-function(m_X,m_W,m_G,main_nonzero,inter_nonzero,both_nonzero,const=c(3),heir=TRUE){
  const<-c(const,-const)
  
  beta_X<-matrix(rnorm(m_X),m_X,1)+as.matrix(sample(const,m_X,replace = T))
  beta_W<-matrix(rnorm(m_W),m_W,1)+as.matrix(sample(const,m_W,replace = T))
  
  
  beta_G<-as.matrix(sample(const,m_G,replace = T))
  beta_I<-as.matrix(sample(const,m_G,replace = T))
  
  
  if(heir){
  beta_G[-c(1:main_nonzero)]<-0   # Only two are nonzero
  beta_I[-c(1:inter_nonzero)]<-0   # Onl three are nonzero
  } else{
    beta_G[-c(1:(both_nonzero+main_nonzero)),]<-0
    beta_I[-c(1:both_nonzero,(both_nonzero+main_nonzero+1):(both_nonzero+main_nonzero+inter_nonzero)),]<-0
  }
  
  beta<-rbind(beta_X,beta_W,beta_G,beta_I)
  return(beta)
}


logistic<-function(X,beta){
  prob<-exp(X%*%beta)/(1+exp(X%*%beta))
  y<-sapply(prob,FUN = function(x) rbinom(1,1,prob=x))
  return(y)
}


# Parameters
tau1<-0.1
lambda<-1
x0<-double(p)
sdErr<-1



