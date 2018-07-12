
# Simulation Setup
t=100 # iteration times

n=100
m_X=2
m_W=1
m_G=30
m_I=m_G
p<-m_X+m_W+m_G+m_I

main_zero=floor(m_G*0.8)
inter_zero=floor(m_G*0.7)


# Design matrix
sim_X2<-function(m_X,m_W,m_G){
  X0<-matrix(rnorm(n*m_X),n,m_X)
  W<-(sample(c(0,1),nrow(X0),replace = T)*2-1)
  G<-matrix(rnorm(n*m_G),n,m_G)
  X<-cbind(X0,W,G,W*G)
  return(X)
}

# Coefficients
sim_beta2<-function(m_X,m_W,m_G,main_zero,inter_zero,bit=FALSE){
  beta_X<-matrix(rnorm(m_X),m_X,1)
  beta_W<-matrix(rnorm(m_W),m_W,1)
  beta_G<-matrix(rnorm(m_G),m_G,1)
  zero_main<-sample(1:m_G,main_zero,replace = F)
  zero<-sample(1:m_G,inter_zero,replace = F)
  beta_G[intersect(zero_main,zero)]<-0   # Only two are nonzero
  beta_I<-matrix(rnorm(m_G),m_G,1)
  beta_I[zero]<-0   # Onl three are nonzero
  
  if(bit){
    bit<-which(abs(beta_I)<0.05)
    beta_I[bit]<-0
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



