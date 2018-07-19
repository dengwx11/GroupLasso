
# Simulation Setup
t=100 # iteration times

n=100
m_X=2
m_W=1
p_G=100
m_G1=p_G
m_G2=p_G
m_I=p_G
p<-m_X+m_W+m_G1+m_G2+m_I

main_zero=floor(p_G*0.7)
inter_zero=floor(p_G*0.8)


# Design matrix
sim_X<-function(m_X,m_W,p_G){
  X0<-matrix(rnorm(n*m_X),n,m_X)
  W<-sample(c(0,1),nrow(X0),replace = T)
  G<-matrix(rnorm(n*p_G),n,p_G)
  X<-cbind(X0,W,G,G,W*G)
  return(X)
}

# Coefficients
sim_beta<-function(m_X,m_W,p_G,main_zero,inter_zero){
  beta_X<-matrix(rnorm(m_X),m_X,1)
  beta_W<-matrix(rnorm(m_W),m_W,1)
  beta_G1<-matrix(rnorm(p_G),p_G,1)
  zero_main<-sample(1:p_G,main_zero,replace = F)
  beta_G1[zero_main]<-0   # Only two are nonzero
  beta_G2<-matrix(rnorm(p_G),p_G,1)
  beta_I<-matrix(rnorm(p_G),p_G,1)
  zero<-sample(1:p_G,inter_zero,replace = F)
  beta_G2[zero]<-0  # Only three are nonzero
  beta_I[zero]<-0   # Onl three are nonzero
  
  bit<-which(abs(beta_I)<0.05)
  beta_G1[bit]<-beta_G2[bit]+beta_G1[bit]
  beta_G2[bit]<-0
  beta_I[bit]<-0
  
  
  beta<-rbind(beta_X,beta_W,beta_G1,beta_G2,beta_I)
  return(beta)
}


logistic<-function(X,beta){
  prob<-exp(X%*%beta)/(1+exp(X%*%beta))
  y<-sapply(prob,FUN = function(x) rbinom(1,1,prob=x))
  return(y)
}


# Parameters
tau1<-10
lambda1<-0.5
lambda2<-0.5
x0<-double(p)




