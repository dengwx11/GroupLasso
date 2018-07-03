### Logistic
f <- function(beta,X,y) { -t(y)%*%(X%*%beta) + sum(log(1+exp(X%*%beta))) } # objective function
gradf <- function(beta,X,y) { -t(X)%*%(y-plogis(X%*%beta)) } # gradient

### Ordinary
f <- function(beta,X,y){ 0.5*norm(X%*%beta - y, "F")^2 }
gradf <- function(beta,X,y){ t(X)%*%(X%*%beta - y) }


### Penalty
split_beta<-function(beta,m_X,m_W,m_G1,m_G2,m_I){
  ma<-rep(1:5,c(m_X,m_W,m_G1,m_G2,m_I))
  beta<-split(beta,ma)
  names(beta)=c("X","W","G1","G2","I")
  return(beta)
}
split_X<-function(X,m_X,m_W,m_G1,m_G2,m_I){
  ma<-rep(1:(2+m_G1+m_G2),c(m_X,m_W,rep(1,m_G1),rep(1,m_G2)))
  ma<-c(ma,rep((2+m_G1+1):(2+m_G1+m_G2),rep(1,m_I)))
  X_split<-split(X,ma)
  X_split<-lapply(X_split,function(x) x=matrix(x,nrow=n))
  return(X_split)
}
# Penalty parameters
group_penalty<-function(X,m_X,m_W,m_G1,m_G2,m_I){
  X_split<-split_X(X,m_X,m_W,m_G1,m_G2,m_I)
  para<-sapply(X_split, function(x) norm(x,'F'))[-1:-2]
  para<-split(para,rep(1:2,c(m_G1,m_G2)))
  names(para)<-c("G1","G2")
  return(para)
}
g <- function(X,beta,m_X,m_W,m_G1,m_G2,m_I,lambda1,lambda2) {
  beta<-split_beta(beta,m_X,m_W,m_G1,m_G2,m_I)
  para<-group_penalty(X,m_X,m_W,m_G1,m_G2,m_I)
  penalty<-lambda1*norm(as.matrix(para$G1*beta$G1),'1')+lambda2*sum(para$G2*sqrt(beta$G2^2+beta$I^2))
  return(penalty)
  }
proxg <- function(X,beta,m_X,m_W,m_G1,m_G2,m_I,tau,lambda1,lambda2) { 
  beta<-split_beta(beta,m_X,m_W,m_G1,m_G2,m_I)
  para<-group_penalty(X,m_X,m_W,m_G1,m_G2,m_I)
  
  beta$G1<-sign(beta$G1)*(sapply(abs(beta$G1) - tau*lambda1*para$G1,FUN=function(x) {max(x,0)})) 
  beta_temp<-cbind(beta$G2,beta$I)
  beta_temp<-split(beta_temp,rep(1:m_G2,rep(1,m_G2)))
  #beta_temp<-apply(beta_temp,1,FUN=function(x) { max(0,1-lambda2*tau/(sum(x^2))^0.5) *x})
  beta_temp<-mapply(FUN=function(x,y) { max(0,1-lambda2*tau*y/(sum(x^2))^0.5) *x},beta_temp,para$G2)
  beta$G2<-beta_temp[1,]
  beta$I<-beta_temp[2,]
  beta<-unlist(beta,use.names = F)
  return(beta)
}
