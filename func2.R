### Logistic
f <- function(beta,X,y) { -t(y)%*%(X%*%beta) + sum(log(1+exp(X%*%beta))) } # objective function
gradf <- function(beta,X,y) { -t(X)%*%(y-plogis(X%*%beta)) } # gradient

### Ordinary
f <- function(beta,X,y){ 0.5*norm(X%*%beta - y, "F")^2 }
gradf <- function(beta,X,y){ t(X)%*%(X%*%beta - y) }


### Penalty
split_beta2<-function(beta,m_X,m_W,m_G,m_I){
  ma<-rep(1:4,c(m_X,m_W,m_G,m_I))
  beta<-split(beta,ma)
  names(beta)=c("X","W","G","I")
  return(beta)
}
split_X2<-function(X,m_X,m_W,m_G,m_I){
  ma<-rep(1:(2+m_G),c(m_X,m_W,rep(1,m_G)))
  ma<-c(ma,rep((2+1):(2+m_G),rep(1,m_I)))
  X_split<-split(X,ma)
  X_split<-lapply(X_split,function(x) x=matrix(x,nrow=n))
  return(X_split)
}
# Penalty parameters
group_penalty2<-function(X,m_X,m_W,m_G,m_I){
  X_split<-split_X2(X,m_X,m_W,m_G,m_I)
  para_I<-sapply(X_split[-1:-2], function(x) norm(as.matrix(x[,2]),'F'))
  para_mI<-sapply(X_split, function(x) norm(x,'F'))[-1:-2]
  para<-list("I"=para_I,"MI"=para_mI)
  return(para)
}
g2 <- function(X,beta,m_X,m_W,m_G,m_I,lambda1,lambda2) {
  beta<-split_beta2(beta,m_X,m_W,m_G,m_I)
  para<-group_penalty2(X,m_X,m_W,m_G,m_I)
  penalty<-lambda1*norm(as.matrix(para$I*beta$I),'1')+lambda2*sum(para$mI*sqrt(beta$G^2+beta$I^2))
  return(penalty)
}
proxg2 <- function(X,beta,m_X,m_W,m_G,m_I,tau,lambda1,lambda2) { 
  beta<-split_beta2(beta,m_X,m_W,m_G,m_I)
  para<-group_penalty2(X,m_X,m_W,m_G,m_I)
  
  #gamma<-mapply(FUN=function(x,y,z1,z2) { if( x!=0 || y!=0 ) { if(tau*lambda1*z1<abs(x)) return(sign(x)*(x-lambda1*tau*z1*sign(x))*(1-lambda2*tau*z2/sqrt((x-lambda1*tau*z1*sign(x))^2+(y)^2))) else return(0)} else return(0) }, beta$I, beta$G, 1, para$MI)
  #alpha<-mapply(FUN=function(x,y,z1,z2) { if( x!=0 || y!=0 ) { if(tau*lambda1*z1<abs(x)) return(y*(1-lambda2*tau*z2/sqrt((x-lambda1*tau*z1*sign(x))^2+(y)^2))) else return(0)} else return(0) }, beta$I, beta$G, 1, para$MI)
  
  beta_temp<-cbind(beta$G,beta$I)
  beta_temp<-split(beta_temp,rep(1:m_G,rep(1,m_G)))
  #beta_temp<-apply(beta_temp,1,FUN=function(x) { max(0,1-lambda2*tau/(sum(x^2))^0.5) *x})
  beta_temp<-mapply(FUN=function(x,y,z) { max(0,1-lambda2*tau*z/(x[1]^2+(sign(x[2])*max(abs(x[2])-lambda1*tau*y,0))^2)^0.5) * c(x[1],sign(x[2])*max(abs(x[2])-lambda1*tau*y,0))},beta_temp,para$I,para$MI)
  beta$G<-beta_temp[1,]
  beta$I<-beta_temp[2,]
  
  #beta$G<-alpha
  #beta$I<-gamma
  
  beta<-unlist(beta,use.names = F)
  return(beta)
}
