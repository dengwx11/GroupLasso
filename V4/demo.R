
sigma<-cov_block(10,.5,3)
X<-sim_X(1,10,sigma,100)
beta<-sim_beta(0,1,10,5,6,T,F)
y<-X%*%beta
