n=200
m_G<-200

sigma<-cov_block(m_G,.5,7)
x<-sim_X(1,m_G,sigma,n)
beta<-sim_beta(0,1,m_G,floor(0.9*m_G),floor(0.9*m_G),T,F)
y0<-x%*%beta

SNR<-10

noise<-rnorm(n,sd=1)
SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
y<-y0+SNRmtl*noise  
colnames(x)<-c(1:dim(x)[2])
truth<-which(beta!=0)
a<-regsubsets(x=x,y=y,method="forward")
summary(a)
