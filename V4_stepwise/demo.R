library(spcov)


n=100
m_G<-20

sigma<-cov_block(m_G,.5,5)
#sigma<-GenerateCliquesCovariance(10,10,0.8)
binprob<-runif(m_G)
x<-sim_X_cate(1,m_G,sigma,n,binprob)
beta<-sim_beta(m_X=0,m_W=1,m_G,main_nonzero=0.1,inter_nonzero=0.1,both_nonzero=0.1,bit=T,heir=F)
y0<-x%*%beta

SNR<-10
noise<-rnorm(n,sd=1)
SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
y<-y0+SNRmtl*noise  
colnames(x)<-c(1:dim(x)[2])
truth<-which(beta!=0)


library(leaps)
a<-regsubsets(x=x,y=y,method="forward",nvmax = 9,force.in = 1)
summary(a)

