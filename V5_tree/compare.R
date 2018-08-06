library(rpart)
library(randomForest)
library(BMA)
library(leaps)
library(BoomSpikeSlab)

n=100
m_G<-50

sigma<-cov_block(m_G,0.3,10)
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

simu<-data.frame(Y=y,X=x)
L<-lm(y~x[,1])
y.res<-L$residuals
simu$Y<-y.res
simu<-simu[,-2]


#### Trees
model <- rpart(Y~ ., data = simu)
par(xpd = NA) # otherwise on some devices the text is clipped
plot(model)
text(model, digits = 3)
print(model, digits = 2)
rsq.rpart(model)


model <- randomForest(Y~.,   data=simu)
#print(model) # view results 
#importance(model)
treerst<-order(importance(model),decreasing = T)[1:length(truth)]

### BMA
bicfit<-bicreg(x[,-1],y.res,strict = T)
bicrst<-bicfit$namesx[order(bicfit$probne0,decreasing = T)][1:length(truth)]

### Stepwise
a<-regsubsets(x=x,y=y,method="forward",nvmax = length(truth),force.in = 1)
steprst<-a$vorder[1:(length(truth)+1)]-1

### Spike-and-Slab
m <- lm.spike(Y~.,niter = 1000,data=simu,bma.method = "SSVS",error.distribution = "gaussian")
summary(m)
plot(m)
