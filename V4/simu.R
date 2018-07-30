library(MASS)
n<-200 #sample size
d<-1000 #number of predictors
x<-mvrnorm(n,rep(0,d),diag(1,d))
sigma<-0.2
truth<-sample(1:d,3)#relevant predictors

#scenario 2.1
#truth<-truth[c(1,2)]
#y<-x[,truth[1]]*x[,truth[2]]+rnorm(n)*sigma

#scenario 2.2
#y<-(1+x[,truth[1]]+x[,truth[2]])*x[,truth[3]]+rnorm(n)*sigma

#scenario 2.3
y<-(x[,truth[1]]+x[,truth[2]])*x[,truth[3]]+rnorm(n)*sigma

#senario 2.4
#y<-x[,truth[1]]*x[,truth[2]]*x[,truth[3]]+rnorm(n)*sigma

#scenario 2.5
#truth<-truth[c(1,2)]
#y<-x[,truth[1]]^2*x[,truth[2]]+rnorm(n)*sigma

#scenario 2.6
#y<-x[,truth[1]]/(x[,truth[2]]+x[,truth[3]])+rnorm(n)*sigma
