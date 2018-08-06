#working directory, need to be changed
dir<-"C:\\Users\\auz5836\\Documents\\GitHub\\GroupLasso\\V4" 

#need to import library(dr) and library(MASS)
source(sprintf("%s\\siri.R",dir))
source(sprintf("%s\\simu.R",dir))
source(sprintf("%s\\siri.fit.R",dir))

############
#SIRI setup#
############
#number of slices
H<-5 

#effective dimension for simple model [1..Q]
Q<-2 

#CV fold
K.fold<-10 

#sample size
n<-dim(x)[1] 

#number of predictors
d<-dim(x)[2] 

#grid of thresholds on chi-square quantiles
#alpha.list<-c(1-1/d,1-.5/d,1-.25/d,1-.05/d,1-.01/d) 
alpha.list<-c(0.9,0.95,0.97,0.99,0.995,0.999)

#sis selection bound
range.sis<-min(d,floor(n/log(n))) 

#simple model selection bound
range.linear<-min(d,floor(0.5*n/H)) 

#augmented model selection bound
range.interact<-min(d,floor(0.25*n/H)) 

#iteration for ISIS
niter<-2 

##########
#Run SIRI#
##########
#Run SIRI and save the result into "results"
results<-siri(x,y,H,Q,K.fold,alpha.list,niter,range.linear,range.interact,range.sis)

######################
#Determine thresholds#
######################
error<-1024 #minimum absolute error
errori<-0 #best threshold index
setA<-NULL
for(i in 1:(length(alpha.list)*(Q+1)))
{
    if(length(results[[i]]$cv$ssr)>0)
    {
      
        if(results[[i]]$cv$ssr<error)
        {
            error<-results[[i]]$cv$ssr
            errori<-i
            setA<-c(results[[i]]$result$linear.set,results[[i]]$result$interact.set)
        } else if(results[[i]]$cv$ssr==error&length(c(results[[i]]$result$linear.set,results[[i]]$result$interact.set))<length(setA)) {
            error<-results[[i]]$cv$ssr
            errori<-i
            setA<-c(results[[i]]$result$linear.set,results[[i]]$result$interact.set)
        }
    }
}
i<-errori
if(i>0)
{
    setA<-c(results[[i]]$result$linear.set,results[[i]]$result$interact.set)
} else {
    setA<-NULL
}
print("Variables selected by SIRI: ")
print(setA) #selected set of variables
truth<-which(beta!=0)
print("True relevant variables: ")
print(truth) #variables used in simulation

for(i in seq_along(results)){
  setA<-c(results[[i]]$result$linear.set,results[[i]]$result$interact.set)
  print(setA)
}
#####################
#Use SIRI to predict#
#####################
#To fit SIRI model given data x, y and 
#obtain predictions on a new set of predictor xnew,
#run the following commands:
alpha = 1-1/d
xold = x[1:floor(n/2),]
yold = y[1:floor(n/2)]
xnew = x[ceiling(n/2):n,]
ynew = y[ceiling(n/2):n]
result = siri.iter(xold,yold,H,Q,alpha,niter,range.linear,range.interact,range.sis)
model = siri.fit(xold,yold,result)
pred = siri.pred(xnew,result,model)
#Prediction results are stored in pred$pred.value and pred$pred.label
