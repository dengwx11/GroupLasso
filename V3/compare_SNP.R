## SNP compare

library(rpart)
library(randomForest)
library(BMA)
library(leaps)
library(BoomSpikeSlab)
library(SIS)

n=100
m_X<-5
m_W<-1
m_G<-50
m_I<-m_G
SNR<-10
tau1<-1

x0<-rep(0,dim(x)[2])
#### Cross Validation finding best lambda for Group Lasso
lamb_candidate<-c(1,1.5,2,2.5,3,4)
lamb_candidate2<-c(0.5,1,1.5,2)
sol_cv_SNP<-opt_lambda(x,y,f, gradf, g, proxg, x0, tau1, max_iters = 100, w = 10, 
                   backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                   eps_n = 1e-15,m_X,m_W,m_G,m_I,K=10,n=100, lamb_candidate, lamb_candidate2,restart=TRUE)
lamb_loc<-which(sol_cv_SNP$mean+sol_cv$var == min(sol_cv_SNP$mean+sol_cv$var), arr.ind = TRUE)
lamb_opt_glasso<-lamb_candidate[lamb_loc[1]]
lamb_opt2_glasso<-lamb_candidate2[lamb_loc[2]]
save(sol_cv_SNP,file="C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//sol_cv_glasso_SNP.RData")

#### Cross Validation finding best lambda for General Lasso
lamb_candidate<-c(15,25,30,35,40,45,50)
lamb_candidate2<-c(1,3,5,7)
sol_cv_SNP<-opt_lambda(x,y,f, gradf, glasso, proxglasso, x0, tau1, max_iters = 100, w = 10, 
                   backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                   eps_n = 1e-15,m_X,m_W,m_G,m_I,K=10,n=100, lamb_candidate, lamb_candidate2,restart=TRUE)
lamb_loc<-which(sol_cv_SNP$mean+sol_cv$var == min(sol_cv_SNP$mean+sol_cv_SNP$var), arr.ind = TRUE)
lamb_opt_lasso<-lamb_candidate[lamb_loc[1]]
lamb_opt2_lasso<-lamb_candidate2[lamb_loc[2]]
save(sol_cv_SNP,file="C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//sol_cv_lasso_SNP.RData")


#### Iteration
treerst_SNP<-list()
bicrst_SNP<-list()
steprst_SNP<-list()
glassorst_SNP<-list()
lassorst_SNP<-list()
sisrst_SNP<-list()

for(i in 1:100){
  print(i)
  #Set Seed
  set.seed(i+1000)
  
  # Generate X and Y
  sigma<-cov_block(m_G,.3,5)
  #sigma<-GenerateCliquesCovariance(10,10,0.8)
  binprob<-runif(m_G)
  x<-sim_X_cate(m_X,m_W,m_G,sigma,n,binprob)
  #x<-sim_X(m_X,m_W,m_G,sigma,n)
  #beta<-sim_beta(m_X=0,m_W=1,m_G,main_nonzero=0.05,inter_nonzero=0.05,both_nonzero=0.1,bit=T,heir=T)
  beta<-sim_beta_const(m_X,m_W=1,m_G,main_nonzero=0.1,inter_nonzero=0.1,both_nonzero=0.01,const=c(3,5),heir=TRUE)
  y0<-x%*%beta
  
  
  noise<-rnorm(n,sd=1)
  SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
  y<-y0+SNRmtl*noise  
  colnames(x)<-c(1:dim(x)[2])
  truth<-which(beta!=0)
  
  simu<-data.frame(X=x,Y=y)
  L<-lm(y~x[,c(1:(m_X+m_W))])
  y.res<-L$residuals
  simu$Y<-y.res
  simu<-simu[,-c(1:(m_X+m_W))]
  
  ### Trees
  model <- randomForest(Y~.,   data=simu)
  #print(model) # view results 
  #importance(model)
  treerst_SNP[[i]]<-order(importance(model),decreasing = T)[1:length(truth)]
  
  ### BMA
  bicfit<-bicreg(x[,-c(1:(m_X+m_W))],y.res,strict = T)
  bicrst_SNP[[i]]<-bicfit$namesx[order(bicfit$probne0,decreasing = T)][1:length(truth)]
  bicrst_SNP[[i]]<-sapply(bicrst_SNP[[i]],function(x) strsplit(x,"X")[[1]][2])
  bicrst_SNP[[i]]<-as.integer(bicrst_SNP[[i]])
  
  ### Stepwise
  a<-regsubsets(x=x,y=y,method="forward",nvmax = 3*length(truth),force.in = c(1:(m_X+m_W)))
  steprst_SNP[[i]]<-a$vorder[1:(3*length(truth))]
  #steprst[[i]]<-steprst[order(steprst[[i]])][[1]]
  
  ### Group Lasso
  
  solg<-FASTA(x,y,f, gradf, g, proxg, x0, tau1, max_iters = 300, w = 10, 
              backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
              eps_n = 1e-15,m_X,m_W,m_G,m_G,lamb_opt_glasso,lamb_opt2_glasso,restart=TRUE)
  glassorst_SNP[[i]]<-which(solg$x!=0)
  
  ### Regular Lasso
  sol<-FASTA(x,y,f, gradf, glasso, proxglasso, x0, tau1, max_iters = 300, w = 10, 
             backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
             eps_n = 1e-15,m_X,m_W,m_G,m_G,lamb_opt_lasso,lamb_opt2_lasso,restart=TRUE)
  lassorst_SNP[[i]]<-which(sol$x!=0)
  
  ### SIS
  #model1<-SIS(x,y,family = "gaussian", penalty = "lasso", tune="bic")
  model2<-SIS(x,y,family = "gaussian", penalty = "lasso", tune="bic",varISIS = "aggr")
  sisrst_SNP[[i]]<-model2$ix
}


#save(truth,"C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//30//truth.RData")
save(bicrst_SNP,file="C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//50//bicrst_SNP.RData")
save(steprst_SNP,file="C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//50//steprst_SNP.RData")
save(glassorst_SNP,file="C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//50//glassorst_SNP.RData")
save(lassorst_SNP,file="C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//50//lassorst_SNP.RData")
save(sisrst_SNP,file="C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//50//sisrst_SNP.RData")


####################### m_G=100

n=100
m_X<-5
m_W<-1
m_G<-100
m_I<-m_G
SNR<-10
tau1<-1


#Set Seed
set.seed(1000)

# Generate X and Y
sigma<-cov_block(m_G,.3,20)
#sigma<-GenerateCliquesCovariance(10,10,0.8)
#binprob<-runif(m_G)
x<-sim_X(m_X,m_W,m_G,sigma,n)
#beta<-sim_beta(m_X=0,m_W=1,m_G,main_nonzero=0.05,inter_nonzero=0.05,both_nonzero=0.1,bit=T,heir=T)
beta<-sim_beta_const(m_X,m_W=1,m_G,main_nonzero=0.1,inter_nonzero=0.1,both_nonzero=0.01,const=c(3,5),heir=TRUE)
y0<-x%*%beta


noise<-rnorm(n,sd=1)
SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
y<-y0+SNRmtl*noise  
colnames(x)<-c(1:dim(x)[2])
truth<-which(beta!=0)

x0<-rep(0,dim(x)[2])
#### Cross Validation finding best lambda for Group Lasso
lamb_candidate<-c(1,1.5,2,2.5,3,4)
lamb_candidate2<-c(0.5,1,1.5,2)
sol_cv_SNP<-opt_lambda(x,y,f, gradf, g, proxg, x0, tau1, max_iters = 100, w = 10, 
                   backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                   eps_n = 1e-15,m_X,m_W,m_G,m_I,K=10,n=100, lamb_candidate, lamb_candidate2,restart=TRUE)
lamb_loc<-which(sol_cv_SNP$mean+sol_cv_SNP$var == min(sol_cv_SNP$mean+sol_cv_SNP$var), arr.ind = TRUE)
lamb_opt_glasso<-lamb_candidate[lamb_loc[1]]
lamb_opt2_glasso<-lamb_candidate2[lamb_loc[2]]
save(sol_cv_SNP,file="C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//sol_cv_glasso_SNP.RData")

#### Cross Validation finding best lambda for General Lasso
lamb_candidate<-c(15,25,30,35,40,45,50)
lamb_candidate2<-c(1,3,5,7)
sol_cv_SNP<-opt_lambda(x,y,f, gradf, glasso, proxglasso, x0, tau1, max_iters = 100, w = 10, 
                   backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                   eps_n = 1e-15,m_X,m_W,m_G,m_I,K=10,n=100, lamb_candidate, lamb_candidate2,restart=TRUE)
lamb_loc<-which(sol_cv_SNP$mean+sol_cv_SNP$var == min(sol_cv_SNP$mean+sol_cv_SNP$var), arr.ind = TRUE)
lamb_opt_lasso<-lamb_candidate[lamb_loc[1]]
lamb_opt2_lasso<-lamb_candidate2[lamb_loc[2]]
save(sol_cv_SNP,file="C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//sol_cv_lasso_SNP.RData")


#### Iteration
treerst_SNP<-list()
bicrst_SNP<-list()
steprst_SNP<-list()
glassorst_SNP<-list()
lassorst_SNP<-list()
sisrst_SNP<-list()

for(portion in c(0.05,0.1,0.15,0.2)){
  for(i in 1:100){
    print(i)
    #Set Seed
    set.seed(i+1000)
    
    # Generate X and Y
    sigma<-cov_block(m_G,.3,20)
    #sigma<-GenerateCliquesCovariance(10,10,0.8)
    #binprob<-runif(m_G)
    x<-sim_X(m_X,m_W,m_G,sigma,n)
    #beta<-sim_beta(m_X=0,m_W=1,m_G,main_nonzero=0.05,inter_nonzero=0.05,both_nonzero=0.1,bit=T,heir=T)
    beta<-sim_beta_const(m_X,m_W=1,m_G,main_nonzero=portion,inter_nonzero=portion,both_nonzero=0.01,const=c(3,5),heir=TRUE)
    y0<-x%*%beta
    
    
    noise<-rnorm(n,sd=1)
    SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
    y<-y0+SNRmtl*noise  
    colnames(x)<-c(1:dim(x)[2])
    truth<-which(beta!=0)
    
    simu<-data.frame(X=x,Y=y)
    L<-lm(y~x[,c(1:(m_X+m_W))])
    y.res<-L$residuals
    simu$Y<-y.res
    simu<-simu[,-c(1:(m_X+m_W))]
    
    ### Trees
    model <- randomForest(Y~.,   data=simu)
    #print(model) # view results 
    #importance(model)
    treerst_SNP[[i]]<-order(importance(model),decreasing = T)[1:length(truth)]
    
    ### BMA
    bicfit<-bicreg(x[,-c(1:(m_X+m_W))],y.res,strict = T)
    bicrst_SNP[[i]]<-bicfit$namesx[order(bicfit$probne0,decreasing = T)][1:length(truth)]
    bicrst_SNP[[i]]<-sapply(bicrst_SNP[[i]],function(x) strsplit(x,"X")[[1]][2])
    bicrst_SNP[[i]]<-as.integer(bicrst_SNP[[i]])
    
    ### Stepwise
    a<-regsubsets(x=x,y=y,method="forward",nvmax = 3*length(truth),force.in = c(1:(m_X+m_W)))
    steprst_SNP[[i]]<-a$vorder[1:(length(3*truth))]
    #steprst[[i]]<-steprst[order(steprst[[i]])][[1]]
    
    ### Group Lasso
    
    solg<-FASTA(x,y,f, gradf, g, proxg, x0, tau1, max_iters = 300, w = 10, 
                backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                eps_n = 1e-15,m_X,m_W,m_G,m_G,lamb_opt_glasso,lamb_opt2_glasso,restart=TRUE)
    glassorst_SNP[[i]]<-which(solg$x!=0)
    
    ### Regular Lasso
    sol<-FASTA(x,y,f, gradf, glasso, proxglasso, x0, tau1, max_iters = 300, w = 10, 
               backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
               eps_n = 1e-15,m_X,m_W,m_G,m_G,lamb_opt_lasso,lamb_opt2_lasso,restart=TRUE)
    lassorst_SNP[[i]]<-which(sol$x!=0)
    
    ### SIS
    #model1<-SIS(x,y,family = "gaussian", penalty = "lasso", tune="bic")
    model2<-SIS(x,y,family = "gaussian", penalty = "lasso", tune="bic",varISIS = "aggr")
    sisrst_SNP[[i]]<-model2$ix
  }
  
  
  #save(truth,"C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//30//truth.RData")
  save(bicrst_SNP,file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//bicrst",portion,"_SNP.RData"))
  save(steprst_SNP,file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//steprst",portion,"_SNP.RData"))
  save(glassorst_SNP,file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//glassorst",portion,"_SNP.RData"))
  save(lassorst_SNP,file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//lassorst",portion,"_SNP.RData"))
  save(sisrst_SNP,file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//sisrst",portion,"_SNP.RData"))
  
}

