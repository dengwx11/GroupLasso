#!/usr/bin/R

suppressMessages( library(magrittr) )
suppressMessages( library(argparser) )

p <- arg_parser("Simulation on Group lasso compared with multiple other methods") %>%
  # add_argument(
  #   arg = "--seed",
  #   help = "Random seed",
  #   type = "integer"
  # ) %>%
  add_argument(
    arg = "--samples",
    help = "Number of samples",
    default = 100,
    type = "integer",
    short = "-N"
  ) %>%
  add_argument(
    arg = "--baseline_dim",
    help = "Number of baseline covariates",
    default = 2,
    type = "integer",
    short = "-m_X"
  ) %>%
  add_argument(
    arg = "--treatment_dim",
    help = "Number of treatment covariates",
    default = 1,
    type = "integer",
    short = "-m_W"
  ) %>%
  add_argument(
    arg = "--gene_dim",
    help = "Number of gene covariates",
    default = 10,
    type = "integer",
    short = "-m_G"
  ) %>%
  add_argument(
    arg = "--out",
    help = "Output file",
    default = ""
  ) %>%
  add_argument(
    arg = "--main_nonzero",
    help = "Proportion of nonzero among all nonzero main effects",
    default = .1,
    type="numeric",
    short = "-MN"
  ) %>%
  add_argument(
    arg = "--interaction_nonzero",
    help = "Proportion of nonzero among all nonzero interaction effects",
    default = .1,
    type="numeric",
    short = "-IN"
  ) %>%
  add_argument(
    arg = "--both_nonzero",
    help = "Proportion of both nonzero among all nonzero interaction effects when hierachical relationship doesn't exist",
    default = .1,
    type="numeric",
    short = "-BN"
  ) %>%
  add_argument(
    arg = "--SNR",
    help = "Signal Noise Ratio",
    type="numeric",
    default = 10
  ) %>%
  # add_argument(
  #   arg = "--Cov_type",
  #   help = "Data type of gene covariates",
  #   default = "SNP",
  #   short="-type"
  # ) %>%
  add_argument(
    arg = "--Hierarchy",
    help = "Whether hierarchical relationship exist",
    default = TRUE,
    type = "logical",
    short="-H"
  ) 
  argv <- parse_args(p)

# if (length(argv) < 3) {
#   stop("Must have at least 12 covariates")
# }
print(argv)

#working directory, need to be changed
dir<-"/home/fas/zhao/wd262/project/predictive" 
#dir<-"/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu"

#need to import library(dr) and library(MASS)
source(sprintf("%s/cv3.R",dir))
source(sprintf("%s/func3.R",dir))
source(sprintf("%s/opt3.R",dir))
source(sprintf("%s/sim4.R",dir))

library(rpart)
library(randomForest)
library(BMA)
library(leaps)
#library(BoomSpikeSlab)
library(SIS)

n=argv$samples
m_X=argv$baseline_dim
m_W=argv$treatment_dim
m_G=argv$gene_dim
m_I=m_G
SNR=argv$SNR
main_nonzero=argv$main_nonzero
inter_nonzero=argv$interaction_nonzero
both_nonzero=argv$both_nonzero
tau1<-1

print(m_X+m_W+m_G+m_I)

print(c(main_nonzero,inter_nonzero))

#Set Seed
set.seed(1000)

# Generate X and Y
sigma<-cov_block(m_G,.3,20)
#sigma<-GenerateCliquesCovariance(10,10,0.8)
#binprob<-runif(m_G)
x<-sim_X(m_X,m_W,m_G,sigma,n)
#beta<-sim_beta(m_X=0,m_W=1,m_G,main_nonzero=0.05,inter_nonzero=0.05,both_nonzero=0.1,bit=T,heir=T)
beta<-sim_beta_const(m_X,m_W=1,m_G,main_nonzero=main_nonzero,inter_nonzero=inter_nonzero,both_nonzero=both_nonzero,const=c(3,5),heir=TRUE)
y0<-x%*%beta


noise<-rnorm(n,sd=1)
SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
y<-y0+SNRmtl*noise  
colnames(x)<-c(1:dim(x)[2])
truth<-which(beta!=0)

x0<-rep(0,dim(x)[2])
#### Cross Validation finding best lambda for Group Lasso
#lamb_candidate<-c(1,1.5,2,2.5,3,4)
#lamb_candidate2<-c(0.5,1,1.5,2)

lamb_candidate<-c(1,1.5)
lamb_candidate2<-c(0.5,1)
sol_cv<-opt_lambda(x,y,f, gradf, g, proxg, x0, tau1, max_iters = 100, w = 10, 
                   backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                   eps_n = 1e-15,m_X,m_W,m_G,m_I,K=10,n=100, lamb_candidate, lamb_candidate2,restart=TRUE)
lamb_loc<-which(sol_cv$mean+sol_cv$var == min(sol_cv$mean+sol_cv$var), arr.ind = TRUE)
lamb_opt_glasso<-lamb_candidate[lamb_loc[1]]
lamb_opt2_glasso<-lamb_candidate2[lamb_loc[2]]
save(sol_cv,file=paste0("/home/fas/zhao/wd262/project/predictive/sol_cv_glasso_",SNR,".RData"))

#### Cross Validation finding best lambda for General Lasso
lamb_candidate<-c(15,25,30,35,40,45,50)
lamb_candidate2<-c(1,3,5,7)

sol_cv<-opt_lambda(x,y,f, gradf, glasso, proxglasso, x0, tau1, max_iters = 100, w = 10, 
                   backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                   eps_n = 1e-15,m_X,m_W,m_G,m_I,K=10,n=100, lamb_candidate, lamb_candidate2,restart=TRUE)
lamb_loc<-which(sol_cv$mean+sol_cv$var == min(sol_cv$mean+sol_cv$var), arr.ind = TRUE)
lamb_opt_lasso<-lamb_candidate[lamb_loc[1]]
lamb_opt2_lasso<-lamb_candidate2[lamb_loc[2]]
save(sol_cv,file=paste0("/home/fas/zhao/wd262/project/predictive/sol_cv_lasso_",SNR,".RData"))


#### Iteration
treerst<-list()
bicrst<-list()
steprst<-list()
glassorst<-list()
lassorst<-list()
sisrst<-list()


  for(i in 1:2){
    print(i)
    #Set Seed
    set.seed(i+1000)
    
    # Generate X and Y
    sigma<-cov_block(m_G,.3,20)
    #sigma<-GenerateCliquesCovariance(10,10,0.8)
    #binprob<-runif(m_G)
    x<-sim_X(m_X,m_W,m_G,sigma,n)
    #beta<-sim_beta(m_X=0,m_W=1,m_G,main_nonzero=0.05,inter_nonzero=0.05,both_nonzero=0.1,bit=T,heir=T)
    beta<-sim_beta_const(m_X,m_W=1,m_G,main_nonzero=main_nonzero,inter_nonzero=inter_nonzero,both_nonzero=both_nonzero,const=c(3,5),heir=TRUE)
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
    treerst[[i]]<-order(importance(model),decreasing = T)[1:length(truth)]
    
    ### BMA
    bicfit<-bicreg(x[,-c(1:(m_X+m_W))],y.res,strict = T)
    bicrst[[i]]<-bicfit$namesx[order(bicfit$probne0,decreasing = T)][1:length(truth)]
    bicrst[[i]]<-sapply(bicrst[[i]],function(x) strsplit(x,"X")[[1]][2])
    bicrst[[i]]<-as.integer(bicrst[[i]])
    
    ### Stepwise
    a<-regsubsets(x=x,y=y,method="forward",nvmax = 3*length(truth),force.in = c(1:(m_X+m_W)))
    steprst[[i]]<-a$vorder[1:(length(truth))]
    steprst[[i]]<-steprst[order(steprst[[i]])][[1]]
    
    ### Group Lasso
    
    solg<-FASTA(x,y,f, gradf, g, proxg, x0, tau1, max_iters = 300, w = 10, 
                backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                eps_n = 1e-15,m_X,m_W,m_G,m_G,lamb_opt_glasso,lamb_opt2_glasso,restart=TRUE)
    glassorst[[i]]<-which(solg$x!=0)
    
    ### Regular Lasso
    sol<-FASTA(x,y,f, gradf, glasso, proxglasso, x0, tau1, max_iters = 300, w = 10, 
               backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
               eps_n = 1e-15,m_X,m_W,m_G,m_G,lamb_opt_lasso,lamb_opt2_lasso,restart=TRUE)
    lassorst[[i]]<-which(sol$x!=0)
    
    ### SIS
    #model1<-SIS(x,y,family = "gaussian", penalty = "lasso", tune="bic")
    model2<-SIS(x,y,family = "gaussian", penalty = "lasso", tune="bic",varISIS = "aggr")
    sisrst[[i]]<-model2$ix
  }
  
  
  #save(truth,"C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//30//truth.RData")
  save(bicrst,file=paste0("/home/fas/zhao/wd262/project/predictive/bicrst_",SNR,".RData"))
  save(steprst,file=paste0("/home/fas/zhao/wd262/project/predictive/steprst_",SNR,".RData"))
  save(glassorst,file=paste0("/home/fas/zhao/wd262/project/predictive/glassorst_",SNR,".RData"))
  save(lassorst,file=paste0("/home/fas/zhao/wd262/project/predictive/lassorst_",SNR,".RData"))
  save(sisrst,file=paste0("/home/fas/zhao/wd262/project/predictive/sisrst_",SNR,".RData"))
  

