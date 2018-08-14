# Data Analysis

##### p=100, nonzero proportion


n=100
m_X<-5
m_W<-1
m_G<-100
m_I<-m_G
SNR<-10
tau1<-1

k=1
size<-list()
aic<-list()
L1<-list()
L2<-list()
TP.all<-list()
TN.all<-list()
TP.prog<-list()
TP.pred<-list()
TN.prog<-list()
TN.pred<-list()

for(portion in c(0.1)){

    
 load(file=paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/100/glassorst",portion,".RData"))
 load(file=paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/100/lassorst",portion,".RData"))
 load(file=paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/100/bicrst",portion,".RData"))
 load(file=paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/100/steprst",portion,".RData"))
 load(file=paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/100/sisrst",portion,".RData"))

size[[k]]<-list()
 size[[k]][[1]]<-sapply(glassorst,length)+m_X+m_W
  size[[k]][[2]]<-sapply(lassorst,length)+m_X+m_W
   size[[k]][[3]]<-sapply(bicrst,length)
    size[[k]][[4]]<-sapply(steprst,length)
     size[[k]][[5]]<-sapply(sisrst,length)



aic[[k]]<-matrix(0,ncol=5,nrow=100)
L1[[k]]<-matrix(0,ncol=5,nrow=100)
L2[[k]]<-matrix(0,ncol=5,nrow=100)
TP.all[[k]]<-matrix(0,ncol=5,nrow=100)
TN.all[[k]]<-matrix(0,ncol=5,nrow=100)
TP.prog[[k]]<-matrix(0,ncol=5,nrow=100)
TP.pred[[k]]<-matrix(0,ncol=5,nrow=100)
TN.prog[[k]]<-matrix(0,ncol=5,nrow=100)
TN.pred[[k]]<-matrix(0,ncol=5,nrow=100)

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
    beta0<-split_beta(beta,m_X,m_W,m_G,m_I)
    
    
    noise<-rnorm(n,sd=1)
    SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
    y<-y0+SNRmtl*noise  
    colnames(x)<-c(1:dim(x)[2])
    truth<-which(beta!=0)
   
   ## Group Lasso
   instance.glasso<-c(c(1:(m_X+m_W)),glassorst[[i]])
   instance.glasso<-instance.glasso[order(instance.glasso, decreasing = FALSE)]
    
    AIC.glasso<-AIC(lm(y~-1+x[,instance.glasso]))

    beta.glasso<-lm(y~-1+x[,instance.glasso])$coef
    temp<-rep(0,length(beta))
    temp[instance.glasso]<-beta.glasso
    beta.glasso<-temp
    L1.glasso<-norm(as.matrix(beta.glasso-beta),"1")
    L2.glasso<-norm(as.matrix(beta.glasso-beta),"2")
    
    TP.all.glasso<-length(which(beta.glasso*beta!=0))/length(instance.glasso)
    TN.all.glasso<-length(which(beta.glasso+beta==0))/length(which(beta.glasso==0))
    beta.glasso0<-split_beta(beta.glasso,m_X,m_W,m_G,m_I)
    TP.prog.glasso<-length(which(beta.glasso0$G*beta0$G!=0))/length(which(beta.glasso0$G!=0))
    TP.pred.glasso<-length(which(beta.glasso0$I*beta0$I!=0))/length(which(beta.glasso0$I!=0))
    TN.prog.glasso<-length(which(beta.glasso0$G+beta0$G==0))/length(which(beta.glasso0$G==0))
    TN.pred.glasso<-length(which(beta.glasso0$I+beta0$I==0))/length(which(beta.glasso0$I==0))
    
    aic[[k]][i,1]<-AIC.glasso
    L1[[k]][i,1]<-L1.glasso
    L2[[k]][i,1]<-L2.glasso
    TP.all[[k]][i,1]<-TP.all.glasso
    TN.all[[k]][i,1]<-TN.all.glasso
    TP.prog[[k]][i,1]<-TP.prog.glasso
    TN.prog[[k]][i,1]<-TN.prog.glasso
    TP.pred[[k]][i,1]<-TP.pred.glasso
    TN.pred[[k]][i,1]<-TN.pred.glasso
    

    ## Lasso
   instance.lasso<-c(c(1:(m_X+m_W)),lassorst[[i]])
   instance.lasso<-instance.lasso[order(instance.lasso, decreasing = FALSE)]
    
    AIC.lasso<-AIC(lm(y~-1+x[,instance.lasso]))

    beta.lasso<-lm(y~-1+x[,instance.lasso])$coef
    temp<-rep(0,length(beta))
    temp[instance.lasso]<-beta.lasso
    beta.lasso<-temp
    L1.lasso<-norm(as.matrix(beta.lasso-beta),"1")
    L2.lasso<-norm(as.matrix(beta.lasso-beta),"2")
    
    TP.all.lasso<-length(which(beta.lasso*beta!=0))/length(instance.lasso)
    TN.all.lasso<-length(which(beta.lasso+beta==0))/length(which(beta.lasso==0))
    beta.lasso0<-split_beta(beta.lasso,m_X,m_W,m_G,m_I)
    TP.prog.lasso<-length(which(beta.lasso0$G*beta0$G!=0))/length(which(beta.lasso0$G!=0))
    TP.pred.lasso<-length(which(beta.lasso0$I*beta0$I!=0))/length(which(beta.lasso0$I!=0))
    TN.prog.lasso<-length(which(beta.lasso0$G+beta0$G==0))/length(which(beta.lasso0$G==0))
    TN.pred.lasso<-length(which(beta.lasso0$I+beta0$I==0))/length(which(beta.lasso0$I==0))


    aic[[k]][i,2]<-AIC.lasso
    L1[[k]][i,2]<-L1.lasso
    L2[[k]][i,2]<-L2.lasso
    TP.all[[k]][i,2]<-TP.all.lasso
    TN.all[[k]][i,2]<-TN.all.lasso
    TP.prog[[k]][i,2]<-TP.prog.lasso
    TN.prog[[k]][i,2]<-TN.prog.lasso
    TP.pred[[k]][i,2]<-TP.pred.lasso
    TN.pred[[k]][i,2]<-TN.pred.lasso

    

    ## BMA
    instance.bic<-bicrst[[i]][which(is.na(bicrst[[i]])==F)]
    instance.bic<-c(c(1:(m_X+m_W)),instance.bic)

    AIC.bic<-AIC(lm(y~-1+x[,instance.bic]))

    beta.bic<-lm(y~-1+x[,instance.bic])$coef
    temp<-rep(0,length(beta))
    temp[instance.bic]<-beta.bic
    beta.bic<-temp
    L1.bic<-norm(as.matrix(beta.bic-beta),"1")
    L2.bic<-norm(as.matrix(beta.bic-beta),"2")
    
    TP.all.bic<-length(which(beta.bic*beta!=0))/length(instance.bic)
    TN.all.bic<-length(which(beta.bic+beta==0))/length(which(beta.bic==0))
    beta.bic0<-split_beta(beta.bic,m_X,m_W,m_G,m_I)
    TP.prog.bic<-length(which(beta.bic0$G*beta0$G!=0))/length(which(beta.bic0$G!=0))
    TP.pred.bic<-length(which(beta.bic0$I*beta0$I!=0))/length(which(beta.bic0$I!=0))
    TN.prog.bic<-length(which(beta.bic0$G+beta0$G==0))/length(which(beta.bic0$G==0))
    TN.pred.bic<-length(which(beta.bic0$I+beta0$I==0))/length(which(beta.bic0$I==0))

    aic[[k]][i,3]<-AIC.bic
    L1[[k]][i,3]<-L1.bic
    L2[[k]][i,3]<-L2.bic
    TP.all[[k]][i,3]<-TP.all.bic
    TN.all[[k]][i,3]<-TN.all.bic
    TP.prog[[k]][i,3]<-TP.prog.bic
    TN.prog[[k]][i,3]<-TN.prog.bic
    TP.pred[[k]][i,3]<-TP.pred.bic
    TN.pred[[k]][i,3]<-TN.pred.bic
    
    

    ## Stepwise
    instance.step<-steprst[[i]]

    AIC.step<-AIC(lm(y~-1+x[,instance.step]))

    beta.step<-lm(y~-1+x[,instance.step])$coef
    temp<-rep(0,length(beta))
    temp[instance.step]<-beta.step
    beta.step<-temp
    L1.step<-norm(as.matrix(beta.step-beta),"1")
    L2.step<-norm(as.matrix(beta.step-beta),"2")
    
    TP.all.step<-length(which(beta.step*beta!=0))/length(instance.step)
    TN.all.step<-length(which(beta.step+beta==0))/length(which(beta.step==0))
    beta.step0<-split_beta(beta.step,m_X,m_W,m_G,m_I)
    TP.prog.step<-length(which(beta.step0$G*beta0$G!=0))/length(which(beta.step0$G!=0))
    TP.pred.step<-length(which(beta.step0$I*beta0$I!=0))/length(which(beta.step0$I!=0))
    TN.prog.step<-length(which(beta.step0$G+beta0$G==0))/length(which(beta.step0$G==0))
    TN.pred.step<-length(which(beta.step0$I+beta0$I==0))/length(which(beta.step0$I==0))

    aic[[k]][i,4]<-AIC.step
    L1[[k]][i,4]<-L1.step
    L2[[k]][i,4]<-L2.step
    TP.all[[k]][i,4]<-TP.all.step
    TN.all[[k]][i,4]<-TN.all.step
    TP.prog[[k]][i,4]<-TP.prog.step
    TN.prog[[k]][i,4]<-TN.prog.step
    TP.pred[[k]][i,4]<-TP.pred.step
    TN.pred[[k]][i,4]<-TN.pred.step

    ## SIS
    instance.sis<-union(c(1:(m_X+m_W)),sisrst[[i]])

    AIC.sis<-AIC(lm(y~-1+x[,instance.sis]))

    beta.sis<-lm(y~-1+x[,instance.sis])$coef
    temp<-rep(0,length(beta))
    temp[instance.sis]<-beta.sis
    beta.sis<-temp
    L1.sis<-norm(as.matrix(beta.sis-beta),"1")
    L2.sis<-norm(as.matrix(beta.sis-beta),"2")
    
    TP.all.sis<-length(which(beta.sis*beta!=0))/length(instance.sis)
    TN.all.sis<-length(which(beta.sis+beta==0))/length(which(beta.sis==0))
    beta.sis0<-split_beta(beta.sis,m_X,m_W,m_G,m_I)
    TP.prog.sis<-length(which(beta.sis0$G*beta0$G!=0))/length(which(beta.sis0$G!=0))
    TP.pred.sis<-length(which(beta.sis0$I*beta0$I!=0))/length(which(beta.sis0$I!=0))
    TN.prog.sis<-length(which(beta.sis0$G+beta0$G==0))/length(which(beta.sis0$G==0))
    TN.pred.sis<-length(which(beta.sis0$I+beta0$I==0))/length(which(beta.sis0$I==0))

    aic[[k]][i,5]<-AIC.sis
    L1[[k]][i,5]<-L1.sis
    L2[[k]][i,5]<-L2.sis
    TP.all[[k]][i,5]<-TP.all.sis
    TN.all[[k]][i,5]<-TN.all.sis
    TP.prog[[k]][i,5]<-TP.prog.sis
    TN.prog[[k]][i,5]<-TN.prog.sis
    TP.pred[[k]][i,5]<-TP.pred.sis
    TN.pred[[k]][i,5]<-TN.pred.sis
 }

  k=k+1
}


rst.summary<-list()
for(k in 1){
  rst.summary[[k]]<-data.frame(L2=apply(L2[[k]], 2, mean))
  rst.summary[[k]]$L2<-apply(L2[[k]], 2, mean)
  rst.summary[[k]]$L1<-apply(L1[[k]], 2, mean)
  rst.summary[[k]]$aic<-apply(aic[[k]], 2, mean)
  rst.summary[[k]]$TP.all<-apply(TP.all[[k]], 2, mean)
  rst.summary[[k]]$TN.all<-apply(TN.all[[k]], 2, mean)
  rst.summary[[k]]$TP.prog<-apply(TP.prog[[k]], 2, mean)
  rst.summary[[k]]$TP.pred<-apply(TP.pred[[k]], 2, mean)
  rst.summary[[k]]$TN.prog<-apply(TN.prog[[k]], 2, mean)
  rst.summary[[k]]$TN.pred<-apply(TN.pred[[k]], 2, mean)
  rst.summary[[k]]$size<-sapply(size[[k]], median)
  rownames(rst.summary[[k]])<-c("glasso","lasso","BMA","stepwise","SIS")
}

library(ggplot2)
library(gridExtra)

portion<-c(0.05,0.1,0.15,0.2)
for(k in 1:4){
png(paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/summary/SNR/",portion[k],".png"), height=200, width=1200)
p<-tableGrob(rst.summary[[k]])
grid.arrange(p)
dev.off()
}


######## SNR, p=100
n=100
m_X<-5
m_W<-1
m_G<-100
m_I<-m_G
SNR<-10
tau1<-1

k=1
# size<-list()
aic<-list()
L1<-list()
L2<-list()
TP.all<-list()
TN.all<-list()
TP.prog<-list()
TP.pred<-list()
TN.prog<-list()
TN.pred<-list()

SNRlist<-c(1,5,10,100)

for(SNR in SNRlist){

    
 load(file=paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/100/glassorst_",SNR,".RData"))
 load(file=paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/100/lassorst_",SNR,".RData"))

# size[[k]]<-list()
#  size[[k]][[1]]<-sapply(glassorst,length)+m_X+m_W
#   size[[k]][[2]]<-sapply(lassorst,length)+m_X+m_W
 


aic[[k]]<-matrix(0,ncol=2,nrow=100)
L1[[k]]<-matrix(0,ncol=2,nrow=100)
L2[[k]]<-matrix(0,ncol=2,nrow=100)
TP.all[[k]]<-matrix(0,ncol=2,nrow=100)
TN.all[[k]]<-matrix(0,ncol=2,nrow=100)
TP.prog[[k]]<-matrix(0,ncol=2,nrow=100)
TP.pred[[k]]<-matrix(0,ncol=2,nrow=100)
TN.prog[[k]]<-matrix(0,ncol=2,nrow=100)
TN.pred[[k]]<-matrix(0,ncol=2,nrow=100)

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
    beta0<-split_beta(beta,m_X,m_W,m_G,m_I)
    
    
    noise<-rnorm(n,sd=1)
    SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
    y<-y0+SNRmtl*noise  
    colnames(x)<-c(1:dim(x)[2])
    truth<-which(beta!=0)
   
   ## Group Lasso
   instance.glasso<-c(c(1:(m_X+m_W)),glassorst[[i]])
   instance.glasso<-instance.glasso[order(instance.glasso, decreasing = FALSE)]
    
    AIC.glasso<-AIC(lm(y~-1+x[,instance.glasso]))

    beta.glasso<-lm(y~-1+x[,instance.glasso])$coef
    temp<-rep(0,length(beta))
    temp[instance.glasso]<-beta.glasso
    beta.glasso<-temp
    L1.glasso<-norm(as.matrix(beta.glasso-beta),"1")
    L2.glasso<-norm(as.matrix(beta.glasso-beta),"2")
    
    TP.all.glasso<-length(which(beta.glasso*beta!=0))/length(instance.glasso)
    TN.all.glasso<-length(which(beta.glasso+beta==0))/length(which(beta.glasso==0))
    beta.glasso0<-split_beta(beta.glasso,m_X,m_W,m_G,m_I)
    TP.prog.glasso<-length(which(beta.glasso0$G*beta0$G!=0))/length(which(beta.glasso0$G!=0))
    TP.pred.glasso<-length(which(beta.glasso0$I*beta0$I!=0))/length(which(beta.glasso0$I!=0))
    TN.prog.glasso<-length(which(beta.glasso0$G+beta0$G==0))/length(which(beta.glasso0$G==0))
    TN.pred.glasso<-length(which(beta.glasso0$I+beta0$I==0))/length(which(beta.glasso0$I==0))
    
    aic[[k]][i,1]<-AIC.glasso
    L1[[k]][i,1]<-L1.glasso
    L2[[k]][i,1]<-L2.glasso
    TP.all[[k]][i,1]<-TP.all.glasso
    TN.all[[k]][i,1]<-TN.all.glasso
    TP.prog[[k]][i,1]<-TP.prog.glasso
    TN.prog[[k]][i,1]<-TN.prog.glasso
    TP.pred[[k]][i,1]<-TP.pred.glasso
    TN.pred[[k]][i,1]<-TN.pred.glasso
    

    ## Lasso
   instance.lasso<-c(c(1:(m_X+m_W)),lassorst[[i]])
   instance.lasso<-instance.lasso[order(instance.lasso, decreasing = FALSE)]
    
    AIC.lasso<-AIC(lm(y~-1+x[,instance.lasso]))

    beta.lasso<-lm(y~-1+x[,instance.lasso])$coef
    temp<-rep(0,length(beta))
    temp[instance.lasso]<-beta.lasso
    beta.lasso<-temp
    L1.lasso<-norm(as.matrix(beta.lasso-beta),"1")
    L2.lasso<-norm(as.matrix(beta.lasso-beta),"2")
    
    TP.all.lasso<-length(which(beta.lasso*beta!=0))/length(instance.lasso)
    TN.all.lasso<-length(which(beta.lasso+beta==0))/length(which(beta.lasso==0))
    beta.lasso0<-split_beta(beta.lasso,m_X,m_W,m_G,m_I)
    TP.prog.lasso<-length(which(beta.lasso0$G*beta0$G!=0))/length(which(beta.lasso0$G!=0))
    TP.pred.lasso<-length(which(beta.lasso0$I*beta0$I!=0))/length(which(beta.lasso0$I!=0))
    TN.prog.lasso<-length(which(beta.lasso0$G+beta0$G==0))/length(which(beta.lasso0$G==0))
    TN.pred.lasso<-length(which(beta.lasso0$I+beta0$I==0))/length(which(beta.lasso0$I==0))


    aic[[k]][i,2]<-AIC.lasso
    L1[[k]][i,2]<-L1.lasso
    L2[[k]][i,2]<-L2.lasso
    TP.all[[k]][i,2]<-TP.all.lasso
    TN.all[[k]][i,2]<-TN.all.lasso
    TP.prog[[k]][i,2]<-TP.prog.lasso
    TN.prog[[k]][i,2]<-TN.prog.lasso
    TP.pred[[k]][i,2]<-TP.pred.lasso
    TN.pred[[k]][i,2]<-TN.pred.lasso

 }

  k=k+1
}

rst.summary<-list()
for(k in 1:4){
  rst.summary[[k]]<-data.frame(L2=apply(L2[[k]], 2, mean))
  rst.summary[[k]]$L2<-apply(L2[[k]], 2, mean)
  rst.summary[[k]]$L1<-apply(L1[[k]], 2, mean)
  rst.summary[[k]]$aic<-apply(aic[[k]], 2, mean)
  rst.summary[[k]]$TP.all<-apply(TP.all[[k]], 2, mean)
  rst.summary[[k]]$TN.all<-apply(TN.all[[k]], 2, mean)
  rst.summary[[k]]$TP.prog<-apply(TP.prog[[k]], 2, mean)
  rst.summary[[k]]$TP.pred<-apply(TP.pred[[k]], 2, mean)
  rst.summary[[k]]$TN.prog<-apply(TN.prog[[k]], 2, mean)
  rst.summary[[k]]$TN.pred<-apply(TN.pred[[k]], 2, mean)
  rownames(rst.summary[[k]])<-c("glasso","lasso")
}



for(k in 1:4){
png(paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/summary/SNR/",SNRlist[k],".png"), height=200, width=1200)
p<-tableGrob(rst.summary[[k]])
grid.arrange(p)
dev.off()
}


### SNP

