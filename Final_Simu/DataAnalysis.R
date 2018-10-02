# Data Analysis
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggridges)
##### p=100, nonzero proportion

setwd("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu")
source('opt3.R')
source("func3.R")
source("sim4.R")
source("cv3.R")

n=100
m_X<-5
m_W<-1
m_G<-100
m_I<-m_G
SNR<-10
tau1<-1

k=1

SSE<-list()
L1<-list()
L2<-list()
TP.all<-list()
FN.all<-list()
TP.prog<-list()
TP.pred<-list()
FN.prog<-list()
FN.pred<-list()
num.pred<-list()

n.method<-6

for(portion in c(0.05,0.1,0.15,0.2)){
  
  
  load(file=paste0("./100/glassorst",portion,".RData"))
  load(file=paste0("./100/lassorst",portion,".RData"))
  load(file=paste0("./100/treerst",portion,".RData"))
  load(file=paste0("./100/steprst",portion,".RData"))
  load(file=paste0("./100/sisrst",portion,".RData"))
  load(file=paste0("./100/bicrst",portion,".RData"))
  

  
  
  
  SSE[[k]]<-matrix(0,ncol=n.method,nrow=100)
  L1[[k]]<-matrix(0,ncol=n.method,nrow=100)
  L2[[k]]<-matrix(0,ncol=n.method,nrow=100)
  TP.all[[k]]<-matrix(0,ncol=n.method,nrow=100)
  FN.all[[k]]<-matrix(0,ncol=n.method,nrow=100)
  TP.prog[[k]]<-matrix(0,ncol=n.method,nrow=100)
  TP.pred[[k]]<-matrix(0,ncol=n.method,nrow=100)
  FN.prog[[k]]<-matrix(0,ncol=n.method,nrow=100)
  FN.pred[[k]]<-matrix(0,ncol=n.method,nrow=100)
  num.pred[[k]]<-matrix(0,ncol=n.method,nrow=100)
  
  for(i in 1:100){
    print(i)
    set.seed(i+1000)
    
    # Generate X and Y
    sigma<-cov_block(m_G,.3,20)
    #sigma<-GenerateCliquesCovariance(10,10,0.8)
    #binprob<-runif(m_G)
    x<-sim_X(m_X,m_W,m_G,sigma,n)
    #beta<-sim_beta(m_X=0,m_W=1,m_G,main_nonzero=0.05,inter_nonzero=0.05,both_nonzero=0.1,bit=T,heir=T)
    beta<-sim_beta_const(m_X,m_W=1,m_G,main_nonzero=.1,inter_nonzero=.1,both_nonzero=0.01,const=c(3,5),heir=TRUE)
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
    
    SSE.glasso<- -logLik(lm(y~-1+x[,instance.glasso]))
    
    beta.glasso<-lm(y~-1+x[,instance.glasso])$coef
    temp<-rep(0,length(beta))
    temp[instance.glasso]<-beta.glasso
    beta.glasso<-temp
    beta.glasso[which(is.na(beta.glasso)==T)]<-0
    L1.glasso<-norm(as.matrix(beta.glasso-beta),"1")
    L2.glasso<-norm(as.matrix(beta.glasso-beta),"2")
    
    
    TP.all.glasso<-(length(which(beta.glasso*beta!=0))-m_X-m_W)/(length(instance.glasso)-m_X-m_W)
    FN.all.glasso<-(length(intersect(which(beta.glasso==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.glasso0<-split_beta(beta.glasso,m_X,m_W,m_G,m_I)
    TP.prog.glasso<-length(which(beta.glasso0$G*beta0$G!=0))/length(which(beta.glasso0$G!=0))
    TP.pred.glasso<-length(which(beta.glasso0$I*beta0$I!=0))/length(which(beta.glasso0$I!=0))
    FN.prog.glasso<-(length(intersect(which(beta.glasso0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.glasso<-(length(intersect(which(beta.glasso0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.glasso<-length(which(beta.glasso0$I!=0))
    
    SSE[[k]][i,1]<-SSE.glasso
    L1[[k]][i,1]<-L1.glasso
    L2[[k]][i,1]<-L2.glasso
    TP.all[[k]][i,1]<-TP.all.glasso
    FN.all[[k]][i,1]<-FN.all.glasso
    TP.prog[[k]][i,1]<-TP.prog.glasso
    FN.prog[[k]][i,1]<-FN.prog.glasso
    TP.pred[[k]][i,1]<-TP.pred.glasso
    FN.pred[[k]][i,1]<-FN.pred.glasso
    num.pred[[k]][i,1]<-num.glasso
    
    ## Lasso
    instance.lasso<-c(c(1:(m_X+m_W)),lassorst[[i]])
    instance.lasso<-instance.lasso[order(instance.lasso, decreasing = FALSE)]
    
    SSE.lasso<- -logLik(lm(y~-1+x[,instance.lasso]))
    
    beta.lasso<-lm(y~-1+x[,instance.lasso])$coef
    temp<-rep(0,length(beta))
    temp[instance.lasso]<-beta.lasso
    beta.lasso<-temp
    beta.lasso[which(is.na(beta.lasso)==T)]<-0
    L1.lasso<-norm(as.matrix(beta.lasso-beta),"1")
    L2.lasso<-norm(as.matrix(beta.lasso-beta),"2")
    
    TP.all.lasso<-(length(which(beta.lasso*beta!=0))-m_X-m_W)/(length(instance.lasso)-m_X-m_W)
    FN.all.lasso<-(length(intersect(which(beta.lasso==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.lasso0<-split_beta(beta.lasso,m_X,m_W,m_G,m_I)
    TP.prog.lasso<-length(which(beta.lasso0$G*beta0$G!=0))/length(which(beta.lasso0$G!=0))
    TP.pred.lasso<-length(which(beta.lasso0$I*beta0$I!=0))/length(which(beta.lasso0$I!=0))
    FN.prog.lasso<-(length(intersect(which(beta.lasso0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.lasso<-(length(intersect(which(beta.lasso0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    
    
    SSE[[k]][i,2]<-SSE.lasso
    L1[[k]][i,2]<-L1.lasso
    L2[[k]][i,2]<-L2.lasso
    TP.all[[k]][i,2]<-TP.all.lasso
    FN.all[[k]][i,2]<-FN.all.lasso
    TP.prog[[k]][i,2]<-TP.prog.lasso
    FN.prog[[k]][i,2]<-FN.prog.lasso
    TP.pred[[k]][i,2]<-TP.pred.lasso
    FN.pred[[k]][i,2]<-FN.pred.lasso
    
    num.lasso<-length(which(beta.lasso0$I!=0))
    num.pred[[k]][i,2]<-num.lasso
    
    
    
    ## Stepwise
    instance.step<-steprst[[i]]
    
    SSE.step<- -logLik(lm(y~-1+x[,instance.step]))
    
    beta.step<-lm(y~-1+x[,instance.step])$coef
    temp<-rep(0,length(beta))
    temp[instance.step]<-beta.step
    beta.step<-temp
    beta.step[which(is.na(beta.step)==T)]<-0
    L1.step<-norm(as.matrix(beta.step-beta),"1")
    L2.step<-norm(as.matrix(beta.step-beta),"2")
    
    TP.all.step<-(length(which(beta.step*beta!=0))-m_X-m_W)/(length(instance.step)-m_X-m_W)
    FN.all.step<-(length(intersect(which(beta.step==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.step0<-split_beta(beta.step,m_X,m_W,m_G,m_I)
    TP.prog.step<-length(which(beta.step0$G*beta0$G!=0))/length(which(beta.step0$G!=0))
    TP.pred.step<-length(which(beta.step0$I*beta0$I!=0))/length(which(beta.step0$I!=0))
    FN.prog.step<-(length(intersect(which(beta.step0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.step<-(length(intersect(which(beta.step0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.step<-length(which(beta.step0$I!=0))
    
    SSE[[k]][i,3]<-SSE.step
    L1[[k]][i,3]<-L1.step
    L2[[k]][i,3]<-L2.step
    TP.all[[k]][i,3]<-TP.all.step
    FN.all[[k]][i,3]<-FN.all.step
    TP.prog[[k]][i,3]<-TP.prog.step
    FN.prog[[k]][i,3]<-FN.prog.step
    TP.pred[[k]][i,3]<-TP.pred.step
    FN.pred[[k]][i,3]<-FN.pred.step
    num.pred[[k]][i,3]<-num.step
    
    ## SIS
    instance.sis<-union(c(1:(m_X+m_W)),sisrst[[i]])
    
    SSE.sis<- -logLik(lm(y~-1+x[,instance.sis]))
    
    beta.sis<-lm(y~-1+x[,instance.sis])$coef
    temp<-rep(0,length(beta))
    temp[instance.sis]<-beta.sis
    beta.sis<-temp
    beta.sis[which(is.na(beta.sis)==T)]<-0
    L1.sis<-norm(as.matrix(beta.sis-beta),"1")
    L2.sis<-norm(as.matrix(beta.sis-beta),"2")
    
    TP.all.sis<-(length(which(beta.sis*beta!=0))-m_X-m_W)/(length(instance.sis)-m_X-m_W)
    FN.all.sis<-(length(intersect(which(beta.sis==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.sis0<-split_beta(beta.sis,m_X,m_W,m_G,m_I)
    TP.prog.sis<-length(which(beta.sis0$G*beta0$G!=0))/length(which(beta.sis0$G!=0))
    TP.pred.sis<-length(which(beta.sis0$I*beta0$I!=0))/length(which(beta.sis0$I!=0))
    FN.prog.sis<-(length(intersect(which(beta.sis0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.sis<-(length(intersect(which(beta.sis0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.sis<-length(which(beta.sis0$I!=0))
    
    SSE[[k]][i,4]<-SSE.sis
    L1[[k]][i,4]<-L1.sis
    L2[[k]][i,4]<-L2.sis
    TP.all[[k]][i,4]<-TP.all.sis
    FN.all[[k]][i,4]<-FN.all.sis
    TP.prog[[k]][i,4]<-TP.prog.sis
    FN.prog[[k]][i,4]<-FN.prog.sis
    TP.pred[[k]][i,4]<-TP.pred.sis
    FN.pred[[k]][i,4]<-FN.pred.sis
    num.pred[[k]][i,4]<-num.sis
    
    
    ## Random Forest
    instance.tree<-union(c(1:(m_X+m_W)),treerst[[i]])
    instance.tree<-instance.tree[order(instance.tree, decreasing = FALSE)]
    
    SSE.tree<- -logLik(lm(y~-1+x[,instance.tree]))
    
    beta.tree<-lm(y~-1+x[,instance.tree])$coef
    temp<-rep(0,length(beta))
    temp[instance.tree]<-beta.tree
    beta.tree<-temp
    beta.tree[which(is.na(beta.tree)==T)]<-0
    L1.tree<-norm(as.matrix(beta.tree-beta),"1")
    L2.tree<-norm(as.matrix(beta.tree-beta),"2")
    
    TP.all.tree<-(length(which(beta.tree*beta!=0))-m_X-m_W)/(length(instance.tree)-m_X-m_W)
    FN.all.tree<-(length(intersect(which(beta.tree==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.tree0<-split_beta(beta.tree,m_X,m_W,m_G,m_I)
    TP.prog.tree<-length(which(beta.tree0$G*beta0$G!=0))/length(which(beta.tree0$G!=0))
    TP.pred.tree<-length(which(beta.tree0$I*beta0$I!=0))/length(which(beta.tree0$I!=0))
    FN.prog.tree<-(length(intersect(which(beta.tree0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.tree<-(length(intersect(which(beta.tree0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.tree<-length(which(beta.tree0$I!=0))
    
    SSE[[k]][i,5]<-SSE.tree
    L1[[k]][i,5]<-L1.tree
    L2[[k]][i,5]<-L2.tree
    TP.all[[k]][i,5]<-TP.all.tree
    FN.all[[k]][i,5]<-FN.all.tree
    TP.prog[[k]][i,5]<-TP.prog.tree
    FN.prog[[k]][i,5]<-FN.prog.tree
    TP.pred[[k]][i,5]<-TP.pred.tree
    FN.pred[[k]][i,5]<-FN.pred.tree
    num.pred[[k]][i,5]<-num.tree

            ## BMA
        instance.bic<-bicrst[[i]][which(is.na(bicrst[[i]])==F)]
        instance.bic<-c(c(1:(m_X+m_W)),instance.bic)
        
        SSE.bic<- -logLik(lm(y~-1+x[,instance.tree]))
    
    beta.bic<-lm(y~-1+x[,instance.bic])$coef
    temp<-rep(0,length(beta))
    temp[instance.bic]<-beta.bic
    beta.bic<-temp
    beta.bic[which(is.na(beta.bic)==T)]<-0
    L1.bic<-norm(as.matrix(beta.bic-beta),"1")
    L2.bic<-norm(as.matrix(beta.bic-beta),"2")
    
    TP.all.bic<-(length(which(beta.bic*beta!=0))-m_X-m_W)/(length(instance.bic)-m_X-m_W)
    FN.all.bic<-(length(intersect(which(beta.bic==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.bic0<-split_beta(beta.bic,m_X,m_W,m_G,m_I)
    TP.prog.bic<-length(which(beta.bic0$G*beta0$G!=0))/length(which(beta.bic0$G!=0))
    TP.pred.bic<-length(which(beta.bic0$I*beta0$I!=0))/length(which(beta.bic0$I!=0))
    FN.prog.bic<-(length(intersect(which(beta.bic0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.bic<-(length(intersect(which(beta.bic0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.bic<-length(which(beta.bic0$I!=0))
    
    SSE[[k]][i,6]<-SSE.bic
    L1[[k]][i,6]<-L1.bic
    L2[[k]][i,6]<-L2.bic
    TP.all[[k]][i,6]<-TP.all.bic
    FN.all[[k]][i,6]<-FN.all.bic
    TP.prog[[k]][i,6]<-TP.prog.bic
    FN.prog[[k]][i,6]<-FN.prog.bic
    TP.pred[[k]][i,6]<-TP.pred.bic
    FN.pred[[k]][i,6]<-FN.pred.bic
    num.pred[[k]][i,6]<-num.bic
    
  }
  
  k=k+1
}


rst.summary<-list()
for(k in 1:4){
  rst.summary[[k]]<-data.frame(L2=apply(L2[[k]], 2, function(x) mean(x,na.rm=T)))
  rst.summary[[k]]$L2<-apply(L2[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$L1<-apply(L1[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$SSE<-apply(SSE[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.all<-apply(TP.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.all<-apply(FN.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.prog<-apply(TP.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.pred<-apply(TP.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.prog<-apply(FN.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.pred<-apply(FN.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$num.pred<-apply(num.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rownames(rst.summary[[k]])<-c("glasso","lasso","Stepwise","SIS","Random Forest","BMA")
}

portion<-c(0.05,0.1,0.15,0.2)
# Ridge Plot
TP.pred.m<-list()
q<-list()
for(k in 1:4){
  TP.pred.m[[k]]<-data.frame("ID"=seq(1,100,1),TP.pred[[k]])
  colnames(TP.pred.m[[k]])<-c("ID","glasso","lasso","Stepwise","SIS","Random Forest","BMA")
  TP.pred.m[[k]]<-melt(TP.pred.m[[k]],id.vars = "ID",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
  q[[k]]<-ggplot(TP.pred.m[[k]],aes(x=value,y=Method))+stat_density_ridges(data=TP.pred.m[[k]], quantile_lines = T,scale = 1, size = 0.25, rel_min_height = 0,fill="skyblue",alpha=.9,color="white")+ggtitle(paste0("Nonzero Interaction Proportion:",portion[k]))+xlab("Predictive Biomarkers PPV")
}
FN.pred.m<-list()
p<-list()
for(k in 1:4){
  FN.pred.m[[k]]<-data.frame("ID"=seq(1,100,1),FN.pred[[k]])
  colnames(FN.pred.m[[k]])<-c("ID","glasso","lasso","Stepwise","SIS","Random Forest","BMA")
  FN.pred.m[[k]]<-melt(FN.pred.m[[k]],id.vars = "ID",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
  p[[k]]<-ggplot(FN.pred.m[[k]],aes(x=value,y=Method))+stat_density_ridges(data=FN.pred.m[[k]], quantile_lines = T,scale = 1, size = 0.25, rel_min_height = 0,fill="skyblue",alpha=.9,color="white")+ggtitle(paste0("Nonzero Interaction Proportion:",portion[k]))+xlab("Predictive Biomarkers FNR")
}

# Bar Plot
TP.pred.bar<-sapply(rst.summary,function(x) x$TP.pred)
TP.pred.bar<-data.frame(portion,t(TP.pred.bar))
colnames(TP.pred.bar)<-c("proportion",c("glasso","lasso","Stepwise","SIS","Random Forest","BMA"))
TP.pred.bar.m<-melt(TP.pred.bar,id.vars = "proportion",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
colnames(TP.pred.bar.m)[3]<-"PPV"
ggplot(data = TP.pred.bar.m, mapping = aes(x = factor(proportion), y = PPV,fill = Method)) + geom_bar(stat = 'identity', position = 'dodge',color="blue")+xlab("Nonzero Interaction Proportion")+ggtitle("PPV vs Nonzero Interaction Proportion")+scale_fill_brewer(palette = "Spectral")

FN.pred.bar<-sapply(rst.summary,function(x) x$FN.pred)
FN.pred.bar<-data.frame(portion,t(FN.pred.bar))
colnames(FN.pred.bar)<-c("proportion",c("glasso","lasso","Stepwise","SIS","Random Forest","BMA"))
FN.pred.bar.m<-melt(FN.pred.bar,id.vars = "proportion",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
colnames(FN.pred.bar.m)[3]<-"FNR"
ggplot(data = FN.pred.bar.m, mapping = aes(x = factor(proportion), y = FNR,fill = Method)) + geom_bar(stat = 'identity', position = 'dodge',color="blue")+xlab("Nonzero Interaction Proportion")+ggtitle("FNR vs Nonzero Interaction Proportion")+scale_fill_brewer(palette = "Spectral")


num.pred.bar<-sapply(rst.summary,function(x) x$num.pred)
num.pred.bar<-rbind(num.pred.bar,c(5,10,15,20))
num.pred.bar<-data.frame(portion,t(num.pred.bar))
colnames(num.pred.bar)<-c("Proportion",c("glasso","lasso","Stepwise","SIS","Random Forest","BMA","Truth"))
num.pred.bar.m<-melt(num.pred.bar,id.vars = "Proportion",measure.vars = c(2:8),variable.name = "Method",na.rm=T)
colnames(num.pred.bar.m)[3]<-"num"
ggplot(data = num.pred.bar.m, mapping = aes(x = factor(Proportion), y = num/m_G,color=Method,group=Method,size=Method)) + geom_line()+xlab("Proportion")+ggtitle("Model Size vs Proportion")+ylab("Estimated Nonzero Predictive Biomarker Proportion")+scale_size_manual(values=c(2,0.8,0.8,.8,.8,.8,2))+theme_bw() +scale_colour_manual(breaks=levels(num.pred.bar.m$Method), values=c("red","green","orange","blue","pink","purple","black"))
theme(plot.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(),
      legend.key = element_blank(), legend.title = element_blank())

rst.summary<-list()
for(k in 1:4){
  rst.summary[[k]]<-data.frame(L2=apply(L2[[k]], 2, function(x) mean(x,na.rm=T)))
  rst.summary[[k]]$L2<-apply(L2[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$L1<-apply(L1[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$SSE<-apply(SSE[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.all<-apply(TP.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.all<-apply(FN.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.prog<-apply(TP.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.pred<-apply(TP.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.prog<-apply(FN.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.pred<-apply(FN.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$num.pred<-apply(num.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rownames(rst.summary[[k]])<-c("glasso","lasso","Stepwise","SIS","Random Forest","BMA")
  colnames(rst.summary[[k]])<-c("L2","L1","SSE","PPV all","FNR all","PPV prog","PPV pred","FNR prog","FNR pred","num pred")
}


portion<-c(0.05,0.1,0.15,0.2)
for(k in 1:4){
  png(paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/summary/proportion/",portion[k],".png"), height=200, width=700)
  p<-tableGrob(round(rst.summary[[k]],digits=3))
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
SSE<-list()
L1<-list()
L2<-list()
TP.all<-list()
FN.all<-list()
TP.prog<-list()
TP.pred<-list()
FN.prog<-list()
FN.pred<-list()
num.pred<-list()


SNRlist<-c(1,5,10,20,100)

n.method=6




for(SNR in SNRlist){
  
  
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//glassorst_",SNR,".RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//lassorst_",SNR,".RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//steprst_",SNR,".RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//sisrst_",SNR,".RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//treerst_",SNR,".RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//bicrst_",SNR,".RData"))
    
  # size[[k]]<-list()
  #  size[[k]][[1]]<-sapply(glassorst,length)+m_X+m_W
  #   size[[k]][[2]]<-sapply(lassorst,length)+m_X+m_W
  
  
  SSE[[k]]<-matrix(0,ncol=n.method,nrow=100)
  L1[[k]]<-matrix(0,ncol=n.method,nrow=100)
  L2[[k]]<-matrix(0,ncol=n.method,nrow=100)
  TP.all[[k]]<-matrix(0,ncol=n.method,nrow=100)
  FN.all[[k]]<-matrix(0,ncol=n.method,nrow=100)
  TP.prog[[k]]<-matrix(0,ncol=n.method,nrow=100)
  TP.pred[[k]]<-matrix(0,ncol=n.method,nrow=100)
  FN.prog[[k]]<-matrix(0,ncol=n.method,nrow=100)
  FN.pred[[k]]<-matrix(0,ncol=n.method,nrow=100)
  num.pred[[k]]<-matrix(0,ncol=n.method,nrow=100)

  
  for(i in 1:100){
    print(i)
    set.seed(i+1000)
    
    # Generate X and Y
    sigma<-cov_block(m_G,.3,20)
    #sigma<-GenerateCliquesCovariance(10,10,0.8)
    #binprob<-runif(m_G)
    x<-sim_X(m_X,m_W,m_G,sigma,n)
    #beta<-sim_beta(m_X=0,m_W=1,m_G,main_nonzero=0.05,inter_nonzero=0.05,both_nonzero=0.1,bit=T,heir=T)
    beta<-sim_beta_const(m_X,m_W=1,m_G,main_nonzero=.1,inter_nonzero=.1,both_nonzero=0.01,const=c(3,5),heir=TRUE)
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
    
    SSE.glasso<- -logLik(lm(y~-1+x[,instance.glasso]))
    
    beta.glasso<-lm(y~-1+x[,instance.glasso])$coef
    temp<-rep(0,length(beta))
    temp[instance.glasso]<-beta.glasso
    beta.glasso<-temp
    beta.glasso[which(is.na(beta.glasso)==T)]<-0
    L1.glasso<-norm(as.matrix(beta.glasso-beta),"1")
    L2.glasso<-norm(as.matrix(beta.glasso-beta),"2")
    
    
    TP.all.glasso<-(length(which(beta.glasso*beta!=0))-m_X-m_W)/(length(instance.glasso)-m_X-m_W)
    FN.all.glasso<-(length(intersect(which(beta.glasso==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.glasso0<-split_beta(beta.glasso,m_X,m_W,m_G,m_I)
    TP.prog.glasso<-length(which(beta.glasso0$G*beta0$G!=0))/length(which(beta.glasso0$G!=0))
    TP.pred.glasso<-length(which(beta.glasso0$I*beta0$I!=0))/length(which(beta.glasso0$I!=0))
    FN.prog.glasso<-(length(intersect(which(beta.glasso0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.glasso<-(length(intersect(which(beta.glasso0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
   
    SSE[[k]][i,1]<-SSE.glasso
    L1[[k]][i,1]<-L1.glasso
    L2[[k]][i,1]<-L2.glasso
    TP.all[[k]][i,1]<-TP.all.glasso
    FN.all[[k]][i,1]<-FN.all.glasso
    TP.prog[[k]][i,1]<-TP.prog.glasso
    FN.prog[[k]][i,1]<-FN.prog.glasso
    TP.pred[[k]][i,1]<-TP.pred.glasso
    FN.pred[[k]][i,1]<-FN.pred.glasso
    
    num.glasso<-length(which(beta.glasso0$I!=0))
    num.pred[[k]][i,1]<-num.glasso
    
    
    ## Lasso
    instance.lasso<-c(c(1:(m_X+m_W)),lassorst[[i]])
    instance.lasso<-instance.lasso[order(instance.lasso, decreasing = FALSE)]
    
    SSE.lasso<- -logLik(lm(y~-1+x[,instance.lasso]))
    
    beta.lasso<-lm(y~-1+x[,instance.lasso])$coef
    temp<-rep(0,length(beta))
    temp[instance.lasso]<-beta.lasso
    beta.lasso<-temp
    beta.lasso[which(is.na(beta.lasso)==T)]<-0
    L1.lasso<-norm(as.matrix(beta.lasso-beta),"1")
    L2.lasso<-norm(as.matrix(beta.lasso-beta),"2")
    
    TP.all.lasso<-(length(which(beta.lasso*beta!=0))-m_X-m_W)/(length(instance.lasso)-m_X-m_W)
    FN.all.lasso<-(length(intersect(which(beta.lasso==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.lasso0<-split_beta(beta.lasso,m_X,m_W,m_G,m_I)
    TP.prog.lasso<-length(which(beta.lasso0$G*beta0$G!=0))/length(which(beta.lasso0$G!=0))
    TP.pred.lasso<-length(which(beta.lasso0$I*beta0$I!=0))/length(which(beta.lasso0$I!=0))
    FN.prog.lasso<-(length(intersect(which(beta.lasso0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.lasso<-(length(intersect(which(beta.lasso0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    
    
    SSE[[k]][i,2]<-SSE.lasso
    L1[[k]][i,2]<-L1.lasso
    L2[[k]][i,2]<-L2.lasso
    TP.all[[k]][i,2]<-TP.all.lasso
    FN.all[[k]][i,2]<-FN.all.lasso
    TP.prog[[k]][i,2]<-TP.prog.lasso
    FN.prog[[k]][i,2]<-FN.prog.lasso
    TP.pred[[k]][i,2]<-TP.pred.lasso
    FN.pred[[k]][i,2]<-FN.pred.lasso
    
    num.lasso<-length(which(beta.lasso0$I!=0))
    num.pred[[k]][i,2]<-num.lasso
    
    
    ## Stepwise
    instance.step<-steprst[[i]]
    
    SSE.step<- -logLik(lm(y~-1+x[,instance.step]))
    
    beta.step<-lm(y~-1+x[,instance.step])$coef
    temp<-rep(0,length(beta))
    temp[instance.step]<-beta.step
    beta.step<-temp
    beta.step[which(is.na(beta.step)==T)]<-0
    L1.step<-norm(as.matrix(beta.step-beta),"1")
    L2.step<-norm(as.matrix(beta.step-beta),"2")
    
    TP.all.step<-(length(which(beta.step*beta!=0))-m_X-m_W)/(length(instance.step)-m_X-m_W)
    FN.all.step<-(length(intersect(which(beta.step==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.step0<-split_beta(beta.step,m_X,m_W,m_G,m_I)
    TP.prog.step<-length(which(beta.step0$G*beta0$G!=0))/length(which(beta.step0$G!=0))
    TP.pred.step<-length(which(beta.step0$I*beta0$I!=0))/length(which(beta.step0$I!=0))
    FN.prog.step<-(length(intersect(which(beta.step0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.step<-(length(intersect(which(beta.step0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    
    SSE[[k]][i,3]<-SSE.step
    L1[[k]][i,3]<-L1.step
    L2[[k]][i,3]<-L2.step
    TP.all[[k]][i,3]<-TP.all.step
    FN.all[[k]][i,3]<-FN.all.step
    TP.prog[[k]][i,3]<-TP.prog.step
    FN.prog[[k]][i,3]<-FN.prog.step
    TP.pred[[k]][i,3]<-TP.pred.step
    FN.pred[[k]][i,3]<-FN.pred.step
    
    num.step<-length(which(beta.step0$I!=0))
    num.pred[[k]][i,3]<-num.step
    
    ## SIS
    instance.sis<-union(c(1:(m_X+m_W)),sisrst[[i]])
    
    SSE.sis<- -logLik(lm(y~-1+x[,instance.sis]))
    
    beta.sis<-lm(y~-1+x[,instance.sis])$coef
    temp<-rep(0,length(beta))
    temp[instance.sis]<-beta.sis
    beta.sis<-temp
    beta.sis[which(is.na(beta.sis)==T)]<-0
    L1.sis<-norm(as.matrix(beta.sis-beta),"1")
    L2.sis<-norm(as.matrix(beta.sis-beta),"2")
    
    TP.all.sis<-(length(which(beta.sis*beta!=0))-m_X-m_W)/(length(instance.sis)-m_X-m_W)
    FN.all.sis<-(length(intersect(which(beta.sis==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.sis0<-split_beta(beta.sis,m_X,m_W,m_G,m_I)
    TP.prog.sis<-length(which(beta.sis0$G*beta0$G!=0))/length(which(beta.sis0$G!=0))
    TP.pred.sis<-length(which(beta.sis0$I*beta0$I!=0))/length(which(beta.sis0$I!=0))
    FN.prog.sis<-(length(intersect(which(beta.sis0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.sis<-(length(intersect(which(beta.sis0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    
    SSE[[k]][i,4]<-SSE.sis
    L1[[k]][i,4]<-L1.sis
    L2[[k]][i,4]<-L2.sis
    TP.all[[k]][i,4]<-TP.all.sis
    FN.all[[k]][i,4]<-FN.all.sis
    TP.prog[[k]][i,4]<-TP.prog.sis
    FN.prog[[k]][i,4]<-FN.prog.sis
    TP.pred[[k]][i,4]<-TP.pred.sis
    FN.pred[[k]][i,4]<-FN.pred.sis
    
    num.sis<-length(which(beta.sis0$I!=0))
    num.pred[[k]][i,4]<-num.sis
    
    
    ## Random Forest
    instance.tree<-union(c(1:(m_X+m_W)),treerst[[i]])
    instance.tree<-instance.tree[order(instance.tree, decreasing = FALSE)]
    
    SSE.tree<- -logLik(lm(y~-1+x[,instance.tree]))
    
    beta.tree<-lm(y~-1+x[,instance.tree])$coef
    temp<-rep(0,length(beta))
    temp[instance.tree]<-beta.tree
    beta.tree<-temp
    beta.tree[which(is.na(beta.tree)==T)]<-0
    L1.tree<-norm(as.matrix(beta.tree-beta),"1")
    L2.tree<-norm(as.matrix(beta.tree-beta),"2")
    
    TP.all.tree<-(length(which(beta.tree*beta!=0))-m_X-m_W)/(length(instance.tree)-m_X-m_W)
    FN.all.tree<-(length(intersect(which(beta.tree==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.tree0<-split_beta(beta.tree,m_X,m_W,m_G,m_I)
    TP.prog.tree<-length(which(beta.tree0$G*beta0$G!=0))/length(which(beta.tree0$G!=0))
    TP.pred.tree<-length(which(beta.tree0$I*beta0$I!=0))/length(which(beta.tree0$I!=0))
    FN.prog.tree<-(length(intersect(which(beta.tree0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.tree<-(length(intersect(which(beta.tree0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    
    SSE[[k]][i,5]<-SSE.tree
    L1[[k]][i,5]<-L1.tree
    L2[[k]][i,5]<-L2.tree
    TP.all[[k]][i,5]<-TP.all.tree
    FN.all[[k]][i,5]<-FN.all.tree
    TP.prog[[k]][i,5]<-TP.prog.tree
    FN.prog[[k]][i,5]<-FN.prog.tree
    TP.pred[[k]][i,5]<-TP.pred.tree
    FN.pred[[k]][i,5]<-FN.pred.tree
    
    num.tree<-length(which(beta.tree0$I!=0))
    num.pred[[k]][i,5]<-num.tree
    
    ## BMA
    instance.bic<-bicrst[[i]][which(is.na(bicrst[[i]])==F)]
    instance.bic<-c(c(1:(m_X+m_W)),instance.bic)
    
    SSE.bic<- -logLik(lm(y~-1+x[,instance.tree]))
    
    beta.bic<-lm(y~-1+x[,instance.bic])$coef
    temp<-rep(0,length(beta))
    temp[instance.bic]<-beta.bic
    beta.bic<-temp
    beta.bic[which(is.na(beta.bic)==T)]<-0
    L1.bic<-norm(as.matrix(beta.bic-beta),"1")
    L2.bic<-norm(as.matrix(beta.bic-beta),"2")
    
    TP.all.bic<-(length(which(beta.bic*beta!=0))-m_X-m_W)/(length(instance.bic)-m_X-m_W)
    FN.all.bic<-(length(intersect(which(beta.bic==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.bic0<-split_beta(beta.bic,m_X,m_W,m_G,m_I)
    TP.prog.bic<-length(which(beta.bic0$G*beta0$G!=0))/length(which(beta.bic0$G!=0))
    TP.pred.bic<-length(which(beta.bic0$I*beta0$I!=0))/length(which(beta.bic0$I!=0))
    FN.prog.bic<-(length(intersect(which(beta.bic0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.bic<-(length(intersect(which(beta.bic0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.bic<-length(which(beta.bic0$I!=0))
    
    SSE[[k]][i,6]<-SSE.bic
    L1[[k]][i,6]<-L1.bic
    L2[[k]][i,6]<-L2.bic
    TP.all[[k]][i,6]<-TP.all.bic
    FN.all[[k]][i,6]<-FN.all.bic
    TP.prog[[k]][i,6]<-TP.prog.bic
    FN.prog[[k]][i,6]<-FN.prog.bic
    TP.pred[[k]][i,6]<-TP.pred.bic
    FN.pred[[k]][i,6]<-FN.pred.bic
    num.pred[[k]][i,6]<-num.bic
    
  }
  
  k=k+1
}

rst.summary<-list()
for(k in 1:5){
  rst.summary[[k]]<-data.frame(L2=apply(L2[[k]], 2, function(x) mean(x,na.rm=T)))
  rst.summary[[k]]$L2<-apply(L2[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$L1<-apply(L1[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$SSE<-apply(SSE[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.all<-apply(TP.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.all<-apply(FN.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.prog<-apply(TP.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.pred<-apply(TP.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.prog<-apply(FN.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.pred<-apply(FN.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$num.pred<-apply(num.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rownames(rst.summary[[k]])<-c("glasso","lasso","Stepwise","SIS","Random Forest","BMA")
}


# Ridge Plot
TP.pred.m<-list()
q<-list()
for(k in 1:5){
  TP.pred.m[[k]]<-data.frame("ID"=seq(1,100,1),TP.pred[[k]])
  colnames(TP.pred.m[[k]])<-c("ID","glasso","lasso","Stepwise","SIS","Random Forest","BMA")
  TP.pred.m[[k]]<-melt(TP.pred.m[[k]],id.vars = "ID",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
  q[[k]]<-ggplot(TP.pred.m[[k]],aes(x=value,y=Method))+stat_density_ridges(data=TP.pred.m[[k]], quantile_lines = T,scale = 1, size = 0.25, rel_min_height = 0,fill="skyblue",alpha=.9,color="white")+ggtitle(paste0("SNR:",SNRlist[k]))+xlab("Predictive Biomarkers PPV")
}
FN.pred.m<-list()
p<-list()
for(k in 1:5){
  FN.pred.m[[k]]<-data.frame("ID"=seq(1,100,1),FN.pred[[k]])
  colnames(FN.pred.m[[k]])<-c("ID","glasso","lasso","Stepwise","SIS","Random Forest","BMA")
  FN.pred.m[[k]]<-melt(FN.pred.m[[k]],id.vars = "ID",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
  p[[k]]<-ggplot(FN.pred.m[[k]],aes(x=value,y=Method))+stat_density_ridges(data=FN.pred.m[[k]], quantile_lines = T,scale = 1, size = 0.25, rel_min_height = 0,fill="skyblue",alpha=.9,color="white")+ggtitle(paste0("SNR:",SNRlist[k]))+xlab("Predictive Biomarkers FNR")
}

# Bar Plot
TP.pred.bar<-sapply(rst.summary,function(x) x$TP.pred)
TP.pred.bar<-data.frame(SNRlist,t(TP.pred.bar))
colnames(TP.pred.bar)<-c("SNR",c("glasso","lasso","Stepwise","SIS","Random Forest","BMA"))
TP.pred.bar.m<-melt(TP.pred.bar,id.vars = "SNR",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
colnames(TP.pred.bar.m)[3]<-"PPV"
ggplot(data = TP.pred.bar.m, mapping = aes(x = factor(SNR), y = PPV,fill = Method)) + geom_bar(stat = 'identity', position = 'dodge',color="blue")+xlab("SNR")+ggtitle("PPV vs SNR")+scale_fill_brewer(palette = "Spectral")

FN.pred.bar<-sapply(rst.summary,function(x) x$FN.pred)FN.pred.bar<-data.frame(SNRlist,t(FN.pred.bar))
colnames(FN.pred.bar)<-c("SNR",c("glasso","lasso","Stepwise","SIS","Random Forest","BMA"))
FN.pred.bar.m<-melt(FN.pred.bar,id.vars = "SNR",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
colnames(FN.pred.bar.m)[3]<-"FNR"
ggplot(data = FN.pred.bar.m, mapping = aes(x = factor(SNR), y = FNR,fill = Method)) + geom_bar(stat = 'identity', position = 'dodge',color="blue")+xlab("SNR")+ggtitle("FNR vs SNR")+scale_fill_brewer(palette = "Spectral")


num.pred.bar<-sapply(rst.summary,function(x) x$num.pred)
num.pred.bar<-rbind(num.pred.bar,c(10,10,10,10,10))
num.pred.bar<-data.frame(SNRlist,t(num.pred.bar))
colnames(num.pred.bar)<-c("SNR",c("glasso","lasso","Stepwise","SIS","Random Forest","BMA","Truth"))
num.pred.bar.m<-melt(num.pred.bar,id.vars = "SNR",measure.vars = c(2:8),variable.name = "Method",na.rm=T)
colnames(num.pred.bar.m)[3]<-"num"
ggplot(data = num.pred.bar.m, mapping = aes(x = factor(SNR), y = num/m_G,color=Method,group=Method,size=Method)) + geom_line()+xlab("SNR")+ggtitle("Model Size vs SNR")+ylab("Estimated Nonzero Predictive Biomarker Proportion")+scale_size_manual(values=c(1.5,0.8,0.8,.8,.8,.8,1.5))+theme_bw() +scale_colour_manual(breaks=levels(num.pred.bar.m$Method), values=c("red","green","orange","blue","pink","purple","black"))
  theme(plot.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(),
        legend.key = element_blank(), legend.title = element_blank())

rst.summary<-list()
for(k in 1:5){
  rst.summary[[k]]<-data.frame(L2=apply(L2[[k]], 2, function(x) mean(x,na.rm=T)))
  rst.summary[[k]]$L2<-apply(L2[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$L1<-apply(L1[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$SSE<-apply(SSE[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.all<-apply(TP.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.all<-apply(FN.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.prog<-apply(TP.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.pred<-apply(TP.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.prog<-apply(FN.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.pred<-apply(FN.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$num.pred<-apply(num.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rownames(rst.summary[[k]])<-c("glasso","lasso","Stepwise","SIS","Random Forest","BMA")
  colnames(rst.summary[[k]])<-c("L2","L1","SSE","PPV all","FNR all","PPV prog","PPV pred","FNR prog","FNR pred","num pred")
}




for(k in 1:5){
  png(paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/summary/SNR/",SNRlist[k],".png"), height=200, width=700)
  p<-tableGrob(round(rst.summary[[k]],digits=3))
  grid.arrange(p)
  dev.off()
}


### SNP, p=100

n=100
m_X<-5
m_W<-1
m_G<-100
m_I<-m_G
SNR<-10
tau1<-1

k=1
# size<-list()
SSE<-list()
L1<-list()
L2<-list()
TP.all<-list()
FN.all<-list()
TP.prog<-list()
TP.pred<-list()
FN.prog<-list()
FN.pred<-list()
num.pred<-list()

n.method=6

for(j in 1){
  
  
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//glassorst_SNP.RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//lassorst_SNP.RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//bicrst_SNP.RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//steprst_SNP.RData"))
  sisrst<-readRDS(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//sisrst_SNP.RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//100//treerst_SNP.RData"))
  
  # size[[k]]<-list()
  #  size[[k]][[1]]<-sapply(glassorst,length)+m_X+m_W
  #   size[[k]][[2]]<-sapply(lassorst,length)+m_X+m_W
  
  glassorst<-glassorst_SNP
  lassorst<- lassorst_SNP
  steprst <- steprst_SNP
  treerst <- treerst_SNP
  bicrst <- bicrst_SNP
  
  SSE[[k]]<-matrix(0,ncol=n.method,nrow=100)
  L1[[k]]<-matrix(0,ncol=n.method,nrow=100)
  L2[[k]]<-matrix(0,ncol=n.method,nrow=100)
  TP.all[[k]]<-matrix(0,ncol=n.method,nrow=100)
  FN.all[[k]]<-matrix(0,ncol=n.method,nrow=100)
  TP.prog[[k]]<-matrix(0,ncol=n.method,nrow=100)
  TP.pred[[k]]<-matrix(0,ncol=n.method,nrow=100)
  FN.prog[[k]]<-matrix(0,ncol=n.method,nrow=100)
  FN.pred[[k]]<-matrix(0,ncol=n.method,nrow=100)
  num.pred[[k]]<-matrix(0,ncol=n.method,nrow=100)
  
  for(i in 1:100){
    print(i)
    set.seed(i+1000)
    
    # Generate X and Y
    sigma<-cov_block(m_G,.3,20)
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

    
    ## Group Lasso
    instance.glasso<-c(c(1:(m_X+m_W)),glassorst[[i]])
    instance.glasso<-instance.glasso[order(instance.glasso, decreasing = FALSE)]
    
    SSE.glasso<- -logLik(lm(y~-1+x[,instance.glasso]))
    
    beta.glasso<-lm(y~-1+x[,instance.glasso])$coef
    temp<-rep(0,length(beta))
    temp[instance.glasso]<-beta.glasso
    beta.glasso<-temp
    beta.glasso[which(is.na(beta.glasso)==T)]<-0
    L1.glasso<-norm(as.matrix(beta.glasso-beta),"1")
    L2.glasso<-norm(as.matrix(beta.glasso-beta),"2")
    
    
    TP.all.glasso<-(length(which(beta.glasso*beta!=0))-m_X-m_W)/(length(instance.glasso)-m_X-m_W)
    FN.all.glasso<-(length(intersect(which(beta.glasso==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.glasso0<-split_beta(beta.glasso,m_X,m_W,m_G,m_I)
    TP.prog.glasso<-length(which(beta.glasso0$G*beta0$G!=0))/length(which(beta.glasso0$G!=0))
    TP.pred.glasso<-length(which(beta.glasso0$I*beta0$I!=0))/length(which(beta.glasso0$I!=0))
    FN.prog.glasso<-(length(intersect(which(beta.glasso0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.glasso<-(length(intersect(which(beta.glasso0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.glasso<-length(which(beta.glasso0$I!=0))
    
    SSE[[k]][i,1]<-SSE.glasso
    L1[[k]][i,1]<-L1.glasso
    L2[[k]][i,1]<-L2.glasso
    TP.all[[k]][i,1]<-TP.all.glasso
    FN.all[[k]][i,1]<-FN.all.glasso
    TP.prog[[k]][i,1]<-TP.prog.glasso
    FN.prog[[k]][i,1]<-FN.prog.glasso
    TP.pred[[k]][i,1]<-TP.pred.glasso
    FN.pred[[k]][i,1]<-FN.pred.glasso
    num.pred[[k]][i,1]<-num.glasso
    
    ## Lasso
    instance.lasso<-c(c(1:(m_X+m_W)),lassorst[[i]])
    instance.lasso<-instance.lasso[order(instance.lasso, decreasing = FALSE)]
    
    SSE.lasso<- -logLik(lm(y~-1+x[,instance.lasso]))
    
    beta.lasso<-lm(y~-1+x[,instance.lasso])$coef
    temp<-rep(0,length(beta))
    temp[instance.lasso]<-beta.lasso
    beta.lasso<-temp
    beta.lasso[which(is.na(beta.lasso)==T)]<-0
    L1.lasso<-norm(as.matrix(beta.lasso-beta),"1")
    L2.lasso<-norm(as.matrix(beta.lasso-beta),"2")
    
    TP.all.lasso<-(length(which(beta.lasso*beta!=0))-m_X-m_W)/(length(instance.lasso)-m_X-m_W)
    FN.all.lasso<-(length(intersect(which(beta.lasso==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.lasso0<-split_beta(beta.lasso,m_X,m_W,m_G,m_I)
    TP.prog.lasso<-length(which(beta.lasso0$G*beta0$G!=0))/length(which(beta.lasso0$G!=0))
    TP.pred.lasso<-length(which(beta.lasso0$I*beta0$I!=0))/length(which(beta.lasso0$I!=0))
    FN.prog.lasso<-(length(intersect(which(beta.lasso0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.lasso<-(length(intersect(which(beta.lasso0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    
    
    SSE[[k]][i,2]<-SSE.lasso
    L1[[k]][i,2]<-L1.lasso
    L2[[k]][i,2]<-L2.lasso
    TP.all[[k]][i,2]<-TP.all.lasso
    FN.all[[k]][i,2]<-FN.all.lasso
    TP.prog[[k]][i,2]<-TP.prog.lasso
    FN.prog[[k]][i,2]<-FN.prog.lasso
    TP.pred[[k]][i,2]<-TP.pred.lasso
    FN.pred[[k]][i,2]<-FN.pred.lasso
    
    num.lasso<-length(which(beta.lasso0$I!=0))
    num.pred[[k]][i,2]<-num.lasso
    
    
    
    ## Stepwise
    instance.step<-steprst[[i]]
    
    SSE.step<- -logLik(lm(y~-1+x[,instance.step]))
    
    beta.step<-lm(y~-1+x[,instance.step])$coef
    temp<-rep(0,length(beta))
    temp[instance.step]<-beta.step
    beta.step<-temp
    beta.step[which(is.na(beta.step)==T)]<-0
    L1.step<-norm(as.matrix(beta.step-beta),"1")
    L2.step<-norm(as.matrix(beta.step-beta),"2")
    
    TP.all.step<-(length(which(beta.step*beta!=0))-m_X-m_W)/(length(instance.step)-m_X-m_W)
    FN.all.step<-(length(intersect(which(beta.step==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.step0<-split_beta(beta.step,m_X,m_W,m_G,m_I)
    TP.prog.step<-length(which(beta.step0$G*beta0$G!=0))/length(which(beta.step0$G!=0))
    TP.pred.step<-length(which(beta.step0$I*beta0$I!=0))/length(which(beta.step0$I!=0))
    FN.prog.step<-(length(intersect(which(beta.step0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.step<-(length(intersect(which(beta.step0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.step<-length(which(beta.step0$I!=0))
    
    SSE[[k]][i,3]<-SSE.step
    L1[[k]][i,3]<-L1.step
    L2[[k]][i,3]<-L2.step
    TP.all[[k]][i,3]<-TP.all.step
    FN.all[[k]][i,3]<-FN.all.step
    TP.prog[[k]][i,3]<-TP.prog.step
    FN.prog[[k]][i,3]<-FN.prog.step
    TP.pred[[k]][i,3]<-TP.pred.step
    FN.pred[[k]][i,3]<-FN.pred.step
    num.pred[[k]][i,3]<-num.step
    
    ## SIS
    instance.sis<-union(c(1:(m_X+m_W)),sisrst[[i]])
    
    SSE.sis<- -logLik(lm(y~-1+x[,instance.sis]))
    
    beta.sis<-lm(y~-1+x[,instance.sis])$coef
    temp<-rep(0,length(beta))
    temp[instance.sis]<-beta.sis
    beta.sis<-temp
    beta.sis[which(is.na(beta.sis)==T)]<-0
    L1.sis<-norm(as.matrix(beta.sis-beta),"1")
    L2.sis<-norm(as.matrix(beta.sis-beta),"2")
    
    TP.all.sis<-(length(which(beta.sis*beta!=0))-m_X-m_W)/(length(instance.sis)-m_X-m_W)
    FN.all.sis<-(length(intersect(which(beta.sis==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.sis0<-split_beta(beta.sis,m_X,m_W,m_G,m_I)
    TP.prog.sis<-length(which(beta.sis0$G*beta0$G!=0))/length(which(beta.sis0$G!=0))
    TP.pred.sis<-length(which(beta.sis0$I*beta0$I!=0))/length(which(beta.sis0$I!=0))
    FN.prog.sis<-(length(intersect(which(beta.sis0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.sis<-(length(intersect(which(beta.sis0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.sis<-length(which(beta.sis0$I!=0))
    
    SSE[[k]][i,4]<-SSE.sis
    L1[[k]][i,4]<-L1.sis
    L2[[k]][i,4]<-L2.sis
    TP.all[[k]][i,4]<-TP.all.sis
    FN.all[[k]][i,4]<-FN.all.sis
    TP.prog[[k]][i,4]<-TP.prog.sis
    FN.prog[[k]][i,4]<-FN.prog.sis
    TP.pred[[k]][i,4]<-TP.pred.sis
    FN.pred[[k]][i,4]<-FN.pred.sis
    num.pred[[k]][i,4]<-num.sis
    
    
    ## Random Forest
    instance.tree<-union(c(1:(m_X+m_W)),treerst[[i]])
    instance.tree<-instance.tree[order(instance.tree, decreasing = FALSE)]
    
    SSE.tree<- -logLik(lm(y~-1+x[,instance.tree]))
    
    beta.tree<-lm(y~-1+x[,instance.tree])$coef
    temp<-rep(0,length(beta))
    temp[instance.tree]<-beta.tree
    beta.tree<-temp
    beta.tree[which(is.na(beta.tree)==T)]<-0
    L1.tree<-norm(as.matrix(beta.tree-beta),"1")
    L2.tree<-norm(as.matrix(beta.tree-beta),"2")
    
    TP.all.tree<-(length(which(beta.tree*beta!=0))-m_X-m_W)/(length(instance.tree)-m_X-m_W)
    FN.all.tree<-(length(intersect(which(beta.tree==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.tree0<-split_beta(beta.tree,m_X,m_W,m_G,m_I)
    TP.prog.tree<-length(which(beta.tree0$G*beta0$G!=0))/length(which(beta.tree0$G!=0))
    TP.pred.tree<-length(which(beta.tree0$I*beta0$I!=0))/length(which(beta.tree0$I!=0))
    FN.prog.tree<-(length(intersect(which(beta.tree0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.tree<-(length(intersect(which(beta.tree0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.tree<-length(which(beta.tree0$I!=0))
    
    SSE[[k]][i,5]<-SSE.tree
    L1[[k]][i,5]<-L1.tree
    L2[[k]][i,5]<-L2.tree
    TP.all[[k]][i,5]<-TP.all.tree
    FN.all[[k]][i,5]<-FN.all.tree
    TP.prog[[k]][i,5]<-TP.prog.tree
    FN.prog[[k]][i,5]<-FN.prog.tree
    TP.pred[[k]][i,5]<-TP.pred.tree
    FN.pred[[k]][i,5]<-FN.pred.tree
    num.pred[[k]][i,5]<-num.tree
    
    ## BMA
    instance.bic<-bicrst[[i]][which(is.na(bicrst[[i]])==F)]
    instance.bic<-c(c(1:(m_X+m_W)),instance.bic)
    
    SSE.bic<- -logLik(lm(y~-1+x[,instance.tree]))
    
    beta.bic<-lm(y~-1+x[,instance.bic])$coef
    temp<-rep(0,length(beta))
    temp[instance.bic]<-beta.bic
    beta.bic<-temp
    beta.bic[which(is.na(beta.bic)==T)]<-0
    L1.bic<-norm(as.matrix(beta.bic-beta),"1")
    L2.bic<-norm(as.matrix(beta.bic-beta),"2")
    
    TP.all.bic<-(length(which(beta.bic*beta!=0))-m_X-m_W)/(length(instance.bic)-m_X-m_W)
    FN.all.bic<-(length(intersect(which(beta.bic==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.bic0<-split_beta(beta.bic,m_X,m_W,m_G,m_I)
    TP.prog.bic<-length(which(beta.bic0$G*beta0$G!=0))/length(which(beta.bic0$G!=0))
    TP.pred.bic<-length(which(beta.bic0$I*beta0$I!=0))/length(which(beta.bic0$I!=0))
    FN.prog.bic<-(length(intersect(which(beta.bic0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.bic<-(length(intersect(which(beta.bic0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.bic<-length(which(beta.bic0$I!=0))
    
    SSE[[k]][i,6]<-SSE.bic
    L1[[k]][i,6]<-L1.bic
    L2[[k]][i,6]<-L2.bic
    TP.all[[k]][i,6]<-TP.all.bic
    FN.all[[k]][i,6]<-FN.all.bic
    TP.prog[[k]][i,6]<-TP.prog.bic
    FN.prog[[k]][i,6]<-FN.prog.bic
    TP.pred[[k]][i,6]<-TP.pred.bic
    FN.pred[[k]][i,6]<-FN.pred.bic
    num.pred[[k]][i,6]<-num.bic
    
  }
  
  k=k+1
}

rst.summary<-list()
for(k in 1){
  rst.summary[[k]]<-data.frame(L2=apply(L2[[k]], 2, function(x) mean(x,na.rm=T)))
  rst.summary[[k]]$L2<-apply(L2[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$L1<-apply(L1[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$SSE<-apply(SSE[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.all<-apply(TP.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.all<-apply(FN.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.prog<-apply(TP.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.pred<-apply(TP.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.prog<-apply(FN.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.pred<-apply(FN.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$num.pred<-apply(num.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rownames(rst.summary[[k]])<-c("glasso","lasso","Stepwise","SIS","Random Forest","BMA")
  rst.summary[[k]][which(is.na(rst.summary[[k]])==T,arr.ind = T)]<-0
}




# Bar Plot
TP.pred.bar<-sapply(rst.summary,function(x) x$TP.pred)
TP.pred.bar<-data.frame("SNP",t(TP.pred.bar))
colnames(TP.pred.bar)<-c("SNP",c("glasso","lasso","Stepwise","SIS","Random Forest","BMA"))
TP.pred.bar.m<-melt(TP.pred.bar,id.vars = "SNP",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
colnames(TP.pred.bar.m)[3]<-"PPV"
ggplot(data = TP.pred.bar.m, mapping = aes(x = factor(SNP), y = PPV,fill = Method)) + geom_bar(stat = 'identity', position = 'dodge',color="blue")+xlab("SNP")+ggtitle("PPV when covariates are SNP")+scale_fill_brewer(palette="Spectral")

FN.pred.bar<-sapply(rst.summary,function(x) x$FN.pred)
FN.pred.bar<-data.frame("SNP",t(FN.pred.bar))
colnames(FN.pred.bar)<-c("SNP",c("glasso","lasso","Stepwise","SIS","Random Forest","BMA"))
FN.pred.bar.m<-melt(FN.pred.bar,id.vars = "SNP",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
colnames(FN.pred.bar.m)[3]<-"FNR"
ggplot(data = FN.pred.bar.m, mapping = aes(x = factor(SNP), y = FNR,fill = Method)) + geom_bar(stat = 'identity', position = 'dodge',color="blue")+xlab("SNP")+ggtitle("FNR when covariates are SNP")+scale_fill_brewer(palette="Spectral")


num.pred.bar<-sapply(rst.summary,function(x) x$num.pred)
num.pred.bar<-rbind(num.pred.bar,c(5,10,20))
num.pred.bar<-data.frame(dim,t(num.pred.bar))
colnames(num.pred.bar)<-c("Dimension",c("glasso","lasso","Stepwise","SIS","Random Forest","BMA","Truth"))
num.pred.bar.m<-melt(num.pred.bar,id.vars = "Dimension",measure.vars = c(2:8),variable.name = "Method",na.rm=T)
colnames(num.pred.bar.m)[3]<-"num"
ggplot(data = num.pred.bar.m, mapping = aes(x = factor(Dimension), y = num/dim,color=Method,group=Method)) + geom_line(size=3)+xlab("Number of Biomarkers")+ggtitle("Model Size vs Number of Biomarkers")+ylab("Estimated Nonzero Predictive Biomarker Proportion")



## Different P, p=200

n=100
m_X<-5
m_W<-1


SNR<-10
tau1<-1

k=1
# size<-list()
SSE<-list()
L1<-list()
L2<-list()
TP.all<-list()
FN.all<-list()
TP.prog<-list()
TP.pred<-list()
FN.prog<-list()
FN.pred<-list()
num.pred<-list()


n.method=6
dim<-c(50,100,200)

for(j in 1:3){
  
  m_G<-dim[j]
  m_I<-m_G
  
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//",dim[j],"//glassorst.RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//",dim[j],"//lassorst.RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//",dim[j],"//bicrst.RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//",dim[j],"//steprst.RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//",dim[j],"//sisrst.RData"))
  load(file=paste0("C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//",dim[j],"//treerst.RData"))
  
  # size[[k]]<-list()
  #  size[[k]][[1]]<-sapply(glassorst,length)+m_X+m_W
  #   size[[k]][[2]]<-sapply(lassorst,length)+m_X+m_W
  
  
  
  SSE[[k]]<-matrix(0,ncol=n.method,nrow=100)
  L1[[k]]<-matrix(0,ncol=n.method,nrow=100)
  L2[[k]]<-matrix(0,ncol=n.method,nrow=100)
  TP.all[[k]]<-matrix(0,ncol=n.method,nrow=100)
  FN.all[[k]]<-matrix(0,ncol=n.method,nrow=100)
  TP.prog[[k]]<-matrix(0,ncol=n.method,nrow=100)
  TP.pred[[k]]<-matrix(0,ncol=n.method,nrow=100)
  FN.prog[[k]]<-matrix(0,ncol=n.method,nrow=100)
  FN.pred[[k]]<-matrix(0,ncol=n.method,nrow=100)
  num.pred[[k]]<-matrix(0,ncol=n.method,nrow=100)
  
  for(i in 1:100){
    print(i)
    set.seed(i+1000)
    
    # Generate X and Y
    sigma<-cov_block(m_G,.3,m_G/5)
    #sigma<-GenerateCliquesCovariance(10,10,0.8)
    #binprob<-runif(m_G)
    x<-sim_X(m_X,m_W,m_G,sigma,n)
    #beta<-sim_beta(m_X=0,m_W=1,m_G,main_nonzero=0.05,inter_nonzero=0.05,both_nonzero=0.1,bit=T,heir=T)
    beta<-sim_beta_const(m_X,m_W=1,m_G,main_nonzero=.1,inter_nonzero=.1,both_nonzero=0.01,const=c(3,5),heir=TRUE)
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
    
    SSE.glasso<- -logLik(lm(y~-1+x[,instance.glasso]))
    
    beta.glasso<-lm(y~-1+x[,instance.glasso])$coef
    temp<-rep(0,length(beta))
    temp[instance.glasso]<-beta.glasso
    beta.glasso<-temp
    beta.glasso[which(is.na(beta.glasso)==T)]<-0
    L1.glasso<-norm(as.matrix(beta.glasso-beta),"1")
    L2.glasso<-norm(as.matrix(beta.glasso-beta),"2")
    
    
    TP.all.glasso<-(length(which(beta.glasso*beta!=0))-m_X-m_W)/(length(instance.glasso)-m_X-m_W)
    FN.all.glasso<-(length(intersect(which(beta.glasso==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.glasso0<-split_beta(beta.glasso,m_X,m_W,m_G,m_I)
    TP.prog.glasso<-length(which(beta.glasso0$G*beta0$G!=0))/length(which(beta.glasso0$G!=0))
    TP.pred.glasso<-length(which(beta.glasso0$I*beta0$I!=0))/length(which(beta.glasso0$I!=0))
    FN.prog.glasso<-(length(intersect(which(beta.glasso0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.glasso<-(length(intersect(which(beta.glasso0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.glasso<-length(which(beta.glasso0$I!=0))
    
    SSE[[k]][i,1]<-SSE.glasso
    L1[[k]][i,1]<-L1.glasso
    L2[[k]][i,1]<-L2.glasso
    TP.all[[k]][i,1]<-TP.all.glasso
    FN.all[[k]][i,1]<-FN.all.glasso
    TP.prog[[k]][i,1]<-TP.prog.glasso
    FN.prog[[k]][i,1]<-FN.prog.glasso
    TP.pred[[k]][i,1]<-TP.pred.glasso
    FN.pred[[k]][i,1]<-FN.pred.glasso
    num.pred[[k]][i,1]<-num.glasso
    
    ## Lasso
    instance.lasso<-c(c(1:(m_X+m_W)),lassorst[[i]])
    instance.lasso<-instance.lasso[order(instance.lasso, decreasing = FALSE)]
    
    SSE.lasso<- -logLik(lm(y~-1+x[,instance.lasso]))
    
    beta.lasso<-lm(y~-1+x[,instance.lasso])$coef
    temp<-rep(0,length(beta))
    temp[instance.lasso]<-beta.lasso
    beta.lasso<-temp
    beta.lasso[which(is.na(beta.lasso)==T)]<-0
    L1.lasso<-norm(as.matrix(beta.lasso-beta),"1")
    L2.lasso<-norm(as.matrix(beta.lasso-beta),"2")
    
    TP.all.lasso<-(length(which(beta.lasso*beta!=0))-m_X-m_W)/(length(instance.lasso)-m_X-m_W)
    FN.all.lasso<-(length(intersect(which(beta.lasso==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.lasso0<-split_beta(beta.lasso,m_X,m_W,m_G,m_I)
    TP.prog.lasso<-length(which(beta.lasso0$G*beta0$G!=0))/length(which(beta.lasso0$G!=0))
    TP.pred.lasso<-length(which(beta.lasso0$I*beta0$I!=0))/length(which(beta.lasso0$I!=0))
    FN.prog.lasso<-(length(intersect(which(beta.lasso0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.lasso<-(length(intersect(which(beta.lasso0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    
    
    SSE[[k]][i,2]<-SSE.lasso
    L1[[k]][i,2]<-L1.lasso
    L2[[k]][i,2]<-L2.lasso
    TP.all[[k]][i,2]<-TP.all.lasso
    FN.all[[k]][i,2]<-FN.all.lasso
    TP.prog[[k]][i,2]<-TP.prog.lasso
    FN.prog[[k]][i,2]<-FN.prog.lasso
    TP.pred[[k]][i,2]<-TP.pred.lasso
    FN.pred[[k]][i,2]<-FN.pred.lasso
    
    num.lasso<-length(which(beta.lasso0$I!=0))
    num.pred[[k]][i,2]<-num.lasso
    
    
    
    ## Stepwise
    instance.step<-steprst[[i]]
    
    SSE.step<- -logLik(lm(y~-1+x[,instance.step]))
    
    beta.step<-lm(y~-1+x[,instance.step])$coef
    temp<-rep(0,length(beta))
    temp[instance.step]<-beta.step
    beta.step<-temp
    beta.step[which(is.na(beta.step)==T)]<-0
    L1.step<-norm(as.matrix(beta.step-beta),"1")
    L2.step<-norm(as.matrix(beta.step-beta),"2")
    
    TP.all.step<-(length(which(beta.step*beta!=0))-m_X-m_W)/(length(instance.step)-m_X-m_W)
    FN.all.step<-(length(intersect(which(beta.step==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.step0<-split_beta(beta.step,m_X,m_W,m_G,m_I)
    TP.prog.step<-length(which(beta.step0$G*beta0$G!=0))/length(which(beta.step0$G!=0))
    TP.pred.step<-length(which(beta.step0$I*beta0$I!=0))/length(which(beta.step0$I!=0))
    FN.prog.step<-(length(intersect(which(beta.step0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.step<-(length(intersect(which(beta.step0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.step<-length(which(beta.step0$I!=0))
    
    SSE[[k]][i,3]<-SSE.step
    L1[[k]][i,3]<-L1.step
    L2[[k]][i,3]<-L2.step
    TP.all[[k]][i,3]<-TP.all.step
    FN.all[[k]][i,3]<-FN.all.step
    TP.prog[[k]][i,3]<-TP.prog.step
    FN.prog[[k]][i,3]<-FN.prog.step
    TP.pred[[k]][i,3]<-TP.pred.step
    FN.pred[[k]][i,3]<-FN.pred.step
    num.pred[[k]][i,3]<-num.step
    
    ## SIS
    instance.sis<-union(c(1:(m_X+m_W)),sisrst[[i]])
    
    SSE.sis<- -logLik(lm(y~-1+x[,instance.sis]))
    
    beta.sis<-lm(y~-1+x[,instance.sis])$coef
    temp<-rep(0,length(beta))
    temp[instance.sis]<-beta.sis
    beta.sis<-temp
    beta.sis[which(is.na(beta.sis)==T)]<-0
    L1.sis<-norm(as.matrix(beta.sis-beta),"1")
    L2.sis<-norm(as.matrix(beta.sis-beta),"2")
    
    TP.all.sis<-(length(which(beta.sis*beta!=0))-m_X-m_W)/(length(instance.sis)-m_X-m_W)
    FN.all.sis<-(length(intersect(which(beta.sis==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.sis0<-split_beta(beta.sis,m_X,m_W,m_G,m_I)
    TP.prog.sis<-length(which(beta.sis0$G*beta0$G!=0))/length(which(beta.sis0$G!=0))
    TP.pred.sis<-length(which(beta.sis0$I*beta0$I!=0))/length(which(beta.sis0$I!=0))
    FN.prog.sis<-(length(intersect(which(beta.sis0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.sis<-(length(intersect(which(beta.sis0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.sis<-length(which(beta.sis0$I!=0))
    
    SSE[[k]][i,4]<-SSE.sis
    L1[[k]][i,4]<-L1.sis
    L2[[k]][i,4]<-L2.sis
    TP.all[[k]][i,4]<-TP.all.sis
    FN.all[[k]][i,4]<-FN.all.sis
    TP.prog[[k]][i,4]<-TP.prog.sis
    FN.prog[[k]][i,4]<-FN.prog.sis
    TP.pred[[k]][i,4]<-TP.pred.sis
    FN.pred[[k]][i,4]<-FN.pred.sis
    num.pred[[k]][i,4]<-num.sis
    
    
    ## Random Forest
    instance.tree<-union(c(1:(m_X+m_W)),treerst[[i]])
    instance.tree<-instance.tree[order(instance.tree, decreasing = FALSE)]
    
    SSE.tree<- -logLik(lm(y~-1+x[,instance.tree]))
    
    beta.tree<-lm(y~-1+x[,instance.tree])$coef
    temp<-rep(0,length(beta))
    temp[instance.tree]<-beta.tree
    beta.tree<-temp
    beta.tree[which(is.na(beta.tree)==T)]<-0
    L1.tree<-norm(as.matrix(beta.tree-beta),"1")
    L2.tree<-norm(as.matrix(beta.tree-beta),"2")
    
    TP.all.tree<-(length(which(beta.tree*beta!=0))-m_X-m_W)/(length(instance.tree)-m_X-m_W)
    FN.all.tree<-(length(intersect(which(beta.tree==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.tree0<-split_beta(beta.tree,m_X,m_W,m_G,m_I)
    TP.prog.tree<-length(which(beta.tree0$G*beta0$G!=0))/length(which(beta.tree0$G!=0))
    TP.pred.tree<-length(which(beta.tree0$I*beta0$I!=0))/length(which(beta.tree0$I!=0))
    FN.prog.tree<-(length(intersect(which(beta.tree0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.tree<-(length(intersect(which(beta.tree0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.tree<-length(which(beta.tree0$I!=0))
    
    SSE[[k]][i,5]<-SSE.tree
    L1[[k]][i,5]<-L1.tree
    L2[[k]][i,5]<-L2.tree
    TP.all[[k]][i,5]<-TP.all.tree
    FN.all[[k]][i,5]<-FN.all.tree
    TP.prog[[k]][i,5]<-TP.prog.tree
    FN.prog[[k]][i,5]<-FN.prog.tree
    TP.pred[[k]][i,5]<-TP.pred.tree
    FN.pred[[k]][i,5]<-FN.pred.tree
    num.pred[[k]][i,5]<-num.tree
    
    ## BMA
    instance.bic<-bicrst[[i]][which(is.na(bicrst[[i]])==F)]
    instance.bic<-c(c(1:(m_X+m_W)),instance.bic)
    
    SSE.bic<- -logLik(lm(y~-1+x[,instance.tree]))
    
    beta.bic<-lm(y~-1+x[,instance.bic])$coef
    temp<-rep(0,length(beta))
    temp[instance.bic]<-beta.bic
    beta.bic<-temp
    beta.bic[which(is.na(beta.bic)==T)]<-0
    L1.bic<-norm(as.matrix(beta.bic-beta),"1")
    L2.bic<-norm(as.matrix(beta.bic-beta),"2")
    
    TP.all.bic<-(length(which(beta.bic*beta!=0))-m_X-m_W)/(length(instance.bic)-m_X-m_W)
    FN.all.bic<-(length(intersect(which(beta.bic==0),which(beta!=0))))/(length(which(beta!=0))-m_X-m_W)
    beta.bic0<-split_beta(beta.bic,m_X,m_W,m_G,m_I)
    TP.prog.bic<-length(which(beta.bic0$G*beta0$G!=0))/length(which(beta.bic0$G!=0))
    TP.pred.bic<-length(which(beta.bic0$I*beta0$I!=0))/length(which(beta.bic0$I!=0))
    FN.prog.bic<-(length(intersect(which(beta.bic0$G==0),which(beta0$G!=0))))/(length(which(beta0$G!=0)))
    FN.pred.bic<-(length(intersect(which(beta.bic0$I==0),which(beta0$I!=0))))/(length(which(beta0$I!=0)))
    num.bic<-length(which(beta.bic0$I!=0))
    
    SSE[[k]][i,6]<-SSE.bic
    L1[[k]][i,6]<-L1.bic
    L2[[k]][i,6]<-L2.bic
    TP.all[[k]][i,6]<-TP.all.bic
    FN.all[[k]][i,6]<-FN.all.bic
    TP.prog[[k]][i,6]<-TP.prog.bic
    FN.prog[[k]][i,6]<-FN.prog.bic
    TP.pred[[k]][i,6]<-TP.pred.bic
    FN.pred[[k]][i,6]<-FN.pred.bic
    num.pred[[k]][i,6]<-num.bic
    
  }
  
  k=k+1
}

rst.summary<-list()
for(k in 1:3){
  rst.summary[[k]]<-data.frame(L2=apply(L2[[k]], 2, function(x) mean(x,na.rm=T)))
  rst.summary[[k]]$L2<-apply(L2[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$L1<-apply(L1[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$SSE<-apply(SSE[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.all<-apply(TP.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.all<-apply(FN.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.prog<-apply(TP.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.pred<-apply(TP.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.prog<-apply(FN.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.pred<-apply(FN.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$num.pred<-apply(num.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rownames(rst.summary[[k]])<-c("glasso","lasso","Stepwise","SIS","Random Forest","BMA")
  rst.summary[[k]][which(is.na(rst.summary[[k]])==T,arr.ind = T)]<-0
}


# Ridge Plot
TP.pred.m<-list()
q<-list()
for(k in 1:4){
  TP.pred.m[[k]]<-data.frame("ID"=seq(1,100,1),TP.pred[[k]])
  colnames(TP.pred.m[[k]])<-c("ID","glasso","lasso","Stepwise","SIS","Random Forest","BMA")
  TP.pred.m[[k]]<-melt(TP.pred.m[[k]],id.vars = "ID",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
  q[[k]]<-ggplot(TP.pred.m[[k]],aes(x=value,y=Method))+stat_density_ridges(data=TP.pred.m[[k]], quantile_lines = T,scale = 1, size = 0.25, rel_min_height = 0,fill="skyblue",alpha=.9,color="white")+ggtitle(paste0("Nonzero Interaction Proportion:",portion[k]))+xlab("Predictive Biomarkers PPV")
}
FN.pred.m<-list()
p<-list()
for(k in 1:4){
  FN.pred.m[[k]]<-data.frame("ID"=seq(1,100,1),FN.pred[[k]])
  colnames(FN.pred.m[[k]])<-c("ID","glasso","lasso","Stepwise","SIS","Random Forest","BMA")
  FN.pred.m[[k]]<-melt(FN.pred.m[[k]],id.vars = "ID",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
  p[[k]]<-ggplot(FN.pred.m[[k]],aes(x=value,y=Method))+stat_density_ridges(data=FN.pred.m[[k]], quantile_lines = T,scale = 1, size = 0.25, rel_min_height = 0,fill="skyblue",alpha=.9,color="white")+ggtitle(paste0("Nonzero Interaction Proportion:",portion[k]))+xlab("Predictive Biomarkers FNR")
}

# Bar Plot
TP.pred.bar<-sapply(rst.summary,function(x) x$TP.pred)
TP.pred.bar<-data.frame(dim,t(TP.pred.bar))
colnames(TP.pred.bar)<-c("Dimension",c("glasso","lasso","Stepwise","SIS","Random Forest","BMA"))
TP.pred.bar.m<-melt(TP.pred.bar,id.vars = "Dimension",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
colnames(TP.pred.bar.m)[3]<-"PPV"
ggplot(data = TP.pred.bar.m, mapping = aes(x = factor(Dimension), y = PPV,fill = Method)) + geom_bar(stat = 'identity', position = 'dodge',color="blue")+xlab("Number of Biomarkers")+ggtitle("PPV vs Number of Biomarkers")+scale_fill_brewer(palette="Spectral")

FN.pred.bar<-sapply(rst.summary,function(x) x$FN.pred)
FN.pred.bar<-data.frame(dim,t(FN.pred.bar))
colnames(FN.pred.bar)<-c("Dimension",c("glasso","lasso","Stepwise","SIS","Random Forest","BMA"))
FN.pred.bar.m<-melt(FN.pred.bar,id.vars = "Dimension",measure.vars = c(2:7),variable.name = "Method",na.rm=T)
colnames(FN.pred.bar.m)[3]<-"FNR"
ggplot(data = FN.pred.bar.m, mapping = aes(x = factor(Dimension), y = FNR,fill = Method)) + geom_bar(stat = 'identity', position = 'dodge',color="black")+xlab("Number of Biomarkers")+ggtitle("FNR vs Number of Biomarkers")+scale_fill_brewer(palette="Set1")


num.pred.bar<-sapply(rst.summary,function(x) x$num.pred)
num.pred.bar<-rbind(num.pred.bar,c(5,10,20))
num.pred.bar<-data.frame(dim,t(num.pred.bar))
colnames(num.pred.bar)<-c("Dimension",c("glasso","lasso","Stepwise","SIS","Random Forest","BMA","Truth"))
num.pred.bar.m<-melt(num.pred.bar,id.vars = "Dimension",measure.vars = c(2:8),variable.name = "Method",na.rm=T)
colnames(num.pred.bar.m)[3]<-"num"
ggplot(data = num.pred.bar.m[which(num.pred.bar.m$Method!="BMA"),], mapping = aes(x = factor(Dimension), y = num/dim,color=Method,group=Method,size=Method)) + geom_line()+xlab("Number of Biomarkers")+ggtitle("Model Size vs Number of Biomarkers")+ylab("Estimated Nonzero Predictive Biomarker Proportion")+scale_size_manual(values=c(1.5,0.8,0.8,.8,.8,1.5))+theme_bw() +scale_colour_manual(breaks=levels(num.pred.bar.m$Method), values=c("red","green","blue","pink","purple","black"))
theme(plot.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(),
      legend.key = element_blank(), legend.title = element_blank())


rst.summary<-list()
for(k in 1:4){
  rst.summary[[k]]<-data.frame(L2=apply(L2[[k]], 2, function(x) mean(x,na.rm=T)))
  rst.summary[[k]]$L2<-apply(L2[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$L1<-apply(L1[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$SSE<-apply(SSE[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.all<-apply(TP.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.all<-apply(FN.all[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.prog<-apply(TP.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$TP.pred<-apply(TP.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.prog<-apply(FN.prog[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$FN.pred<-apply(FN.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rst.summary[[k]]$num.pred<-apply(num.pred[[k]], 2, function(x) mean(x,na.rm=T))
  rownames(rst.summary[[k]])<-c("glasso","lasso","Stepwise","SIS","Random Forest","BMA")
  colnames(rst.summary[[k]])<-c("L2","L1","SSE","PPV all","FNR all","PPV prog","PPV pred","FNR prog","FNR pred","num pred")
}


portion<-c(0.05,0.1,0.15,0.2)
for(k in 1:4){
  png(paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/summary/proportion/",portion[k],".png"), height=200, width=700)
  p<-tableGrob(round(rst.summary[[k]],digits=3))
  grid.arrange(p)
  dev.off()
}


##### Cross Validation

dir<-"C:\\Users\\auz5836\\Documents\\GitHub\\GroupLasso\\Final_Simu\\50\\"
sol_cv_glasso<-list()
for(i in 1:20){
  sol_cv_glasso[[i]]<-readRDS(paste0(dir,"sol_cv_glasso_",i,".RData"))
}
sol_cv_mean_pred<-matrix(0,nrow=10,ncol=5)
sol_cv_MSE_beta<-matrix(0,nrow=10,ncol=5)
sol_cv_SD_pred<-matrix(0,nrow=10,ncol=5)
sol_cv_var_num<-matrix(0,nrow=10,ncol=5)
sol_cv_mean_num<-matrix(0,nrow=10,ncol=5)
for(i in 1:20){
  sol_cv_mean_pred<-sol_cv_mean_pred+sol_cv_glasso[[i]]$mean_pred
  sol_cv_MSE_beta<-sol_cv_MSE_beta+sol_cv_glasso[[i]]$MSE_beta
  sol_cv_SD_pred<-sol_cv_SD_pred+sol_cv_glasso[[i]]$SD_pred
  sol_cv_var_num<-sol_cv_var_num+sol_cv_glasso[[i]]$Var_num
  sol_cv_mean_num<-sol_cv_mean_num+sol_cv_glasso[[i]]$mean_num
}

sol_cv_mean_pred<- sol_cv_mean_pred/20
sol_cv_MSE_beta<- sol_cv_MSE_beta/20
sol_cv_SD_pred<- sol_cv_SD_pred/20
sol_cv_var_num<-sol_cv_var_num/20
sol_cv_mean_num<-sol_cv_mean_num/20



sol_cv_mean_pred.m<-data.frame("lamb1"=seq(1,19,2),sol_cv_mean_pred)
colnames(sol_cv_mean_pred.m)<-c("lamb1",c(1:5))
sol_cv_mean_pred.m<-melt(sol_cv_mean_pred.m,id.vars = "lamb1",measure.vars = c(2:6),variable.name = "lamb2")
p <- ggplot(sol_cv_mean_pred.m, aes(lamb1, lamb2)) + geom_tile(aes(fill = value),
colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")+ ggtitle("Heatmap for Mean Values of SSE on Outcome Predictions in Cross Validation")

sol_cv_MSE_beta.m<-data.frame("lamb1"=seq(1,19,2),sol_cv_MSE_beta)
colnames(sol_cv_MSE_beta.m)<-c("lamb1",c(1:5))
sol_cv_MSE_beta.m<-melt(sol_cv_MSE_beta.m,id.vars = "lamb1",measure.vars = c(2:6),variable.name = "lamb2")
p <- ggplot(sol_cv_MSE_beta.m, aes(lamb1, lamb2)) + geom_tile(aes(fill = value),
                                                               colour = "white") + scale_fill_gradient(low = "white",high = "steelblue") + ggtitle("Heatmap for Mean Values of MSE on Coefficient Estimations in Cross Validation")
                                               

sol_cv_var_num.m<-data.frame("lamb1"=seq(1,19,2),sol_cv_var_num)
colnames(sol_cv_var_num.m)<-c("lamb1",c(1:5))
sol_cv_var_num.m<-melt(sol_cv_var_num.m,id.vars = "lamb1",measure.vars = c(2:6),variable.name = "lamb2")
p <- ggplot(sol_cv_var_num.m, aes(lamb1, lamb2)) + geom_tile(aes(fill = value),
                                                               colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")+ ggtitle("Heatmap for Mean Values of Model Size Variance in Cross Validation")

sol_cv_mean_num.m<-data.frame("lamb1"=seq(1,19,2),sol_cv_mean_num)
colnames(sol_cv_mean_num.m)<-c("lamb1",c(1:5))
sol_cv_mean_num.m<-melt(sol_cv_mean_num.m,id.vars = "lamb1",measure.vars = c(2:6),variable.name = "lamb2")
p <- ggplot(sol_cv_mean_num.m, aes(lamb1, lamb2)) + geom_tile(aes(fill = value),
                                                             colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")+ ggtitle("Heatmap for Mean Values of Model Size in Cross Validation")

