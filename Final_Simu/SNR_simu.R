n=100
m_X<-5
m_W<-1
m_G<-100
m_I<-m_G
SNR<-10
tau1<-1


SNRlist<-c(1,100)


for(SNR in SNRlist){
  
  #### Iteration
#   treerst<-list()
#   bicrst<-list()
#   steprst<-list()
  glassorst<-list()
  lassorst<-list()
  # sisrst<-list()
  
  for(i in 1:100){
    print(i)
    #Set Seed
    set.seed(i+1000)
    
    # Generate X and Y
    sigma<-cov_block(m_G,.3,5)
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
    
#     simu<-data.frame(X=x,Y=y)
#     L<-lm(y~x[,c(1:(m_X+m_W))])
#     y.res<-L$residuals
#     simu$Y<-y.res
#     simu<-simu[,-c(1:(m_X+m_W))]
#     
#     ### Trees
#     model <- randomForest(Y~.,   data=simu)
#     #print(model) # view results 
#     #importance(model)
#     treerst[[i]]<-order(importance(model),decreasing = T)[1:length(truth)]
#     
#     ### BMA
#     bicfit<-bicreg(x[,-c(1:(m_X+m_W))],y.res,strict = T)
#     bicrst[[i]]<-bicfit$namesx[order(bicfit$probne0,decreasing = T)][1:length(truth)]
#     bicrst[[i]]<-sapply(bicrst[[i]],function(x) strsplit(x,"X")[[1]][2])
#     bicrst[[i]]<-as.integer(bicrst[[i]])
#     
#     ### Stepwise
#     a<-regsubsets(x=x,y=y,method="forward",nvmax = 3*length(truth),force.in = c(1:(m_X+m_W)))
#     steprst[[i]]<-a$vorder[1:(length(truth))]
#     steprst[[i]]<-steprst[order(steprst[[i]])][[1]]
    
    x0<-rep(0,dim(x)[2])
    ### Group Lasso
    lamb_opt_glasso<-20
    lamb_opt2_glasso<-1
    solg<-FASTA(x,y,f, gradf, g, proxg, x0, tau1, max_iters = 300, w = 10, 
                backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                eps_n = 1e-15,m_X,m_W,m_G,m_G,lamb_opt_glasso,lamb_opt2_glasso,restart=TRUE)
    glassorst[[i]]<-order(abs(solg$x[-c(1:(m_W+m_X))]),decreasing = T)[1:(length(truth)-m_X-m_W)]+m_X+m_W
   
    
    ### Regular Lasso
    lamb_opt_lasso<-250
    lamb_opt2_lasso<-5
    sol<-FASTA(x,y,f, gradf, glasso, proxglasso, x0, tau1, max_iters = 300, w = 10, 
               backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
               eps_n = 1e-15,m_X,m_W,m_G,m_G,lamb_opt_lasso,lamb_opt2_lasso,restart=TRUE)
    lassorst[[i]]<-order(abs(sol$x[-c(1:(m_W+m_X))]),decreasing = T)[1:(length(truth)-m_X-m_W)]+m_X+m_W
   
#     ### SIS
#     #model1<-SIS(x,y,family = "gaussian", penalty = "lasso", tune="bic")
#     model2<-SIS(x,y,family = "gaussian", penalty = "lasso", tune="bic",varISIS = "aggr")
#     sisrst[[i]]<-model2$ix
  }
  
  
  #save(truth,"C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//30//truth.RData")
#   save(bicrst,file="C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//50//bicrst_",SNR,".RData")
#   save(steprst,file="C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//50//steprst_",SNR,".RData")
  save(glassorst,file=paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/100/glassorst_",SNR,".RData"))
  save(lassorst,file=paste0("/Users/wenxuandeng/GoogleDrive/sucksalt/group_lasso/code/GroupLasso/Final_Simu/100/lassorst_",SNR,".RData"))
#   save(sisrst,file="C://Users//auz5836//Documents//GitHub//GroupLasso//Final_Simu//50//sisrst_",SNR,".RData")
  
}