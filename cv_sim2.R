# Cross Validation CI plots
require(plotrix)
library(gridExtra)
library(grid)
CI<-plotCI(seq(from=0.5,to=5,by=0.25), sol_cv$mean, ui=sol_cv$mean+sqrt(sol_cv$var),
       li=sol_cv$mean-sqrt(sol_cv$var), xlab="lambda", ylab="Prediction error")
png("./G.png", height=600, width=1200)
p<-tableGrob(CI)
grid.arrange(p)
dev.off()


# Simulation
K=10

TP_G1<-rep(0,t)
TP_G2<-rep(0,t)
TP_I<-rep(0,t)

cnt_G1<-rep(0,t)
cnt_G2<-rep(0,t)
cnt_I<-rep(0,t)

TP_G<-rep(0,t)

FDR<-rep(0,t) # False Discovery Rate
FNR<-rep(0,t) # False Omission Rate


for(i in 1:t){
  print(i)
  
  true_beta<-sim_beta2(m_X,m_W,m_G,main_zero,inter_zero,bit=T)
  
  
  X<-sim_X2(m_X,m_W,m_G)
  y<-X%*%true_beta+rnorm(n,sd=sdErr)                          
  
  true_beta<-split_beta2(true_beta,m_X,m_W,m_G,m_I)

  #y<-logistic(X,true_beta)
  sol_cv<-opt_lambda2(X,y,f, gradf, g2, proxg2, x0, tau1, max_iters = 500, w = 10, 
                     backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                     eps_n = 1e-15,m_X,m_W,m_G,m_I,K=10,n=100,from=0.5,to=5,by=0.25,restart=TRUE)
  lamb_candidate<-seq(from=0.5,to=2,by=0.25)
  lamb_opt<-lamb_candidate[which.min(sol_cv$mean+sol_cv$var)]
  
  
  lamb_opt<-0.75
  sol<-FASTA2(X,y,f, gradf, g2, proxg2, x0, tau1, max_iters = 1000, w = 10, 
               backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
               eps_n = 1e-15,m_X,m_W,m_G,m_I,lamb_opt,restart=TRUE)
  
  
  estbeta<-split_beta2(sol$x,m_X,m_W,m_G,m_I)
  estbeta$G1<-estbeta$G1*(abs(estbeta$G1)>0.05)
  estbeta$G2<-estbeta$G2*(abs(estbeta$G2)>0.05)
  #estbeta$I<-estbeta$I*(abs(estbeta$I)>0.05)
  num_nonzero<-length(which(true_beta$I!=0))
  estbeta$I<-estbeta$I*(abs(estbeta$I)>=abs(estbeta$I[order(abs(estbeta$I),decreasing = T)][num_nonzero]))
  estbeta$G<-estbeta$G1+estbeta$G2
  num_nonzero<-length(which(true_beta$G!=0))
  estbeta$G<-estbeta$G*(abs(estbeta$G)>=abs(estbeta$G[order(abs(estbeta$G),decreasing = T)][num_nonzero]))
  
  
  
  TP_G1[i]<-sum(1*(sign(true_beta$G1)==sign(estbeta$G1)))
  cnt_G1[i]<-sum(1*(estbeta$G1!=0))
  
  TP_G2[i]<-sum(1*(sign(true_beta$G2)==sign(estbeta$G2)))
  cnt_G2[i]<-sum(1*(estbeta$G2!=0))
  
  TP_G[i]<-sum(1*(sign(true_beta$G)==sign(estbeta$G)))
  
  TP_I[i]<-sum(1*(sign(true_beta$I)==sign(estbeta$I)))
  cnt_I[i]<-sum(1*(estbeta$I!=0))
  
  FDR[i]<-length(setdiff(which(estbeta$I!=0),which(true_beta$I!=0)))/length(which(estbeta$I!=0))
  FNR[i]<-length(setdiff(which(true_beta$I!=0),which(estbeta$I!=0)))/length(which(true_beta$I!=0))
}

FDR<-FDR[!is.na(FDR)]
FNR<-FNR[!is.na(FNR)]

rst_X<-data.frame("True Base"=true_beta$X,"Est Base"=estbeta$X)
rst_W<-data.frame("True Treatment"=true_beta$W,"Est Treatment"=estbeta$W)
rst_G<-data.frame("True Inter"=true_beta$I,"Est Inter"=estbeta$I, 
                  "True G1"=true_beta$G, "Est G1"=estbeta$G)
rst_G

png("C:\\Users\\auz5836\\Desktop\\papers\\simulation\\Wenxuan\\GroupLasso\\G.png", height=600, width=1200)
p<-tableGrob(rst_G)
grid.arrange(p)
dev.off()

