# Cross Validation CI plots
# require(plotrix)
# library(gridExtra)
# library(grid)
# CI<-plotCI(seq(from=0.5,to=5,by=0.25), sol_cv$mean, ui=sol_cv$mean+sqrt(sol_cv$var),
#            li=sol_cv$mean-sqrt(sol_cv$var), xlab="lambda", ylab="Prediction error")
# png("./G.png", height=600, width=1200)
# p<-tableGrob(CI)
# grid.arrange(p)
# dev.off()


# Simulation
K=10


# Signal/Noise Ratio
signal_to_ratio<-c(1,10,100)
#signal_to_ratio<-c(1,10)

# FDR
TP_G<-matrix(0,nrow = t, ncol = length(signal_to_ratio))
TP_I<-matrix(0,nrow = t, ncol = length(signal_to_ratio))

cnt_G<-matrix(0,nrow = t, ncol = length(signal_to_ratio))
cnt_I<-matrix(0,nrow = t, ncol = length(signal_to_ratio))

TP_G<-matrix(0,nrow = t, ncol = length(signal_to_ratio))

FDR_I<-matrix(0,nrow = t, ncol = length(signal_to_ratio)) # False Discovery Rate
FNR_I<-matrix(0,nrow = t, ncol = length(signal_to_ratio)) # False Negative Rate
FDR_G<-matrix(0,nrow = t, ncol = length(signal_to_ratio)) # False Discovery Rate
FNR_G<-matrix(0,nrow = t, ncol = length(signal_to_ratio)) # False Negative Rate

MSE_G<-matrix(0,nrow = t, ncol = length(signal_to_ratio))
MSE_base<-matrix(0,nrow = t, ncol = length(signal_to_ratio))


# Compare FNR/FDR with TP

FNR_record_I<-vector("list", length = length(signal_to_ratio))
FNR_record_G<-vector("list", length = length(signal_to_ratio))
FDR_record_I<-vector("list", length = length(signal_to_ratio))
FDR_record_G<-vector("list", length = length(signal_to_ratio))

TPR_record_I<-vector("list", length = length(signal_to_ratio))
TPR_record_G<-vector("list", length = length(signal_to_ratio))
TDR_record_I<-vector("list", length = length(signal_to_ratio))
TDR_record_G<-vector("list", length = length(signal_to_ratio))

# Strong/Weak Comparison
compare_sign_G_total<-list(NULL,NULL,NULL,NULL)
compare_sign_I_total<-list(NULL,NULL,NULL,NULL)
names(compare_sign_G_total)<-c("0","1","2","3")
names(compare_sign_I_total)<-c("0","1","2","3")

bk<-c(0.1,0.4,0.85)
split_rst<-function(rst,bk,type='G'){
  bk<-bk[order(bk)]
  if(type=='G'){
    rst<-rst[which(rst[,3]!=0),c(3:4)]
    cls<-sapply(bk,function(x) 1*(abs(rst[,1])>x))
    cls<-apply(cls,1,sum)
    rst_split<-split(rst,cls)
  }else if(type=='I'){
    rst<-rst[which(rst[,1]!=0),c(1:2)]
    cls<-sapply(bk,function(x) 1*(abs(rst[,1])>x))
    cls<-apply(cls,1,sum)
    rst_split<-split(rst,cls)
  }
  
  compare_sign<-sapply(rst_split,function(x) 1*(sign(x[,1])==sign(x[,2])))
  if(is.matrix(compare_sign)){
    compare_sign<-split(compare_sign, rep(1:ncol(compare_sign), each = nrow(compare_sign)))
  }
  
  
  
  return(compare_sign)
}

# Record best lambda
lamb_opt<-matrix(0,nrow=t,ncol = length(signal_to_ratio))

# Main/Interaction Ratio
main_inter_ratio<-data.frame("ratio"=NULL,"I"=NULL,"G"=NULL)
MIratio<-function(rst){
  rst<-rst[which(rst[,1]!=0),]
  rst$ratio<-log(abs(rst$True.G1/rst$True.Inter))
  rst$I<-1*(sign(rst$True.Inter)==sign(rst$Est.Inter))
  rst$G<-1*(sign(rst$True.G1)==sign(rst$Est.G1))
  
  return(rst[,c(5:7)])
}


for(i in 1:t){
  
  for(j in seq_along(signal_to_ratio)){
    print(c(i,j))
    
    SNR<-signal_to_ratio[j]
    
    true_beta<-sim_beta2(m_X,m_W,m_G,main_zero,inter_zero,bit=T)
    
    
    X<-sim_X2(m_X,m_W,m_G)
    y0<-X%*%true_beta
    noise<-rnorm(n,sd=1)
    SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
    y<-y0+SNRmtl*noise                         
    
    true_beta<-split_beta2(true_beta,m_X,m_W,m_G,m_I)
    
    #y<-logistic(X,true_beta)
    sol_cv<-opt_lambda2(X,y,f, gradf, g2, proxg2, x0, tau1, max_iters = 500, w = 10,
                        backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5,
                        eps_n = 1e-15,m_X,m_W,m_G,m_I,K=10,n=100,from=0.1,to=1.5,by=0.2,restart=TRUE)
    lamb_candidate<-seq(from=0.1,to=1.5,by=0.2)
    lamb_opt[i,j]<-lamb_candidate[which.min(sol_cv$mean+sol_cv$var)]
    
    
    #lamb_opt<-2
    sol<-FASTA2(X,y,f, gradf, g2, proxg2, x0, tau1, max_iters = 1000, w = 10, 
                backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                eps_n = 1e-15,m_X,m_W,m_G,m_I,lamb_opt[i,j],restart=TRUE)
    # sol<-FASTA2(X,y,f, gradf, g2, proxg2, x0, tau1, max_iters = 1000, w = 10, 
    #             backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
    #             eps_n = 1e-15,m_X,m_W,m_G,m_I,lamb_opt[i,j],restart=TRUE)
    
    
    estbeta<-split_beta2(sol$x,m_X,m_W,m_G,m_I)
    estbeta$G<-estbeta$G*(abs(estbeta$G)>0.05)
    estbeta$I<-estbeta$I*(abs(estbeta$I)>0.05)
    #estbeta$I<-estbeta$I*(abs(estbeta$I)>0.05)
    #num_nonzero<-length(which(true_beta$I!=0))
    #estbeta$I<-estbeta$I*(abs(estbeta$I)>=abs(estbeta$I[order(abs(estbeta$I),decreasing = T)][num_nonzero]))
    #num_nonzero<-length(which(true_beta$G!=0))
    #estbeta$G<-estbeta$G*(abs(estbeta$G)>=abs(estbeta$G[order(abs(estbeta$G),decreasing = T)][num_nonzero]))
    
    
    TP_G[i,j]<-sum(1*(sign(true_beta$G)==sign(estbeta$G)))
    
    TP_I[i,j]<-sum(1*(sign(true_beta$I)==sign(estbeta$I)))
    cnt_I[i,j]<-sum(1*(estbeta$I!=0))
    
    FDR_I[i,j]<-length(setdiff(which(estbeta$I!=0),which(true_beta$I!=0)))/length(which(estbeta$I!=0))
    FNR_I[i,j]<-length(setdiff(which(true_beta$I!=0),which(estbeta$I!=0)))/length(which(true_beta$I!=0))
    FDR_G[i,j]<-length(setdiff(which(estbeta$G!=0),which(true_beta$G!=0)))/length(which(estbeta$G!=0))
    FNR_G[i,j]<-length(setdiff(which(true_beta$G!=0),which(estbeta$G!=0)))/length(which(true_beta$G!=0))
    
    if(FDR_I[i]!=0){
      FDR_record_I[[j]]<-c(FDR_record_I[[j]], estbeta$I[setdiff(which(estbeta$I!=0),which(true_beta$I!=0))])
    }
    if(FDR_G[i]!=0){
      FDR_record_G[[j]]<-c(FDR_record_G[[j]], estbeta$G[setdiff(which(estbeta$G!=0),which(true_beta$G!=0))])
    }
    if(FNR_I[i]!=0){
      FNR_record_I[[j]]<-c(FNR_record_I[[j]], true_beta$I[setdiff(which(true_beta$I!=0),which(estbeta$I!=0))])
    }
    if(FDR_G[i]!=0){
      FNR_record_G[[j]]<-c(FNR_record_G[[j]], true_beta$G[setdiff(which(true_beta$G!=0),which(estbeta$G!=0))])
    }
    if(TP_I[i]!=0){
      TPR_record_I[[j]]<-c(TPR_record_I[[j]], true_beta$I[intersect(which(true_beta$I!=0),which(estbeta$I!=0))])
      TDR_record_I[[j]]<-c(TDR_record_I[[j]], estbeta$I[intersect(which(true_beta$I!=0),which(estbeta$I!=0))])
    }
    if(TP_G[i]!=0){
      TPR_record_G[[j]]<-c(TPR_record_G[[j]], true_beta$G[intersect(which(true_beta$G!=0),which(estbeta$G!=0))])
      TDR_record_G[[j]]<-c(TDR_record_G[[j]], estbeta$G[intersect(which(true_beta$G!=0),which(estbeta$G!=0))])
    }
    
    MSE_G[i,j]<-norm(sol$x[-c(1:(m_X+m_W))],'2')^2/m_G
    MSE_base[i,j]<-norm(sol$x[c(1:(m_X+m_W))],'2')^2/(m_X+m_W)
    
    
    rst_G<-data.frame("True Inter"=true_beta$I,"Est Inter"=estbeta$I, 
                      "True G1"=true_beta$G, "Est G1"=estbeta$G)
    compare_sign_G<-split_rst(rst_G,bk,type='G')
    compare_sign_I<-split_rst(rst_G,bk,type='I')
    
    compare_sign_G_total<-sapply(c("0","1","2","3"),function(x) compare_sign_G_total[[x]]=c(compare_sign_G_total[[x]], compare_sign_G[[x]]))
    compare_sign_I_total<-sapply(c("0","1","2","3"),function(x) compare_sign_I_total[[x]]=c(compare_sign_I_total[[x]], compare_sign_I[[x]]))
    
    main_inter_ratio_temp<-MIratio(rst_G)
    main_inter_ratio<-rbind(main_inter_ratio,main_inter_ratio_temp)
  }
}

# Save results to files
write.table(lamb_opt,"./lamb_opt_200", col.names = F, row.names = F,quote = F)
write.table(FDR_I,"./FDR_I_200", col.names = F, row.names = F,quote = F)
write.table(FDR_G,"./FDR_G_200", col.names = F, row.names = F,quote = F)
write.table(FNR_I,"./FNR_I_200", col.names = F, row.names = F,quote = F)
write.table(FNR_G,"./FNR_G_200", col.names = F, row.names = F,quote = F)
write.table(TP_I,"./TP_I_200", col.names = F, row.names = F,quote = F)
write.table(TP_G,"./TP_G_200", col.names = F, row.names = F,quote = F)
write.table(cnt_I,"./cnt_I_200", col.names = F, row.names = F,quote = F)
write.table(MSE_G,"./MSE_G_200", col.names = F, row.names = F,quote = F)
write.table(MSE_base,"./MSE_base_200", col.names = F, row.names = F,quote = F)
write.table(main_inter_ratio,"./main_inter_ratio_200",row.names = F,quote = F)

save(FDR_record_I,file="./FDR_record_I_200.RData")
save(FDR_record_G,file="./FDR_record_G_200.RData")
save(FNR_record_I,file="./FNR_record_I_200.RData")
save(FNR_record_G,file="./FNR_record_G_200.RData")
save(TDR_record_I,file="./TDR_record_I_200.RData")
save(TDR_record_G,file="./TDR_record_G_200.RData")
save(TPR_record_I,file="./TPR_record_I_200.RData")
save(TPR_record_G,file="./TPR_record_G_200.RData")
save(compare_sign_G_total, file="./compare_sign_G_200.RData")
save(compare_sign_I_total, file="./compare_sign_I_200.RData")



# Load results from files
#load("C:\\Users\\auz5836\\Desktop\\papers\\simulation\\Wenxuan\\GroupLasso\\compare_sign_G.RData")












# FDR_I<-FDR_I[!is.na(FDR_I)]
# FNR_I<-FNR_I[!is.na(FNR_I)]
# FDR_G<-FDR_I[!is.na(FDR_G)]
# FNR_G<-FNR_I[!is.na(FNR_G)]
# 
# rst_X<-data.frame("True Base"=true_beta$X,"Est Base"=estbeta$X)
# rst_W<-data.frame("True Treatment"=true_beta$W,"Est Treatment"=estbeta$W)
# rst_G<-data.frame("True Inter"=true_beta$I,"Est Inter"=estbeta$I, 
#                   "True G1"=true_beta$G, "Est G1"=estbeta$G)
# rst_G
# 
# png("C:\\Users\\auz5836\\Desktop\\papers\\simulation\\Wenxuan\\GroupLasso\\G.png", height=600, width=1200)
# p<-tableGrob(rst_G)
# grid.arrange(p)
# dev.off()

