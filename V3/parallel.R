args=commandArgs(TRUE)
ncores=1

### Simulation
SNR<-10
true_beta<-sim_beta_const(m_X,m_W,m_G,main_nonzero,inter_nonzero,both_nonzero,c(1,3),heir = F)


X<-sim_X(m_X,m_W,m_G,n)
y0<-X%*%true_beta
noise<-rnorm(n,sd=1)
SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
y<-y0+SNRmtl*noise                         

true_beta<-split_beta(true_beta,m_X,m_W,m_G,m_I)


lamb_candidate<-c(1,2,3,4,5,7,10,13,17)
lamb_candidate2<-c(0.5,1,1.5,2,3)

### CV on group lasso
sol_cv<-opt_lambda(X,y,f, gradf, g, proxg, x0, tau1, max_iters = 100, w = 10, 
                   backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                   eps_n = 1e-15,m_X,m_W,m_G,m_I,K=10,n=100, lamb_candidate, lamb_candidate2,restart=TRUE)
lamb_loc<-which(sol_cv$mean+sol_cv$var == min(sol_cv$mean+sol_cv$var), arr.ind = TRUE)


lamb_opt<-lamb_candidate[lamb_loc[1]]
lamb_opt2<-lamb_candidate2[lamb_loc[2]]
sol_group<-FASTA(X,y,f, gradf, g, proxg, x0, tau1, max_iters = 1000, w = 10, 
           backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
           eps_n = 1e-15,m_X,m_W,m_G,m_I,lamb_opt,lamb_opt2,restart=TRUE)
estbeta_group<-split_beta(sol_group$x,m_X,m_W,m_G,m_I)

### CV on standard Lasso

lamb_candidate<-c(5,10,15,20,25,30,35,40)
lamb_candidate2<-c(2:6)

sol_cv<-opt_lambda(X,y,f, gradf, glasso, proxglasso, x0, tau1, max_iters = 100, w = 10, 
                   backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                   eps_n = 1e-15,m_X,m_W,m_G,m_I,K=10,n=100, lamb_candidate, lamb_candidate2,restart=TRUE)
lamb_loc<-which(sol_cv$mean+sol_cv$var == min(sol_cv$mean+sol_cv$var), arr.ind = TRUE)


lamb_opt<-lamb_candidate[lamb_loc[1]]
lamb_opt2<-lamb_candidate2[lamb_loc[2]]
sol_std<-FASTA(X,y,f, gradf, glasso, proxglasso, x0, tau1, max_iters = 1000, w = 10, 
                 backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                 eps_n = 1e-15,m_X,m_W,m_G,m_I,lamb_opt,lamb_opt2,restart=TRUE)
estbeta_std<-split_beta(sol_std$x,m_X,m_W,m_G,m_I)

### Save

saveRDS(estbeta_group,file=paste("./v3/rst/group_est",args,".rds"))
saveRDS(estbeta_std,file=paste("./v3/rst/std_est",args,".rds"))
saveRDS(true_beta,file=paste("./v3/rst/true",args,".rds"))