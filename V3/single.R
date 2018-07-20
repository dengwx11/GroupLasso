
SNR<-100
true_beta<-sim_beta(m_X,m_W,m_G,main_zero,inter_zero,bit=T)


X<-sim_X(m_X,m_W,m_G)
y0<-X%*%true_beta
noise<-rnorm(n,sd=1)
SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
y<-y0+SNRmtl*noise                         

true_beta<-split_beta(true_beta,m_X,m_W,m_G,m_I)

lamb_opt<-5
lamb_opt2<-2
sol<-FASTA(X,y,f, gradf, g, proxg, x0, tau1, max_iters = 500, w = 10, 
            backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
            eps_n = 1e-15,m_X,m_W,m_G,m_I,lamb_opt,lamb_opt2,restart=TRUE)
estbeta<-split_beta(sol$x,m_X,m_W,m_G,m_I)

f0(unlist(true_beta),X,y)
f0(unlist(estbeta),X,y)

estbeta$G<-(1+lamb_opt2)*estbeta$G
estbeta$I<-(1+lamb_opt2)*estbeta$I
estbeta$G<-estbeta$G*(abs(estbeta$G)>0.1)
estbeta$I<-estbeta$I*(abs(estbeta$I)>0.1)

rst_X<-data.frame("True Base"=true_beta$X,"Est Base"=estbeta$X)
rst_W<-data.frame("True Treatment"=true_beta$W,"Est Treatment"=estbeta$W)
rst_G<-data.frame("True Inter"=true_beta$I,"Est Inter"=estbeta$I,
                  "True G1"=true_beta$G, "Est G1"=estbeta$G)
rst_G

cor<-cor(X[,-c(1:(m_X+m_W))])
pair<-data.frame("row"=which(cor(X[,-c(1:(m_X+m_W))])>0.2)%%dim(cor)[2],"col"=ceiling(which(cor(X[,-c(1:(m_X+m_W))])>0.2)/dim(cor)[2]))

