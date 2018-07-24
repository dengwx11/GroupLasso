
#SNR<-10
true_beta<-sim_beta_const(m_X,m_W,m_G,main_nonzero,inter_nonzero,c(10))


X<-sim_X(m_X,m_W,m_G)
y<-logistic(X,true_beta)
#y0<-X%*%true_beta
#noise<-rnorm(n,sd=1)
#SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
#y<-y0+SNRmtl*noise                         

true_beta<-split_beta(true_beta,m_X,m_W,m_G,m_I)

lamb_opt<-.7
lamb_opt2<-.2
sol<-FASTA(X,y,f, gradf, g, proxg, x0, tau1, max_iters = 800, w = 10, 
            backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
            eps_n = 1e-15,m_X,m_W,m_G,m_I,lamb_opt,lamb_opt2,restart=TRUE)
estbeta<-split_beta(sol$x,m_X,m_W,m_G,m_I)

f0(unlist(true_beta),X,y)
f0(unlist(estbeta),X,y)

estbeta$G<-10*estbeta$G
estbeta$I<-10*estbeta$I
#estbeta$G<-estbeta$G*(abs(estbeta$G)>1)
#estbeta$I<-estbeta$I*(abs(estbeta$I)>1)

rst_X<-data.frame("True Base"=true_beta$X,"Est Base"=estbeta$X)
rst_W<-data.frame("True Treatment"=true_beta$W,"Est Treatment"=estbeta$W)
rst_G<-data.frame("True Inter"=true_beta$I,"Est Inter"=estbeta$I,
                  "True G"=true_beta$G, "Est G"=estbeta$G)
rst_G

cor<-cor(X[,-c(1:(m_X+m_W))])
pair<-data.frame("row"=which(cor(X[,-c(1:(m_X+m_W))])>0.2)%%dim(cor)[2],"col"=ceiling(which(cor(X[,-c(1:(m_X+m_W))])>0.2)/dim(cor)[2]))

sum(1*(rst_G$True.G!=0))
sum(1*(rst_G$Est.G!=0))
length(intersect(which(rst_G$True.G1!=0),which(rst_G$Est.G1!=0)))
sum(1*(rst_G$True.Inter!=0))
sum(1*(rst_G$Est.Inter!=0))
length(intersect(which(rst_G$True.Inter!=0),which(rst_G$Est.Inter!=0)))

# Testing
prob<-exp(X_test%*%sol$x)/(1+exp(X_test%*%sol$x))
y_predicted<-1*(prob>0.5)
cbind(y_test,y_predicted)
sum(1*(y_test==y_predicted))

prob<-exp(X_train%*%sol$x)/(1+exp(X_train%*%sol$x))
y_predicted<-1*(prob>0.5)
cbind(y_train,y_predicted)
sum(1*(y_train==y_predicted))