load("C://Users//auz5836//Desktop//papers//Data//ge.pinfo4Wenxuan07242018.RData")
response0<-read.csv("C://Users//auz5836//Desktop//papers//Data//pinfo4Wenxuan07242018.csv")
response<-response0
response<-response[which(response$r1!=""),]
response$r1<-droplevels(response$r1)
response$r1<-as.numeric(response$r1)-1
response$Trt01p<-(as.numeric(response$Trt01p)-1.5)*2
response$ISSSCRN<-(as.numeric(response$ISSSCRN)-1.5)*2
response$STEMCELL<-(as.numeric(response$STEMCELL)-1.5)*2
response$X<-as.character(response$X)
ge<-ge[,match(response$X,colnames(ge))]
ge<-t(ge)
ge<-ge[,-322]
ge<-scale(ge)
X<-cbind(response$ISSSCRN,response$PRLINED-2,response$STEMCELL,response$Trt01p,ge,response$Trt01p*ge)
m_X<-3
m_W<-1
m_G<-dim(ge)[2]
m_I<-m_G
p<-m_X+m_W+m_G+m_I
y<-response$r1


# Validation
X_train<-X[1:250,]
y_train<-y[1:250]
X_test<-X[251:392,]
y_test<-y[251:392]

# Parameters
tau1<-.1
lambda<-1
x0<-double(p)
sdErr<-1
lamb_opt<-.3
lamb_opt2<-.2

# Training
sol<-FASTA(X,y,f, gradf, g, proxg, x0, tau1, max_iters = 800, w = 10, 
           backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
           eps_n = 1e-15,m_X,m_W,m_G,m_I,lamb_opt,lamb_opt2,restart=TRUE)
sol<-FASTA(X_train,y_train,f, gradf, g, proxg, x0, tau1, max_iters = 800, w = 10, 
           backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
           eps_n = 1e-15,m_X,m_W,m_G,m_I,lamb_opt,lamb_opt2,restart=TRUE)
estbeta<-split_beta(sol$x,m_X,m_W,m_G,m_I)


#estbeta$G<-10*estbeta$G
#estbeta$I<-10*estbeta$I


rst_G<-data.frame("Est Inter"=estbeta$I, "Est G"=estbeta$G)
rownames(rst_G)<-colnames(ge)
rst_G<-rst_G[order(rst_G$Est.Inter),]

# Testing
prob<-exp(X_test%*%sol$x)/(1+exp(X_test%*%sol$x))
y_predicted<-1*(prob>0.5)
cbind(y_test,y_predicted)
sum(1*(y_test==y_predicted))

prob<-exp(X_train%*%sol$x)/(1+exp(X_train%*%sol$x))
y_predicted<-1*(prob>0.5)
cbind(y_train,y_predicted)
sum(1*(y_train==y_predicted))
