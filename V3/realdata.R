load("C://Users//auz5836//Desktop//papers//Data//ge.pinfo4Wenxuan07242018.RData")
response<-read.csv("C://Users//auz5836//Desktop//papers//Data//pinfo4Wenxuan07242018.csv")
response<-response[which(response$r1!=""),]
response$r1<-droplevels(response$r1)
response$r1<-as.numeric(response$r1)-1
response$Trt01p<-(as.numeric(response$Trt01p)-1.5)*2
response$ISSSCRN<-(as.numeric(response$ISSSCRN)-1.5)*2
response$STEMCELL<-(as.numeric(response$STEMCELL)-1.5)*2
response$X<-as.character(response$X)
ge<-ge[,match(response$X,colnames(ge))]
ge<-t(ge)
X<-cbind(response$ISSSCRN,response$PRLINED,response$STEMCELL,response$Trt01p,ge,response$Trt01p*ge)
m_X<-3
m_W<-1
m_G<-dim(ge)[2]
y<-response$r1


sol<-FASTA(X,y,f, gradf, g, proxg, x0, tau1, max_iters = 1500, w = 10, 
           backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
           eps_n = 1e-15,m_X,m_W,m_G,m_I,lamb_opt,lamb_opt2,restart=TRUE)