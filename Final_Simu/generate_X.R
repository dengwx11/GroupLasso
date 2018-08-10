
#### different dimensions, p=30
n=100
m_G<-30

for(i in 1:100){
  sigma<-cov_block(m_G,.5,5)
  #sigma<-GenerateCliquesCovariance(10,10,0.8)
  binprob<-runif(m_G)
  x<-sim_X_cate(5,1,m_G,sigma,n,binprob)
  write.table(x,paste0("C:\\Users\\auz5836\\Documents\\GitHub\\GroupLasso\\Final_Simu\\x\\30\\x_",m_G,"_",i), col.names = F, row.names = F,quote = F)
}

#### different dimensions, p=100
n=100
m_G<-100

for(i in 1:100){
  sigma<-cov_block(m_G,.5,10)
  #sigma<-GenerateCliquesCovariance(10,10,0.8)
  binprob<-runif(m_G)
  x<-sim_X_cate(5,1,m_G,sigma,n,binprob)
  write.table(x,paste0("C:\\Users\\auz5836\\Documents\\GitHub\\GroupLasso\\Final_Simu\\x\\",m_G,"\\x_",m_G,"_",i), col.names = F, row.names = F,quote = F)
}

#### different dimensions, p=200
n=100
m_G<-200

for(i in 1:100){
  sigma<-cov_block(m_G,.5,30)
  #sigma<-GenerateCliquesCovariance(10,10,0.8)
  binprob<-runif(m_G)
  x<-sim_X_cate(5,1,m_G,sigma,n,binprob)
  write.table(x,paste0("C:\\Users\\auz5836\\Documents\\GitHub\\GroupLasso\\Final_Simu\\x\\",m_G,"\\x_",m_G,"_",i), col.names = F, row.names = F,quote = F)
}

#### different dimensions, p=1000
n=100
m_G<-1000

for(i in 1:100){
  sigma<-cov_block(m_G,.5,10)
  #sigma<-GenerateCliquesCovariance(10,10,0.8)
  binprob<-runif(m_G)
  x<-sim_X_cate(5,1,m_G,sigma,n,binprob)
  write.table(x,paste0("C:\\Users\\auz5836\\Documents\\GitHub\\GroupLasso\\Final_Simu\\x\\",m_G,"\\x_",m_G,"_",i), col.names = F, row.names = F,quote = F)
}