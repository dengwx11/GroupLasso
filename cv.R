# Cross Validation
divide<-function(K,n){
  i.mix = sample(1:n)
  folds = vector(mode="list",length=K)
  temp<-split(c(1:n),1:K)
  for(i in 1:K){
    folds[[i]]<-i.mix[temp[[i]]]
  }
  return(folds)
}


cv.FASTA<-function(X,y,f, gradf, g, proxg, x0, tau1, max_iters = 500, w = 10, 
                   backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                   eps_n = 1e-15,m_X,m_W,m_G1,m_G2,m_I,lambda1,lambda2,K=10,n=100,restart){
  
  folds<-divide(K,n)
  
  sol_cv<-list()
  TestErr<-rep(0,K)
  
  for(k in 1:K){
    
    # Generating training and test datasets for cross validation
    test_X<-X[folds[[k]],]
    train_X<-X[setdiff(c(1:n),folds[[k]]),]
    test_y<-y[folds[[k]]]
    train_y<-y[setdiff(c(1:n),folds[[k]])]
    
    # Training model on training dataset
    sol_cv[[k]]<-FASTA(train_X,train_y,f, gradf, g, proxg, x0, tau1, max_iters = 500, w = 10, 
               backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
               eps_n = 1e-15,m_X,m_W,m_G1,m_G2,m_I,lambda1,lambda2,restart)
    x0<-sol_cv[[k]]$x
    
    # Test on the validation group
    TestErr[k]<-f(sol_cv[[k]]$x,test_X,test_y)/length(test_y)

    
  }
  return(list(Err=TestErr,start=x0))
}


opt_lambda<-function(X,y,f, gradf, g, proxg, x0, tau1, max_iters = 500, w = 10, 
                     backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                     eps_n = 1e-15,m_X,m_W,m_G1,m_G2,m_I,K=10,n=100,from=0,to=2,by=0.5,restart){
  
  lamb_candidate<-seq(to,from,-by)
  
  TestErr<-rep(0,length(lamb_candidate))
  VarErr<-rep(0,length(lamb_candidate))
  
  for(i in seq_along(lamb_candidate)){
    rst<-cv.FASTA(X,y,f, gradf, g, proxg, x0, tau1, max_iters = 500, w = 10, 
                      backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                      eps_n = 1e-15,m_X,m_W,m_G1,m_G2,m_I,lamb_candidate[i],lamb_candidate[i],K,n,restart=TRUE)
    cv.Err<-rst$Err
    x0<-rst$start
    TestErr[i]<-mean(cv.Err)
    VarErr[i]<-var(cv.Err)
    print(paste("lambda=",lamb_candidate[i]))
    print(cv.Err)
  }
  
  return(list(mean=rev(TestErr),var=rev(VarErr)))
}



