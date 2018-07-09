FASTA2<-function (X,y,f, gradf, g2, proxg2, x0, tau1, max_iters = 100, w = 10, 
                 backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                 eps_n = 1e-15,m_X,m_W,m_G,m_I,lambda1,lambda2,restart=TRUE) 
{
  residual <- double(max_iters)
  normalizedResid <- double(max_iters)
  taus <- double(max_iters)
  fVals <- double(max_iters)
  objective <- double(max_iters + 1)
  totalBacktracks <- 0
  backtrackCount <- 0
  x1 <- x0
  d1 <- x1
  f1 <- f(d1,X,y)
  fVals[1] <- f1
  gradf1 <- gradf(d1,X,y)
  if (recordIterates) {
    iterates <- matrix(0, length(x0), max_iters + 1)
    iterates[, 1] <- x1
  }else {
    iterates <- NULL
  }
  maxResidual <- -Inf
  minObjectiveValue <- Inf
  objective[1] <- f1 + g2(X,x0,m_X,m_W,m_G,m_I,lambda1,lambda2)
  t0=1
  y0<-x0
  for (i in 1:max_iters) {
    x0 <- x1
    gradf0 <- matrix(gradf1)
    tau0 <- tau1
    x1hat <- x0 - tau0 * c(gradf0)
    x1 <- proxg2(X,x1hat,m_X,m_W,m_G,m_I, tau0, lambda1,lambda2)
    Dx <- matrix(x1 - x0)
    d1 <- x1
    f1 <- f(d1,X,y)
    if (backtrack) {
      M <- max(fVals[max(i - w, 1):max(i - 1, 1)])
      backtrackCount <- 0
      prop <- (f1 - 1e-12 > M + t(Dx) %*% gradf0 + 0.5 * 
                 (norm(Dx, "f")^2)/tau0) && (backtrackCount < 
                                               20)
      while (prop) {
        tau0 <- tau0 * stepsizeShrink
        x1hat <- x0 - tau0 * c(gradf0)
        x1 <- proxg2(X,x1hat,m_X,m_W,m_G,m_I, tau0, lambda1, lambda2)
        d1 <- x1
        f1 <- f(d1,X,y)
        Dx <- matrix(x1 - x0)
        backtrackCount <- backtrackCount + 1
        prop <- (f1 - 1e-12 > M + t(Dx) %*% gradf0 + 
                   0.5 * (norm(Dx, "f")^2)/tau0) && (backtrackCount < 
                                                       20)
      }
      totalBacktracks <- totalBacktracks + backtrackCount
    }
    
    if(restart){
      y1<-x1
      check<-t(x1-y1)%*%(y1-y0)
      if(check>0){
        t1=1
      }else{
        t1<-(1+(1+4*t0^2)^0.5)/2
      }
      x1<-x1+(t0-1)/t1*(y1-y0)
      y0<-y1
    }
    
    taus[i] <- tau0
    residual[i] <- norm(Dx, "f")/tau0
    maxResidual <- max(maxResidual, residual[i])
    normalizer <- max(norm(gradf0, "f"), norm(as.matrix(x1 - 
                                                          x1hat), "f")/tau0) + eps_n
    normalizedResid[i] <- residual[i]/normalizer
    fVals[i] <- f1
    objective[i + 1] <- f1 + g2(X,x1,m_X,m_W,m_G,m_I,lambda1,lambda2)
    newObjectiveValue <- objective[i + 1]
    if (recordIterates) {
      iterates[, i + 1] <- x1
    }
    if (newObjectiveValue < minObjectiveValue) {
      bestObjectiveIterate <- x1
      minObjectiveValue <- min(minObjectiveValue, newObjectiveValue)
    }
    
    gradf1 <- gradf(d1,X,y)
    Dg <- matrix(gradf1 + (x1hat - x0)/tau0)
    dotprod <- t(Dx) %*% Dg
    tau_s <- norm(Dx, "f")^2/dotprod
    tau_m <- dotprod/norm(Dg, "f")^2
    tau_m <- max(tau_m, 0)
    if (abs(dotprod) < 1e-15) 
      break
    if (2 * tau_m > tau_s) {
      tau1 <- tau_m
    }else {
      tau1 <- tau_s - 0.5 * tau_m
    }
    if ((tau1 <= 0) || is.infinite(tau1) || is.nan(tau1)) {
      tau1 <- tau0 * 1.5
    }
    
    
  }
  if (recordIterates) {
    iterates <- iterates[, 1:(i + 1), drop = FALSE]
  }
  return(list(x = bestObjectiveIterate, objective = objective[1:(i + 
                                                                   1)], fVals = fVals[1:i], totalBacktracks = totalBacktracks, 
              residual = residual[1:i], taus = taus[1:i], iterates = iterates))
}


