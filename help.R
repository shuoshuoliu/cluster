source("/storage/work/s/sql5911/Liu-Columbia/sc-clustering/source-single.R")


## For running the algorithm of choosing K
fast_k <- function(q,K,U,maxIter=30,zeta,verbose=FALSE) {
  #q: lower dim
  #K: number of clusters
  #U: the original data without spatial
  
  int_values=fast_init(K,U,q)
  z=int_values$z #initial label
  z_int=int_values$z
  W=int_values$W
  Lambda=int_values$Lambda
  mu=int_values$mu #mu, Sigma: a list containing matrices, mu[[s]][,k]
  Sigma=int_values$Sigma
  
  iter <- 1
  loglik <- rep(NA,maxIter)
  loglik[1] <- -1000000
  p=ncol(U)
  n=nrow(U)
  
  Ux=matrix(NA,n,K)
  Cki_ara=array(NA, dim = c(q, q, K))
  Ex_temp=array(NA, dim = c(n, q, K))
  
  for (iter in 2:maxIter) { #iteration
    #print(iter)
    lambda_inverse=pinv(diag(Lambda)) 
    WtLW <- t(W) %*% lambda_inverse%*%W #q*q
    XLW <- t(W) %*%lambda_inverse%*%t(U) #W'*Lambda*U: q*n
    for (k in 1:K){
      C_k=WtLW+solve(Sigma[,,k]) #q*q
      Ci=pinv(C_k)#inverse of C_k: q*q
      Cki_ara[,,k]=Ci
      
      res=multi_det_SkCpp2(X=U, Lam_vec0=Lambda, W0=W, Ck=C_k, Muk=mu[k,], Sigmak=Sigma[,,k])
      #Ux[, k] <- const-0.5 * res$logdSk + 0.5 * res$mSk - log_multi*(coordinate==TRUE) # negative loglike
      Ux[, k] <- -0.5 * res$logdSk + 0.5 * res$mSk # negative loglike
      Ex_temp[,,k] <- t(Ci %*% (XLW + matrix(rep(solve(Sigma[,,k])%*%mu[k,], n), ncol = n, byrow = FALSE))) #<xi>: n*q
    }
    maxA1 <- apply(-Ux, 1, max)
    Ux <- (-Ux - matrix(rep(maxA1, each = K), nrow = nrow(Ux), ncol = K,byrow=TRUE))
    loglik_more_vec <- rowSums(exp(Ux))
    LL <- sum(log(loglik_more_vec) + maxA1)
    Rik <- exp(Ux) / matrix(rep(loglik_more_vec, each = K), nrow = nrow(Ux), ncol = K,byrow=TRUE)
    z <- apply(Ux, 1, which.max)  # label
    
    # z <- apply(Ux, 1, which.min)  # label
    # LL=-sum(Ux)
    # Rik=t(apply(exp(-Ux), 1, function(row) row / sum(row))) #!!! too large numbers-return NaN
    N <- colSums(Rik)
    
    #update mu
    for (k in 1:K) {mu[k,] <- t(Rik[, k])%*%Ex_temp[,,k]/N[k]} #mu:K*q matrix
    #update Sigma
    Sigma = update_Sigma0(R=Rik, Ez=Ex_temp, Ci=Cki_ara, Mu=mu, N=N)
    #update W
    W = update_W0(X=U, R=Rik, Ez=Ex_temp, Ci=Cki_ara, N=N)
    #update Lambda
    Lambda = as.vector(update_Lam(R=Rik, X=U, W=W, Ez=Ex_temp, Ci=Cki_ara))
    loglik[iter]=LL
    
    if (verbose) {cat("iter =", iter, ", loglik =", loglik[iter], "\n")}
    if (abs(loglik[iter])>10000000){break}
    if (abs((loglik[iter] - loglik[iter - 1]) / loglik[iter - 1]) <zeta) {break}
  }
  
  if(n<1000){
    const=1
  }else(const=2)
  pp=p+p*q + K*(q+q*(q+1)/2) #number of parameters
  bic <-  -2*LL + const*pp* log(log(n)+pp)# HBIC
  #bic <-  -2*LL + pp* log(log(n)+pp)*2.5# HBIC for DLPFC
  #bic=-2*LL +pp*log(n+log(pp+n))#mBIC
  #bic=-2*LL + pp* log(n)* log(log(p+n)) #in DR.SC
  return(list(z=z,z_int=z_int,bic=bic))
}


