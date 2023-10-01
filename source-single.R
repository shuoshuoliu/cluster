library(dbscan)
#library(irlba)
library(mclust)
library(pracma)
library(Rcpp)
library(parallel)
library(doParallel)
library(foreach)
library(bigstatsr)
library(knn.covertree)

sourceCpp('/storage/work/s/sql5911/Liu-Columbia/sc-clustering/utilSimulDRcluster.cpp')
sourceCpp('/storage/work/s/sql5911/Liu-Columbia/sc-clustering/wpca.cpp')

#svd issue appreas in approxPCA

approxPCA <- function(X, q){ ## speed the computation for initial values.
  n <- nrow(X)
  svdX  <- svd(X, nu=q,nv = q)
  PCs <- svdX$u %*% diag(svdX$d[1:q])
  loadings <- svdX$v
  dX <- PCs %*% t(loadings) - X
  Lam_vec <- colSums(dX^2)/n
  return(list(PCs = PCs, loadings = loadings, Lam_vec = Lam_vec))
}

wpca <- function(X, q, weighted=TRUE){
  if(!is.matrix(X)) stop("wpca: X must be a matrix!")
  if(q< 1) stop("wpca: q must be a positive integer!")
  X <- scale(X, scale=F) # centralize
  out <- wpcaCpp(X, q, weighted)
  return(out)
}

mycluster <- function(Z, G,init="GMM", int.model='EEE'){ #Mclust package: for each sample
  set.seed(111)
  mclus2 <- Mclust(Z, G=G, modelNames=int.model,verbose=FALSE)
  return(mclus2)
}

parfun_int <- function(k, Z, int.model='EEE'){ 
  #for each sample
  #Z: low-dimensional representation
  #k: # of clusters
  mclus2 <- mycluster(Z, k, int.model)
  yveck <- mclus2$classification
  Mu0k <- mclus2$parameters$mean #q*K matrix
  Sigma0k <- mclus2$parameters$variance$sigma #array[,,k]
  return(list(yveck=yveck, Mu0k = Mu0k, Sigma0k=Sigma0k,K=mclus2$G))
}

fast_init=function(K,U,q){
  # needs to determine q first
  mu <- matrix(NA, K,q)
  Sigma <- array(NA, dim = c(q, q, K))
  princ <- wpca(U, q, weighted = TRUE) #approxPCA(U,q) #use PCA for dim redu
  #princ <- approxPCA(U,q)
  W0 <- princ$loadings
  hZ <- princ$PCs
  res=parfun_int(K,hZ)
  
  z=res$yveck #labels
  W=W0
  P=matrix(1/K,K,K)
  Lambda=runif(ncol(U),0.5,2) #is a vector
  for (k in 1:K){
    mu[k,]=res$Mu0k[,k]
    Sigma[,,k]=res$Sigma0k[,,k] 
  }
  return(list(z=z,W=W,P=P,Lambda=Lambda,mu=mu,Sigma=Sigma,K=res$K))
}

fast <- function(q=NULL,K,U,platform="ST",maxIter=100,standardize=FALSE,nG=5,zeta=0.1,
                 distance="rankcor",verbose=FALSE) {
  #Note: q=NULL means choosing q, K=NULL means choosing K,platform means whether spatial data
  #q: lower dim
  #K: number of clusters
  #U: the original data: n*(p+2), 2 is for spatial coordinates
  message("Running Algorithm FAST!")
  
  if(platform=='ST'){
    U=U[,1:(ncol(U)-2)] #if spatial, then take the first p
  }
  
  if (standardize==TRUE){
    U <- scale(U, scale=TRUE)
  }else if (standardize==FALSE){
    U <- scale(U, scale=FALSE)
  }
  
  if(is.null(q)){ #If K is not given, estimate K first
    q=selectq(U)$q
  }
  
  if(length(K)>1){ #If K is a vector, estimate K first
    if (standardize==FALSE){U=scale(U, center = FALSE, scale = TRUE)}
    K_set=K
    K_temp=selectK(U, q,K_set=K_set,zeta=zeta)
    #print(K_temp)
    K=which.min(K_temp)+min(K_set)-1 #plus 1 since k_set starts from 2
    cat("FAST K is", K, sep="\n")
  }
  
  int_values=fast_init(K,U,q)
  z=int_values$z #initial label
  #z_int=int_values$z
  W=int_values$W
  P=int_values$P #z, W, P, Lambda: list
  Lambda=int_values$Lambda
  mu=int_values$mu #mu, Sigma: a list containing matrices, mu[[s]][,k]
  Sigma=int_values$Sigma
  
  iter <- 1
  loglik <- rep(NA,maxIter)
  loglik[1] <- -1000000
  p=ncol(U)
  n=nrow(U)
  
  log_multi=rep(0,n)
  Ux=matrix(NA,n,K)
  Cki_ara=array(NA, dim = c(q, q, K))
  Ex_temp=array(NA, dim = c(n, q, K))
  
  if (platform=='ST'){
    if (distance=="l2"){
      id_matrix=kNN(U[,(p-1):p], k = nG)$id #run kNN by euclidean 
    }else if (distance=="rankcor"){
      id_matrix=find_knn(U[,(p-1):p], k = nG,distance="rankcor")$index
    }
    Y <- t(apply(id_matrix, 1, function(row_indices) {#calculate y in the paper, #count matrix n*K
      table(factor(z[row_indices], levels = 1:K)) #the last two columns are coordinates
    }))
  }
  
  for (iter in 2:maxIter) { #iteration
    lambda_inverse=pinv(diag(Lambda)) 
    WtLW <- t(W) %*% lambda_inverse%*%W #q*q
    XLW <- t(W) %*%lambda_inverse%*%t(U) #W'*Lambda*U: q*n
    for (k in 1:K){
      C_k=WtLW+solve(Sigma[,,k]) #q*q
      Ci=pinv(C_k)#inverse of C_k: q*q
      Cki_ara[,,k]=Ci
      res=multi_det_SkCpp2(X=U, Lam_vec0=Lambda, W0=W, Ck=C_k, Muk=mu[k,], Sigmak=Sigma[,,k])
      if(platform=="ST"){
        for (i in 1:n){
          log_multi[i]=log(dmultinom(x=Y[i,], prob=P[k,]))
        }
      }
      Ux[, k] <- -0.5 * res$logdSk + 0.5 * res$mSk - log_multi # negative loglike
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
    #Update P
    if (platform=="ST"){
      P_temp <- (t(Rik) %*% Y)# * (t(Rik) %*% Y > 0) + 0.01 * (t(Rik) %*% Y <= 0)
      P=t(apply(P_temp,1, function(x) x/sum(x)))
    }
    loglik[iter]=LL
    
    if (verbose) {cat("iter =", iter, ", loglik =", loglik[iter], "\n")}
    if (abs(loglik[iter])>10000000){break}
    if (abs((loglik[iter] - loglik[iter - 1]) / loglik[iter - 1]) < zeta) {break} #used 0.001 in the simulation
  }
  a=loglik[-1]
  a=a[!is.na(a)]
  
  Ezz=calculateX(Rik,Ex_temp)
  #X.low=U%*%W%*%pinv(t(W)%*%W)
  
  return(list(z=z,loglike=a,X.low=Ezz,q=q,K=K))
}


selectK <- function(U, q,K_set,zeta){
  num_cores  <- parallel::detectCores()
  cl  <- parallel::makeCluster(num_cores-1)
  registerDoParallel(cl)
  
  message("FAST: starting parallel computing...",sep="\n")
  #Run forloop in Parallel
  x<- foreach(i=K_set, .combine=c) %dopar% {
    source("/storage/work/s/sql5911/Liu-Columbia/sc-clustering/help.R")
    res=fast_k(q=q,K=i,U=U,zeta=zeta)
    c(res$bic)
  }
  stopCluster(cl) #Stop cluster
  return(x)
}


## select factor number q --------------------------------------------------
selectq <- function(X, qmax=16){
  mnlamjFun <- function(eigvals, j){
    p <- length(eigvals)
    lamj <- eigvals[j]
    Sum <- 0
    for(l in (j+1):p){
      Sum <- Sum + 1/(eigvals[l] - lamj)
    }
    res <- Sum + 1/ ((3*lamj + eigvals[j+1])/4 - lamj)
    return(res/(p-j))
  }
  mtnlamjFun <- function(n, eigvals, j){
    p <- length(eigvals)
    rhojn <-  (p-j)/(n-1)
    res <- -(1-rhojn)/ eigvals[j] + rhojn * mnlamjFun(eigvals, j)
    return(res)
  }
  p <- ncol(X)
  n=nrow(X)
  corMat <- cor(X)
  evalues <- eigen(corMat)$values
  hq1 <- sum(evalues>1+sqrt(p/(n-1)))
  if(hq1 < qmax){
    hq <- hq1
  }else{ # ajdust the eigvalues
    adj.eigvals <- sapply(1:(p-1), function(j) -1/mtnlamjFun(n, evalues, j))
    hq <- sum(adj.eigvals >1) # overselect
  }
  if(hq > qmax || hq < 5) hq <- qmax
  
  propvar <- sum(evalues[1:hq]) / sum(evalues)
  res <- list()
  res$q <- hq
  res$propvar <- sum(evalues[1:hq]) / sum(evalues)
  
  return(res)
}


#### Generate Spatial data with ST platform
gendata_ST <- function(height=30, width=30, platform="ST", p =100, q=10, K=7, umax=15,beta,G=4,homo=FALSE, tau=8, error=5){
  if(q <2) stop("error:gendata_sp::q must be greater than 2!")
  
  n <- height * width # # of cell in each indviduals 
  
  if(platform=="ST"){
    beta= beta
  }else if(platform=='scRNAseq'){
    beta = 0
  }
  ## generate deterministic parameters, fixed after generation
  if (homo==TRUE){
    Lambda <- rep(2,p)
  }else{
    Lambda <- runif(p,0.1,error)
  }
  
  W <- matrix(rt(p*q,2), p, q)
  W <- qr.Q(qr(W))
  mu <- matrix(0, q,  K)
  diagmat = array(0, dim = c(q, q, K))
  if(q > K){
    q1 <- floor(K/2)
    for(j in 1:q1){
      if(j <= (q1/2)) mu[j,j] <- tau
      if(j > (q1/2)) mu[j,j] <- -tau
    }
    mu[(q1+1):q, K] <- -tau
    
  }else if(q <= K){
    for(k in 1:K)
      mu[,k] <- rep(tau/8 *k, q) #
  }
  
  for(k in 1:K){
    tmp2  <- rep(1, q)
    if(k <= K/2){
      tmp2=runif(q,2,umax)#[q] <- 10
    }
    diag(diagmat[,,k]) <- tmp2
  }
  
  Mu <- t(mu)
  Sigma <- diagmat
  
  y <- sampler.mrf(iter = n, sampler = "Gibbs", h = height, w = width, ncolors = K, nei = G, param = beta,initialise = FALSE)
  y <- c(y) + 1
  
  Z <- matrix(0, n, q)
  for(k in 1:K){
    nk <- sum(y==k)
    Z[y==k, ] <- MASS::mvrnorm(nk, Mu[k,], Sigma[,,k])
  }
  Ez <- colMeans(Z)
  Mu <- Mu - matrix(Ez, K, q, byrow=T) # center Z
  X <- Z %*% t(W) + MASS::mvrnorm(n, mu=rep(0,p), Sigma=diag(Lambda))
  
  # make position
  pos <- cbind(rep(1:height, width), rep(1:height, each=width))
  
  counts <- t(X) - min(X)
  p <- ncol(X); n <- nrow(X)
  rownames(counts) <- paste0("gene", seq_len(p))
  colnames(counts) <- paste0("spot", seq_len(n))
  counts <- as.matrix(exp(counts)-1)
  ## Make array coordinates - filled rectangle
  
  if(platform=="ST"){
    cdata <- list()
    cdata$row <- pos[,1]
    cdata$col <- pos[,2]
    cdata <- as.data.frame(do.call(cbind, cdata))
    cdata$imagerow <- cdata$row
    cdata$imagecol <- cdata$col 
    cdata$Z=Z
    row.names(cdata) <- colnames(counts)
    
    seu <-  CreateSeuratObject(counts= counts, meta.data=cdata) #
  }else if(platform=='scRNAseq'){
    seu <-  CreateSeuratObject(counts= counts)
  }else{
    stop("gendata_RNAExp: Unsupported platform \"", platform, "\".")
  }
  
  seu$true_clusters <- y
  return(seu)
}


calculate_top_variances <- function(matrix_data, num_top) {
  # Calculate variances for each column
  variances <- apply(matrix_data, 2, var)
  
  # Sort variances in descending order and get the top 'num_top' indices
  top_indices <- order(variances, decreasing = TRUE)[1:min(num_top, length(variances))]
  
  # Extract the top columns from the original matrix
  top_matrix <- matrix_data[, top_indices]
  
  return(top_matrix)
}