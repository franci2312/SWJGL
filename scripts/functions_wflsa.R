## Functions ##
soft <- JGL:::soft

## Cross-validation for WFLSA ##
wflsa_cv <- function(data_list , lam1_vec , lam2_vec , nfold = 5, lambda_soft = 0.05,
                     seed = 1 , cov = F , weight_list = NULL , penalize_diagonal = F){
  
  k <- length(data_list)
  n <- nobs_fold <- vector(mode = 'list' , length = k)
  for(kk in 1:k){
    n[[kk]] <- dim(data_list[[kk]])[1]
    nobs_fold[[kk]] <- n[[kk]]/nfold
  }
  
  lk <- array(dim = c(nfold , length(lam1_vec) , length(lam2_vec)))
  
  set.seed(seed)
  
  for(kk in 1:k){
    shuffle_idx <- sample(1:n[[kk]], n[[kk]] , replace = F)
    data_list[[kk]] <- data_list[[kk]][shuffle_idx,]
  }
  
  for(hh in 1:nfold){
    
    cat("\n",hh,"fold")
    
    idx <- vector(mode = 'list' , length = k)
    obs_test_hh <- obs_train_hh <- mean_train_hh <-obs_test_hh_cen <-obs_train_hh_cen <- list()
    
    for(kk in 1:k){
      idx[[kk]] <- ((hh-1)*nobs_fold[[kk]]+1):(hh*nobs_fold[[kk]])
      obs_test_hh[[kk]] <- data_list[[kk]][idx[[kk]],]
      obs_train_hh[[kk]] <- data_list[[kk]][-idx[[kk]],]
      mean_train_hh[[kk]] <- colMeans(obs_train_hh[[kk]])
      obs_test_hh_cen[[kk]] <-  scale(obs_test_hh[[kk]]  , center = mean_train_hh[[kk]]  , scale = F)
      obs_train_hh_cen[[kk]] <-  scale(obs_train_hh[[kk]] , center = T , scale = F)
    }
    
    A <- vector(mode = 'list' , length = k)
    for(kk in 1:k){
      cc <- cov(obs_train_hh[[kk]])
      if(cov){
        A[[kk]] <- cc
      }else{
        A[[kk]] <- solve(soft(cc , lam = lambda_soft , penalize.diagonal =  penalize_diagonal))  
      }
    }
    
    models <- vector(mode = 'list' , length = length(lam1_vec))
    for(i in 1:length(lam1_vec)){
      models[[i]] <-vector(mode = 'list' , length = length(lam2_vec))
    }
    
    for(i in 1:length(lam1_vec)){
      for(j in 1:length(lam2_vec)){
        lam1 <- lam1_vec[i]
        lam2 <- lam2_vec[j]
        models[[i]][[j]] <- wflsa(A = A , L = 2 , lam1 = lam1 , lam2 = lam2 ,
                                  penalize_diagonal = penalize_diagonal, weight_list = weight_list)
      }
    }
    
    S_test <- vector(mode = 'list' , length = k)
    for(kk in 1:k){
      S_test[[kk]] <- cov(obs_test_hh_cen[[kk]])
    }
    
    for(i in 1:length(lam1_vec)){
      for(j in 1:length(lam2_vec)){
        
        ll <- c()
        
        for(kk in 1:k){
          if(cov){
            part1 <- -log(det(solve(models[[i]][[j]][[kk]])))
            part2 <- sum(diag(S_test[[kk]] %*% solve(models[[i]][[j]][[kk]])))
            ll[kk] <- dim(obs_test_hh)[1]*(part1+part2)
          }else{
            part1 <- -log(det(models[[i]][[j]][[kk]]))
            part2 <- sum(diag(S_test[[kk]] %*% models[[i]][[j]][[kk]]))
            ll[kk] <- dim(obs_test_hh)[1]*(part1+part2) 
          }
        }
        lk[hh , i , j] <- sum(ll)
      }
    }
    
  }
  lk_mean <- apply(lk , c(2,3),mean , na.rm = T)
  # sd <- apply(lk , 2 , sd)/sqrt(nfold)
  rownames(lk_mean) <- as.character(round(lam1_vec,2))
  colnames(lk_mean) <- as.character(round(lam2_vec,2))
  
  return(list(lk = lk_mean,
              lam1_vec = lam1_vec,
              lam2_vec = lam2_vec,
              best_lambdas = c(lam1_vec[which(lk_mean == min(lk_mean),arr.ind = T)[1,1]],
                               lam2_vec[which(lk_mean == min(lk_mean),arr.ind = T)[1,2]])))  
  
  
}

## Weighted Fused Lasso Signal Approximator (WFLSA) ##

wflsa <- function (data_list, lam1, lam2,lambda_soft = 0.01, penalize_diagonal  = F, weight_list = NULL, cov = F) 
{
  L <- 2
  k <- length(data_list)
  A <-cc <-  vector(mode = 'list' , length = k)
  for(kk in 1:k){
    cc[[kk]] <- cov(data_list[[kk]])
    if(cov){
      A[[kk]] <- cc[[kk]]
    }else{
      A[[kk]] <- solve(soft(cc[[kk]] , lam = lambda_soft , penalize.diagonal =  penalize_diagonal))  
    }
  }
  
  
  trueA = A
  if (is.matrix(A[[1]])) {
    p = dim(A[[1]])[1]
  }
  if (is.vector(A[[1]])) {
    p = length(A[[1]])
  }
  K = length(A)
  
  if(is.null(weight_list)){
    weight_list <- list()
    for(k in 1:K){
      weight_list[[k]] <- rep(1,K)
    }
  }
  
  X = list()
  for (k in 1:K) {
    X[[k]] = A[[1]] * NA
  }
  if (is.matrix(A[[1]])) {
    fusions = array(FALSE, dim = c(K, K, p, p))
  }
  if (is.vector(A[[1]])) {
    fusions = array(FALSE, dim = c(K, K, p, 1))
  }
  newc = list()
  for (k in 1:K) {
    others = setdiff(1:K, k)
    others.smaller.k = 1:(k - 1)
    newc[[k]] = A[[1]] * 0
    for (o in others) {
      newc[[k]] = newc[[k]] + weight_list[[k]][o]*((A[[o]] - A[[k]] < -1e-04) - 
                                                     (A[[o]] - A[[k]] > 1e-04))
    }
  }
  for (iter in 1:(K - 1)) {
    ordermats = list()
    for (k in 1:K) {
      others = setdiff(1:K, k)
      others.smaller.k = 1:(k - 1)
      ordermats[[k]] = A[[1]] * 0
      for (o in others) {
        ordermats[[k]] = ordermats[[k]] + (A[[k]] - 
                                             A[[o]] > 1e-04)
      }
      if (k > 1) {
        for (o in others.smaller.k) {
          ordermats[[k]] = ordermats[[k]] + (abs(A[[o]] - 
                                                   A[[k]]) < 1e-04)
        }
      }
      ordermats[[k]] = ordermats[[k]] + 1
    }
    betas.g = list()
    for (k in 1:K) {
      betas.g[[k]] = A[[k]] - lam2/L * newc[[k]]
    }
    new.ordermats = list()
    for (k in 1:K) {
      others = setdiff(1:K, k)
      others.smaller.k = 1:(k - 1)
      new.ordermats[[k]] = A[[1]] * 0
      for (o in others) {
        new.ordermats[[k]] = new.ordermats[[k]] + (betas.g[[k]] - 
                                                     betas.g[[o]] > 1e-04)
      }
      if (k > 1) {
        for (o in others.smaller.k) {
          new.ordermats[[k]] = new.ordermats[[k]] + 
            (abs(betas.g[[o]] - betas.g[[k]]) < 1e-04)
        }
      }
      new.ordermats[[k]] = new.ordermats[[k]] + 1
    }
    for (k in 1:K) {
      for (kp in 1:K) {
        fusions[k, kp, , ] = fusions[k, kp, , ] + ((ordermats[[k]] - 
                                                      1 == ordermats[[kp]]) & (new.ordermats[[k]] < 
                                                                                 new.ordermats[[kp]])) + ((ordermats[[k]] + 
                                                                                                             1 == ordermats[[kp]]) & (new.ordermats[[k]] > 
                                                                                                                                        new.ordermats[[kp]])) + (abs(A[[k]] - A[[kp]]) < 
                                                                                                                                                                   1e-04)
        fusions = (fusions > 0) * 1
      }
    }
    for (k in 1:K) {
      for (kp in 1:K) {
        others = setdiff(1:K, c(k, kp))
        for (o in others) {
          bothfused = fusions[k, o, , ] & fusions[kp, 
                                                  o, , ]
          fusions[k, kp, , ] = fusions[k, kp, , ] | 
            bothfused
        }
      }
    }
    for (k in 1:K) {
      others = setdiff(1:K, k)
      fusemean = trueA[[k]]
      denom = A[[1]] * 0 + 1
      for (o in others) {
        fusemean = fusemean + fusions[k, o, , ] * trueA[[o]]
        denom = denom + fusions[k, o, , ]
      }
      A[[k]] = fusemean/denom
    }
    newc = list()
    for (k in 1:K) {
      others = setdiff(1:K, k)
      others.smaller.k = 1:(k - 1)
      newc[[k]] = A[[1]] * 0
      for (o in others) {
        newc[[k]] = newc[[k]] + weight_list[[k]][o]*((A[[o]] - A[[k]] < -1e-04) - 
                                                       (A[[o]] - A[[k]] > 1e-04))
      }
    }
  }
  for (k in 1:K) {
    betas.g[[k]] = A[[k]] - lam2/L * newc[[k]]
  }
  for (k in 1:K) {
    X[[k]] = soft(betas.g[[k]], lam = lam1/L, penalize.diagonal = penalize_diagonal)
  }
  return(X)
}

## BIC and EBIC for WFLSA
wflsa_bic <-  function(data_list , lam1_vec , lam2_vec , lambda_soft = 0.05,
                       cov = F , weight_list = NULL , penalize_diagonal = F, parallel = F){
  
  k <- length(data_list)
  n  <- vector(mode = 'list' , length = k)
  for(kk in 1:k){
    n[[kk]] <- dim(data_list[[kk]])[1]
  }
  
  bic <- ebic <- matrix(10000,nrow = length(lam1_vec) , ncol = length(lam2_vec))
  
  A <-cc <-  vector(mode = 'list' , length = k)
  for(kk in 1:k){
    cc[[kk]] <- cov(data_list[[kk]])
    if(cov){
      A[[kk]] <- cc[[kk]]
    }else{
      A[[kk]] <- solve(soft(cc[[kk]] , lam = lambda_soft , penalize.diagonal =  penalize_diagonal))  
    }
  }
  
  
  if(!parallel){
    models <- vector(mode = 'list' , length = length(lam1_vec))
    for(i in 1:length(lam1_vec)){
      models[[i]] <-vector(mode = 'list' , length = length(lam2_vec))
    }
    
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = length(lam1_vec), # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar
    
    for(i in 1:length(lam1_vec)){
      setTxtProgressBar(pb, i)
      
      for(j in 1:length(lam2_vec)){
        # cat("\n","lam1 =",lam1_vec[i] , "\t", "lam2 =",lam2_vec[j])
        lam1 <- lam1_vec[i]
        lam2 <- lam2_vec[j]
        models[[i]][[j]] <- wflsa(data_list = data_list ,lambda_soft = lambda_soft , lam1 = lam1 , lam2 = lam2 ,
                                  penalize_diagonal = penalize_diagonal, weight_list = weight_list,cov = cov)
        
        for(kk in 1:k){
          if(cov){
            # part1 <- -log(det(solve(models[[i]][[j]][[kk]])))
            # part2 <- sum(diag(S_test[[kk]] %*% solve(models[[i]][[j]][[kk]])))
            # # bic[i,j] <- dim(obs_test_hh)[1]*(part1+part2)
            # ebic[i,j] <-  extendedBIC(gamma = 0 , omegahat = solve(models[[i]][[j]][[kk]]) , S = cc[[kk]], n = n[[kk]])
            # 
          }else{
            part1 <- -log(det(models[[i]][[j]][[kk]]))
            part2 <- sum(diag(cc[[kk]] %*% models[[i]][[j]][[kk]]))
            ebic[i,j] <-  extendedBIC(gamma = 0 , omegahat = models[[i]][[j]][[kk]] , S = cc[[kk]], n = n[[kk]])
          }
        }
      }
    }
    close(pb)
  }  else{
    calculate_models <- function(i){
      models <- vector(mode = 'list', length = length(lam2_vec))
      for(j in 1:length(lam2_vec)){
        # cat("\n","lam1 =",lam1_vec[i] , "\t", "lam2 =",lam2_vec[j])
        lam1 <- lam1_vec[i]
        lam2 <- lam2_vec[j]
        models[[j]] <- wflsa(data_list = data_list ,lambda_soft = lambda_soft , lam1 = lam1 , lam2 = lam2 ,
                             penalize_diagonal = penalize_diagonal, weight_list = weight_list,cov = cov)
        
      }
      return(models)
    }
    require(parallel)
    cores <- detectCores() -1
    models <- mclapply(1:length(lam1_vec) , calculate_models)
    for(i in 1:length(lam1_vec)){
      for(j in 1:length(lam2_vec)){
        for(kk in 1:k){
          if(cov){
            # part1 <- -log(det(solve(models[[i]][[j]][[kk]])))
            # part2 <- sum(diag(S_test[[kk]] %*% solve(models[[i]][[j]][[kk]])))
            # # bic[i,j] <- dim(obs_test_hh)[1]*(part1+part2)
            # ebic[i,j] <-  extendedBIC(gamma = 0 , omegahat = solve(models[[i]][[j]][[kk]]) , S = cc[[kk]], n = n[[kk]])
            # 
          }else{
            part1 <- -log(det(models[[i]][[j]][[kk]]))
            part2 <- sum(diag(cc[[kk]] %*% models[[i]][[j]][[kk]]))
            ebic[i,j] <-  extendedBIC(gamma = 0 , omegahat = models[[i]][[j]][[kk]] , S = cc[[kk]], n = n[[kk]])
          }
        }
      }
    }
    
  }
  
  rownames(ebic) <- as.character(round(lam1_vec,2))
  colnames(ebic) <- as.character(round(lam2_vec,2))
  
  return(list(ebic = ebic,
              bic = bic,
              lam1_vec = lam1_vec,
              lam2_vec = lam2_vec,
              best_lambdas_ebic = c(lam1_vec[which(ebic == min(ebic),arr.ind = T)[1]],
                                    lam2_vec[which(ebic == min(ebic),arr.ind = T)[2]]),
              best_lambdas_bic = c(lam1_vec[which(bic == min(bic),arr.ind = T)[1]],
                                   lam2_vec[which(bic == min(bic),arr.ind = T)[2]])))  
  
}

## BIC and EBIC for WFLSA
jgl_bic <-  function(data_list , lam1_vec , lam2_vec  ,
                     weight_list = NULL , penalize_diagonal = F, parallel = F){
  
  k <- length(data_list)
  n  <- vector(mode = 'list' , length = k)
  for(kk in 1:k){
    n[[kk]] <- dim(data_list[[kk]])[1]
  }
  cc <-  vector(mode = 'list' , length = k)
  for(kk in 1:k){
    cc[[kk]] <- cov(data_list[[kk]])
  }
  
  bic <- ebic <- matrix(10000,nrow = length(lam1_vec) , ncol = length(lam2_vec))
  
  if(!parallel){
    models <- vector(mode = 'list' , length = length(lam1_vec))
    for(i in 1:length(lam1_vec)){
      models[[i]] <-vector(mode = 'list' , length = length(lam2_vec))
    }
    
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = length(lam1_vec), # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar
    
    for(i in 1:length(lam1_vec)){
      setTxtProgressBar(pb, i)
      
      for(j in 1:length(lam2_vec)){
        # cat("\n","lam1 =",lam1_vec[i] , "\t", "lam2 =",lam2_vec[j])
        lam1 <- lam1_vec[i]
        lam2 <- lam2_vec[j]
        models[[i]][[j]] <- JGL::JGL(Y = data_list, penalty = 'fused', lambda1 = lam1 ,
                                     lambda2 = lam2 , penalize.diagonal = F)$theta
        
        for(kk in 1:k){
          part1 <- -log(det(models[[i]][[j]][[kk]]))
          part2 <- sum(diag(cc[[kk]] %*% models[[i]][[j]][[kk]]))
          ebic[i,j] <-  extendedBIC(gamma = 0 , omegahat = models[[i]][[j]][[kk]] , S = cc[[kk]], n = n[[kk]])
          
        }
      }
    }
    close(pb)
  }  else{
    calculate_models <- function(i){
      models <- vector(mode = 'list', length = length(lam2_vec))
      for(j in 1:length(lam2_vec)){
        # cat("\n","lam1 =",lam1_vec[i] , "\t", "lam2 =",lam2_vec[j])
        lam1 <- lam1_vec[i]
        lam2 <- lam2_vec[j]
        models[[j]] <- JGL::JGL(Y = data_list, penalty = 'fused', lambda1 = lam1 ,
                                lambda2 = lam2 , penalize.diagonal = F)$theta
        
      }
      return(models)
    }
    require(parallel)
    cores <- detectCores() -1
    models <- mclapply(1:length(lam1_vec) , calculate_models)
    for(i in 1:length(lam1_vec)){
      for(j in 1:length(lam2_vec)){
        for(kk in 1:k){
          part1 <- -log(det(models[[i]][[j]][[kk]]))
          part2 <- sum(diag(cc[[kk]] %*% models[[i]][[j]][[kk]]))
          ebic[i,j] <-  extendedBIC(gamma = 0 , omegahat = models[[i]][[j]][[kk]] , S = cc[[kk]], n = n[[kk]])
          
        }
      }
    }
    
  }
  
  rownames(ebic) <- as.character(round(lam1_vec,2))
  colnames(ebic) <- as.character(round(lam2_vec,2))
  
  return(list(ebic = ebic,
              bic = bic,
              lam1_vec = lam1_vec,
              lam2_vec = lam2_vec,
              best_lambdas_ebic = c(lam1_vec[which(ebic == min(ebic),arr.ind = T)[1]],
                                    lam2_vec[which(ebic == min(ebic),arr.ind = T)[2]]),
              best_lambdas_bic = c(lam1_vec[which(bic == min(bic),arr.ind = T)[1]],
                                   lam2_vec[which(bic == min(bic),arr.ind = T)[2]])))  
  
}

