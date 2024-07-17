wjgl <- function (Y, penalty = "fused", lambda1, lambda2, rho = 1, weights = "equal", 
                  penalize.diagonal = FALSE, maxiter = 500, tol = 1e-05, warm = NULL,
                  screening = "fast", truncate = 1e-05, weight_list = NULL) 
{
  
  return.whole.theta = TRUE
  p = dim(Y[[1]])[2]
  K = length(Y)
  n = rep(0, K)
  for (k in 1:K) {
    n[k] = dim(Y[[k]])[1]
  }
  if (length(dimnames(Y[[1]])[[2]]) == 0) {
    for (k in 1:K) {
      dimnames(Y[[k]])[[2]] = paste("V", 1:p, sep = "")
    }
  }
  for (k in 1:K) {
    for (j in 1:p) {
      Y[[k]][, j] = Y[[k]][, j] - mean(Y[[k]][, j])
    }
  }
  if (length(weights) == 1) {
    if (weights == "equal") {
      weights = rep(1, K)
    }
  }
  if (length(weights) == 1) {
    if (weights == "sample.size") {
      weights = n/sum(n)
    }
  }
  if (screening == "memory.efficient") {
    if (penalty == "fused") {
      connected = JGL:::screen.fgl(Y, lambda1, lambda2, weights)
    }
    if (penalty == "group") {
      connected = JGL:::screen.ggl(Y, lambda1, lambda2, weights)
    }
  }
  if (screening == "fast") {
    connected = rep(TRUE, p)
  }
  if (!((screening == "memory.efficient") | (screening == 
                                             "fast"))) {
    stop("screening must equal \"fast\" or \"memory.efficient\".")
  }
  S = vector("list", length = K)
  for (k in 1:K) {
    ntemp = dim(Y[[k]])[1]
    S[[k]] = cov(Y[[k]][, connected]) * (ntemp - 1)/ntemp
  }
  lam1 <-  lambda1
  lam2 <-  lambda2
  if (length(lam1) > 1) {
    lam1 = lam1[connected, connected]
  }
  if (length(lam2) > 1) {
    lam2 = lam2[connected, connected]
  }
  if (penalty == "fused") {
    if (K == 2) {
      crit1 = list()
      for (k in 1:K) {
        crit1[[k]] = abs(S[[k]]) * weights[k] > lam1 + 
          lam2
      }
      S.sum = matrix(0, sum(connected), sum(connected))
      for (k in 1:K) {
        S.sum = S.sum + weights[k] * S[[k]]
      }
      S.sum = abs(S.sum)
      crit2 = S.sum > 2 * lam1
    }
    if (K > 2) {
      crit1 = list()
      for (k in 1:K) {
        crit1[[k]] = abs(S[[k]]) * weights[k] > lam1
      }
      crit2 = matrix(0, sum(connected), sum(connected))
    }
    critboth = crit2
    for (k in 1:K) {
      critboth = critboth + crit1[[k]]
    }
    critboth = (critboth != 0)
    diag(critboth) = 1
  }
  if (penalty == "group") {
    tempsum = matrix(0, sum(connected), sum(connected))
    for (k in 1:K) {
      tempsum = tempsum + (pmax(weights[k] * abs(S[[k]]) - 
                                  lam1, 0))^2
    }
    critboth = tempsum > lam2^2
    diag(critboth) = 1
  }
  g1 <- igraph:::graph.adjacency(critboth)
  cout = igraph:::clusters(g1)
  blocklist = list()
  unconnected = c()
  if (min(cout$membership) == 0) {
    cout$membership = cout$membership + 1
  }
  for (i in 1:(cout$no)) {
    if (sum(cout$membership == i) == 1) {
      unconnected <- c(unconnected, which(cout$membership == 
                                            i))
    }
    if (sum(cout$membership == i) > 1) {
      blocklist[[length(blocklist) + 1]] <- which(cout$membership == 
                                                    i)
    }
  }
  connected[unconnected] = FALSE
  connected.index = rep(0, p)
  connected.index[connected] = 1:sum(connected)
  unconnected = !connected
  theta = list()
  for (k in 1:K) {
    theta[[k]] = matrix(0, sum(connected), sum(connected))
    if (sum(connected) > 0) {
      dimnames(theta[[k]])[[1]] = dimnames(theta[[k]])[[2]] = dimnames(Y[[k]])[[2]][connected]
    }
  }
  Yu = list()
  for (k in 1:K) {
    Yu[[k]] = Y[[k]][, unconnected]
  }
  if (length(lambda1) == 1) {
    lam1.unconnected = lambda1
  }
  if (length(lambda1) > 1) {
    lam1.unconnected = diag(lambda1)[unconnected]
  }
  if (length(lambda2) == 1) {
    lam2.unconnected = lambda2
  }
  if (length(lambda2) > 1) {
    lam2.unconnected = diag(lambda2)[unconnected]
  }
  if (!penalize.diagonal) {
    lam1.unconnected = lam1.unconnected * 0
    if (penalty == "group") {
      lam2.unconnected = lam2.unconnected * 0
    }
  }
  if (sum(unconnected) > 0) {
    theta.unconnected = admm.iters.unconnected.weighted(Yu, lambda1 = lam1.unconnected, 
                                               lambda2 = lam2.unconnected, penalty = penalty, rho = rho, 
                                               weights = weights, maxiter = maxiter, tol = tol,weighted_list = weight_list)$Z
    for (k in 1:K) {
      names(theta.unconnected[[k]]) = dimnames(Y[[k]])[[2]][!connected]
    }
  }
  if (sum(unconnected) == 0) {
    theta.unconnected = NULL
  }
  if (length(blocklist) > 0) {
    for (i in 1:length(blocklist)) {
      bl <- blocklist[[i]]
      Ybl = list()
      for (k in 1:K) {
        Ybl[[k]] = Y[[k]][, bl]
      }
      if (length(lambda1) == 1) {
        lam1.bl = lambda1
      }
      if (length(lambda1) > 1) {
        lam1.bl = lambda1[bl, bl]
      }
      if (length(lambda2) == 1) {
        lam2.bl = lambda2
      }
      if (length(lambda2) > 1) {
        lam2.bl = lambda2[bl, bl]
      }
      lam1.bl = JGL:::penalty.as.matrix(lam1.bl, dim(Ybl[[1]])[2], 
                                  penalize.diagonal = penalize.diagonal)
      if (penalty == "fused") {
        lam2.bl = JGL:::penalty.as.matrix(lam2.bl, dim(Ybl[[1]])[2], 
                                    penalize.diagonal = TRUE)
      }
      if (penalty == "group") {
        lam2.bl = JGL:::penalty.as.matrix(lam2.bl, dim(Ybl[[1]])[2], 
                                    penalize.diagonal = penalize.diagonal)
      }
      if (length(warm) == 0) {
        warm.bl = NULL
      }
      if (length(warm) > 0) {
        warm.bl = list()
        for (k in 1:K) {
          warm.bl[[k]] = warm[[k]][bl, bl]
        }
      }
      Thetabl = admm.iters.weighted(Ybl, lam1.bl, lam2.bl, penalty = penalty, 
                           rho = rho, weights = weights, penalize.diagonal = TRUE, 
                           maxiter = maxiter, tol = tol, warm = warm.bl,weighted_list = weight_list)
      for (k in 1:K) {
        theta[[k]][connected.index[bl], connected.index[bl]] = Thetabl$Z[[k]]
      }
    }
  }
  if (dim(theta[[1]])[1] > 0) {
    for (k in 1:K) {
      rounddown = abs(theta[[k]]) < truncate
      diag(rounddown) = FALSE
      theta[[k]] = theta[[k]] * (1 - rounddown)
    }
  }
  if (!return.whole.theta) {
    out = list(theta = theta, theta.unconnected = theta.unconnected, 
               connected = connected)
  }
  if (return.whole.theta) {
    whole.theta = list()
    for (k in 1:K) {
      whole.theta[[k]] = matrix(0, p, p)
      diag(whole.theta[[k]])[unconnected] = theta.unconnected[[k]]
      whole.theta[[k]][connected, connected] = theta[[k]]
      dimnames(whole.theta[[k]])[[1]] = dimnames(whole.theta[[k]])[[2]] = dimnames(Y[[k]])[[2]]
    }
    out = list(theta = whole.theta, connected = connected)
  }
  class(out) = "jgl"
  return(out)
}

admm.iters.unconnected.weighted <- function (Y, lambda1, lambda2, penalty = "fused", rho = 1, rho.increment = 1, 
                                             weights, maxiter = 1000, tol = 1e-05, warm = NULL,weighted_list=NULL) 
{
  K = length(Y)
  for (k in 1:K) {
    Y[[k]] = as.matrix(Y[[k]])
  }
  p = dim(Y[[1]])[2]
  n = weights
  S = list()
  for (k in 1:K) {
    ntemp = dim(Y[[k]])[1]
    S[[k]] = rep(0, p)
    for (j in 1:p) {
      S[[k]][j] = var(Y[[k]][, j]) * (ntemp - 1)/ntemp
    }
  }
  theta = list()
  for (k in 1:K) {
    theta[[k]] = 1/S[[k]]
  }
  Z = list()
  for (k in 1:K) {
    Z[[k]] = rep(0, p)
  }
  W = list()
  for (k in 1:K) {
    W[[k]] = rep(0, p)
  }
  lam1 = lambda1
  lam2 = lambda2
  iter = 0
  diff_value = 10
  while ((iter == 0) || (iter < maxiter && diff_value > tol)) {
    theta.prev = theta
    for (k in 1:K) {
      B = n[k] * S[[k]] - rho * (Z[[k]] - W[[k]])
      theta[[k]] = 1/(2 * rho) * (-B + sqrt(B^2 + 4 * 
                                              rho * n[k]))
    }
    A = list()
    for (k in 1:K) {
      A[[k]] = theta[[k]] + W[[k]]
    }
    if (penalty == "fused") {
      if (K == 2) {
        Z = JGL:::flsa2(A, rho, lam1, lam2, penalize.diagonal = TRUE)
      }
      if (K > 2) {
        Z = wflsa_A(A,L =  rho,lam1 =  lam1,lam2 =  lam2,penalize_diagonal =  TRUE,weight_list = weighted_list)
      }
    }
    if (penalty == "group") {
      Z = JGL:::dsgl(A, rho, lam1, lam2, penalize.diagonal = TRUE)
    }
    for (k in 1:K) {
      W[[k]] = W[[k]] + (theta[[k]] - Z[[k]])
    }
    iter = iter + 1
    diff_value = 0
    for (k in 1:K) {
      diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]]))/sum(abs(theta.prev[[k]]))
    }
    rho = rho * rho.increment
  }
  diff = 0
  for (k in 1:K) {
    diff = diff + sum(abs(theta[[k]] - Z[[k]]))
  }
  out = list(theta = theta, Z = Z, diff = diff, iters = iter)
  return(out)
}

admm.iters.weighted <- function (Y, lambda1, lambda2, penalty = "fused", rho = 1, rho.increment = 1, 
                                 weights, penalize.diagonal, maxiter = 1000, tol = 1e-05, 
                                 warm = NULL, weighted_list=NULL) 
{
  K = length(Y)
  p = dim(Y[[1]])[2]
  n = weights
  ns = c()
  for (k in 1:K) {
    ns[k] = dim(Y[[k]])[1]
  }
  S = list()
  for (k in 1:K) {
    S[[k]] = cov(Y[[k]]) * (ns[k] - 1)/ns[k]
  }
  theta = list()
  for (k in 1:K) {
    theta[[k]] = diag(1/diag(S[[k]]))
  }
  Z = list()
  for (k in 1:K) {
    Z[[k]] = matrix(0, p, p)
  }
  W = list()
  for (k in 1:K) {
    W[[k]] = matrix(0, p, p)
  }
  lam1 = JGL:::penalty.as.matrix(lambda1, p, penalize.diagonal = penalize.diagonal)
  if (penalty == "fused") {
    lam2 = JGL:::penalty.as.matrix(lambda2, p, penalize.diagonal = TRUE)
  }
  if (penalty == "group") {
    lam2 = JGL:::penalty.as.matrix(lambda2, p, penalize.diagonal = penalize.diagonal)
  }
  iter = 0
  diff_value = 10
  while ((iter == 0) || (iter < maxiter && diff_value > tol)) {
    if (FALSE) {
      print(paste("iter=", iter))
      if (penalty == "fused") {
        # print(paste("crit=", crit(theta, S, n = rep(1, 
        #                                             K), lam1, lam2, penalize.diagonal = penalize.diagonal)))
        # print(paste("crit=", crit(Z, S, n = rep(1, K), 
        #                           lam1, lam2, penalize.diagonal = penalize.diagonal)))
      }
      if (penalty == "group") {
        # print(paste("crit=", gcrit(theta, S, n = rep(1, 
        #                                              K), lam1, lam2, penalize.diagonal = penalize.diagonal)))
      }
    }
    theta.prev = theta
    for (k in 1:K) {
      edecomp = eigen(S[[k]] - rho * Z[[k]]/n[k] + rho * 
                        W[[k]]/n[k])
      D = edecomp$values
      V = edecomp$vectors
      D2 = n[k]/(2 * rho) * (-D + sqrt(D^2 + 4 * rho/n[k]))
      theta[[k]] = V %*% diag(D2) %*% t(V)
    }
    A = list()
    for (k in 1:K) {
      A[[k]] = theta[[k]] + W[[k]]
    }
    if (penalty == "fused") {
      if (K == 2) {
        Z = JGL:::flsa2(A, rho, lam1, lam2, penalize.diagonal = TRUE)
      }
      if (K > 2) {
        Z = wflsa_A(A = A, L = rho,lam1 =  lam1,lam2 =  lam2,penalize_diagonal = TRUE, weight_list = weighted_list)
      }
    }
    if (penalty == "group") {
      Z = JGL:::dsgl(A, rho, lam1, lam2, penalize.diagonal = TRUE)
    }
    for (k in 1:K) {
      W[[k]] = W[[k]] + (theta[[k]] - Z[[k]])
    }
    iter = iter + 1
    diff_value = 0
    for (k in 1:K) {
      diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]]))/sum(abs(theta.prev[[k]]))
    }
    rho = rho * rho.increment
  }
  diff = 0
  for (k in 1:K) {
    diff = diff + sum(abs(theta[[k]] - Z[[k]]))
  }
  out = list(theta = theta, Z = Z, diff = diff, iters = iter)
  return(out)
}

## BIC and EBIC for WFLSA
wjgl_bic <-  function(data_list , lam1_vec , lam2_vec  ,
                     weight_list = NULL , penalize_diagonal = F, parallel = F, tol = 1e-5){
  
  lam1_vec <- sort(lam1_vec,decreasing = T)
  lam2_vec <- sort(lam2_vec,decreasing = T)

  k <- length(data_list)
  n  <- vector(mode = 'list' , length = k)
  for(kk in 1:k){
    n[[kk]] <- dim(data_list[[kk]])[1]
  }
  cc <-  vector(mode = 'list' , length = k)
  for(kk in 1:k){
    cc[[kk]] <- cov(data_list[[kk]])
  }
  
  aic <- ebic <- matrix(10000,nrow = length(lam1_vec) , ncol = length(lam2_vec))
  
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
      lam1 <- lam1_vec[i]
      
      models[[i]][[1]] <- wjgl(Y = data_list, weight_list = weight_list, penalty = 'fused', lambda1 = lam1 ,
                               lambda2 =  lam2_vec[1] , penalize.diagonal = F,tol = tol)$theta
      
      for(j in 2:length(lam2_vec)){
        # cat("\n","lam1 =",lam1_vec[i] , "\t", "lam2 =",lam2_vec[j])
        lam2 <- lam2_vec[j]
        # models[[i]][[j]] <- JGL::JGL(Y = data_list, penalty = 'fused', lambda1 = lam1 ,
        #                              lambda2 = lam2 , penalize.diagonal = F)$theta
        models[[i]][[j]] <- wjgl(Y = data_list, weight_list = weight_list, penalty = 'fused', lambda1 = lam1 ,
                                     lambda2 = lam2 , penalize.diagonal = F,tol = tol, warm = models[[i]][[j-1]])$theta
        # models[[i]][[j]] <- wjgl(Y = data_list, weight_list = weight_list, penalty = 'fused', lambda1 = lam1 ,
        #                          lambda2 = lam2 , penalize.diagonal = F,tol = tol)$theta
        
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
      lam1 <- lam1_vec[i]
      
      models[[1]] <- wjgl(Y = data_list, weight_list = weight_list, penalty = 'fused', lambda1 = lam1 ,
                          lambda2 = lam2_vec[1] , penalize.diagonal = F,tol = tol)$theta
      for(j in 2:length(lam2_vec)){
        # cat("\n","lam1 =",lam1_vec[i] , "\t", "lam2 =",lam2_vec[j])
        lam2 <- lam2_vec[j]
        # models[[j]] <- JGL::JGL(Y = data_list, penalty = 'fused', lambda1 = lam1 ,
        #                         lambda2 = lam2 , penalize.diagonal = F)$theta
        models[[j]] <- wjgl(Y = data_list, weight_list = weight_list, penalty = 'fused', lambda1 = lam1 ,
                                 lambda2 = lam2 , penalize.diagonal = F,tol = tol,warm = models[[j-1]])$theta
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
          aic[i,j] <- aic(omegahat = models[[i]][[j]][[kk]],S = cc[[kk]],n = n[[kk]])
        }
      }
    }
    
  }
  
  rownames(ebic) <- as.character(round(lam1_vec,2))
  colnames(ebic) <- as.character(round(lam2_vec,2))
  
  rownames(aic) <- as.character(round(lam1_vec,2))
  colnames(aic) <- as.character(round(lam2_vec,2))
  
  ebic <- ebic[c(length(lam1_vec):1),c(length(lam2_vec):1)]
  aic <- aic[c(length(lam1_vec):1),c(length(lam2_vec):1)]
  
  lam1_vec <- sort(lam1_vec,decreasing = F)
  lam2_vec <- sort(lam2_vec,decreasing = F)
  
  return(list(ebic = ebic,
              aic = aic,
              lam1_vec = lam1_vec,
              lam2_vec = lam2_vec,
              best_lambdas_ebic = c(lam1_vec[which(ebic == min(ebic),arr.ind = T)[1]],
                                    lam2_vec[which(ebic == min(ebic),arr.ind = T)[2]]),
              best_lambdas_aic = c(lam1_vec[which(aic == min(aic),arr.ind = T)[1]],
                                   lam2_vec[which(aic == min(aic),arr.ind = T)[2]])))  
  
}


wflsa_A <- function (A,L, lam1, lam2, penalize_diagonal  = F, weight_list = NULL) 
{
  k <- length(A)
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

aic <- function (omegahat, S, n) 
{
  p = nrow(omegahat)
  es = sum(omegahat[upper.tri(omegahat)] != 0)
  return(-log(det(omegahat)) + sum(diag(omegahat %*% S)) + 
           es * (2/n))
}
