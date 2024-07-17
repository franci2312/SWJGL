pacman::p_load("parallel","iDINGO", "JGL", "ggokabeito", "igraph")

# euclidean distance ------------------------------------------------------

# @param pixel_coords is a p x 2 matrix of spatial coordinates
distance_mat <- function(pixel_coords){ 
  x1 <- pixel_coords[,1]
  x2 <- pixel_coords[,2]
  diff_x1 <- outer(x1, x1, "-")
  diff_x2 <- outer(x2, x2, "-")
  dist_mat <- sqrt(diff_x1^2 + diff_x2^2)
  return(dist_mat)
}


# Hausdorff distance ------------------------------------------------------

# This function computes the Hausdorff distance between two sets 

# @param A coordinates of the first set
# @param B coordinates of the second set
pairwise.hdist <- function(A, B){
  dist.all <- distance_mat(rbind(A, B))
  out.diag <- dist.all[1:nrow(A),(nrow(A)+1):nrow(dist.all)]
  minA <- apply(out.diag, 1, min) 
  minB <- apply(out.diag, 2, min)
  hausd <- max(max(minA), max(minB))
  return(hausd)
}

# e.g. A = grouped.x$Coords[[1]], B = grouped.x$Coords[[2]]


# DW ----------------------------------------------------------------------

# This function computes the Hausdorff distance between a list of coordinates 
# matrices and the corresponding weights

# @param coords is a list of length K (number of clusters) in which each element
# is a matrix of pixel coordinates 

DW <- function(coords){
  D <- W <- matrix(rep(0, length(coords)^2), length(coords))
  diag(D) <- diag(W) <- 0
  require(parallel)
  which.pair <- lower.tri(D)
  for(i in 1:ncol(D)){
    cat("layer", i, "\n")
    tmp <- coords[which.pair[,i]]
    tmp.vec <- numeric(length(which.pair[,i]))
    if(length(tmp) > 0){
      dd <- unlist(mclapply(tmp, function(x) pairwise.hdist(coords[[i]], x)))
      tmp.vec[which.pair[,i]] <- dd
      D[,i] <- tmp.vec
      tmp.vec[which.pair[,i]] <- 1/dd
      W[,i] <- tmp.vec
    }
  }
  return(list(D = D, W = W))
}


# create.weights ----------------------------------------------------------

# @param coords is a list of length K (number of clusters) in which each element
# is a matrix of pixel coordinates 

create.weights <- function(coords){
  W.tmp <- DW(coords)$W
  weights <- W.tmp/max(W.tmp) # kind of normalization
  weights <- weights + t(weights)
  diag(weights) <- rep(1, nrow(weights))
  W <- list()
  K <- length(coords)
  for(k in 1:K){
    W[[k]] <- weights[k,]
  }
  return(W)
}

# in Layers  --------------------------------------------------------------

which.cluster <- function(x = x.rid, cluster = colData(spe)$spatialLIBD,
                          coords = spatialCoords(spe)){
  lev <- levels(cluster)
  Res <- Coords <- list()
  for(i in 1:length(lev)){
    Res[[i]] <- x[,cluster %in% lev[i]]
    Coords[[i]] <- coords[cluster %in% lev[i],]
  }
  return(list(Mat = Res, Coords = Coords))
}


# extract weights and data ------------------------------------------------

# @param obj spatial object
# @param x logcounts matrix
# @param idx.sample idx of the image (we have 12 different images)

extract.wd <- function(obj = spe, x = x.rid, idx.sample){
  lab <- unique(colData(obj)$sample_id)[idx.sample]
  tmp.cols <- (colData(obj)$sample_id == lab)
  x.idx <- x[,tmp.cols]
  x.group <- which.cluster(x.idx, colData(obj)$spatialLIBD[tmp.cols],
                           spatialCoords(obj)[tmp.cols,])
  W <- create.weights(x.group$Coords)
  data <- list()
  for(k in 1:length(W)){
    data[[k]] <- t(as.matrix(x.group$Mat[[k]]))
  }
  return(list(W = W, data = data))
}


# Comparison wflsa wjgl ---------------------------------------------------

metriche <- function(previsti, osservati){
  n <-  table(previsti,osservati)
  precision <- n[2,2]/(n[2,2] + n[2,1])
  recall <- n[2,2]/(n[2,2] + n[1,2])
  F1 <- (2*precision*recall)/(precision + recall)
  FDR <- n[2,1]/(n[2,1] + n[2,2])
  return(list(precision = precision, recall = recall, F1 = F1, FDR = FDR))
}

falsi = function(previsti, osservati){
  n <-  table(previsti,osservati)
  ec <- 1-sum(diag(n)/sum(n))
  fn <- n[1,2]/(n[1,2]+n[2,2])
  fp <- n[2,1]/(n[1,1]+n[2,1])
  return(list(EC = ec, FP = fp, FN = fn))
}





