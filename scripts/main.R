pacman::p_load("parallel","iDINGO", "JGL", "ggokabeito", "igraph")

# load data ---------------------------------------------------------------

#BiocManager::install("spatialLIBD")
library("spatialLIBD")
spe <- fetch_data(type = "spe")

x <- logcounts(spe)

# selected genes from Maynard (Suppl. Table S5)
imp.genes <- c("CTGF", "LMO4", "ETV1", "PCP4", "CXCL14", "NTNG2", "NR4A2", "CUX2", "RELN", "CALB1", "B3GALT2", 
               "CBLN2", "RXFP1", "CDH24", "GABRA5", "CARTPT", "OPRK1", "RPRM", "CACNA1E", "FEZF2", "HTR2C", 
               "PDE1A", "VAT1L", "GSG1L", "C1QL2", "MEF2C", "ETV1", "BHLHE22", "SEMA3E", "TRIB2", "KCNIP2", 
               "KITLG", "SLITRK1", "SCN3B", "KCNK2", "SYNPR", "CCK", "FAM3C", "ANXA1", "FXYD6", "RASGRF2", 
               "DKK3", "RORB", "TLE1", "GPR6", "TOX", "PCDH17", "POU3F1", "BEND5", "UNC5D", "ATP2B4", "CNR1", 
               "PCP4", "TLE4", "CRIM1", "CPNE7", "TLE4", "PENK", "SYT6", "TH", "BCL11B", "SATB2", "FOXP2", 
               "TBR1", "CRYM", "RAC3", "PDYN", "KCNA1", "SNCG", "SYT17", "CARTPT", "KCNH4", "CHRNA7", "MFGE8", 
               "LXN", "ADRA2A", "AKR1C2", "CHRNA3", "MARCKSL1", "CYR61", "AKR1C3", "NEFH", "LMO3", "TMEM163", 
               "LGALS1", "PRSS12", "SCN4B", "PLXND1", "LHX2", "COL6A1", "NR4A3", "GRIK4", "SYT2", "WFS1", 
               "INPP4B", "OPN3", "ID2", "OMA1", "IGFBP4", "COL24A1", "CUX1", "FOXO1", "S100A10", "KCNIP2", 
               "KCNIP1", "TLE3", "SEMA3C", "DTX4", "MDGA1", "SOX5", "OTX1", "SV2C", "LIX1", "PCDH20", "TPBG", 
               "RGS8", "IGSF11", "LDB2", "SYT10", "NPY2R", "CYP39A1", "CNTN6", "DIAPH3", "PPP1R1B", "CRYM", "SYT9")

which.genes <- rowData(spe)@rownames[rowData(spe)$gene_name %in% imp.genes]
x.rid <- x[rownames(x) %in% which.genes,]
# x.rid contains on the rows the important genes and on the cols all the spots

# retained genes
genes <- rowData(spe)$gene_name[rowData(spe)$gene_name %in% imp.genes]

source("weights.R")

# extract weights and data ------------------------------------------------

x.1 <- extract.wd(idx.sample = 9) # 9th frame 

# Precision matrices estimation ------------------------------------------

K <- length(x.1$data)
cc <- c()
for(kk in 1:K){
  cc <- c(cc,which(diag(cov(x.1$data[[kk]]))==0) )
}

for(kk in 1:K){
  x.1$data[[kk]] <- x.1$data[[kk]][,-unique(cc)]
}

source("functions_wflsa.R")
source("functions_wjgl.R")
source("graph.R")


# JGL ---------------------------------------------------------------------

# scale data
x <- list()
for(k in 1:K) {
  x[[k]]<-scale(x.1$data[[k]],center = T,scale = T)
}

gene_names <- genes[-c(unique(cc))]

# mod_grid <- wjgl_bic(data_list = x, lam1_vec = seq(1e-3,1, length.out = 10),
#                     lam2_vec = seq(1e-3,1, length.out = 10), weight_list = x.1$W,
#                     penalize_diagonal = F,parallel = T, tol = 1e-3)
# mod_grid$best_lambdas_ebic[1] # 0.2563333
# mod_grid$best_lambdas_ebic[2] # 0.012

mod_jgl <- wjgl(Y = x, penalty = "fused", lambda1 = 0.2563333, lambda2 = 0.012,
                penalize.diagonal = F, weight_list = x.1$W, tol = 1e-4)



par(mfrow = c(3,3))
for(kk in 1:K){
  theta0 <- mod_jgl$theta[[kk]]
  dimnames(theta0) <- NULL
  adj <- (theta0!=0) + 0
  graph.hat <- igraph::graph_from_adjacency_matrix(adjmatrix = adj, 
                               mode = 'undirected',diag = F,weighted = F)
  plot.igraph(x = graph.hat, layout = layout_in_circle(graph = graph.hat),
              vertex.size = 10)
}
par(mfrow = c(1,1))

# save(mod_jgl, file = "mod_jgl.RData")

# save graphs 
pdf(file = "layer4.pdf")
graph4 <- circ.graph(mod = mod_jgl, idx = 4, weighted.edges = F)
dev.off()

# flsa --------------------------------------------------------------------

# mod_grid <- wflsa_bic(x, lam1_vec = seq(1e-3,1, length.out = 10), 
#                       lam2_vec = seq(1e-3,1, length.out = 10), 
#                       lambda_soft = 0.0001,cov = F,weight_list = x.1$W,
#                       penalize_diagonal = F, parallel = T)
# par(mfrow = c(1,1))
# image(mod_grid$lam1_vec , mod_grid$lam2_vec,mod_grid$ebic)
# points(mod_grid$best_lambdas_ebic[1],mod_grid$best_lambdas_ebic[2],pch = 16)

# mod_grid <- wflsa_bic(x, lam1_vec = seq(0.2,0.4, length.out = 10), 
#                      lam2_vec = seq(0.1,0.2, length.out = 10), 
#                       lambda_soft = 0.0001, cov = F, weight_list = x.1$W,
#                      penalize_diagonal = F, parallel = T)

# mod_grid$best_lambdas_ebic[1] # 0.2666667
# mod_grid$best_lambdas_ebic[2] # 0.1222222


mod_flsa <- wflsa(x, lam1 = 0.2666667, lam2 = 0.1222222, lambda_soft = 0.0001,
                  penalize_diagonal = F, weight_list = x.1$W, cov = F)


par(mfrow = c(3,3))
for(kk in 1:K){
  theta0 <- mod_flsa[[kk]]
  dimnames(theta0) <- NULL
  adj <- (theta0!=0) + 0
  graph.hat <- igraph::graph_from_adjacency_matrix(adjmatrix = adj, 
                              mode = 'undirected',diag = F,weighted = F)
  plot.igraph(x = graph.hat, layout = layout_in_circle(graph = graph.hat),
              vertex.size = 10)
}
par(mfrow = c(1,1))

# save(mod_flsa, file = "mod_flsa.RData")

# Comparison wjgl wflsa ----------------------------------------------------

adj <- lapply(mod_flsa, function(k) (k != 0) + 0)
prev <- unlist(lapply(adj, as.vector))

adj_jgl <- lapply(mod_jgl$theta, function(k) (k != 0) + 0)
obs <- unlist(lapply(adj_jgl, as.vector))

M <- metriche(prev, obs)
FR <- falsi(prev, obs)

cbind(knitr::kable(c(unlist(FR),0)), knitr::kable(unlist(M)))

# Edge density ------------------------------------------------------------

mod <- mod_jgl

for(k in 1:length(mod$theta)){
  adj <- (mod$theta[[k]] != 0) + 0
  graph.hat <- igraph::graph_from_adjacency_matrix(adjmatrix = adj, 
                           mode = 'undirected', diag = F, weighted = F)  
  print(edge_density(graph.hat))
}





