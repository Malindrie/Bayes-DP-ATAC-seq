library(Seurat)
library(Signac)
set.seed(1234)

load("RObjects/Basal/regions/promoter_peaks/TAC_1_atac_exprs.rds")

# Normalize data
TAC_1_atac_norm <- apply(TAC_1_atac_exprs, 2, function(x)
  log1p((x / sum(x))*10000))

# Scale data
TAC_1_atac_norm <- scale(TAC_1_atac_norm)

gexpr <- TAC_1_atac_norm


#fiting Bayesian non-parametric mixture models
library(dirichletprocess)

Dirichlet_Normal <- function(y){
  
  dp <- dirichletprocess:: DirichletProcessGaussian(y)
  dp <- dirichletprocess::Fit(dp, 1000, progressBar = FALSE)
  burned_dp <- Burn(dp, 100)
  
  return(burned_dp)
}


DP_KNN_10000 <- apply(gexpr, 1, FUN = Dirichlet_Normal)
save(DP_KNN_10000, file="./RObjects/Basal/KNN/DP_KNN_10000.rds")


# # with parallelization
# library(future.apply)
# plan(multisession, workers = 12)
# 
# DP_KNN <- future_apply(gexpr, MARGIN = 1L, FUN = Dirichlet_Normal, future.seed = 0xBEEF)

####################################################
# Running Dirichlet process changing the initial number of clusters
require(dirichletprocess)
require(ggplot2)
require(dplyr)
require(tidyr)
library(Seurat)

set.seed(110010101)

load("RObjects/Basal/KNN/TAC_1_atac_exprs.rds")

TAC_1_atac_norm <- apply(TAC_1_atac_exprs, 2, function(x)
  log1p((x / sum(x))*10000))

TAC_1_atac_norm <- scale(TAC_1_atac_norm)

gexpr <- TAC_1_atac_norm


#fiting Bayesian non-parametric mixture models
library(dirichletprocess)

Dirichlet_Normal <- function(y){
  
  dp <- dirichletprocess::DirichletProcessGaussian(y)
  dp <- dirichletprocess::Initialise(dp, numInitialClusters = length(y))
  dp <- dirichletprocess::Fit(dp, 1000)
  
  return(dp)
}


DP_KNN_all_10000 <- apply(gexpr, 1, FUN = Dirichlet_Normal)
save(DP_KNN_all_10000, file="./RObjects/Basal/KNN/DP_KNN_all_10000.rds")


# conduct Gelman and Rubin diagnostic

require(coda)

Gelman_diag <- NULL

for (i in 1:length(DP_KNN_10000)){
  
  numClusters <- vapply(DP_KNN_10000[[i]]$weightsChain, function(x) length(x), numeric(1))
  numClusters2 <- vapply(DP_KNN_all_10000[[i]]$weightsChain, function(x) length(x), numeric(1))
  
  chains <- mcmc.list(mcmc(cbind(Alpha = DP_KNN_10000[[i]]$alphaChain, 
                                 NumClusters = numClusters, 
                                 Likelihood = DP_KNN_10000[[i]]$likelihoodChain)),
                      mcmc(cbind(Alpha= DP_KNN_all_10000[[i]]$alphaChain, 
                                 NumClusters = numClusters2,
                                 Likelihood = DP_KNN_all_10000[[i]]$likelihoodChain)))
  
  Gelman_diag[i] <- gelman.diag(chains)[[2]]
  
}

Gelman_mpsrf <- data.frame(mpsrf=Gelman_diag, row.names = names(DP_KNN_10000))

save(Gelman_mpsrf, file="./RObjects/aggre/Gelman_mpsrf.rds")


# select regions with Gelman_mpsrf < 1.1
Gelman_mpsrf <- as.matrix(Gelman_mpsrf)
selec_regions <- as.matrix(Gelman_mpsrf[Gelman_mpsrf < 1.1,])
