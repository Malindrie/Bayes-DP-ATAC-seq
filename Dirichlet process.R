#fiting Bayesian non-parametric mixture models
library(dirichletprocess)

Likelihood.poisson <- function(mdobj, x, theta){
  return(as.numeric(dpois(x, theta[[1]])))
}


PriorDraw.poisson <- function(mdobj, n){
  draws <- rgamma(n, mdobj$priorParameters[1], mdobj$priorParameters[2])
  theta <- list(array(draws, dim=c(1,1,n)))
  return(theta)
}


PosteriorDraw.poisson <- function(mdobj, x, n=1){
  priorParameters <- mdobj$priorParameters
  lambda <- rgamma(n, priorParameters[1] + sum(x), priorParameters[2] + nrow(x))
  return(list(array(lambda, dim=c(1,1,n))))
}


Predictive.poisson <- function(mdobj, x){
  priorParameters <- mdobj$priorParameters
  pred <- numeric(length(x))
  for(i in seq_along(x)){
    alphaPost <- priorParameters[1] + x[i]
    betaPost <- priorParameters[2] + 1
    pred[i] <- (priorParameters[2] ^ priorParameters[1]) / gamma(priorParameters[1])
    pred[i] <- pred[i] * gamma(alphaPost) / (betaPost^alphaPost)
    pred[i] <- pred[i] * (1 / prod(factorial(x[i])))
  }
  return(pred)
}


poisMd <- MixingDistribution(distribution="poisson",
                             priorParameters = c(1, 1),
                             conjugate="conjugate")


Dirichlet_Poisson <- function(y, poisMd){
  
  dp <- dirichletprocess::DirichletProcessCreate(y, poisMd)
  dp <- dirichletprocess::Initialise(dp)
  dp <- dirichletprocess::Fit(dp, 1000)
  
  return(dp)
}


# running dirichletprocess changing intnumclusters
set.seed(110010101)

Dirichlet_Poisson_2 <- function(y, poisMd){
  
  dp <- dirichletprocess::DirichletProcessCreate(y, poisMd)
  dp <- dirichletprocess::Initialise(dp, numInitialClusters = length(y))
  dp <- dirichletprocess::Fit(dp, 1000)
  
  return(dp)
}


# dirichletprocess diagnostics
require(coda)

Gelman_diag <- NULL

for (i in 1:length(DP_var_5000)){
  
  numClusters <- vapply(DP_var_5000[[i]]$weightsChain, function(x) length(x), numeric(1))
  numClusters2 <- vapply(DP_var_all5000_2[[i]]$weightsChain, function(x) length(x), numeric(1))
  
  chains <- mcmc.list(mcmc(cbind(Alpha = DP_var_5000[[i]]$alphaChain, 
                                 NumClusters = numClusters, 
                                 Likelihood = DP_var_5000[[i]]$likelihoodChain)),
                      mcmc(cbind(Alpha= DP_var_all5000_2[[i]]$alphaChain, 
                                 NumClusters = numClusters2,
                                 Likelihood = DP_var_all5000_2[[i]]$likelihoodChain)))
  
  Gelman_diag[i] <- gelman.diag(chains)[[2]]
  
}

Gelman_mpsrf <- data.frame(mpsrf=Gelman_diag, row.names = names(DP_var_5000))


