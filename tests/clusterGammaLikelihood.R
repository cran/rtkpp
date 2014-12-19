# TODO: Add comment
#
# Author: iovleff
#-----------------------------------------------------------------------
library(rtkpp)
data(iris)

gamma_model <- clusterGamma( iris[1:4], nbCluster = 2:4, modelNames = clusterGammaNames(shapeBetweenCluster = "all")
                           , strategy = clusterStrategy(nbTry = 1, nbInit = 2))

data<-gamma_model@component@data
nbSample <- nrow(data)
nbVariable <- ncol(data)
nbCluster <- gamma_model@nbCluster
prop <- gamma_model@pk
shape <- gamma_model@component@shape
scale <- gamma_model@component@scale

f <-vector(length=nbSample)
lnComp <- vector(length=nbCluster)

for (i in 1:nbSample)
{
  for (k in 1:nbCluster)
  { lnComp[k] = log(prop[k]) + sum(dgamma(data[i,], shape=shape[k,], scale=scale[k,],log=TRUE)); }
  lmax <- max(lnComp)

  for (k in 1:nbCluster)
  { lnComp[k] =  lnComp[k] -lmax;}

  f[i] = log(sum(exp(lnComp))) + lmax;
}

cat("Computed log-likelihood: ", sum(f), "\n")
cat("Model log-likelihood: ", gamma_model@lnLikelihood, "\n")
