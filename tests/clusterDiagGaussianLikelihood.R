# TODO: Add comment
#
# Author: iovleff
#-----------------------------------------------------------------------
library(rtkpp)
data(iris)

gauss_model <- clusterDiagGaussian(iris[1:4], nbCluster = 2:8, modelNames = diagGaussianNames()
              , strategy = clusterStrategy(nbTry = 3, nbInit = 5), criterion = "ICL")

data<-gauss_model@data
nbSample <- nrow(data)
nbVariable <- ncol(data)
nbCluster <- gauss_model@nbCluster
prop <- gauss_model@pk
mean <- gauss_model@mean
sigma <- gauss_model@sigma

f <-vector(length=nbSample)
lnComp <- vector(length=nbCluster)

for (i in 1:nbSample)
{
  for (k in 1:nbCluster)
  { lnComp[k] = log(prop[k]) + sum(dnorm(data[i,], mean[k,], sigma[k,],log=TRUE)); }
  lmax <- max(lnComp)

  for (k in 1:nbCluster)
  { lnComp[k] =  lnComp[k] -lmax;}

  f[i] = log(sum(exp(lnComp))) + lmax;
}

cat("Computed log-likelihood: ", sum(f), "\n")
cat("Model log-likelihood: ", gauss_model@lnLikelihood, "\n")
