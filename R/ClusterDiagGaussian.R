#-----------------------------------------------------------------------
#     Copyright (C) 2012-2014  Serge Iovleff, University Lille 1, Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
#' @include IClusterModel.R
NULL

#-----------------------------------------------------------------------
#' Create an instance of the [\code{\linkS4class{ClusterDiagGaussian}}] class
#'
#' This function computes the optimal diagonal Gaussian mixture model according
#' to the [\code{criterion}] among the list of model given in [\code{modelNames}]
#' and the number of clusters given in [\code{nbCluster}], using the strategy specified in [\code{strategy}].
#'
#' @param data frame or matrix containing the data. Rows correspond to observations
#' and columns correspond to variables. If the data set contains NA values, they
#' will be estimated during the estimation process.
#' @param nbCluster  [\code{\link{vector}}] listing the number of clusters to test.
#' @param modelNames [\code{\link{vector}}] of model names to run. By default all diagonal
#' Gaussian models are estimated. All the model names are given by the method
#' [\code{\link{diagGaussianNames}}].
#' @param strategy a [\code{\linkS4class{ClusterStrategy}}] object containing
#' the strategy to run. clusterStrategy() method by default.
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL". Default is "ICL".
#'
#' @examples
#' ## A quantitative example with the famous geyser data set
#' data(geyser)
#' ## add 10 missing values
#' x = geyser;
#' x[round(runif(5,1,nrow(geyser))), 1] <- NA
#' x[round(runif(5,1,nrow(geyser))), 2] <- NA
#' ## with default values
#' model <- clusterDiagGaussian(data=x, nbCluster=2:3, strategy = clusterFastStrategy())
#'
#' ## use graphics functions
#' \dontrun{
#' plot(model)
#' }
#'
#' ## get summary
#' summary(model)
#' ## print model
#' print(model)
#' ## get estimated missing values
#' model["data"][model["missings"]]
#'
#' @return An instance of the [\code{\linkS4class{ClusterDiagGaussian}}] class.
#' @author Serge Iovleff
#' @export
#'
clusterDiagGaussian <- function(data, nbCluster=2, modelNames=diagGaussianNames(), strategy=clusterStrategy(), criterion="ICL")
{
  # check nbCluster
  nbClusterModel = length(nbCluster);
  nbClusterMin = min(nbCluster);
  nbClusterMax = max(nbCluster);
  if (nbClusterMin < 1) { stop("The number of clusters must be greater or equal to 1")}

  # check criterion
  if(sum(criterion %in% c("BIC","AIC", "ICL")) != 1)
  { stop("criterion is not valid. See ?clusterDiagGaussian for the list of valid criterion")}

  # check data
  data = as.matrix(data)
  if (nrow(data) <= 3*nbClusterMax) {stop("There is not enough individuals (rows) in the data set")}
  if (ncol(data) < 1) {stop("Error: empty data set")}

  # check modelNames
  if (!validDiagGaussianNames(modelNames))
  { stop("modelNames is not valid. See ?diagGaussianNames for the list of valid model names")}

  # check strategy
  if(class(strategy)[1] != "ClusterStrategy")
  {stop("strategy is not a Cluster Stategy class (must be an instance of the class ClusterStrategy).")}
  validObject(strategy);

  # Create model
  model = new("ClusterDiagGaussian", data)
  model@missings = which(is.na(data), arr.ind=TRUE);
  model@strategy = strategy;
  # start estimation of the models
  ResFlag <- FALSE;
  res1Flag <- FALSE;
  if (nbClusterMin == 1)
  {
    model = .diagGaussianNoCluster(model);
    res1Flag <- TRUE;
    ind <- which(nbCluster == 1, arr.ind = TRUE);
    nbCluster <- nbCluster[-ind];
  }
  if (length(nbCluster) >0)
  {
    resFlag = .Call("clusterMixture", model, nbCluster, modelNames, strategy, criterion, PACKAGE="rtkpp");
  }
  # set names
  if (resFlag != TRUE && res1Flag != TRUE) {cat("WARNING: An error occur during the clustering process");}
  else
  {
    colnames(model@mean)  <- colnames(model@data);
    colnames(model@sigma) <- colnames(model@data);
  }
  model
}

#' Definition of the [\code{\linkS4class{ClusterDiagGaussian}}] class
#'
#' This class defines a diagonal Gaussian mixture Model.
#'
#' This class inherits from the [\code{\linkS4class{IClusterModel}}] class.
#' A diagonal gaussian model is a mixture model of the form:
#' \deqn{
#'   f({x}|\boldsymbol{\theta})
#'   =\sum_{k=1}^K p_k \prod_{j=1}^d \phi(x_j;\mu_{jk},\sigma^2_{jk})
#'    \quad x \in {R}^d.
#' }
#'
#' @slot mean  Matrix with the mean of the jth variable in the kth cluster.
#' @slot sigma  Matrix with the standard deviation of the jth variable in the kth cluster.
#' @seealso [\code{\linkS4class{IClusterModel}}] class
#'
#' @examples
#' getSlots("ClusterDiagGaussian")
#' data(geyser)
#' new("ClusterDiagGaussian", data=geyser)
#'
#' @author Serge Iovleff
#'
#' @name ClusterDiagGaussian
#' @rdname ClusterDiagGaussian-class
#' @aliases ClusterDiagGaussian-class
#' @exportClass ClusterDiagGaussian
#'
setClass(
    Class="ClusterDiagGaussian",
    representation( mean = "matrix", sigma = "matrix"),
    contains=c("IClusterModel"),
    prototype=list( mean   = matrix(nrow=0,ncol=0), sigma = matrix(nrow=0,ncol=0) ),
    validity=function(object)
    {
      if (nrow(object@mean)!=object@nbCluster)
      {stop("mean must have nbCluster rows.")}
      if (ncol(object@mean)!=ncol(object@data))
      {stop("mean must have nbVariable columns.")}
      if (nrow(object@sigma)!=object@nbCluster)
      {stop("sigma must have nbCluster rows.")}
      if (ncol(object@sigma)!=ncol(object@data))
      {stop("sigma must have nbVariable columns.")}
      if (!validDiagGaussianNames(object@modelName))
      {stop("Invalid Gaussian model name.")}
      return(TRUE)
    }
)

#-----------------------------------------------------------------------
#' Initialize an instance of a rtkpp class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterDiagGaussian}}] class.
#' Used internally in the `rtkpp' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterDiagGaussian"),
    definition=function(.Object, data, nbCluster=2, modelName="gaussian_pk_sjk")
    {
      .Object <- callNextMethod(.Object, data, nbCluster, modelName)
      # resize
      nbVariable <- ncol(.Object@data)
      # for modelName
      if(missing(modelName)) {.Object@modelName<-"gaussian_pk_sjk"}
      else  {.Object@modelName<-modelName}
      .Object@mean <- matrix(0, nbCluster, nbVariable)
      .Object@sigma <- matrix(1, nbCluster, nbVariable)
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname print-methods
#' @aliases print print,ClusterDiagGaussian-method
#'
setMethod(
  f="print",
  signature=c("ClusterDiagGaussian"),
  function(x,...){
    cat("****************************************\n")
    callNextMethod()
    cat("****************************************\n")
    for(k in 1:length(x@pk))
    {
      cat("*** Cluster: ",k,"\n")
      cat("* Proportion = ", format(x@pk[k]), "\n")
      cat("* Mean(s)    = ", format(x@mean[k,]), "\n")
      cat("* Sd(s)      = ", format(x@sigma[k,]), "\n")
      cat("****************************************\n")
    }
  }
)

#' @rdname show-methods
#' @aliases show-ClusterDiagGaussian,ClusterDiagGaussian,ClusterDiagGaussian-method
setMethod(
    f="show",
    signature=c("ClusterDiagGaussian"),
    function(object)
    {
      cat("****************************************\n")
      callNextMethod()
      cat("****************************************\n")
      for(k in 1:length(object@pk))
      {
        cat("*** Cluster: ",k,"\n")
        cat("* Proportion = ", format(object@pk[k]), "\n")
        cat("* Means      = ", format(object@mean[k,]), "\n")
        cat("* Variances  = ", format(object@sigma[k,]), "\n")
        cat("****************************************\n")
      }
    }
)

#' @rdname summary-methods
#' @aliases summary summary,ClusterDiagGaussian-method
#'
setMethod(
    f="summary",
    signature=c("ClusterDiagGaussian"),
    function(object, ...)
    {
      cat("**************************************************************\n")
      callNextMethod()
      cat("**************************************************************\n")
    }
)

#' Plotting of a class [\code{\linkS4class{ClusterDiagGaussian}}]
#'
#' Plotting data from a [\code{\linkS4class{ClusterDiagGaussian}}] object
#' using the estimated parameters and partition.
#'
#' @param x an object of class [\code{\linkS4class{ClusterDiagGaussian}}]
#' @param y a list of variables to plot (subset). Variables names or indices.
#' If missing all the variables are represented.
#' @param ... further arguments passed to or from other methods
#'
#' @importFrom graphics plot
#' @aliases plot-ClusterDiagGaussian
#' @docType methods
#' @rdname plot-ClusterDiagGaussian-method
#' @export
#'
#' @seealso \code{\link{plot}}
#' @examples
#'   ## the famous iris data set
#'   data(iris)
#'   model <- clusterDiagGaussian(iris[1:4], 3, strategy = clusterFastStrategy())
#' \dontrun{
#'   plot(model)
#'   plot(model, c(1,3))
#'   plot(model, c("Sepal.Length","Sepal.Width"))
#'   }
#'
setMethod(
    f="plot",
    signature=c("ClusterDiagGaussian"),
    function(x, y, ...)
    { # use generic plot
      .clusterPlot(x, y, .dGauss,...);
    }
)

# wrapper of dnorm
# x a vector with the point
.dGauss <- function(x, j, k, model)
{ dnorm(x, (model@mean)[k, j] , (model@sigma)[k, j])}

# wrapper of dnorm
# x a vector with the point
.diagGaussianNoCluster <- function(model)
{
  nbSample   <- nrow(model@data);
  nbVariable <- ncol(model@data);
  model@mean  <- matrix(colMeans(model@data, na.rm = TRUE), nrow=1, ncol = nbVariable);
  model@sigma <- matrix(apply(model@data, 2, na.rm = TRUE, sd), nrow=1, ncol = nbVariable);
  model@nbCluster <- 1;
  model@pk <- c(1);
  model@tik <- matrix(1, nrow= nrow(model@data), ncol =1);
  model@lnFi <- rowSums( t(apply(model@data, 1, mean=as.vector(model@mean), sd = as.vector(model@sigma), log = TRUE, dnorm)) );
  model@zi <- as.integer(rep(1, nbVariable));
  model@lnLikelihood <- sum(model@lnFi);
  model@nbFreeParameter <- 2 * nbVariable;
  model@criterion <- -2 * model@lnLikelihood + model@nbFreeParameter * log(nbSample);
  model@modelName <- c("gaussian_pk_sjk");
  model
}
