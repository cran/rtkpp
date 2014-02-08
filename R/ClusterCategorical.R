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
#' @include ClusterModelNames.R IClusterModel.R global.R
NULL

#-----------------------------------------------------------------------
#' Create an instance of the [\code{\linkS4class{ClusterCategorical}}] class
#'
#' This function computes the optimal Categorical mixture model according
#' to the [\code{criterion}] among the list of model given in [\code{modelNames}]
#' and the number of clusters given in [\code{nbCluster}], using the strategy
#' specified in [\code{strategy}].
#'
#' @param data frame or matrix containing the data. Rows correspond to observations
#' and columns correspond to variables. If the data set contains NA values, they
#' will be estimated during the estimation process.
#' @param nbCluster  [\code{\link{vector}}] listing the number of clusters to test.
#' @param modelNames [\code{\link{vector}}] of model names to run. By default
#' the categorical models "categorical_pk_pjk" and "categorical_p_pjk" are estimated.
#' @param strategy a [\code{\linkS4class{ClusterStrategy}}] object containing
#' the strategy to run. clusterStrategy() method by default.
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL". Default is "ICL".
#'
#' @examples
#' ## A quantitative example with the birds data set
#' data(birds)
#' ## add 10 missing values
#' x = birds;
#' x[round(runif(5,1,nrow(birds))), 2] <- NA
#' x[round(runif(5,1,nrow(birds))), 4] <- NA
#' ## with default values
#' model <- clusterCategorical(data=x, nbCluster=2:3, strategy = clusterFastStrategy())
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
#' missingValues(model)
#'
#' @return An instance of the [\code{\linkS4class{ClusterCategorical}}] class.
#' @author Serge Iovleff
#' @export
#'
clusterCategorical <- function(data, nbCluster=2, modelNames=c( "categorical_pk_pjk", "categorical_p_pjk"), strategy=clusterFastStrategy(), criterion="ICL")
{
  # check nbCluster
  nbClusterModel = length(nbCluster);
  nbClusterMin   = min(nbCluster);
  nbClusterMax   = max(nbCluster);
  if (nbClusterMin < 2) { stop("The number of clusters must be greater or equal to 2\n")}

  # check criterion
  if(sum(criterion %in% c("BIC","AIC", "ICL")) != 1)
  { stop("criterion is not valid. See ?clusterCategorical for the list of valid criterion\n")}

  # get data
  data <- as.matrix(data);
  if (nrow(data) <= 3*nbClusterMax) {stop("There is not enough individuals (rows) in the data set\n")}
  if (ncol(data) <= 1) {stop("Error: empty data set or not enough columns (must be greater than 1 for Categorical variables)\n")}

  # check modelNames
  if (is.null(modelNames)) { modelNames = c( "categorical_pk_pjk", "categorical_p_pjk")}
  if (!clusterValidCategoricalNames(modelNames))
  { stop("modelNames is not valid. See ?CategoricalNames for the list of valid model names\n")}

  # check strategy
  if(class(strategy)[1] != "ClusterStrategy")
  {stop("strategy is not a Cluster Strategy class (must be an instance of the class ClusterStrategy).\n")}
  validObject(strategy);

  # Create model
  model = new("ClusterCategorical", data)
  model@strategy = strategy;

  # start estimation of the models
  resFlag = .Call("clusterMixture", model, nbCluster, modelNames, strategy, criterion, PACKAGE="rtkpp")
  # set names
  # dimnames(model@plkj) <- list(NULL, colnames(model@component@data), NULL) # not working
  # should be done on the C++ side
  if (resFlag != 1) {cat("WARNING: An error occur during the clustering process\n");}
  dim(model@component@plkj) <- c(model@component@nbModalities, model@nbCluster, ncol(data));
  model
}

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterCategoricalComponent}}] class
#'
#' This class defines a categorical component of a mixture Model.
#'
#' @slot plkj Array with the probability of the jth variable in the kth cluster to be l.
#' @slot nbModalities Integer with the (maximal) number of modalities of the categorical  data.
#' @seealso [\code{\linkS4class{IClusterComponent}}] class
#'
#' @examples
#' getSlots("ClusterCategoricalComponent")
#'
#' @author Serge Iovleff
#'
#' @name ClusterCategoricalComponent
#' @rdname IClusterComponent-class
#' @aliases ClusterCategoricalComponent-class
#' @exportClass ClusterCategoricalComponent
#'
setClass(
  Class="ClusterCategoricalComponent",
  representation( plkj = "array", nbModalities = "numeric"),
  contains=c("IClusterComponent"),
  validity=function(object)
  {
    dims <- dim(object@plkj)

    if (round(object@nbModalities)!=object@nbModalities)
    {stop("nbModalities must be an integer.")}
    if (dims[1]!=object@nbModalities)
    {stop("First dimension in plkj must be nbModalities.")}
    if (!clusterValidCategoricalNames(object@modelName))
    {stop("Invalid Categorical model name.")}
    return(TRUE)
  }
)
#' Initialize an instance of a rtkpp class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterCategoricalComponent}}] class.
#' Used internally in the `rtkpp' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterCategoricalComponent"),
    definition=function(.Object, data=matrix(nrow=0, ncol=0), nbCluster=2, modelName="categorical_pk_pjk")
    {
      # for data
      if(missing(data)) {stop("data is mandatory in ClusterCategoricalComponent.")}
      # check model name
      if (!clusterValidCategoricalNames(modelName)) { stop("modelName is invalid");}
      # prepare data and compute number of modalities
      data = as.data.frame(data)
      nbModalities = 0;
      for ( j in 1:ncol(data) )
      {
        nbModalities <- max(nbModalities, nlevels(factor(data[,j])))
        data[,j] <- as.integer(factor(data[,j]))
      }
      # initialize base class
      .Object <- callNextMethod(.Object, data, modelName);
      .Object@nbModalities = nbModalities;
      # create slots
      nbVariable = ncol(.Object@data);
      .Object@plkj <- array(data = 1/nbModalities, dim=c(nbModalities,nbCluster,nbVariable))
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname extract-methods
#' @aliases [,ClusterCategoricalComponent-method
setMethod(
    f="[",
    signature(x = "ClusterCategoricalComponent"),
    definition=function(x, i, j, drop)
    {
      if ( missing(j) )
      {
        switch(EXPR=i,
            "plkj"={return(x@plkj)},
            stop("This attribute doesn't exist !")
        )
      }
      else
      {
        if (!is.numeric(j)) {stop("j must be an integer.")}
        if (round(j)!=j)    {stop("j must be an integer.")}
        switch(EXPR=i,
            "plkj"={return(x@plkj[,j,])},
            stop("This attribute doesn't exist !")
        )
      }
    }
)

#' @param k the number of the cluster to print
#' @rdname print-methods
#' @aliases print print,ClusterCategoricalComponent-method
setMethod(
    f="print",
    signature=c("ClusterCategoricalComponent"),
    function(x,k,...)
    {
      cat("* probabilities = \n");
      print(format(x@plkj[,k,]), quote=FALSE);
    }
)

#' @rdname show-methods
#' @aliases show-ClusterCategoricalComponent,ClusterCategoricalComponent,ClusterCategoricalComponent-method
setMethod(
    f="show",
    signature=c("ClusterCategoricalComponent"),
    function(object)
    {
      cat("* probabilities = \n");
      print(format(object@plkj), quote=FALSE);
    }
)

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterCategorical}}] class
#'
#' This class defines a categorical mixture Model. Inherits from the
#'[\code{\linkS4class{IClusterModelBase}}] class. A categorical mixture model is
#' a mixture model of the form
#'
#' \deqn{
#'   f({x}|\boldsymbol{\theta}) \\
#'   =\sum_{k=1}^K p_k \prod_{j=1}^d \mathcal{M}(x_j;p_{jk},1) \\
#'    \quad {x} \in \{1,\ldots,L\}^d.
#' }
#'
#' @slot component  A [\code{\linkS4class{ClusterCategoricalComponent}}] with the
#' probabilities of the categorical component
#'
#' @examples
#'   getSlots("ClusterCategorical")
#'   data(birds)
#'   new("ClusterCategorical", data=birds)
#'
#' @author Serge Iovleff
#'
#' @name ClusterCategorical-class
#' @rdname ClusterCategorical-class
#' @aliases ClusterCategorical-class
#' @exportClass ClusterCategorical
setClass(
  Class="ClusterCategorical",
  representation( component = "ClusterCategoricalComponent"),
  contains=c("IClusterModelBase"),
  validity=function(object)
  {
    dims <- dim(object@component@plkj)

    if (round(object@component@nbModalities)!=object@component@nbModalities)
    {stop("nbModalities must be an integer.")}
    if (dims[1]!=object@component@nbModalities)
    {stop("First dimension in plkj must be nbModalities.")}
    if (dims[2]!=object@nbCluster)
    {stop("Second dimension in plkj must be nbCluster.")}
    if (dims[3]!=ncol(object@component@data))
    {stop("Third dimension in plkj must be nbCluster.")}
    if (!clusterValidCategoricalNames(object@component@modelName))
    {stop("Invalid Categorical model name.")}
    return(TRUE)
  }
)

#' Initialize an instance of a rtkpp class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterCategorical}}] class.
#' Used internally in the `rtkpp' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterCategorical"),
    definition=function(.Object, data=matrix(nrow=0, ncol=0), nbCluster=2, modelName="categorical_pk_pjk")
    {
      # for data
      if(missing(data)) {stop("data is mandatory in ClusterCategorical.")}
      .Object@component = new("ClusterCategoricalComponent", data, nbCluster, modelName);
      .Object <- callNextMethod(.Object, nrow(.Object@component@data), nbCluster);
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname print-methods
#' @aliases print print,ClusterCategorical-method
#'
setMethod(
  f="print",
  signature=c("ClusterCategorical"),
  function(x,...)
  {
    cat("****************************************\n")
    callNextMethod();
    cat("****************************************\n")
    for(k in 1:length(x@pk))
    {
      cat("*** Cluster: ",k,"\n")
      cat("* Proportion = ", format(x@pk[k]), "\n");
      print(x@component,k);
      cat("* probabilities = \n"); print(format(x@component@plkj[,k,]), quote=FALSE)
      cat("****************************************\n")
    }
  }
)

#' @rdname show-methods
#' @aliases show-ClusterCategorical,ClusterCategorical,ClusterCategorical-method
setMethod(
  f="show",
  signature=c("ClusterCategorical"),
  function(object)
  {
    cat("****************************************\n")
    callNextMethod();
    show(object@component);
    cat("****************************************\n")
    for(k in 1:length(object@pk))
    {
      cat("*** Cluster: ",k,"\n")
      cat("* Proportion = ", format(object@pk[k]), "\n")
      print(object@component,k);
      cat("****************************************\n")
    }
  }
)

#' @rdname summary-methods
#' @aliases summary summary,ClusterCategorical-method
setMethod(
  f="summary",
  signature=c("ClusterCategorical"),
  function(object, ...)
  {
    cat("**************************************************************\n")
    callNextMethod()
    summary(object@component);
    cat("* nbModalities   = ", format(object@component@nbModalities), "\n")
    cat("**************************************************************\n")
  }
)

#' Plotting of a class [\code{\linkS4class{ClusterCategorical}}]
#'
#' Plotting data from a [\code{\linkS4class{ClusterCategorical}}] object
#' using the estimated parameters and partition.
#'
#' @param x an object of class [\code{\linkS4class{ClusterCategorical}}]
#' @param y a number between 1 and K-1.
#' @param ... further arguments passed to or from other methods
#'
#' @importFrom graphics plot
#' @aliases plot-ClusterCategorical
#' @docType methods
#' @rdname plot-ClusterCategorical-method
#' @export
#'
#' @seealso \code{\link{plot}}
#' @examples
#'   ## the car data set
#'   data(car)
#'   model <- clusterCategorical(car, 3, strategy = clusterFastStrategy())
#' \dontrun{
#'   plot(model)
#'   }
#'
setMethod(
    f="plot",
    signature=c("ClusterCategorical"),
    function(x, y, ...)
    {
      # total number of cluster in the data set
      nbCluster = ncol(x@tik);
      # check y, no y => display all dimensions
      if (missing(y)) { y=1:(nbCluster-1); }
      else
      { if (round(y)!=y) {stop("y must be an integer.")}
        if (y>nbCluster-1)
        stop("y should not be greater than K-1")
        y <- 1:y
      }
      # get representation
      Y=.visut(x@tik, nbCluster);
      if (nbCluster == 2) { ndim = 1;}
      else { ndim = ncol(Y);}
      # Compute gaussian statistics
      mean  <- matrix(0, nrow = nbCluster, ncol =ndim)
      sigma <- matrix(1, nrow = nbCluster, ncol =ndim)
      for (k in 1:nbCluster)
      {
        wcov = cov.wt(as.matrix(Y), x@tik[,k], method = "ML");
        mean[k,]  = wcov$center;
        sigma[k,] = sqrt(diag(wcov$cov))
      }
      # create gaussian model
      gauss<-new("ClusterDiagGaussian", Y, nbCluster = x@nbCluster)
      gauss@component@mean  = mean
      gauss@component@sigma = sigma
      gauss@pk   = x@pk
      gauss@tik  = x@tik
      gauss@lnFi = x@lnFi
      gauss@zi   = x@zi
      gauss@component@missing = x@component@missing
      gauss@lnLikelihood    = x@lnLikelihood
      gauss@criterion       = x@criterion
      gauss@nbFreeParameter = x@nbFreeParameter
      gauss@strategy        = x@strategy
      .clusterPlot(gauss, y, .dGauss,...);
    }
)

# get logisitic representation
.visut <- function(t, gp)
{ m <- min(t[,gp]);
  if (m==0) t[,gp] = t[,gp] + 1e-30;
  return(scale(log(sweep(t,1,t[,gp],FUN="/")+ 1e-30), center=TRUE, scale=FALSE)[,-gp])
}

