#-----------------------------------------------------------------------
#     Copyright (C) 2012-2014  Inria
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
#' @include ClusterStrategy.R
NULL

#' Interface Class [\code{\linkS4class{IClusterModel}}] for Cluster models.
#'
#' This class encapsulate the common parameters of all the Cluster models.
#'
#' A Cluster model is a model of the form
#' \deqn{
#'   f({x}|\boldsymbol{\theta})
#'     \sum_{k=1}^K p_k h({x};\boldsymbol{\lambda}_k,\boldsymbol{\alpha})
#'    \quad {x} \in J.
#' }
#'where h can be either a pdf or a discrete probability.
#'
#' @slot data      \code{\link{matrix}} of size \eqn{n \times p} with the data set to cluster.
#' If the original data set had NA values, they have been estimated in the estimation process.
#' @slot nbCluster Integer with the number of cluster of the model.
#' @slot pk        Vector of size K with the proportions of each mixture.
#' @slot tik       Matrix of size \eqn{n \times K} with the posterior probability of
#' the ith individual to belong to kth cluster.
#' @slot lnFi        Vector of size n with the log-likelihood of the ith individuals.
#' @slot zi        Vector of integer of size n  with the attributed class label of the individuals.
#' @slot missings   \code{\link{matrix}} of two columns with the indexes (i,j) of the missing values.
#' @slot lnLikelihood Real given the ln-liklihood of the Cluster model.
#' @slot criterion Real given the value of the AIC, BIC or ICL criterion.
#' @slot nbFreeParameter Integer given the number of free parameters of the model.
#' @slot modelName mixture model name.
#' @slot strategy  the instance of the [\code{\linkS4class{ClusterStrategy}}] used in the
#' estimation process of the mixture.
#'
#' @examples
#'   getSlots("IClusterModel")
#'
#' @author Serge Iovleff
#'
#' @name IClusterModel
#' @rdname ClusterModels-class
#' @aliases IClusterModel-class
#' @exportClass IClusterModel
setClass(
  Class="IClusterModel",
  representation( data = "matrix"
                , nbCluster = "numeric"
                , pk = "numeric"
                , tik = "matrix"
                , lnFi = "numeric"
                , zi = "integer"
                , missings = "matrix"
                , lnLikelihood = "numeric"
                , criterion = "numeric"
                , nbFreeParameter = "numeric"
                , modelName = "character"
                , strategy = "ClusterStrategy"
                , "VIRTUAL"
                ),
  prototype=list( data = matrix(nrow=0,ncol=0)
                , nbCluster = 0
                , pk = vector("numeric")
                , tik = matrix(nrow=0, ncol=0)
                , lnFi  = vector("numeric")
                , zi  = vector("integer")
                , missings = matrix(nrow=0, ncol=2)
                , lnLikelihood = -Inf
                , criterion = -Inf
                , nbFreeParameter = 0
                , modelName = character(1)
                , strategy = clusterStrategy()
                ),
  # validity function
  validity=function(object)
  {
    nbSample  = nrow(object@data)
    nbCluster = object@nbCluster
    # check nbCluster
    if (round(nbCluster)!=object@nbCluster)
    {stop("nbCluster must be an integer.")}
    if( nbCluster < 2 )
    {  stop("nbCluster must be greater than 1.")}
    # check pk
    if (length(object@pk) != nbCluster)
    {stop("pk must have length nbCluster.")}
    # check tik
    if (ncol(object@tik) != nbCluster)
    {stop("tik must have nbCluster columns.")}
    if (nrow(object@tik) != nbSample)
    {stop("tik must have nbSample rows.")}
    # check lnFi
    if (length(object@lnFi) != nbSample)
    {stop("fi must have nbSample size.")}
    # check zi
    if (length(object@zi) != nbSample)
    {stop("zi must have nbSample size.")}
    # check nbFreeParameter
    if (round(object@nbFreeParameter)!=object@nbFreeParameter)
    {stop("nbFreeParameter must be an integer.")}
    return(TRUE)
  }
)

#-----------------------------------------------------------------------
#' Initialize an instance of a rtkpp class.
#'
#' Initialization method of the [\code{\linkS4class{IClusterModel}}] class.
#' Used internally in the `rtkpp' package.
#'
#' @rdname initialize-methods
#' @keywords internal
#'
setMethod(
    f="initialize",
    signature=c("IClusterModel"),
    definition=function(.Object, data, nbCluster, modelName)
    {
      # for data
      if(missing(data)) {stop("data is mandatory.")}
      .Object@data<-as.matrix(data)
      # for nbCluster
      if(missing(nbCluster)) {stop("nbCluster is mandatory.")}
      .Object@nbCluster<-nbCluster
      # for nbCluster
      if(missing(modelName)) {stop("modelName is mandatory.")}
      .Object@modelName<-modelName
      # resize
      nbSample    <- nrow(.Object@data)
      nbVariable  <- ncol(.Object@data)
      .Object@pk  <- rep(1/.Object@nbCluster, nbCluster)
      .Object@tik <- matrix(1/nbCluster, nbSample, nbCluster)
      .Object@lnFi  <- rep(0, nbSample)
      .Object@zi  <- as.integer(rep(1, nbSample))
      .Object@strategy <- clusterStrategy()
      # validObject(.Object) will be called at the end of the initialization process
      # in the derived classes
      return(.Object)
    }
)

###################################################################################
# @name [
# @docType methods
#' @rdname extract-methods
#' @aliases [,IClusterModel-method
setMethod(
    f="[",
    signature(x = "IClusterModel"),
    definition=function(x,i,j,drop){
      if ( missing(j) ){
        switch(EXPR=i,
            "data"={return(x@data)},
            "nbCluster"={return(x@nbCluster)},
            "pk"={return(x@pk)},
            "tik"={return(x@tik)},
            "lnFi"={return(x@lnFi)},
            "zi"={return(x@zi)},
            "missing"={return(x@missing)},
            "lnLikelihood"={return(x@lnLikelihood)},
            "criterion"={return(x@criterion)},
            "modelName"={return(x@modelName)},
            stop("This attribute doesn't exist !")
        )
      }else{
        stop("This attribute is not a list !")
      }
    }
)

#-----------------------------------------------------------------------
#' @rdname print-methods
#' @aliases print print,IClusterModel-method
#'
setMethod(
  f="print",
  signature=c("IClusterModel"),
  function(x,...)
  {
    cat("* nbCluster    = ", x@nbCluster, "\n")
    cat("* lnLikelihood = ", x@lnLikelihood,"\n")
    cat("* criterion    = ", x@criterion, "\n")
    cat("* model name   = ", x@modelName, "\n")
  }
)

#-----------------------------------------------------------------------
#' @rdname show-methods
#' @aliases show show,IClusterModel-method
setMethod(
    f="show",
    signature=c("IClusterModel"),
    function(object)
    {
      if(length(object@data)!=0)
      {
        nrowShow <- min(10,nrow(object@data))
        ncolShow <- min(10,ncol(object@data))
        cat("* data (limited to 10 samples and 10 variables) =\n")
        print(format(object@data[1:nrowShow,1:ncolShow]),quote=FALSE)
      }
      cat("* ... ...\n")
      cat("* nbCluster    = ", object@nbCluster, "\n")
      cat("* lnLikelihood = ", object@lnLikelihood,"\n")
      cat("* criterion    = ", object@criterion, "\n")
      cat("* model name   = ", object@modelName, "\n")
    }
)

#-----------------------------------------------------------------------
#' @rdname summary-methods
#' @aliases summary summary,IClusterModel-method
setMethod(
    f="summary",
    signature=c("IClusterModel"),
    function(object,...)
    {
      cat("* nbSample       = ", nrow(object@data), "\n")
      cat("* nbVariable     = ", ncol(object@data), "\n")
      cat("* nbCluster      = ", object@nbCluster, "\n")
      cat("* lnLikelihood   = ", object@lnLikelihood,"\n")
      cat("* criterion      = ", object@criterion, "\n")
      cat("* nbFreeParameter= ", object@nbFreeParameter, "\n")
      cat("* model name     = ", object@modelName, "\n")
    }
)
