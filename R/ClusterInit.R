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
#' Create an instance of [\code{\linkS4class{ClusterInit}}] class
#'
#' The initialization step is a two stages process: the proper initialization step
#' and some (optionnals) iterations of an algorithm [\code{\link{clusterAlgo}}].
#'
#' @details
#' There is three ways to initialize the parameters:
#' \itemize{
#'   \item \code{random} {The initial parameters of the mixture are chosen randomly.}
#'   \item \code{class}  {The initial membership of individuals are sampled randomly.}
#'   \item \code{fuzzy}  {The initial probabilities of membership of individuals are
#'                        sampled randomly.}
#' }
#' A few iteration of an algorithm [\code{\link{clusterAlgo}}] are then performed.
#' It is strongly recommended to use a few number of iterations of the \code{SEM}
#' or \code{CEM} algorithms after initialization. This allow to detect "bad"
#' initialization starting point of the estimation algorithm.
#'
#' These two stages are repeated until \code{nbInit} is reached. The initial
#' point with the best log-likelihood is conserved as the initial starting point.
#'
#' @param method Character string with the initialisation method.
#' Possible values: "random", "class", "fuzzy". Default value is "class".
#' @param algo Character string with the initialisation algorithm.
#' Possible values: "EM", "CEM", "SEM", "SemiSEM". Default value is "SEM".
#' @param nbInit integer defining the number of initialization point to test. Default value is 5.
#' @param nbIteration Integer defining the number of iteration in \code{algo}.
#' nbIteration must be a positive integer. Default values is 20. Not used if  \code{algo} = NULL.
#' @param epsilon Real defining the epsilon value for the algorithm. Default value: 0.01.
#'
#' @examples
#'  clusterInit(method = "class", nbInit=1, algo="CEM",nbIteration=50, epsilon=0.00001)
#'  clusterInit(nbIteration=0) # no algorithm
#'
#' @return a [\code{\linkS4class{ClusterInit}}] object
#' @author Serge Iovleff
#' @export
clusterInit <- function( method="class", nbInit=5,  algo = "SEM", nbIteration=20, epsilon=0.01)
{ return(new("ClusterInit", nbInit=nbInit, algo=new("ClusterAlgo", algo, nbIteration, epsilon)))}


#' Constructor of the [\code{\linkS4class{ClusterInit}}] class
#'
#' This class encapsulates the parameters of initialization methods of the
#' rtkpp Cluster estimation method.
#'
#' @slot method Character string with the initialization method to use. Default value: "random"
#' @slot nbInit Integer defining the number of initialization to perform. Default value: 5.
#' @slot algo An instance of \code{\linkS4class{ClusterAlgo}} class.
#' Default value: \code{clusterAlgo("SEM", 20, 0)}.
#'
#' @examples
#'   getSlots("ClusterInit")
#'   new("ClusterInit")
#'   new("ClusterInit", nbInit=1)
#'
#' @author Serge Iovleff
#'
#' @name ClusterInit
#' @rdname ClusterInit-class
#' @aliases ClusterInit-class
#' @exportClass ClusterInit
#'
setClass(
  Class="ClusterInit",
  slots=c(method="character", nbInit = "numeric", algo = "ClusterAlgo"),
  prototype=list(method="class", nbInit = 5, algo = clusterAlgo("SEM", 20, 0)),
  # validity function
  validity=function(object)
  {
    # for method
    if ( sum(object@method %in% c("random","class","fuzzy")) != 1 )
    {stop("Initialization method is not valid. See ?clusterInit for the list of available initialization method.")}
    # for nbInit
    if (round(object@nbInit)!=object@nbInit)
    {stop("nbIInit must be an integer.")}
    if( object@nbInit < 1 ) # can't be zero
    {stop("nbInit must be strictly greater than 0.");}
    # for algo
    if (!is.null(object@algo))
    {
      if(class(object@algo)[1] != "ClusterAlgo")
      {stop("algo is not of a Cluster algorithm (must be an instance of the class ClusterAlgo).")}
      if (!validObject(object@algo))
      {stop("algo is not of a valid algorithm. See ?clusterAlgo).")}
    }
    return(TRUE)
  }
)


#' Initialize an instance of a rtkpp class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterInit}}] class.
#' Used internally in the `rtkpp' package.
#'
#' @rdname initialize-methods
#' @keywords internal
#'
setMethod(
  f="initialize",
  signature=c("ClusterInit"),
  definition=function(.Object,method="class",nbInit = 5,algo= clusterAlgo("SEM", 20, 0))
  {
    # for method
    if(missing(method)) {.Object@method<-"class"}
    else  {.Object@method<-method}
    # for nbIteration
    if( missing(nbInit) ){ .Object@nbInit<-5 }
    else{.Object@nbInit<-nbInit}
    # for algo
    if(missing(algo)){ .Object@algo<-clusterAlgo("SEM", 20, 0.1) }
    else{.Object@algo<-algo}
    # validate
    validObject(.Object)
    return(.Object)
  }
)

#' @aliases print-init,ClusterInit,ClusterInit-method
#' @rdname print-methods
setMethod(
  f="print",
  signature=c("ClusterInit"),
  function(x,...){
    function(object){
      cat("****************************************\n")
      cat("*** Cluster init:\n")
      cat("* method               = ", object@method, "\n")
      cat("* number of init       = ", object@nbInit, "\n")
      cat("* algorithm            = ", object@algo@algo, "\n")
      cat("* number of iterations = ", object@algo@nbIteration, "\n")
      cat("* epsilon              = ", object@algo@epsilon, "\n")
      cat("****************************************\n")
    }
  }
)

#' @rdname show-methods
#' @aliases show-init,ClusterInit,ClusterInit-method
setMethod(
  f="show",
  signature=c("ClusterInit"),
  function(object){
    cat("****************************************\n")
    cat("*** Cluster init:\n")
    cat("* method              = ", object@method, "\n")
    cat("* number of init      = ", object@nbInit, "\n")
    cat("* algorithm            = ", object@algo@algo, "\n")
    cat("* number of iterations = ", object@algo@nbIteration, "\n")
    cat("* epsilon              = ", object@algo@epsilon, "\n")
    cat("****************************************\n")
  }
)

#' @rdname extract-methods
#' @aliases [,ClusterInit-method
setMethod(
  f="[",
  signature(x = "ClusterInit"),
  definition=function(x,i,j,drop){
    if ( missing(j) ){
      switch(EXPR=i,
        "method"={return(x@method)},
        "nbInit"={return(x@nbInit)},
        "algo"={return(x@algo)},
        stop("This attribute doesn't exist !")
        )
      }
    else
    {stop("This attribute is not a list !")}
  }
)

#' @name [
#' @rdname extract-methods
#' @aliases [<-,ClusterInit-method
setReplaceMethod(
  f="[",
  signature(x = "ClusterInit"),
  definition=function(x,i,j,value){
    if ( missing(j) )
    {
      switch(EXPR=i,
        "method"={x@method<-value},
        "nbInit"={x@nbInit<-value},
        "algo"={x@algo<-value},
        stop("This attribute doesn't exist !")
      )
    }
    else
    { stop("This attribute is not a list !")}
    validObject(x)
    return(x)
  }
)
