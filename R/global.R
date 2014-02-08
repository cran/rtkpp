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
#' Extract parts of a rtkpp S4 class
#'
#' @param x object from which to extract element(s) or in which to replace element(s).
#' @param i the name of the element we want to extract or replace.
#' @param j if the element designing by i is complex, j specifying elements to extract or replace.
#' @param drop For matrices and arrays.  If TRUE the result is coerced to the lowest
#' possible dimension (see the examples).  This only works for extracting elements,
#' not for the replacement.  See drop for further details.
#' @param value	typically an array-like R object of a similar class as the element
#' of x we want to replace.
#' @name [
#' @aliases Extract
#' @docType methods
#' @rdname extract-methods
#'
NULL

#-----------------------------------------------------------------------
#' Print a rtkpp class to standard output.
#'
#' @param x a rtkpp object: a \code{\linkS4class{ClusterStrategy}},
#' a \code{\linkS4class{ClusterInit}} or a \code{\linkS4class{ClusterAlgo}}.
#' @param ... further arguments passed to or from other methods
#'
#' @return NULL. Prints to standard out.
#'
#' @name print
#' @rdname print-methods
#' @docType methods
#' @exportMethod print
#'
#' @seealso \code{\link{print}}
#' @examples
#'   ## for cluster strategy
#'   strategy <- clusterStrategy()
#'   print(strategy)
#'   ## for cluster init
#'   init <- clusterInit()
#'   print(init)
#'   ## for cluster algo
#'   algo <- clusterAlgo()
#'   print(algo)
#'
NULL

#-----------------------------------------------------------------------
#' Show description of a rtkpp class to standard output.
#'
#' @param object a rtkpp object: a \code{\linkS4class{ClusterStrategy}},
#' a \code{\linkS4class{ClusterInit}} or a \code{\linkS4class{ClusterAlgo}}.
#'
#' @return NULL. Prints to standard out.
#'
#' @importFrom methods show
#' @name show
#' @docType methods
#' @rdname show-methods
#' @exportMethod show
#'
#' @seealso \code{\link{show}}
#' @examples
#'   ## for strategy
#'   strategy <- clusterStrategy()
#'   show(strategy)
#'   ## for cluster init
#'   init <- clusterInit()
#'   show(init)
#'   ## for cluster algo
#'   algo <- clusterAlgo()
#'   show(algo)
NULL

#-----------------------------------------------------------------------
#' Produce summary of a rtkpp class. 
#'
#' @param object any cluster model deriving from a \code{\linkS4class{IClusterModelBase}} object.
#' @param ... further arguments passed to or from other methods
#'
#' @return NULL. Summaries to standard out.
#'
#' @name summary
#' @docType methods
#' @rdname summary-methods
#' @exportMethod summary
#'
#' @seealso \code{\link{summary}}
#' @examples
#'   data(geyser)
#'   model <- clusterDiagGaussian(geyser,3)
#'   summary(model)
#'
NULL

#' Return the missing values of a component or a cluster class.
#'
#' The missing methods allow the user to get the extrapolated mssing
#' values from a mixture model.
#'
#' @param x an object that can return the extrapolated missing values
#'
#' @return A matrix with three columns (row index, column index, value)
#'
#' @name missingValues
#' @docType methods
#' @rdname missingValues-methods
#' @exportMethod missingValues
#'
#' @examples
#'   data(geyser)
#'   model <- clusterDiagGaussian(geyser,3)
#'   missingValues(model)
setGeneric(
  name = "missingValues",
  function(x)
  { standardGeneric("missingValues")}
)


