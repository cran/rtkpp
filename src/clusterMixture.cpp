/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  rtkpp
 * created on: 28 juil. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file clusterMixture.cpp
 *  @brief In this file we launch the computation for estimating a mixture model.
 **/


#include "RTKpp.h"

using namespace STK;

/** @param model ClusterDiagModel S4 class
 *  @param nbCluster a vector with the number of clusters to test
 */
RcppExport SEXP clusterMixture( SEXP model, SEXP nbCluster, SEXP modelNames, SEXP strategy, SEXP critName )
{
  BEGIN_RCPP
  // create a launcher
  ClusterLauncher launcher(model, nbCluster, modelNames, strategy, critName);
  // return result
  return Rcpp::wrap(launcher.run());

  END_RCPP
}

/** @param model ClusterDiagModel S4 class
 *  @param nbCluster a vector with the number of clusters to test
 */
RcppExport SEXP clusterMixtureHeterogene( SEXP model, SEXP nbCluster, SEXP strategy, SEXP critName )
{
  BEGIN_RCPP
  // create a launcher
  ClusterLauncher launcher(model, nbCluster, strategy, critName);
  // return result
  return Rcpp::wrap(launcher.run());

  END_RCPP
}
