/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff, University Lille 1, Inria

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 4 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_MixtureFacade.cpp
 *  @brief In this file we implement the ClusterFacade class.
 **/

#include "RTKpp.h"
#include "Clustering/include/STK_MixtureInit.h"
#include "Clustering/include/STK_MixtureAlgo.h"
#include "Clustering/include/STK_MixtureStrategy.h"

namespace STK
{

ClusterFacade::~ClusterFacade() { if (p_strategy_) delete p_strategy_;}


/* create a FullStrategy */
void ClusterFacade::createFullStrategy(Rcpp::S4 R_strategy)
{
  // get fields of the strategies
  int nbTry = R_strategy.slot("nbTry");
  int nbShortRun = R_strategy.slot("nbShortRun");
  Rcpp::S4 R_initMethod = R_strategy.slot("initMethod");
  Rcpp::S4 R_shortAlgo = R_strategy.slot("shortAlgo");
  Rcpp::S4 R_longAlgo = R_strategy.slot("longAlgo");

  // get fields of the initMethod
  std::string method = R_initMethod.slot("method");
  Clust::initType init = Clust::stringToInit(method);
  int nbInitRun = R_initMethod.slot("nbInitRun");
  Rcpp::S4 R_initAlgo = R_initMethod.slot("algo");

  // get fields of the initAlgo
  std::string initAlgoName = R_initAlgo.slot("algo");
  Clust::algoType initAlgo = Clust::stringToAlgo(initAlgoName);
  int nbInitIter = R_initAlgo.slot("nbIteration");
  Real initEpsilon = R_initAlgo.slot("epsilon");

  // get fields of the shortAlgo
  std::string shortAlgoName = R_shortAlgo.slot("algo");
  Clust::algoType shortAlgo = Clust::stringToAlgo(shortAlgoName);
  int nbShortIter = R_shortAlgo.slot("nbIteration");
  Real shortEpsilon = R_shortAlgo.slot("epsilon");

  // get fields of the longAlgo
  std::string longAlgoName = R_longAlgo.slot("algo");
  Clust::algoType longAlgo = Clust::stringToAlgo(longAlgoName);
  int nbLongIter = R_longAlgo.slot("nbIteration");
  Real longEpsilon = R_longAlgo.slot("epsilon");

  // create STK objects
  IMixtureInit* p_init = Clust::createInit(init, 1, initAlgo, nbInitIter, initEpsilon);
  IMixtureAlgo* p_shortAlgo = Clust::createAlgo(shortAlgo, nbShortIter, shortEpsilon);
  IMixtureAlgo* p_longAlgo = Clust::createAlgo(longAlgo, nbLongIter, longEpsilon);
  p_strategy_ = Clust::createFullStrategy(p_model_, nbTry, nbInitRun, p_init, nbShortRun, p_shortAlgo, p_longAlgo);
}

bool ClusterFacade::run()
{
  bool flag = false;
  if (p_strategy_)
  {
    // just check if the model is fresh or has been used
    if (p_strategy_->run()) { flag = true;}
    else { msg_error_ = p_strategy_->error();}
    p_model_->imputationStep();
    p_model_->finalizeStep();
  }
  else
  { msg_error_ = STKERROR_NO_ARG(MixtureFacade::run(),strategy is not set);}
  return flag;
}


}  // namespace STK




