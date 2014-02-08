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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  rtkpp
 * created on: 4 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_ClusterLauncher.cpp
 *  @brief In this file we implement the ClusterLauncher which
 *  construct properly a mixture model.
 **/


#include "RTKpp.h"

using namespace Rcpp;

namespace STK
{

/** facade design pattern.
 * The ClusterLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
ClusterLauncher::ClusterLauncher( SEXP model, SEXP nbCluster, SEXP modelNames, SEXP strategy, SEXP r_critName )
                                : s4_model_(model)
                                , s4_strategy_(strategy)
                                , v_nbCluster_(nbCluster)
                                , v_modelNames_(modelNames)
                                , critName_(Rcpp::as<std::string>(r_critName))
                                , handler_()
                                , manager_(handler_)
                                , p_composer_(0)
                                , idData_()
                                , idMixtureModel_(Clust::unknown_mixture_)
                                , isFreeProp_(false)
    {}
/* destructor. */
ClusterLauncher::~ClusterLauncher() { if (p_composer_) delete p_composer_;}

/* run the estimation */
bool ClusterLauncher::run()
{
  Real criter;
  criter = selectBestModel();
  if (criter == Arithmetic<Real>::max() || !Arithmetic<Real>::isFinite(criter))
    return false;
  // copy common part
  std::string idModelName;
  handler_.getIdModel( idData_, idModelName);
  idMixtureModel_ = Clust::stringToMixture(idModelName);
  s4_model_.slot("criterion")    = criter;
  s4_model_.slot("modelName")    = mixtureToString(idMixtureModel_, isFreeProp_);
  s4_model_.slot("nbCluster")    = p_composer_->nbCluster();
  s4_model_.slot("lnLikelihood") = p_composer_->lnLikelihood();
  s4_model_.slot("nbFreeParameter")= p_composer_->nbFreeParameter();
  s4_model_.slot("pk")           = wrap(p_composer_->pk());
  s4_model_.slot("tik")          = wrap(p_composer_->tik());
  s4_model_.slot("zi")           = wrap(p_composer_->zi());
  NumericVector fi = s4_model_.slot("lnFi");
  NumericVector zi = s4_model_.slot("zi");
  for (int i=0; i< fi.length(); ++i)
  {
    fi[i] = p_composer_->computeLnLikelihood(i);
    zi[i] += (1 - baseIdx);  // set base 1 for the class labels
  }
  // get specific parameters
  getParameters();
  return true;
}

/* get the parameters */
Real ClusterLauncher::selectBestModel()
{
  // wrap data matrix with Rcpp and wrap Rcpp matrix with STK++ matrix
  NumericMatrix m_data = s4_model_.slot("data");
  Real criter = s4_model_.slot("criterion");
  RcppMatrix<double> data(m_data);

  int nbSample   = m_data.rows();
  IMixtureComposer* p_current =0;
  IMixtureCriterion* p_criterion =0;

  try
  {
    // check if a model is free prop or not and add data to handler (no copy)
    Array1D<bool> v_free;
    for (int l= 0; l <v_modelNames_.size() ; ++l)
    {
      std::string idData = "model" + typeToString<int>(l);
      std::string idModel(as<std::string>(v_modelNames_[l]));
      // transform R model names to STK++ model names
      // check have been done on the R side so.... Let's go
      bool freeProp;
      Clust::Mixture model = Clust::stringToMixture(idModel, freeProp);
      handler_.addData(m_data, idData, Clust::mixtureToString(model));
      v_free.push_back(freeProp);
    }
    // create criterion
    if (critName_ == "BIC") { p_criterion = new BICMixtureCriterion();}
    if (critName_ == "AIC") { p_criterion = new AICMixtureCriterion();}
    if (critName_ == "ICL") { p_criterion = new ICLMixtureCriterion();}

    // start the estimation process, should end with the best model according to
    // the criteria
    p_composer_ = 0;
    for (int k=0; k <v_nbCluster_.length(); ++k)
    {
      int K = v_nbCluster_[k];
      for (int l=0; l <v_modelNames_.size(); ++l)
      {
        // create composer
        if (v_free[l]) { p_current = new MixtureComposer(nbSample, K);}
        else           { p_current = new MixtureComposerFixedProp(nbSample, K);}
        // create current mixture and register it
        std::string idData = "model" + typeToString<int>(l);
        static_cast<MixtureComposer*>(p_current)->createMixture(manager_, idData);

        // create facade and strategy
        ClusterFacade facade(p_current);
        facade.createFullStrategy(s4_strategy_);
        // run estimation and get results if possible
        if (facade.run())
        {
          // compute criterion and update model if necessary
          p_criterion->setModel(p_current);
          if (!p_criterion->run()) { delete p_current; p_current = 0;}
          else if (criter > p_criterion->value())
          {
            if (p_composer_) { std::swap(p_current, p_composer_);}
            else             { p_composer_ = p_current; p_current = 0;}
            idData_     = idData;
            isFreeProp_ = v_free[l];
            criter = p_criterion->value();
          }
        }
        // release current composer
        if (p_current) { delete p_current; p_current = 0;}
      }
    }
    // release
    delete p_criterion;
    return criter;
  } catch (Exception const& e)
  {
    if (p_current) delete p_current;
    if (p_criterion) delete p_criterion;
    ::Rf_error(e.error().c_str()) ;
  }
  // failed
  return Arithmetic<Real>::max();
}

void ClusterLauncher::getParameters()
{
  Clust::MixtureClass mixtClass = Clust::MixtureToMixtureClass(idMixtureModel_);
  switch (mixtClass)
  {
    case Clust::Gaussian_:
      getDiagGaussianParameters();
      break;
    case Clust::Gamma_:
      getGammaParameters();
      break;
    case Clust::Categorical_:
      getCategoricalParameters();
      break;
    default:
      break;
  }
}

/* get the diagonal Gaussian parameters */
void ClusterLauncher::getDiagGaussianParameters()
{
  // get parameters
  Array2D<Real> params;
  static_cast<MixtureComposer*>(p_composer_)->getParameters(manager_,idData_, params);
  // get dimensions
  int K = params.sizeRows()/2, nbVariable = params.sizeCols();
  // get results
  Array2D<Real> mean(K, nbVariable), sigma(K, nbVariable);
  for (int k=0; k<K; ++k)
  {
    mean.row(k)   = params.row(2*k);
    sigma.row(k) = params.row(2*k+1);
  }
  // save results in s4_model
  s4_model_.slot("mean")   = wrap(mean);
  s4_model_.slot("sigma") = wrap(sigma);
  // get data
  RcppMatrix<double> m_data;
  static_cast<MixtureComposer*>(p_composer_)->getData(manager_, idData_, m_data);
  s4_model_.slot("data") = (Rcpp::Matrix< RcppMatrix<double>::Rtype_>)m_data;
}

/* get the gamma parameters */
void ClusterLauncher::getGammaParameters()
{
  // get parameters
  Array2D<Real> params;
  static_cast<MixtureComposer*>(p_composer_)->getParameters(manager_,idData_, params);
  // get dimensions
  int K = params.sizeRows()/2, nbVariable = params.sizeCols();
  // get results
  Array2D<Real> shape(K, nbVariable), scale(K, nbVariable);
  for (int k=0; k<K; ++k)
  {
    shape.row(k) = params.row(2*k);
    scale.row(k) = params.row(2*k+1);
  }
  // save results in s4_model
  s4_model_.slot("shape") = wrap(shape);
  s4_model_.slot("scale") = wrap(scale);
  // get data
  RcppMatrix<double> m_data;
  static_cast<MixtureComposer*>(p_composer_)->getData(manager_, idData_, m_data);
  s4_model_.slot("data") = (Rcpp::Matrix< RcppMatrix<double>::Rtype_>)m_data;
}

/* get the diagonal Categorical parameters */
void ClusterLauncher::getCategoricalParameters()
{
  // get parameters
  Array2D<Real> params;
  static_cast<MixtureComposer*>(p_composer_)->getParameters(manager_,idData_, params);
  params.shift(0,0);
  // save results in s4_model
  s4_model_.slot("plkj") = wrap(params);
  // get data
  RcppMatrix<int> m_data;
  static_cast<MixtureComposer*>(p_composer_)->getData(manager_, idData_, m_data);
  s4_model_.slot("data") = (Rcpp::Matrix< RcppMatrix<int>::Rtype_>) m_data;
}

/* select best model*/
}  // namespace STK

