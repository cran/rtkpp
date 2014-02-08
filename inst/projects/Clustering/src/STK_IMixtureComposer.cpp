/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2012  Serge Iovleff

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_IMixtureComposer.h
 *  @brief In this file we implement the abstract base class for mixture models.
 **/

#include <cmath>
#ifdef STK_MIXTURE_DEBUG
#include "Arrays/include/STK_Display.h"
#endif
#include "../include/STK_IMixtureComposer.h"
#include "STatistiK/include/STK_Law_Categorical.h"
#include "STatistiK/include/STK_Stat_Functors.h"

namespace STK
{
IMixtureComposer::IMixtureComposer( int nbSample, int nbVariable, int nbCluster)
                                  : IStatModelBase(nbSample, nbVariable)
                                  , nbCluster_(nbCluster)
                                  , prop_(nbCluster), tik_(nbSample, nbCluster), tk_(nbCluster), zi_(nbSample)
                                  , state_(Clust::modelCreated_)
{ initializeMixtureParameters(); }

/* copy constructor */
IMixtureComposer::IMixtureComposer( IMixtureComposer const& model)
                                  : IStatModelBase(model)
                                  , nbCluster_(model.nbCluster_)
                                  , prop_(model.prop_)
                                  , tik_(model.tik_)
                                  , zi_(model.zi_)
                                  , state_(model.state_)
{}
/* destructor */
IMixtureComposer::~IMixtureComposer() {}

/* @brief Initialize the model before at its first use.
 *  This function can be overloaded in derived class for initialization of
 *  the specific model parameters. It should be called prior to any used of
 *  the class.
 *  @sa IMixture,MixtureBridge,MixtureComposer
 **/
void IMixtureComposer::initializeStep()
{
  // (re)initialize the parameters tik, zi etc.
  initializeMixtureParameters();
  // (re)initialize the likelihood (others parameters are not modified)
  setLnLikelihood(-Arithmetic<Real>::infinity());
  // compute proportions
  pStep();
  setState(Clust::modelInitialized_);
}

/* initialize randomly the labels zi of the model */
void IMixtureComposer::randomClassInit()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering IMixtureComposer::randomClassInit()\n");
#endif
  if (state() < 2) { initializeStep();}
  Law::Categorical law(prop_);
  for (int i = zi_.begin(); i< zi_.end(); ++i)
  { zi_.elt(i) = law.rand();}
  cStep();
  eStep();
  // model intialized
  setState(Clust::modelParamInitialized_);
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("IMixtureComposer::randomClassInit() done\n");
#endif
}

/* initialize randomly the posterior probabilities tik of the model */
void IMixtureComposer::randomFuzzyInit()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("Entering IMixtureComposer::randomFuzzyInit()\n");
#endif
  if (state() < 2) { initializeStep();}
  randomFuzzyTik();
  eStep();
  // model intialized
  setState(Clust::modelParamInitialized_);
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("IMixtureComposer::randomFuzzyInit() done\n");
#endif
}

/* cStep */
int IMixtureComposer::cStep()
{
  tik_ = 0.;
  for (int i=tik_.beginRows(); i < tik_.endRows(); i++)
  { tik_.elt(i, zi_[i]) = 1.;}
  // count the number of individuals in each class
  tk_.move(Stat::sum(tik_));
  return tk_.minElt();
}

/* simulate zi  */
int IMixtureComposer::sStep()
{
  // simulate zi
  for (int i = zi_.begin(); i< zi_.end(); ++i)
  { zi_.elt(i) = Law::Categorical::rand(tik_.row(i));}
  return cStep();
}
/* compute tik, default implementation. */
Real IMixtureComposer::eStep()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering IMixtureComposer::eStep()\n");
#endif
  Real sum = 0.; tk_ =0.;
  for (int i = tik_.beginRows(); i < tik_.endRows(); ++i)
  {
    // get maximal element of ln(x_i,\theta_k) + ln(p_k)
    Array2DPoint<Real> lnComp(tik_.cols());
    for (int k=tik_.beginCols(); k< tik_.endCols(); k++)
    { lnComp[k] = std::log(prop_[k])+lnComponentProbability(i,k);}
    int kmax;
    Real max = lnComp.maxElt(kmax);
    // set zi_
    zi_[i] = kmax;
    // compute sum_k p_k exp{lnCom_k - lnComp_kmax}
    tik_.row(i) = (lnComp - max).exp();
    Real sum2 = tik_.row(i).sum();
    tik_.row(i) /= sum2;
    // count the expected number of individuals in each class
    tk_ += tik_.row(i);
    // compute lnLikelihood
    sum += max + std::log(sum2);
  }
  setLnLikelihood(sum);
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("IMixtureComposer::eStep() done\n");
#endif
  return tk_.minElt();
}

/* @return the computed likelihood of the i-th sample.
 *  @param i index of the sample
 **/
Real IMixtureComposer::computeLnLikelihood(int i) const
{
  // get maximal value
  CArrayPoint<Real> lnComp(prop_.size());
  for (int k = prop_.begin(); k< prop_.end(); ++k)
  { lnComp[k] = std::log(prop_[k]) + lnComponentProbability(i, k);}
  int kmax;
  Real lnCompMax = lnComp.maxElt(kmax);
  // compute result
  return std::log((lnComp-lnCompMax).exp().sum())+lnCompMax;
}

/* @return the computed log-likelihood. */
Real IMixtureComposer::computeLnLikelihood() const
{
  Real res = 0.0;
  for (int i = tik().beginRows(); i< tik().endRows(); ++i)
  { res += computeLnLikelihood(i);}
  return res;
}

/* @return the computed ICL criteria. */
Real IMixtureComposer::computeICL() const
{
  Real res = 0.0;
  for (int j = tik().beginCols(); j< tik().endCols(); ++j)
  { res += (tik_.col(j) * (tik_.col(j)+1e-15).log()).sum();}

  return (- 2. * lnLikelihood() + nbFreeParameter() * lnNbSample() - 2. * res);
}

/* estimate the proportions and the parameters of the components of the
 *  model given the current tik/zi mixture parameters values.
 **/
void IMixtureComposer::mStep()
{ pStep();
  /* implement specific parameters estimation in concrete class. */
}

/* Compute prop using the ML estimate, default implementation. */
void IMixtureComposer::pStep()
{ prop_ = Stat::mean(tik_);}

/* Compute Zi using the Map estimate, default implementation. */
void IMixtureComposer::mapStep()
{
  for (int i = zi_.begin(); i< zi_.end(); ++i)
  {
    int k;
    tik_.row(i).maxElt(k);
    zi_.elt(i) = k;
  }
}

/* Create the parameters of the  mixture model. */
void IMixtureComposer::initializeMixtureParameters()
{
  prop_ = 1./nbCluster_;
  tik_  = 1./nbCluster_;
  tk_   = Real(nbSample())/nbCluster_;
  zi_   = baseIdx;
}

/* generate random tik_ */
void IMixtureComposer::randomFuzzyTik()
{
  RandBase generator;
  tk_ = 0.;
  for (int i = tik_.beginRows(); i < tik_.endRows(); ++i)
  {
    // create a reference on the i-th row
    Array2DPoint<Real> tikRowi(tik_.row(i), true);
    generator.randUnif(tikRowi);
    tikRowi = tikRowi * prop_;
    tikRowi /= tikRowi.sum();
    tk_ += tikRowi;
  }
}

} // namespace STK

