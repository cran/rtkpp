/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015 Serge Iovleff

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: Oct 24, 2013
 * Author:   Serge Iovleff
 **/

/** @file STK_Poisson_lk.h
 *  @brief In this file we implement the Poisson_lk class
 **/

#ifndef STK_POISSON_LK_H
#define STK_POISSON_LK_H

#include "STK_PoissonBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class Poisson_lk;

namespace Clust
{
/** @ingroup Clustering
 *  Traits class for the Poisson_lk traits policy. */
template<class _Array>
struct MixtureTraits< Poisson_lk<_Array> >
{
  typedef _Array                  Array;
  typedef typename Array::Type    Type;
  typedef ParametersHandler<Clust::Poisson_lk_> ParamHandler;
};

} // namespace hidden

/** Specialization of the ParametersHandler struct for Poisson_lk model */
template <>
struct ParametersHandler<Clust::Poisson_lk_>: public PoissonHandlerBase<  ParametersHandler<Clust::Poisson_lk_> >
{
  /** Array of the rates */
  MixtureParametersSet<Real> lambda_;
  /** @return the value of lambda of the kth cluster and jth variable */
  inline Real const& lambdaImpl(int k, int j) const { return lambda_[k];}
  /** copy operator */
  inline ParametersHandler& operator=( ParametersHandler const& other)
  { lambda_ = other.lambda_; return *this; }
  /** copy operator using an array/expression storing the values */
  template<class Array>
  inline ParametersHandler& operator=( ExprBase<Array> const& param)
  {
    for (int k= param.beginRows(); k < param.endRows(); ++k)
    { lambda_[k] = param.row(k).mean();}
    return *this;
  }

  /** default constructor. All lambdas are initialized to 1. */
  inline ParametersHandler( int nbCluster): lambda_(nbCluster){}
  /** copy constructor.
   * @param param the parameters to copy.
   **/
  inline ParametersHandler( ParametersHandler const& param): lambda_(param.lambda_){}
  /** Initialize the parameters with an array/expression of value */
  template<class Array>
  inline ParametersHandler(int nbCluster, ExprBase<Array> const& param): lambda_(nbCluster)
  {
    for (int k= param.beginRows(); k < param.endRows(); ++k)
    { lambda_[k] = param.row(k).mean();}
  }
  /** destructor */
  inline ~ParametersHandler() {}
  /** Initialize the parameters of the model.
   *  This function initialize the parameters and the statistics.
   **/
  inline void resize(Range const& range)
  { lambda_.initialize(1.);}
  /** Store the intermediate results of the Mixture.
   *  @param iteration Provides the iteration number beginning after the burn-in period.
   **/
  inline void storeIntermediateResults(int iteration)
  { lambda_.storeIntermediateResults(iteration);}
  /** Release the stored results. This is usually used if the estimation
   *  process failed.
   **/
  inline void releaseIntermediateResults()
  { lambda_.releaseIntermediateResults();}
  /** set the parameters stored in stat_proba_ and release stat_proba_. */
  inline void setParameters() { lambda_.setParameters();}
};

/** @ingroup Clustering
 *  The Poisson mixture model @c Poisson_lk has a probability function of the form
 * \f[
 *    P(\mathbf{x}=(n_1,\ldots,n_d)|\theta)
 *     = \sum_{k=1}^K p_k \prod_{j=1}^d e^{-\lambda_{k}} \frac{\lambda_{k}^{n_j}}{n_j!}.
 * \f]
 **/
template<class Array>
class Poisson_lk : public PoissonBase<Poisson_lk<Array> >
{
  public:
    typedef PoissonBase<Poisson_lk<Array> > Base;
    using Base::p_tik;
    using Base::p_nk;
    using Base::param_;
    using Base::p_data;
    using Base::nbVariable;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline Poisson_lk( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline Poisson_lk( Poisson_lk const& model) : Base(model) {}
    /** destructor */
    inline ~Poisson_lk() {}
    /** @return the value of lambda of the kth cluster and jth variable */
    inline Real lambdaImpl(int k, int j) const { return param_.lambda_[k];}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0., lambda = param_.lambda_.param_[k];
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      { sum += Law::Poisson::lpdf(p_data()->elt(i,j), lambda);}
      return sum;
    }
    /** Initialize randomly the parameters of the Poisson mixture. */
    void randomInit();
    /** Compute the weighted probabilities. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const { return this->nbCluster();}
};

/* Initialize randomly the parameters of the Poisson mixture. */
template<class Array>
void Poisson_lk<Array>::randomInit()
{
  Real m = p_data()->template cast<Real>().mean();
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  { param_.lambda_[k] = Law::Exponential::rand(m);}
}


/* Compute the modalities probabilities */
template<class Array>
bool Poisson_lk<Array>::mStep()
{
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    param_.lambda_[k]= (p_data()->transpose() * p_tik()->col(k)).sum()
                      /(p_data()->sizeCols()*p_nk()->elt(k));
  }
  return true;
}

} // namespace STK

#endif /* STK_POISSON_LK_H */
