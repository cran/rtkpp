/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013 Serge Iovleff

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

/** @file STK_Categorical_pjk.h
 *  @brief In this file we implement the Categorical_pjk class
 **/

#ifndef STK_CATEGORICAL_PJK_H
#define STK_CATEGORICAL_PJK_H

#include "STK_CategoricalBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class Categorical_pjk;

namespace Clust
{
/** @ingroup Clustering
 *  Traits class for the Categorical_pjk traits policy. */
template<class Array_>
struct MixtureTraits< Categorical_pjk<Array_> >
{
  typedef Array_ Array;
  typedef ParametersHandler<Clust::Categorical_pjk_> ParamHandler;
};

} // namespace hidden


/** Specialization of the ParametersHandler struct for Categorical_pk model */
template <>
struct ParametersHandler<Clust::Categorical_pjk_>: public CategoricalHandlerBase<  ParametersHandler<Clust::Categorical_pjk_> >
{
  /** Vector and statistics of the probabilities */
  MixtureParametersSet<ArrayXX> proba_;
  /** @return the probability of the kth cluster, jth variable, lth modality */
  inline Real const& probaImpl(int k, int j, int l) const { return proba_[k](l,j);}
  /** @return the probability law of the kth cluster for the jth variable */
  inline VectorX probaImpl(int k, int j) const { return proba_[k].col(j);}
  /** copy operator */
  inline ParametersHandler& operator=( ParametersHandler const& other)
  { proba_ = other.proba_; return *this; }
  /** copy operator using an array/expression storing the values */
  template<class Array>
  inline ParametersHandler& operator=( ExprBase<Array> const& param)
  {
    int nbModalities = param.sizeRows()/proba_().size();
    for (int k1= proba_().begin(), k2= param.beginRows(); k1 < proba_().end(); k1++, k2+=nbModalities)
    {
      Range rangeModalities = proba_[k1].rows();
      for (int j= param.beginCols();  j< param.endCols(); ++j)
      {
        for (int l1 = rangeModalities.begin(), l2= 0; l1 < rangeModalities.end(); ++l1, l2++)
        { proba_[k1](l1, j) = param(k2 + l2, j);}
      }
    }
    return *this;
  }

  /** default constructor */
  ParametersHandler(int nbCluster): proba_(nbCluster) {}
  /** copy constructor */
  ParametersHandler(ParametersHandler const& model):proba_(model.proba_) {}
  /** Initialize the parameters with an array/expression of value */
  template<class Array>
  ParametersHandler( int nbCluster, ExprBase<Array> const& param)
                   : proba_(nbCluster)
  {
    int nbModalities = param.sizeRows()/nbCluster;
    proba_.resize(nbModalities, param.cols());
    for (int k2= param.beginRows(), k= param.beginRows(); k2 < param.endRows(); k2+=nbModalities, k++)
    {
      for (int j= param.beginCols();  j< param.endCols(); ++j)
      {
        for (int l= 0; l < nbModalities; ++l)
        { proba_[k](baseIdx+l, j) = param(k2 + l, j);}
      }
    }
  }

  /** Initialize the parameters of the model.
   *  This function initialize the parameter proba and the statistics.
   **/
  void resize(Range const& rangeModalities, Range const& rangeData )
  {
    proba_.resize(rangeModalities, rangeData);
    proba_.initialize(1./rangeModalities.size());
  }
  /** Store the intermediate results of the Mixture.
   *  @param iteration Provides the iteration number beginning after the burn-in period.
   **/
  inline void storeIntermediateResults(int iteration)
  { proba_.storeIntermediateResults(iteration);}
  /** Release the stored results. This is usually used if the estimation
   *  process failed.
   **/
  inline void releaseIntermediateResults()
  { proba_.releaseIntermediateResults();}
  /** set the parameters stored in stat_proba_ and release stat_proba_. */
  inline void setParameters() { proba_.setParameters();}
};

/** @ingroup Clustering
 *  The diagonal Categorical mixture model @c Categorical_pjk is
 *  the most general diagonal Categorical model and have a probability
 *  function of the form
 * \f[
 *    P(\mathbf{x}=(l_1,\ldots,l_d)|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d p_{kl_j}^j.
 * \f]
 **/
template<class Array>
class Categorical_pjk : public CategoricalBase<Categorical_pjk<Array> >
{
  public:
    typedef CategoricalBase<Categorical_pjk<Array> > Base;
    using Base::p_tik; using Base::param_;
    using Base::p_data;
    using Base::modalities_;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Categorical_pjk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Categorical_pjk( Categorical_pjk const& model): Base(model) {}
    /** destructor */
    inline ~Categorical_pjk() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0., prob;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      {
        if ( (prob = param_.proba_[k](p_data()->elt(i,j), j)) <= 0.) return -Arithmetic<Real>::infinity();
        sum += std::log(prob);
       }
      return sum;
    }
    /** Initialize randomly the parameters of the Categorical mixture. */
    void randomInit();
    /** Compute the weighted probabilities. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*((this->nbModalities_-1).sum());}
};

/* Initialize randomly the parameters of the Categorical mixture. */
template<class Array>
void Categorical_pjk<Array>::randomInit()
{
  for (int k = p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    param_.proba_[k].randUnif();
    for (int j=param_.proba_[k].beginCols(); j< param_.proba_[k].endCols(); ++j)
    {
      param_.proba_[k].col(j) /= param_.proba_[k].col(j).sum();
    }
  }
}


/* Compute the modalities probabilities */
template<class Array>
bool Categorical_pjk<Array>::mStep()
{
  for (int k = p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    param_.proba_[k] = 0.;
    for (int j = p_data()->beginCols(); j < p_data()->endCols(); ++j)
    {
      // count the number of modalities weighted by the tik
      for (int i = p_data()->beginRows(); i < p_data()->endRows(); ++i)
      { param_.proba_[k].col(j)[(*p_data())(i, j)] += (*p_tik())(i, k);}
      // normalize the probabilities
      param_.proba_[k].col(j) /= param_.proba_[k].col(j).sum();
    }
  }
  return true;
}

} // namespace STK

#endif /* STK_CATEGORICAL_PJK_H */
