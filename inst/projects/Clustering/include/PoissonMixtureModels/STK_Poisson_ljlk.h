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

/** @file STK_Poisson_ljlk.h
 *  @brief In this file we implement the Poisson_ljlk class
 **/

#ifndef STK_POISSON_LJLK_H
#define STK_POISSON_LJLK_H

#include "STK_PoissonBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class Poisson_ljlk;

namespace Clust
{
/** @ingroup Clustering
 *  Traits class for the Poisson_ljlk traits policy. */
template<class _Array>
struct MixtureTraits< Poisson_ljlk<_Array> >
{
  typedef _Array                  Array;
  typedef typename Array::Type    Type;
  typedef Poisson_ljlk_Parameters  Parameters;
  typedef Array2D<Real>           Param;
};

} // namespace hidden

/** @ingroup Clustering
 *  The Poisson mixture model @c Poisson_ljlk is a Poisson model
 *  with a probability function of the form
 * \f[
 *    P(\mathbf{x}=(n_1,\ldots,n_d)|\theta)
 *     = \sum_{k=1}^K p_k \prod_{j=1}^d e^{-\lambda_{j}\lambda_{k}} \frac{(\lambda_{j}\lambda_{k})^{n_j}}{n_j!}.
 * \f]
 **/
template<class Array>
class Poisson_ljlk : public PoissonBase<Poisson_ljlk<Array> >
{
  public:
    typedef PoissonBase<Poisson_ljlk<Array> > Base;
    typedef typename Clust::MixtureTraits< Poisson_ljlk<Array> >::Parameters Parameters;

    using Base::p_tik;
    using Base::p_nk;
    using Base::components;
    using Base::p_data;
    using Base::param;


    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Poisson_ljlk( int nbCluster) : Base(nbCluster), stat_lambdaj_() {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Poisson_ljlk( Poisson_ljlk const& model)
                : Base(model), stat_lambdaj_(model.stat_lambdaj_) {}
    /** destructor */
    ~Poisson_ljlk() {}
    /** Initialize the component of the model.
     *  This function initialize the shared parameter sigma_  for all the
     *  components.
     **/
    void initializeModelImpl()
    {
      lambdaj_.resize(p_data()->cols());
      lambdaj_ = 1./this->nbVariable();
      for (int k= baseIdx; k < components().end(); ++k)
      { param(k).p_lambdaj_ = &lambdaj_;}
      stat_lambdaj_.initialize(p_data()->cols());
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResultsImpl(int iteration)
    { stat_lambdaj_.update(lambdaj_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResultsImpl()
    { stat_lambdaj_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParametersImpl()
    {
      lambdaj_ = stat_lambdaj_.param_;
      stat_lambdaj_.release();
    }
    /** Initialize randomly the parameters of the Poisson mixture. */
    void randomInit();
    /** Compute the weighted probabilities. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()+this->nbVariable();}

  protected:
    /** Vector of scaling */
    PointX lambdaj_;
    /** Common standard deviation */
    MixtureStatVector stat_lambdaj_;
};

/* Initialize randomly the parameters of the Poisson mixture. */
template<class Array>
void Poisson_ljlk<Array>::randomInit()
{
  for (int k = baseIdx; k < components().end(); ++k)
  {
    for (int j=p_data()->beginCols(); j< p_data()->endCols(); ++j)
    {
      Real m = (Real)p_data()->col(j).sum() / this->nbSample();
      param(k).lambdak_ = Law::Exponential::rand(m)/lambdaj_[j];
    }
  }
}


/* Compute the lambdas */
template<class Array>
bool Poisson_ljlk<Array>::mStep()
{
  lambdaj_  =  (Stat::sumByRow(*p_tik()).transpose() * (*p_data()))
             / (Stat::sumByRow(*p_tik()) * Stat::sumByRow(*p_data())).sum();
  PointX lk = Stat::sumByRow(*p_data()).transpose() * (*p_tik())/(*p_nk());

  for (int k = baseIdx; k < components().end(); ++k)
  { param(k).lambdak_ = lk[k];}
  return true;
}

} // namespace STK

#endif /* STK_POISSON_LJLK_H */
