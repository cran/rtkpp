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

/** @file STK_Poisson_ljk.h
 *  @brief In this file we implement the Poisson_ljk class
 **/

#ifndef STK_POISSON_LJK_H
#define STK_POISSON_LJK_H

#include "STK_PoissonBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class Poisson_ljk;

namespace Clust
{
/** @ingroup Clustering
 *  Traits class for the Poisson_ljk traits policy. */
template<class _Array>
struct MixtureTraits< Poisson_ljk<_Array> >
{
  typedef _Array                  Array;
  typedef typename Array::Type    Type;
  typedef Poisson_ljk_Parameters  Parameters;
  typedef Array2D<Real>           Param;
};

} // namespace hidden

/** @ingroup Clustering
 *  The Poisson mixture model @c Poisson_ljk is the most general Poisson model
 *  and have a probability function of the form
 * \f[
 *    P(\mathbf{x}=(n_1,\ldots,n_d)|\theta)
 *     = \sum_{k=1}^K p_k \prod_{j=1}^d e^{-\lambda_{jk}} \frac{\lambda_{jk}^{n_j}}{n_j!}.
 * \f]
 **/
template<class Array>
class Poisson_ljk : public PoissonBase<Poisson_ljk<Array> >
{
  public:
    typedef PoissonBase<Poisson_ljk<Array> > Base;
    typedef typename Clust::MixtureTraits< Poisson_ljk<Array> >::Parameters Parameters;

    using Base::p_tik;
    using Base::components;
    using Base::p_data;
    using Base::param;


    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Poisson_ljk( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Poisson_ljk( Poisson_ljk const& model) : Base(model) {}
    /** destructor */
    ~Poisson_ljk() {}
    /** Initialize randomly the parameters of the Poisson mixture. */
    void randomInit();
    /** Compute the weighted probabilities. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*this->nbVariable();}
};

/* Initialize randomly the parameters of the Poisson mixture. */
template<class Array>
void Poisson_ljk<Array>::randomInit()
{
  for (int k = baseIdx; k < components().end(); ++k)
  {
    for (int j=p_data()->beginCols(); j< p_data()->endCols(); ++j)
    {
      Real m = (Real)p_data()->col(j).sum() / this->nbSample();
      param(k).lambda_[j] = Law::Exponential::rand(m);
    }
  }
}


/* Compute the modalities probabilities */
template<class Array>
bool Poisson_ljk<Array>::mStep()
{
  for (int k = baseIdx; k < components().end(); ++k)
  {
    for (int j = p_data()->beginCols(); j < p_data()->endCols(); ++j)
    {
      param(k).lambda_[j] = p_data()->col(j).template cast<Real>().wmean(p_tik()->col(k));
    }
  }
  return true;
}

} // namespace STK

#endif /* STK_POISSON_LJK_H */
