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
  typedef Poisson_lk_Parameters  Parameters;
  typedef Array2D<Real>           Param;
};

} // namespace hidden

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
    typedef typename Clust::MixtureTraits< Poisson_lk<Array> >::Parameters Parameters;

    using Base::p_tik;
    using Base::p_nk;
    using Base::components;
    using Base::p_data;
    using Base::p_param;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Poisson_lk( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Poisson_lk( Poisson_lk const& model) : Base(model) {}
    /** destructor */
    ~Poisson_lk() {}
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
  Real m = p_data()->mean();
  for (int k = baseIdx; k < components().end(); ++k)
  {
    p_param(k)->lambda_ = Law::Exponential::rand(m);
  }
}


/* Compute the modalities probabilities */
template<class Array>
bool Poisson_lk<Array>::mStep()
{
  for (int k = baseIdx; k < components().end(); ++k)
  { p_param(k)->lambda_ = (p_data()->transpose() * p_tik()->col(k)).sum()
                        / (p_data()->sizeCols()*p_nk()->elt(k));
  }
  return true;
}

} // namespace STK

#endif /* STK_POISSON_LK_H */
