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

/* Project:  stkpp::Clustering
 * created on: Dec 9, 2014
 * Authors: Serge Iovleff
 **/

/** @file STK_PoissonBase.h
 *  @brief In this file we implement the base class for the exponential models
 **/

#ifndef STK_POISSONBASE_H
#define STK_POISSONBASE_H

#include "../STK_IMixtureModel.h"
#include "../STK_MixtureParameters.h"
#include <STatistiK/include/STK_Law_Poisson.h>

namespace STK
{
/** @ingroup Clustering
 *  Base class for the Poisson models Parameter Handler
 **/
template<class Derived>
struct PoissonHandlerBase: public IRecursiveTemplate<Derived>
{
  /** @return the value of lambda of the kth cluster and jth variable */
  inline Real lambda(int k, int j) const { return this->asDerived().lambdaImpl(k,j);}
};

/** @ingroup Clustering
 *  Base class for the Poisson models
 **/
template<class Derived>
class PoissonBase: public IMixtureModel<Derived >
{
  public:
    typedef IMixtureModel<Derived > Base;
    using Base::p_tik; using Base::param_;
    using Base::p_data;
    using Base::nbCluster;

  protected:
    /** default constructor
     *  @param nbCluster number of cluster in the model
     **/
    inline PoissonBase( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline PoissonBase( PoissonBase const& model) : Base(model) {}
    /** destructor */
    inline ~PoissonBase() {}

  public:
    /** @return the value of lambda of the kth cluster and jth variable */
    inline Real lambda(int k, int j) const { return param_.lambda(k,j);}
    /** Initialize the parameters of the model. */
    void initializeModelImpl() { param_.resize(p_data()->cols());}
    /** @return a value to impute for the jth variable of the ith sample*/
    Real impute(int i, int j) const
    {
      Real sum = 0.;
      for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
      { sum += p_tik()->elt(i,k) * lambda(k,j);}
      return sum;
    }
    /** @return a simulated value for the jth variable of the ith sample
     *  in the kth cluster.
     *  @param i,j,k indexes of the data to simulate
     **/
    inline Real rand(int i, int j, int k) const
    { return Law::Poisson::rand(lambda(k,j));}
};

} // namespace STK

#endif /* STK_POISSONBASE_H */
