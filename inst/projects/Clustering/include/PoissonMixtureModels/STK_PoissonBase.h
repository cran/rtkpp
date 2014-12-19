/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014 Serge IOVLEFF

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
#include "STK_PoissonParameters.h"

#include "STatistiK/include/STK_Law_Poisson.h"

namespace STK
{

/** @ingroup Clustering
 *  Base class for the Poisson models
 **/
template<class Derived>
class PoissonBase : public IMixtureModel<Derived >
{
  public:
    typedef IMixtureModel<Derived > Base;

    using Base::p_tik;
    using Base::components;
    using Base::p_data;
    using Base::p_param;


  protected:
    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline PoissonBase( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline PoissonBase( PoissonBase const& model) : Base(model) {}
    /** destructor */
    inline ~PoissonBase() {}

  public:
    /** @return a value to impute for the jth variable of the ith sample*/
    Real impute(int i, int j) const
    {
      Real sum = 0.;
      for (int k= p_tik()->beginCols(); k < components().end(); ++k)
      { sum += p_tik()->elt(i,k) * p_param(k)->lambda(j);}
      return sum;
    }
    /** @return a simulated value for the jth variable of the ith sample
     *  @param i,j indexes of the value to sample
     **/
    Real sample(int i, int j) const
    {
      int k = Law::Categorical::rand(p_tik()->row(i));
      return Law::Poisson::rand(p_param(k)->lambda(j));
    }
    /** get the parameters of the model
     *  @param params the array to fill with the parameters of the model
     **/
    void getParameters(Array2D<Real>& params) const;
    /** @return the parameters of the model in an array of size (K * d). */
    ArrayXX getParametersImpl() const;
    /** Write the parameters on the output stream os */
    void writeParameters(ostream& os) const;
};

/* Write the parameters on the output stream os */
template<class Derived>
void PoissonBase<Derived>::writeParameters(ostream& os) const
{
  PointX lambda(p_data()->cols());
  for (int k= baseIdx; k < components().end(); ++k)
  {
    // store shape and scale values in an array for a nice output
    for (int j= p_data()->beginCols();  j < p_data()->endCols(); ++j)
    { lambda[j] = p_param(k)->lambda(j);}
    os << _T("---> Component ") << k << _T("\n");
    os << _T("lambda = ") << lambda;
  }
}

/*get the parameters of the model*/
template<class Derived>
void PoissonBase<Derived>::getParameters(Array2D<Real>& params) const
{
  params.resize(this->nbCluster(), p_data()->cols());
  for (int k= params.beginRows(); k < params.endRows(); ++k)
  {
    for (int j= p_data()->beginCols();  j < p_data()->endCols(); ++j)
    { params(k, j) = p_param(k)->lambda(j);}
  }
}
/* get the parameters of the model in an array of size (K * d). */
template<class Derived>
ArrayXX PoissonBase<Derived>::getParametersImpl() const
{
  ArrayXX params;
  params.resize(this->nbCluster(), p_data()->cols());
  for (int k= params.beginRows(); k < params.endRows(); ++k)
  {
    for (int j= p_data()->beginCols();  j < p_data()->endCols(); ++j)
    { params(k, j) = p_param(k)->lambda(j);}
  }
  return params;
}


} // namespace STK

#endif /* STK_POISSONBASE_H */
