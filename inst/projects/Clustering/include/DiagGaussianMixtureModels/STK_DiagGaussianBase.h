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
 * created on: Dec 4, 2013
 * Authors: Serge Iovleff
 **/

/** @file STK_DiagGaussianBase.h
 *  @brief In this file we implement the base class for the Gaussian diagonal models
 **/

#ifndef STK_DIAGGAUSSIANBASE_H
#define STK_DIAGGAUSSIANBASE_H

#include "../STK_IMixtureModel.h"
#include "../STK_MixtureParameters.h"
#include <STatistiK/include/STK_Law_Normal.h>
#include <STatistiK/include/STK_Law_Uniform.h>

namespace STK
{

/** @ingroup Clustering
 *  Base class for the Diagonal Gaussian model Parameter Handler
 **/
template<class Derived>
struct DiagGaussianHandlerBase: public IRecursiveTemplate<Derived>
{
  /** @return the value of the mean of the kth cluster and jth variable */
  inline Real const& mean(int k, int j) const { return this->asDerived().meanImpl(k,j);}
  /** @return the value of the mean of the kth cluster and jth variable */
  inline Real const& sigma(int k, int j) const { return this->asDerived().sigmaImpl(k,j);}
};

/** @ingroup Clustering
 *  Base class for the diagonal Gaussian models
 **/
template<class Derived>
class DiagGaussianBase : public IMixtureModel<Derived >
{
  public:
    typedef IMixtureModel<Derived > Base;
    typedef typename Clust::MixtureTraits<Derived>::ParamHandler ParamHandler;
    using Base::p_tik; using Base::param_;
    using Base::p_data;

  protected:
    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline DiagGaussianBase( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline DiagGaussianBase( DiagGaussianBase const& model) : Base(model) {}
    /** destructor */
    inline ~DiagGaussianBase() {}

  public:
    /** @return mean of the kth cluster and jth variable */
    inline Real const& mean(int k, int j) const { return param_.meanImpl(k,j);}
    /** @return sigma of the kth cluster and jth variable */
    inline Real const& sigma(int k, int j) const { return param_.sigmaImpl(k,j);}
    /** Initialize the parameters of the model. */
    inline void initializeModelImpl() { param_.resize(p_data()->cols());}
    /** @return an imputation value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute */
    Real impute(int i, int j) const
    {
      Real sum = 0.;
      for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
      { sum += p_tik()->elt(i,k) * mean(k,j);}
      return sum;
    }
    /** @return a simulated value for the jth variable of the ith sample
     * in the kth cluster
     * @param i,j,k indexes of the data to simulate */
    inline Real rand(int i, int j, int k) const
    { return Law::Normal::rand(mean(k, j), sigma(k,j));}

  protected:
    PointX& mean(int k) { return param_.mean_[k];}
    /** sample randomly the mean of each component by sampling randomly a row
     *  of the data set.
     **/
    void randomMean();
    /** compute the weighted mean of a Gaussian mixture. */
    bool updateMean();
};

template<class Derived>
void DiagGaussianBase<Derived>::randomMean()
{
  // indexes array
  VectorXi indexes(p_data()->rows());
  for(int i=p_data()->beginRows(); i< p_data()->endRows(); ++i) { indexes[i] = i;}
  Range rind = p_data()->rows();
  // sample without repetition
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    // random number in [0, end-k[
    int i = (int)Law::Uniform::rand(rind.begin(), rind.end());
    // get ith individuals
    mean(k).copy(p_data()->row(indexes[i]));
    // exchange it with nth
    indexes.swap(i, rind.lastIdx());
    // decrease
    rind.decLast(1);
  }
}

template<class Derived>
bool DiagGaussianBase<Derived>::updateMean()
{
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    for (int j=p_data()->beginCols(); j< p_data()->endCols(); ++j)
    { mean(k)[j] = p_data()->col(j).wmean(p_tik()->col(k));}
  }
  return true;
}

} // namespace STK

#endif /* STK_DIAGGAUSSIANBASE_H */
