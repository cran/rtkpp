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

/** @file STK_CategoricalBase.h
 *  @brief In this file we implement the base class for the Categorical
 *  diagonal models
 **/

#ifndef STK_CATEGORICALBASE_H
#define STK_CATEGORICALBASE_H

#include "../STK_IMixtureModel.h"
#include "../STK_MixtureParameters.h"
#include <STatistiK/include/STK_Law_Categorical.h>

namespace STK
{

/** @ingroup Clustering
 *  Base class for the categorical models Parameter Handler
 **/
template<class Derived>
struct CategoricalHandlerBase: public IRecursiveTemplate<Derived>
{
  /** @return the probability of the kth cluster, jth variable, lth modality */
  inline Real const& proba(int k, int j, int l) const { return this->asDerived().probaImpl(k,j,l);}
  /** @return the probability of the kth cluster for the jth variable */
  inline VectorX proba(int k, int j) const { return this->asDerived().probaImpl(k,j);}
};

/** @ingroup Clustering
 *  Base class for the Categorical models
 **/
template<class Derived>
class CategoricalBase : public IMixtureModel<Derived >
{
  protected:
    typedef IMixtureModel<Derived> Base;
    using Base::p_tik;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     *  @param nbCluster number of cluster in the model
     **/
    inline CategoricalBase( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline CategoricalBase( CategoricalBase const& model)
                          : Base(model), modalities_(model.modalities_) {}
    /** destructor */
    inline ~CategoricalBase() {}

  public:
    /** @return the array with the number of modalities of each columns in data set */
    inline PointXi const& nbModalities() const { return nbModalities_;}
    /** @return the range of the modalities */
    inline Range const& modalities() const { return modalities_;}
    /** @return the probability of the kth cluster, jth variable, lth modality */
    inline Real proba(int k, int j, int l) const { return param_.probaImpl(k,j,l);}
    /** @return the probability of the kth cluster for the jth variable */
    inline VectorX proba(int k, int j) const { return param_.probaImpl(k,j);}
    /** Initialize the model. Resize the probability arrays of each component.*/
    void initializeModelImpl()
    {
      // compute the maximal number of modalities
      nbModalities_.resize(p_data()->cols());
      int amin = Arithmetic<int>::max(), amax = Arithmetic<int>::min();
      for (int j= p_data()->beginCols(); j < p_data()->endCols(); ++j)
      {
        int min = p_data()->col(j).minElt(), max = p_data()->col(j).maxElt();
        amin = std::min(amin, min); amax = std::max(amax, max);
        nbModalities_[j] = max-min+1;
      }
      // set range of the modalities
      modalities_ = _R(amin, amax);
      // resize vectors of probabilities
      param_.resize(modalities_,p_data()->cols());
    }
    /** @return an imputation value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute */
    int impute(int i, int j) const;
    /** @return a simulated value for the jth variable of the ith sample
     * in the kth cluster
     * @param i,j,k indexes of the data to simulate */
    inline Real rand(int i, int j, int k) const
    { return Law::Categorical::rand(proba(k,j));}

  protected:
    /** Array with the number of modalities of each columns of the data set */
    PointXi nbModalities_;
    /** range of the modalities */
    Range modalities_;
};

/* Implementation  */
template<class Derived>
int CategoricalBase<Derived>::impute(int i, int j) const
{
  int lmax = modalities_.begin();
  Real pmax = -Arithmetic<Real>::max();
  // compute for each modality the pondered probability of occurrence
  for (int l=modalities_.begin(); l< modalities_.end(); ++l)
  {
    Real p = 0.;
    for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
    { p += p_tik()->elt(i,k) * proba(k, j, l);}

    if (pmax < p) { pmax = p; lmax = l;}
  }
  return lmax;
}


} // namespace STK

#endif /* STK_CATEGORICALBASE_H */
