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

/** @file STK_GammaBase.h
 *  @brief In this file we implement the base class for the gamma models
 **/

#ifndef STK_GAMMABASE_H
#define STK_GAMMABASE_H

#include "../STK_IMixtureModel.h"
#include "../STK_MixtureParameters.h"
//#include "STK_GammaParameters.h"

#include <STatistiK/include/STK_Law_Categorical.h>
#include <STatistiK/include/STK_Law_Gamma.h>
#include <STatistiK/include/STK_Stat_Functors.h>

#include <Analysis/include/STK_Algo_FindZero.h>
#include <Analysis/include/STK_Funct_raw.h>
#include <Analysis/include/STK_Funct_gamma.h>

namespace STK
{

namespace hidden
{
/** @ingroup hidden
 *  Functor computing the derivative of the lnLikelihood of a gamma_ajk_bjk model */
class invPsiMLog : public IFunction<invPsiMLog >
{
  public:
    inline invPsiMLog( Real const& y): y_(y)  {}
    /** @return the value of the function at a
     * @param a a positive real value
     **/
    inline Real fImpl(Real const& a) const
    { return (y_ + std::log(a) - Funct::psi_raw(a));}
    /** @return the minimal value of the function at x */
    inline Real xminImpl() const { return 0;}

  private:
    Real y_;
};

/** @ingroup hidden
 *  Functor computing the difference between the psi function  and a fixed value
 **/
class invPsi : public IFunction<invPsi >
{
  public:
    /** initialize y_ */
    inline invPsi( Real const& y): y_(y) {}
    /** @return the value of the function at a
     *  @param x a positive real value
     **/
    inline Real fImpl(Real const& x) const { return (y_ - Funct::psi_raw(x));}
    /** @return the minimal value of the function at x */
    inline Real xminImpl() const { return 0;}
  private:
    Real y_;
};

} // namespace hidden

/** base class of the Gamma models Parameter Handler */
template<class Derived>
struct ParametersHandlerGammaBase: public IRecursiveTemplate<Derived>
{
  /** @return the mean of the kth cluster and jth variable */
  inline Real const& shape(int k, int j) const { return this->asDerived().shapeImpl(k,j);}
  /** @return the mean of the kth cluster and jth variable */
  inline Real const& scale(int k, int j) const { return this->asDerived().scaleImpl(k,j);}

  /** mean and statistics of the means */
  MixtureParametersSet<PointX> mean_;
  /** mean and statistics of the log-means */
  MixtureParametersSet<PointX> meanLog_;
  /** standard deviation and statistics */
  MixtureParametersSet<PointX> variance_;
  /** copy operator */
  inline ParametersHandlerGammaBase& operator=( Derived const& other)
  { mean_ = other.mean_; meanLog_ = other.meanLog_; variance_ = other.variance_; return *this; }

  /** default constructor */
  ParametersHandlerGammaBase( int nbCluster)
                            : mean_(nbCluster), meanLog_(nbCluster), variance_(nbCluster) {}
  /** copy constructor */
  ParametersHandlerGammaBase( ParametersHandlerGammaBase const& model)
                            : mean_(model.mean_), meanLog_(model.meanLog_), variance_(model.variance_) {}
  /** destructor */
  inline ~ParametersHandlerGammaBase() {}
  /** Initialize the parameters of the model.
   *  This function initialize the parameters and the statistics.
   **/
  void resize(Range const& range)
  {
    mean_.resize(range);
    mean_.initialize(1.);
    meanLog_.resize(range);
    meanLog_.initialize(0.);
    variance_.resize(range);
    variance_.initialize(1.);
  }
};

/** @ingroup Clustering
 *  Base class for the gamma models
 **/
template<class Derived>
class GammaBase : public IMixtureModel<Derived >
{
  public:
    typedef IMixtureModel<Derived > Base;
    using Base::p_tik; using Base::param_;
    using Base::p_nk;
    using Base::p_data;

  protected:
    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline GammaBase( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline GammaBase( GammaBase const& model): Base(model) {}
    /** destructor */
    inline ~GammaBase() {}

  public:
    /** @return the shape of the kth cluster and jth variable */
    inline Real shape(int k, int j) const { return param_.shape(k,j);}
    /** @return the scale of the kth cluster and jth variable */
    inline Real scale(int k, int j) const { return param_.scale(k,j);}
    /** Initialize the parameters of the model. */
    void initializeModelImpl() { param_.resize(p_data()->cols());}
    /** @return a value to impute for the jth variable of the ith sample*/
    Real impute(int i, int j) const
    {
      Real sum = 0.;
      for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
      { sum += p_tik()->elt(i,k) * shape(k,j) * scale(k,j);}
      return sum;
    }
    /** @return a simulated value for the jth variable of the ith sample
     *  in the kth cluster.
     *  @param i,j,k indexes of the data to simulate
     **/
    inline Real rand(int i, int j, int k) const
    { return Law::Gamma::rand(shape(k,j), scale(k,j));}

  protected:
    /** compute the Q(theta) value. */
    Real qValue() const;
    /** compute the weighted moments of a gamma mixture. */
    bool moments();
    /** get the weighted mean of the jth variable of the kth cluster. */
    inline Real meanjk( int j, int k) { return param_.mean_[k][j];}
    /** get the weighted variance of the jth variable of the kth cluster. */
    inline Real variancejk( int j, int k) { return param_.variance_[k][j];}
    /** get the mean of the weighted means of the kth cluster. */
    inline Real meank( int k) { return param_.mean_[k].mean();}
    /** get the mean of the weighted variances of the kth cluster. */
    inline Real variancek( int k) { return param_.variance_[k].mean();}
};

/* compute safely the weighted moments of a gamma law. */
template<class Derived>
bool GammaBase<Derived>::moments()
{
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    CVectorX tikColk(p_tik()->col(k), true); // create a reference
    for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
    {
      // mean
      Real mean =  p_data()->col(j).wmean(tikColk);
      if ( (mean<=0) || isNA(mean) ) { return false;}
      param_.mean_[k][j] = mean;
      // mean log
      Real meanLog =  p_data()->col(j).log().wmean(tikColk);
      if (isNA(meanLog)) { return false;}
      param_.meanLog_[k][j] = meanLog;
      // variance
      Real variance =  p_data()->col(j).wvariance(mean, tikColk);
      if ((variance<=0)||isNA(variance)){ return false;}
      param_.variance_[k][j] = variance;
    }
  }
  return true;
}

/*get the parameters of the model*/
template<class Derived>
Real GammaBase<Derived>::qValue() const
{
  Real value = 0.;
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    Real sumjk=0.0;
    for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
    {
      Real a = shape(k,j), b = scale(k,j);
      sumjk += a * (param_.meanLog_[k][j]-std::log(b))
             - param_.mean_[k][j]/b - STK::Funct::gammaLn(a);
    }
    value += p_nk()->elt(k) * sumjk;
  }
  return value;
}

} // namespace STK

#endif /* STK_GAMMABASE_H */
