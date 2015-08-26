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
 * created on: 21 mai 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_MixtureParameters.h
 *  @brief In this file we define the structures enclosing the different kind of
 *  parameters.
 **/

#ifndef STK_MIXTUREPARAMETERS_H
#define STK_MIXTUREPARAMETERS_H

#include <STatistiK/include/STK_Stat_Online.h>
#include <Arrays/include/STK_Array1D.h>

namespace STK
{
namespace hidden
{
  template<bool, class> struct TypeDispatch;
  template<class TypeParam> struct TypeDispatch<true, TypeParam>
  {
    typedef TypeParam Result;
  };
  template<class TypeParam> struct TypeDispatch<false, TypeParam>
  {
    typedef typename TypeParam::Type Result;
  };
} // namespace hidden
/** @ingroup Clustering
 *  @brief utility class used in order to manage parameters stored in an
 *  arbitrary structure.
 *
 *  This class can also handle single
 */
template<class TypeParam>
struct MixtureParameters
{
  enum
  {
    isAritmetic_ = hidden::isArithmetic<TypeParam>::yes
  };
  typedef typename hidden::TypeDispatch< isAritmetic_ , TypeParam>::Result Type;
  /** Array of the probabilities */
  TypeParam param_;
  /** Array of the probabilities statistics */
  Stat::Online<TypeParam, Type> stat_param_;

  /** default constructor */
  inline MixtureParameters(): param_(), stat_param_() {}
  /** copy constructor.
   * @param param the parameters to copy.
   **/
  inline MixtureParameters( MixtureParameters const& param)
                          : param_(param.param_)
                          , stat_param_(param.stat_param_)
  {}
  /** copy operator */
  inline MixtureParameters& operator=( MixtureParameters const& other)
  { param_ = other.param_; return *this; }
  /** get all parameters */
  inline TypeParam const&  operator()() const { return param_;}
  /** get parameters */
  inline TypeParam&  operator()() { return param_;}
  /** resize 1D containers */
  inline void resize(Range const& range)
  { param_.resize(range); stat_param_.resize(range); }
  /** resize 2D containers */
  inline void resize(Range const& rows, Range const& cols)
  { param_.resize(rows, cols); stat_param_.resize(rows, cols); }
  /** initialize the parameters with a value and rlease the statistics */
  inline void initialize(Type const& value)
  { param_ = value; stat_param_.release();}
  /** destructor */
  inline ~MixtureParameters() {}
  /** Store the intermediate results of the Mixture.
   *  @param iteration Provides the iteration number beginning after the burn-in period.
   **/
  inline void storeIntermediateResults(int iteration)
  { stat_param_.update(param_);}
  /** Release the stored results. This is usually used if the estimation
   *  process failed.
   **/
  inline void releaseIntermediateResults()
  { stat_param_.release();}
  /** set the parameters stored in stat_param and release stat_param. */
  inline void setParameters()
  { param_ = stat_param_.mean_; stat_param_.release();}
};

/** @ingroup Clustering
 *  @brief utility class used in order to manage a set of parameters stored in
 *  an arbitrary structure.
 *  @tparam TypeParam the type of the structure storing the parameters
 */
template<class TypeParam>
struct MixtureParametersSet
{
  enum
  {
    isAritmetic_ = hidden::isArithmetic<TypeParam>::yes
  };
  typedef typename hidden::TypeDispatch< isAritmetic_ , TypeParam>::Result Type;
  /** Array of the probabilities */
  Array1D<TypeParam> param_;
  /** Array of the probabilities statistics */
  Array1D<Stat::Online<TypeParam, Type> > stat_param_;
  /** default constructor */
  inline MixtureParametersSet( int nbCluster)
                             : param_(nbCluster), stat_param_(nbCluster)
  {}
  /** copy constructor.
   * @param param the parameters to copy.
   **/
  inline MixtureParametersSet( MixtureParametersSet const& param)
                             : param_(param.param_)
                             , stat_param_(param.stat_param_)
  {}
  /** destructor */
  inline ~MixtureParametersSet() {}
  /** get all parameters */
  inline Array1D<TypeParam> const&  operator()() const { return param_;}
  /** get all parameters */
  inline Array1D<TypeParam>&  operator()() { return param_;}
  /** get parameters of the kth component */
  inline TypeParam const&  operator[](int k) const { return param_[k];}
  /** get parameters of the kth component */
  inline TypeParam&  operator[](int k) { return param_[k];}
  /** resize 1D containers */
  inline void resize(Range const& range)
  {
    STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(TypeParam);
    for (int k= stat_param_.begin(); k < stat_param_.end(); ++k)
    {
      param_[k].resize(range);
      stat_param_[k].resize(range);
    }
  }
  /** resize 2D containers */
  inline void resize(Range const& rows, Range const& cols)
  {
    for (int k= param_.begin(); k < param_.end(); ++k)
    {
      param_[k].resize(rows, cols);
      stat_param_[k].resize(rows, cols);
    }
  }
  /** initialize parameters with a value and release statistics */
  inline void initialize(Type const& value)
  {
    for (int k= param_.begin(); k < param_.end(); ++k)
    { param_[k] = value; stat_param_[k].release();}
  }
  /** Store the intermediate results of the Mixture.
   *  @param iteration Provides the iteration number beginning after the burn-in period.
   **/
  inline void storeIntermediateResults(int iteration)
  {
    for (int k= stat_param_.begin(); k < stat_param_.end(); ++k)
    { stat_param_[k].update(param_[k]);}
  }
  /** Release the stored statistics. This is usually used if the estimation
   *  process failed.
   **/
  inline void releaseIntermediateResults()
  {
    for (int k= stat_param_.begin(); k < stat_param_.end(); ++k)
    { stat_param_[k].release();}
  }
  /** set the parameters stored in stat_param and release stat_param. */
  inline void setParameters()
  {
    for (int k= stat_param_.begin(); k < stat_param_.end(); ++k)
    { param_[k] = stat_param_[k].mean_; stat_param_[k].release();}
  }
};

} // namespace STK

#endif /* STK_MIXTUREPARAMSTAT_H*/
