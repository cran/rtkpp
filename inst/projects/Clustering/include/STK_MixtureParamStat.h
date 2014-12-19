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
 * created on: 11 dec. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_MixtureParamStat.h
 *  @brief In this file we define the structures that will be used in order to
 *  compute the statistics on the parameters.
 **/

#ifndef STK_MIXTUREPARAMSTAT_H
#define STK_MIXTUREPARAMSTAT_H

#include "Arrays/include/STK_Array2DVector.h"

namespace STK
{
/** @ingroup Clustering
 *  @brief utility class used in order to store statistics about the
 *  values of two vectors of parameters during the estimation process.
 */
struct MixtureStatReal
{
  MixtureStatReal(): param_(), iter_(0) {}
  MixtureStatReal( MixtureStatReal const& stat)
                 : param_(stat.param_), iter_(stat.iter_)
  {}
  /** initialize the structure */
  inline void initialize() { param_ = 0.; iter_ =0;}
  /** release the computed parameters */
  inline void release() { param_ = 0.; iter_ = 0;}
  /** update the parameters using the current estimated parameters
   *  @param param the current value of the parameter
   **/
  inline void update(Real const& param)
  { iter_++; param_ += (param - param_)/iter_;}
  /** First Vector of the parameters */
  Real param_;
  /** number of stored values */
  int iter_;
};

/** @ingroup Clustering
 *  @brief utility class used in order to store statistics about the
 *  values of a vector of parameters during the estimation process.
 */
struct MixtureStatVector
{
  MixtureStatVector(): param_(), iter_(0) {}
  MixtureStatVector( MixtureStatVector const& stat)
                   : param_(stat.param_), iter_(stat.iter_)
  {}
  /** initialize the structure */
  inline void initialize(Range const& range)
  { param_.resize(range); param_ = 0; iter_ =0;}
  /** release the computed parameters */
  void release() { param_ = 0; iter_ = 0;}
  /** update the parameters using the current estimate param
   * @param param the current value of the parameters
   **/
  template<class OtherArray>
  void update(OtherArray const& param) { iter_++; param_ += (param - param_)/iter_;}
  /** Vector of the parameters */
  VectorX param_;
  /** number of stored values */
  int iter_;
};


} // namespace STK

#endif /* STK_MIXTUREPARAMSTAT_H*/
