/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015  Serge Iovleff

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

 Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 */

/*
 * Project:  stkpp::StatModel
 * created on: 13 août 2011
 * Purpose:  Create a gaussian statistical model.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_GaussianModel.h
 *  @brief In this file we define the GassianModel class.
 **/

#ifndef STK_GAUSSIANMODEL_H
#define STK_GAUSSIANMODEL_H

#include "STK_IGaussianModel.h"
#include <STatistiK/include/STK_Stat_MultivariateReal.h>
#include <STatistiK/include/STK_MultiLaw_Normal.h>

namespace STK
{

/** @ingroup StatModels
 *  @brief Compute the maximum likelihood estimates of a complete Gaussian
 *  statistical model.
 **/
template <class Array>
class GaussianModel : public IGaussianModel<Array>
{
  public:
    typedef IGaussianModel<Array> Base;
    using Base::p_data_;
    using Base::nbVariable;
    using Base::mean_;
    using Base::p_law_;

    typedef typename Array::Col ColVector;
    typedef typename Array::Row RowVector;
    /** constructor.
     * @param p_data pointer on the data set
     */
    GaussianModel( Array const* p_data);
    /** constructor.
     * @param data reference on the data set
     */
    GaussianModel( Array const& data);
    /** destructor. */
    ~GaussianModel();
    /** implementation of the Gaussian statistical model
     * @return @c true if no error occur and @c false otherwise.
     */
    bool run();
    /** implementation of the weighted Gaussian statistical model
     * @param weights the weights of the samples
     * @return @c true if no error occur and @c false otherwise.
     */
    bool run(ColVector const& weights);
    /** get the empirical covariance
     * @return the empirical covariance
     */
    inline ArraySquareX const& covariance() const { return cov_;}
  protected:
    /** ArrayXX of the empirical covaiance */
    ArraySquareX cov_;
    /** compute the empirical covariance matrix. */
    void compCovariance();
    /** compute the empirical weighted covariance matrix.
     * @param weights the weights of the samples
     **/
    void compWeightedCovariance(ColVector const& weights);
};

/* constructor */
template <class Array>
GaussianModel<Array>::GaussianModel( Array const* p_data)
                            : Base(p_data)
                            , cov_(p_data_->cols())
{
  this->setNbFreeParameter(nbVariable() + (nbVariable()* (nbVariable()-1))/2);
}

/* constructor */
template <class Array>
GaussianModel<Array>::GaussianModel( Array const& data)
                            : Base(data)
                            , cov_(data.cols())
{
  setNbFreeParameter(nbVariable() + (nbVariable()* (nbVariable()-1))/2);
}

/* destructor */
template <class Array>
GaussianModel<Array>::~GaussianModel()
{ if (p_law_) delete p_law_;}

/* implementation of the Gaussian statistical model
 * @return @c true if no error occur and @c false otherwise.
 */
template <class Array>
bool GaussianModel<Array>::run()
{
  // compute the mean
  this->compMean();
  // compute the covariance matrix
  compCovariance();
  // create p_law_ (will be deleted in base class)
  // update gaussian law (will be deleted in base class)
  if (!p_law_) p_law_ = new MultiLaw::Normal<RowVector>(mean_, cov_);
  else static_cast<MultiLaw::Normal<RowVector>*>(p_law_)->setParameters(mean_, cov_);
  // compute log likelihood of the gaussian law
  this->setLnLikelihood(static_cast<MultiLaw::Normal<RowVector>* >(p_law_)->lnLikelihood(*p_data_ ));
  // everything ok
  return true;
}

/* implementation of the weighted Gaussian statistical model
 * @param weights the weights of the samples
 * @return @c true if no error occur and @c false otherwise.
 */
template <class Array>
bool GaussianModel<Array>::run(ColVector const& weights)
{
  // compute the mean
  this->compWeightedMean(weights);
  // compute the covariance matrix
  compWeightedCovariance(weights);
  // create p_law_ (will be deleted in base class)
  // update gaussian law (will be deleted in base class)
  if (!p_law_) p_law_ = new MultiLaw::Normal<RowVector>(mean_, cov_);
  else static_cast<MultiLaw::Normal<RowVector>*>(p_law_)->setParameters(mean_, cov_);
  // compute log likelihood of the gaussian law
  this->setLnLikelihood(static_cast<MultiLaw::Normal<RowVector>* >(p_law_)->lnLikelihood(*p_data_ ));
  // everything ok
  return true;
}

/** compute the empirical covariance matrix. */
template <class Array>
void GaussianModel<Array>::compCovariance()
{ Stat::covariance(*p_data_,cov_);}
/** compute the empirical weighted covariance matrix.
 * @param weights the weights of the samples
 **/
template <class Array>
void GaussianModel<Array>::compWeightedCovariance(ColVector const& weights)
{ Stat::covariance(*p_data_, weights, cov_);}

} // namespace STK

#endif /* STK_GAUSSIANMODEL_H */
