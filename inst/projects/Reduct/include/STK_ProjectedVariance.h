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
 * Project:  stkpp::AAModels
 * Purpose:  Implementation of the ILinearReduct interface using the total variance.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/


#ifndef STK_PROJECTEDVARIANCE_H
#define STK_PROJECTEDVARIANCE_H

#include <STatistiK/include/STK_Stat_MultivariateReal.h>
#include <Algebra/include/STK_SymEigen.h>
#include "STK_ILinearReduct.h"

namespace STK
{

/** @ingroup Reduct
 *  @brief A ProjectedVariance is an implementation of the abstract
 *  @c ILinearReduct interface.
 *
 *  ProjectedVariance (PCA) reduce the dimension of data by maximizing the
 *  projected varaince on an affine subspace of dimension d.
**/
template<class Array>
class ProjectedVariance : public ILinearReduct<Array, Vector>
{
  public:
    typedef ILinearReduct<Array, Vector> Base;
    using Base::p_data_;
    using Base::p_reduced_;
    using Base::axis_;
    using Base::idx_values_;
    /** default constructor */
    ProjectedVariance();
    /** Constructor.
     *  @param p_data a pointer on the constant data set to reduce.
     **/
    ProjectedVariance(Array const* p_data);
    /** Constructor.
     *  @param data a constant reference on the data set to reduce.
     **/
    ProjectedVariance(Array const& data);
    /** Copy constructor.
     * @param reducer the reducer to copy
     **/
    ProjectedVariance(ProjectedVariance const& reducer);
    /** Destructor */
    virtual ~ProjectedVariance();
    /** clone pattern
     *  @return a pointer on the clone of this
     **/
    inline virtual ProjectedVariance* clone() const
    { return new ProjectedVariance(*this);}

  protected:
    /** the covariance Array */
    ArraySquareX covariance_;

  private:
    /** Find the axis by maximizing the Index. */
    virtual void maximizeStep();

    /** Find the axis by maximizing the weighed Index.
     * @param weights the weights of the samples.
     **/
    virtual void maximizeStep(Vector const& weights);
    /** compute axis and index. */
    void computeAxis();
    /**  update the class if a new data set is set. */
    virtual void update();
};

/* default constructor */
template<class Array>
ProjectedVariance<Array>::ProjectedVariance(): Base() {}
/* Constructor.
 *  @param p_data a pointer on the constant data set to reduce.
 **/
template<class Array>
ProjectedVariance<Array>::ProjectedVariance( Array const* p_data)
                                           : Base(p_data)
{}
/* Constructor.
 *  @param data a constatn reference on the data set to reduce.
 **/
template<class Array>
ProjectedVariance<Array>::ProjectedVariance( Array const& data)
                                           : Base(data)
{}
/* Copy constructor.
 * @param reducer the reducer to copy
 **/
template<class Array>
ProjectedVariance<Array>::ProjectedVariance( ProjectedVariance const& reducer)
                                           : Base(reducer)
{}

 /* Destructor */
template<class Array>
ProjectedVariance<Array>::~ProjectedVariance()
{}

/* Find the axis by maximizing the Index. */
template<class Array>
void ProjectedVariance<Array>::maximizeStep()
{
#ifdef STK_REDUCT_DEBUG
  if (!p_data_)
  { STKRUNTIME_ERROR_NO_ARG(ProjectedVariance::maximizeStep,data is not set);}
#endif
  Stat::covariance(*p_data_, covariance_, true);
  computeAxis();
}

/* Find the axis by maximizing the weighed Index.
 * @param weights the weights of the samples.
 **/
template<class Array>
void ProjectedVariance<Array>::maximizeStep(Vector const& weights)
{
#ifdef STK_REDUCT_DEBUG
  if (!p_data_)
  { STKRUNTIME_ERROR_NO_ARG(ProjectedVariance::maximizeStep,data is not set);}
#endif
  Stat::covariance(*p_data_, weights, covariance_, true);
  computeAxis();
}

/* compute axis and index. */
template<class Array>
void ProjectedVariance<Array>::computeAxis()
{
  SymEigen<ArraySquareX> eigen(covariance_);
  eigen.run();

  // compute the range of the axis
  Range range(p_data_->beginCols(), std::min(this->dim_, p_data_->sizeCols()));
  // copy axis and index values
  axis_.resize(p_data_->cols(), range);
  idx_values_.resize(range);
  axis_       = eigen.rotation().col(range);
  idx_values_ = eigen.eigenValues().sub(range);
}

/* update the class if a new data set is set */
template<class Array>
void ProjectedVariance<Array>::update()
{
  covariance_.clear();
}


} // namespace STK

#endif /* STK_PROJECTEDVARIANCE_H_ */
