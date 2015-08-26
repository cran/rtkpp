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
 * created on: 17 avr. 2010
 * Purpose:  Abstract class for the computation of the Index in the SLAAM.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_ILinearReduct.h In this file we define the interface base class
 *  ILinearReduct.
 **/

#ifndef STK_ILINEAREDUCT_H
#define STK_ILINEAREDUCT_H

#include "STK_IReducer.h"
#include <Arrays/include/STK_Array2DVector.h>

namespace STK
{
/** @ingroup Reduct
 *  @brief A ILinearReduct is an interface base class for reduction method
 *  using linear reduction. In order to find the linear reduction,
 *  the derived classes maximize an index.
 *
 * The class computes the optimal axis (stored in the @c axis_ )
 * attribute and the @c projected data set (stored in the @c p_reduct_ attribute
 * of the base class IReducer) when user uses the virtual method @c run()
 * (for not weighted observations) or @c run(weights) (for weighted observations).
 *
 * The Array axis_ is computed by maximizing some criteria defined in
 * derived classes. It is constructed using the pure virtual functions:
 * @code
 *  virtual void maximizeStep() =0;
 *  virtual void maximizeStep(weights) =0;
 * @endcode
 */
template<class Array, class Weights>
class ILinearReduct : public IReducer<Array, Weights>
{
  public:
    typedef IReducer<Array, Weights> Base;
    using Base::p_data_;
    using Base::p_reduced_;
    /** Default Constructor. */
    ILinearReduct();
    /** Constructor.
     *  @param p_data a pointer on the constant data set to reduce.
     **/
    ILinearReduct( Array const* p_data);
    /** Constructor.
     *  @param data a constant reference on the data set to reduce.
     **/
    ILinearReduct( Array const& data);
    /** copy Constructor.
     *  @param reducer the reducer to copy.
     **/
    ILinearReduct( ILinearReduct const& reducer);
    /** virtual destructor  */
    virtual ~ILinearReduct();
    /** @return a constant reference Vector of the value of the Index */
    inline VectorX const& criteriaValues() const { return idx_values_; }
    /** @return a constant reference Array of the axis */
    inline Array const& axis() const { return axis_; }
    /** Compute the projection matrix @c axis_ by maximizing the criteria and
     *  project the data set in order to obtain @c p_projected_.
     **/
    virtual bool run();
    /** Compute the projection matrix set by maximizing the weighted criteria
     *  and project the data set in order to obtain @c p_projected_.
     *  @param weights the weights to used in the index.
     */
    virtual bool run(Weights const& weights);

  protected:
    /** The values of the index for each axis. */
    VectorX idx_values_;
    /** The computed axis. */
    Array axis_;

  private:
    /** Find the axis by maximizing the Index. */
    virtual void maximizeStep() =0;
    /** Find the axis by maximizing the weighted Index.
     *  @param weights the weights to used
     **/
    virtual void maximizeStep( Weights const& weights) =0;
    /** Compute the projection of the data set on the Axis. */
    void projectionStep();
};

/*
 * Constructor.
 * @param data the input data set
 */
template<class Array, class Weights>
ILinearReduct<Array, Weights>::ILinearReduct(): Base() {}

/*
 * Constructor.
 * @param data the input data set
 */
template<class Array, class Weights>
ILinearReduct<Array, Weights>::ILinearReduct( Array const* p_data): Base(p_data) {}
/*
 * Constructor.
 * @param data the input data set
 */
template<class Array, class Weights>
ILinearReduct<Array, Weights>::ILinearReduct( Array const& data): Base(data) {}
/* copy Constructor.
 *  @param reducer the reducer to copy.
 **/
template<class Array, class Weights>
ILinearReduct<Array, Weights>::ILinearReduct( ILinearReduct const& reducer)
                                            : Base(reducer)
                                            , idx_values_(reducer.idx_values_)
                                            , axis_(reducer.axis_)
{}
/* Destructor */
template<class Array, class Weights>
ILinearReduct<Array, Weights>::~ILinearReduct() {}

/* Compute the Index.
 *  @param nbAxis number of Axis to compute
 */
template<class Array, class Weights>
bool ILinearReduct<Array, Weights>::run()
{
  if (!this->p_data_)
  { this->msg_error_ = STKERROR_NO_ARG(MultivariateArray::run(),data is not set);
    return false;
  }
  try
  {
    // maximize the Index and compute the axis
    maximizeStep();
    // project data
    projectionStep();

  } catch (Exception const& e)
  {
    this->msg_error_ = e.error();
    return false;
  }
  return true;
}

/*
 * Compute the weighted index.
 * @param weights the weights to used
 * @param nbAxis number of Axis to compute
 */
template<class Array, class Weights>
bool ILinearReduct<Array, Weights>::run( Weights const& weights)
{
  if (!this->p_data_)
  { this->msg_error_ = STKERROR_NO_ARG(MultivariateArray::run(),data is not set);
    return false;
  }
  try
  {
    // maximize the Index and compute the axis
    maximizeStep(weights);
    // project data
    projectionStep();

  } catch (Exception const& e)
  {
    this->msg_error_ = e.error();
    return false;
  }
  return true;
}

/* Compute the reduction of the data set on the Axis. */
template<class Array, class Weights>
void ILinearReduct<Array, Weights>::projectionStep()
{
  // check if p_reduced exists
  if (!p_reduced_) p_reduced_ = new Array;
  // compute matrix multiplication
  *p_reduced_ =   (*p_data_) * axis_;
}


} // namespace STK

#endif /* STK_ILINEARREDUCT_H */
