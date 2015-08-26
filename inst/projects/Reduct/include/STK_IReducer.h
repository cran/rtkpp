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
 * Project:  stkpp::Reduct
 * created on: 17 avr. 2010
 * Purpose:  Abstract class for the computation of the Index in the SLAAM.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_IReducer.h In this file we define the interface base
 *  class IReducer.
 **/

#ifndef STK_IREDUCER_H
#define STK_IREDUCER_H

#include "Sdk/include/STK_IRunner.h"

namespace STK
{

/** @ingroup Reduct
 *  @brief Interface base class for reducing methods.
 *  A reducer is a class implementing a dimension reduction techniques.
 *
 * The class receives a matrix in input of size (n,p).
 * - n is the number of samples,
 * - p is the number of variables.
 * The observations can be weighted.
 *
 * The derived class will compute a @em reduced data set of dimension (n,d).
 * At this level there is no assumption about the type of the data that will
 * be reduced. However, it is expected they can be part of a linear space.
 */
template<class Array, class Weights>
class IReducer : public IRunnerUnsupervised<Array, Weights>
{
  protected:
    typedef IRunnerUnsupervised<Array, Weights> Runner;
    /** Default constructor. */
    IReducer();
    /** Constructor with a pointer on the constant data set.
     *  @param p_data a pointer on the data set to reduce.
     * */
    IReducer( Array const* p_data);
    /** Constructor with a constant reference on the data set.
     *  @param data the data set to reduce.
     * */
    IReducer( Array const& data);
    /** Copy constructor. Use the operator = of the clas ArrayXX
     *  @param reducer The reducer to copy.
     * */
    IReducer( IReducer const& reducer);

  public:
    /** virtual destructor. Delete reduced data set if allocated. */
    virtual ~IReducer();
    /** get the number of dimension.
     *  @return The number of dimension computed
     **/
    inline int dim() const { return dim_;}
    /** get a pointer on the reduced data set
     *  @return a constant pointer on the data set
     **/
    inline Array* p_reduced() const { return p_reduced_; }
    /** set the number of dimension.
     *  @param dim the number of dimension to set
     **/
    inline void setDimension( const int& dim) { dim_ = dim;}
    /** clear allocated memory */
    inline void clear()
    { if (p_reduced_) { delete p_reduced_; p_reduced_ = 0;} }

  protected:
    /** dimension of the reduced data set */
    int dim_;
    /** The reduced data set. */
    Array* p_reduced_;
};

/* Default constructor. */
template<class Array, class Weights>
IReducer<Array, Weights>::IReducer(): Runner(), dim_(0), p_reduced_(0) {}

/* Constructor with a pointer on the constant data set.
 *  @param p_data the data set to reduce.
 **/
template<class Array, class Weights>
IReducer<Array, Weights>::IReducer( Array const* p_data)
                                  : Runner(p_data), dim_(0), p_reduced_(0)
{}
/* Constructor with a constant reference on the data set.
 *  @param data the data set to reduce.
 * */
template<class Array, class Weights>
IReducer<Array, Weights>::IReducer( Array const& data)
                                  : Runner(data), dim_(0), p_reduced_(0)
{}
/* Copy constructor.
 *  @param data The data set to reduce.
 * */
template<class Array, class Weights>
IReducer<Array, Weights>::IReducer( IReducer const& reducer)
                                  : Runner(reducer), dim_(reducer.dim_), p_reduced_(0)
{
  if (reducer.p_reduced_)
  { p_reduced_ = reducer.p_reduced_->clone();}
}

/*
 * Destructor
 */
template<class Array, class Weights>
IReducer<Array, Weights>::~IReducer()
{ if (p_reduced_) delete p_reduced_;}


} // namespace STK

#endif /* STK_IREDUCER_H */
