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
 * Project:  stkpp::
 * created on: 5 avr. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Kernel_RationalQuadratic.h
 *  @brief In this file we define the class and methods for computing a RationalQuadratic Kernel.
 **/


#ifndef STK_KERNEL_RATIONALQUADRATIC_H
#define STK_KERNEL_RATIONALQUADRATIC_H

#include "STK_Kernel_IKernelBase.h"

namespace STK
{

namespace Kernel
{
/** @ingroup Kernel
 * The RationalQuadratic Kernel is a kernel of the form
 * \f[
 * k(x,y) = 1 - \left( \frac{\|x-y\|^2}{\|x-y\|^2+h} \right)
 * \f]
 * where @e h represents the bandwidth of the kernel.
 */
template<class Array>
class RationalQuadratic : public IKernelBase<Array>
{
  public:
    typedef IKernelBase<Array> Base;
    using Base::p_data_;
    using Base::gram_;
    using Base::symmetrize;
    /** constructor with a constant pointer on the data set
     *  @param p_data a pointer on a data set that will be "kernelized"
     *  @param shift the shift to use in the kernel
     **/
    RationalQuadratic( Array const* p_data, Real const& shift= 1.)
                     : Base(p_data), shift_(shift)
    { if (shift_ == 0.)
        STKDOMAIN_ERROR_1ARG(RationalQuadratic::RationalQuadratic,shift,shift must be!=0);
    }
    /** constructor with a constant pointer on the data set
     *  @param data a reference on a data set that will be "kernelized"
     *  @param shift the size of the windows to use in the kernel
     **/
    RationalQuadratic( Array const& data, Real const& shift= 1.)
                     : Base(data),shift_(shift)
    { if (shift_ == 0.)
        STKDOMAIN_ERROR_1ARG(RationalQuadratic::RationalQuadratic,shift,shift must be!=0);
    }
    /** destructor */
    virtual ~RationalQuadratic() {}
    /** @return the shift of the kernel */
    Real const& shift() const {return shift_;}
    /** set the shift of the kernel */
    void setWidth(Real const& shift) { shift_ = shift;}
    /** compute the kernel */
    virtual bool run();

  private:
    /** shift of the kernel */
    Real shift_;
};

template<class Array>
bool RationalQuadratic<Array>::run()
{
  typedef typename Array::Row RowVector;
  gram_.resize(p_data_->sizeRows());
  gram_.shift(p_data_->beginRows());
  for (int j= gram_.begin(); j < gram_.end(); ++j)
  {
    // create a reference on the current row
    RowVector row_j(p_data_->row(j), true);
    for (int i= gram_.begin(); i < j; ++i)
    { Real aux = (row_j - p_data_->row(i)).norm2();
       gram_(i,j) = 1 - aux/(aux + shift_);}
  }
  // symmetrize
  symmetrize();
  return true;
}
} // namespace Kernel

} // namespace STK

#endif /* STK_KERNEL_RATIONALQUADRATIC_H */
