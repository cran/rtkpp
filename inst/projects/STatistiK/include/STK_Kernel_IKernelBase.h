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
 * Project:  stkpp::Stat::Kernel
 * created on: 5 avr. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Kernel_IKernelBase.h
 *  @brief In this file we define the Interface base class for computing a Kernels.
 **/


#ifndef STK_KERNEL_IKERNELBASE_H
#define STK_KERNEL_IKERNELBASE_H

#include "Arrays/include/STK_CArraySquare.h"

namespace STK
{
namespace Kernel
{
/** @ingroup Kernel
 *  Interface Base class for the kernels classes.
 */
template<class Array>
class IKernelBase : public IRunnerBase
{
  protected:
    /** constructor with a constant pointer on the data set
     *  @param p_data a pointer on a data set that will be "kernelized"
     **/
    inline IKernelBase(Array const* p_data) : p_data_(p_data), gram_() {}
    /** constructor with a constant reference on the data set
     *  @param data a reference on a data set that will be "kernelized"
     **/
    inline IKernelBase(Array const& data) : p_data_(&data), gram_() {}

  public:
    /** destructor */
    inline virtual ~IKernelBase() {}
    /** @return the gram matrix*/
    inline CSquareX const& k() const { return gram_;}
    /** @return the gram matrix (bis) */
    inline CSquareX const& gram() const { return gram_;}
    /** set the current data set
     *  @param data the data set to "kernelized"
     **/
    inline void setData(Array const& data) { p_data_ = &data;}

  protected:
    /** A constant pointer on the data set */
    Array const* p_data_;
    /** the resulting gram_ matrix */
    CSquareX gram_;
    /** symmetrize the gram_ matrix using the upper part */
    void symmetrize()
    {
      // lower part
      for (int j= gram_.begin(); j < gram_.end(); ++j)
      {
        gram_(j,j) = 1.;
        for (int i= gram_.begin(); i < j; ++i)
        { gram_(j,i) = gram_(i,j);}
      }
    }
};

} // namespace Kernel

} // namespace STK

#endif /* STK_KERNEL_IKERNELBASE_H */
