/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2012  Serge Iovleff

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
 * Project:  stkpp::Arrays
 * created on: 26 nov. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ExprBaseDot.h
 *  @brief In this file we define the dot product and its particular cases.
 **/


#ifndef STK_EXPRBASEDOT_H
#define STK_EXPRBASEDOT_H

#include "Sdk/include/STK_StaticAssert.h"

namespace STK
{

/*  @returns the dot product of *this with other. */
template<class Derived>
template<class Rhs>
typename STK::hidden::Promote<typename STK::hidden::Traits<Derived>::Type, typename Rhs::Type>::result_type const
inline ExprBase<Derived>::dot(ExprBase<Rhs> const& other) const
{
  STK_STATIC_ASSERT_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  return BinaryOperator< ProductOp<Type, typename hidden::Traits<Rhs>::Type>
                       , Derived
                       , Rhs>(this->asDerived(), other.asDerived()).sum();
}

/*  @returns the dot product of *this with other. */
template<class Derived>
template<class Rhs>
typename STK::hidden::Promote<typename STK::hidden::Traits<Derived>::Type, typename Rhs::Type>::result_type const

inline ExprBase<Derived>::dotSafe(ExprBase<Rhs> const& other) const
{
  STK_STATIC_ASSERT_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  return BinaryOperator< ProductOp<Type, typename hidden::Traits<Rhs>::Type>
                       , Derived
                       , Rhs>(this->asDerived(), other.asDerived()).sumSafe();
}


} // namespace STK

#endif /* STK_ARRAYBASEDOT_H */
