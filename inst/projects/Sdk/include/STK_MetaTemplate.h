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
 * Project:  stkpp::
 * created on: Sep 23, 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_MetaTemplate.h
 *  @brief This file contains generic meta-programming classes which are not
 *  specifically related to STK++.
 **/


#ifndef STK_METATEMPLATE_H
#define STK_METATEMPLATE_H

namespace STK
{
namespace hidden
{

struct TrueType {  enum { value = 1 }; };
struct FalseType { enum { value = 0 }; };

template<bool Condition, typename Then, typename Else>
struct If;

template<typename Then, typename Else>
struct If <true, Then, Else> { typedef Then Result; };

template<typename Then, typename Else>
struct If <false, Then, Else> { typedef Else Result; };

/* C/C++ fundamental Types */
template<typename T>
struct isArithmetic { enum { yes = false, no = true}; };

/* specializations */
template<> struct isArithmetic<float>         { enum { yes = true, no = false }; };
template<> struct isArithmetic<double>        { enum { yes = true, no = false }; };
template<> struct isArithmetic<long double>   { enum { yes = true, no = false }; };

template<> struct isArithmetic<bool>          { enum { yes = true, no = false }; };
template<> struct isArithmetic<char>          { enum { yes = true, no = false }; };
template<> struct isArithmetic<signed char>   { enum { yes = true, no = false }; };
template<> struct isArithmetic<unsigned char> { enum { yes = true, no = false }; };
template<> struct isArithmetic<signed short>  { enum { yes = true, no = false }; };
template<> struct isArithmetic<unsigned short>{ enum { yes = true, no = false }; };
template<> struct isArithmetic<signed int>    { enum { yes = true, no = false }; };
template<> struct isArithmetic<unsigned int>  { enum { yes = true, no = false }; };
template<> struct isArithmetic<signed long>   { enum { yes = true, no = false }; };
template<> struct isArithmetic<unsigned long> { enum { yes = true, no = false }; };
//template<> struct isArithmetic<signed long long>   { enum { yes = true, no = false }; };
//template<> struct isArithmetic<unsigned long long> { enum { yes = true, no = false }; };

/* C/C++ fundamental Types */
template<typename T>
struct isInt { enum { yes = false, no = true}; };

template<> struct isInt<bool>          { enum { yes = true, no = false }; };
template<> struct isInt<char>          { enum { yes = true, no = false }; };
template<> struct isInt<signed char>   { enum { yes = true, no = false }; };
template<> struct isInt<unsigned char> { enum { yes = true, no = false }; };
template<> struct isInt<signed short>  { enum { yes = true, no = false }; };
template<> struct isInt<unsigned short>{ enum { yes = true, no = false }; };
template<> struct isInt<signed int>    { enum { yes = true, no = false }; };
template<> struct isInt<unsigned int>  { enum { yes = true, no = false }; };
template<> struct isInt<signed long>   { enum { yes = true, no = false }; };
template<> struct isInt<unsigned long> { enum { yes = true, no = false }; };
//template<> struct isInt<signed long long>   { enum { yes = true, no = false }; };
//template<> struct isInt<unsigned long long> { enum { yes = true, no = false }; };

/** check if T and U are of the same type. */
template<typename T, typename U> struct isSame { enum { value = 0 }; };
template<typename T> struct isSame<T,T> { enum { value = 1 }; };

/** @ingroup hidden
  * Convenient struct to Promote the result Type of some binary functors.
  */
template<typename Type1, typename Type2>
struct Promote
{ typedef typename If<(sizeof(Type1) > sizeof(Type2)), Type1, Type2>::Result result_type;};
/** @ingroup STKernel
  * Specialization when we have the same type.
  */
template<typename Type>
struct Promote<Type, Type>
{ typedef Type result_type;};


}// namespace hidden

} // namespace STK

#endif /* STK_METATEMPLATE_H_ */
