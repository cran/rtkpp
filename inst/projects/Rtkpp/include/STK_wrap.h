/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  rtkpp
 * created on: 22 août 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_wrap.h
 *  @brief In this file we define a stk++ wrapper converter for expressions.
 **/


#ifndef STK_WRAP_H
#define STK_WRAP_H

namespace STK
{
namespace hidden
{
// forward declaration
template<class Derived, int structure = Traits<Derived>::structure_ >
struct WrapHelper;


/** @ingroup hidden
 * Specialization of WrapHelpher for vector_
 */
template<class Derived>
struct WrapHelper<Derived, Arrays::vector_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Vector<Rtype_> Result;

  static SEXP wrapImpl(Derived const& vec)
  {
    Result res(vec.size());
    for(int i=vec.begin(), iRes=0; i< vec.end(); i++, iRes++)
    { res(iRes) = vec.elt(i);}
    return Rcpp::wrap(res);
  }
};

// specialization for vector_
template<class Derived>
struct WrapHelper<Derived, Arrays::point_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Vector<Rtype_> Result;

  static SEXP wrapImpl(Derived const& vec)
  {
    Result res(vec.size());
    for(int i=vec.begin(), iRes=0; i< vec.end(); i++, iRes++)
    { res(iRes) = vec.elt(i);}
    return Rcpp::wrap(res);
  }
};

// specialization for array2D_
template<class Derived>
struct WrapHelper<Derived, Arrays::array2D_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Matrix<Rtype_> Result;

  static SEXP wrapImpl(Derived const& matrix)
  {
    Result res(matrix.sizeRows(), matrix.sizeCols());
    for(int j=matrix.beginCols(), jRes=0; j< matrix.endCols(); j++, jRes++)
    {
      for(int i=matrix.beginRows(), iRes=0; i< matrix.endRows(); i++, iRes++)
      { res(iRes, jRes) = matrix.elt(i,j);}
    }
    return Rcpp::wrap(res);
  }
};

// specialization for array2D_
template<class Derived>
struct WrapHelper<Derived, Arrays::square_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Matrix<Rtype_> Result;

  static SEXP wrapImpl(Derived const& matrix)
  {
    Result res(matrix.sizeRows(), matrix.sizeCols());
    for(int j=matrix.beginCols(), jRes=0; j< matrix.endCols(); j++, jRes++)
    {
      for(int i=matrix.beginRows(), iRes=0; i< matrix.endRows(); i++, iRes++)
      { res(iRes, jRes) = matrix.elt(i,j);}
    }
    return Rcpp::wrap(res);
  }
};

} // namespace hidden

/** template method allowing to wrap expressions in a RMatrix and a SEXP
 *  @param expr the expression to wrap
 *  @return a copy of expr in a SEXP
 **/
template<typename Derived>
SEXP wrap( ExprBase<Derived> const& expr)
{ return hidden::WrapHelper<Derived>::wrapImpl(expr.asDerived());}

} // namespace STK


#endif /* STK_WRAP_H */
