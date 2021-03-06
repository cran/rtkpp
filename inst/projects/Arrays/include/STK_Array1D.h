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
 * Project:  stkpp::Arrays
 * created on: 26 nov. 2007
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Array1D.h
  * @brief Implementation of the final class Array1D
 **/

#ifndef STK_ARRAY1D_H
#define STK_ARRAY1D_H

#include <Arrays/include/STK_ExprBase.h>
#include "STK_IArray1D.h"

namespace STK
{

template<class Type, int Size_ = UnknownSize > class Array1D;


namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for Array1D class.
 **/
template<class Type_, int Size_>
struct Traits< Array1D<Type_, Size_> >
{
  typedef Array1D<Type_, 1> Row;
  typedef Array1D<Type_, Size_> Col;
  typedef Array1D<Type_, UnknownSize> SubRow;
  typedef Array1D<Type_, UnknownSize> SubCol;
  typedef Array1D<Type_, UnknownSize> SubArray;
  typedef Array1D<Type_, UnknownSize> SubVector;

  typedef Type_ Type;
  typedef typename RemoveConst<Type_>::Type const& ReturnType;

  enum
  {
    structure_ = Arrays::vector_,
    orient_    = Arrays::by_col_,
    size_      = Size_,
    sizeCols_  = 1,
    sizeRows_  = Size_,
    storage_   = Arrays::dense_ // always dense
  };
};

} // namespace hidden

/** @ingroup Arrays
 *  @brief Templated one dimensional Arrays.
 * 
 * An Array1D is a templated non-oriented container (even if the @c Traits
 * struct define it as column oriented) implementing the interface
 * class @c IArray1D.
 *
 * @note It is a final class for the curious recursive paradigm.
 *
 * @tparam the Type of the objects stored in the @c Array1D
 **/
template<class Type, int Size_ >
class Array1D : public IArray1D< Array1D<Type, Size_> >
{
  public:
    typedef typename hidden::Traits<Array1D<Type, Size_> >::Row Row;
    typedef typename hidden::Traits<Array1D<Type, Size_> >::Col Col;
    typedef typename hidden::Traits<Array1D<Type, Size_> >::SubRow SubRow;
    typedef typename hidden::Traits<Array1D<Type, Size_> >::SubCol SubCol;
    typedef typename hidden::Traits<Array1D<Type, Size_> >::SubVector SubVector;
    typedef typename hidden::Traits<Array1D<Type, Size_> >::SubArray SubArray;
    /** Type for the Array1DBase Class. */
    typedef IArray1D< Array1D<Type, Size_> > Base;
    /** Default constructor. */
    Array1D() : Base(){}
    /** constructor with a specified Range
     *  @param I range of the container
     **/
    Array1D( Range const& I) : Base(I) {}
    /** Misc constructor with beg and end, initialization with a constant.
     *  @param I range of the container
     *  @param v initial value of the container
     **/
    Array1D( Range const& I, Type const& v) : Base(I, v) {}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    Array1D( const Array1D &T, bool ref =false) : Base(T, ref) {}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I range of the data to wrap
     **/
    Array1D( Array1D const& T, Range const& I) : Base(T, I) {}
    /** Wrapper constructor : the container is a reference of a C-Array.
     *  @param q pointer on data
     *  @param I range of the data
     **/
    Array1D( Type* q, Range const& I) : Base(q, I) {}
    /** destructor: allocated memory is liberated by AllocatorBase base class. */
    ~Array1D() {}
    /** operator = : overwrite the Array1D with T.
     *  We resize the object if this and T does not have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
    Array1D& operator=(Array1D const& T)
    {
      // check size
      if (this->range()!=T.range()) this->resize(T.range());
      // copy without ovelapping.
      if (this->begin() < T.begin())
      { for (int i=this->begin(), j=T.begin(); i<this->end(); i++, j++) this->elt(i) = T.elt(j);}
      else
      { for (int i=this->lastIdx(), j=T.lastIdx(); i>=this->begin(); i--, j--) this->elt(i) = T.elt(j);}
      return *this;
    }
    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    Array1D& operator=(Type const& v)
    {
      for (int i=this->begin(); i<this->end(); i++) this->elt(i)= v;
      return *this;
    }
    /** Copy an other type of 1D array in an Array1D.
     *  @param T the array to copy
     **/
    template<class OtherArray>
    Array1D& operator=(ITContainer1D<OtherArray> const& T)
    {
      // check size
      if (this->range()!=T.range()) this->resize(T.range());
      for (int i=this->begin(), j=T.begin(); i<this->end(); i++, j++) this->elt(i) = T.elt(j);
      return *this;
    }
    /** Copy an other type of array/expression in an Array1D.
     *  @param T the array/expression to copy
     **/
    template<class OtherArray>
    Array1D& operator=(ExprBase<OtherArray> const& T)
    {
      // check size
      if (this->size()!=T.size()) this->resize(T.range());
      for (int i=this->begin(); i<this->end(); i++) this->elt(i)= T[i];
      return *this;
    }
};

} // namespace STK

#endif // STK_ARRAY1D_H
