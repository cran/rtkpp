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
 * Project:  stkpp::Array
 * created on: 10 août 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ICArray.h
 *  @brief Interface base class for the CArray, this is an internal header file,
 *  included by other Containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like CArray, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_DENSEARRAYBASE_H
#define STK_DENSEARRAYBASE_H

#include "STK_ArrayBase.h"

#include "STK_ExprBaseVisitor.h"
#include "STK_ExprBaseDot.h"
#include "STK_ExprBaseProduct.h"
#include "STK_ArrayBaseApplier.h"
#include "STK_ArrayBaseAssign.h"
#include "STK_ArrayBaseInitializer.h"

namespace STK
{

namespace hidden
{
template<class Derived, int Structure_ = Derived::structure_>
struct CheckArray;

template<class Derived>
struct CheckArray<Derived, Arrays::array2D_>
{
  static bool resize(Derived const& array, Range const& I, Range const& J)
  { return (array.rows() == I && array.cols() == J) ? false : true;}
  static bool shift(Derived const& array, int beginRow, int beginCol)
  { return (array.beginRows() == beginRow && array.beginCols() == beginCol) ? false : true;}
  static bool resize(Derived const& array, Range const& I)
  { return (array.rows() == I && array.cols() == I) ? false : true;}
  static bool shift(Derived const& array, int begin)
  { return (array.beginRows() == begin && array.beginCols() == begin) ? false : true;}
};
template<class Derived>
struct CheckArray<Derived, Arrays::upper_triangular_>
{
  static bool resize(Derived const& array, Range const& I, Range const& J)
  { return (array.rows() == I && array.cols() == J) ? false : true;}
  static bool shift(Derived const& array, int beginRow, int beginCol)
  { return (array.beginRows() == beginRow && array.beginCols() == beginCol) ? false : true;}
  static bool resize(Derived const& array, Range const& I)
  { return (array.rows() == I && array.cols() == I) ? false : true;}
  static bool shift(Derived const& array, int begin)
  { return (array.beginRows() == begin && array.beginCols() == begin) ? false : true;}
};
template<class Derived>
struct CheckArray<Derived, Arrays::lower_triangular_>
{
  static bool resize(Derived const& array, Range const& I, Range const& J)
  { return (array.rows() == I && array.cols() == J) ? false : true;}
  static bool shift(Derived const& array, int beginRow, int beginCol)
  { return (array.beginRows() == beginRow && array.beginCols() == beginCol) ? false : true;}
  static bool resize(Derived const& array, Range const& I)
  { return (array.rows() == I && array.cols() == I) ? false : true;}
  static bool shift(Derived const& array, int begin)
  { return (array.beginRows() == begin && array.beginCols() == begin) ? false : true;}
};

template<class Derived, int Structure_>
struct CheckArray
{
  static bool resize(Derived const& array, Range const& I)
  { return (array.range() == I) ? false : true;}
  static bool shift(Derived const& array, int begin)
  { return (array.begin() == begin) ? false : true;}
};
}
/** @class ICArray
  * @ingroup Arrays
  *
  * @brief Interface class for all CArray, CArrayPoint, CArrayVector and CArraySquare.
  *
  * This class is the base that is inherited by all objects (matrix, vector,
  * point) which are not expression and stored as CArrays. The common API for
  * these objects is contained in this class.
  *
  * This is essentially a wrapper of a CAllocator
  *
  * @tparam Derived is the derived type, e.g., a matrix type.
  */
template<class Derived>
class ICArray : public ArrayBase<Derived>
{
  public:
    typedef ArrayBase<Derived> Base;

    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Row  Row;
    typedef typename hidden::Traits<Derived>::Col  Col;
    typedef typename hidden::Traits<Derived>::SubRow SubRow;
    typedef typename hidden::Traits<Derived>::SubCol SubCol;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;
    typedef typename hidden::Traits<Derived>::SubArray SubArray;

    typedef typename hidden::Traits<Derived>::Allocator Allocator;
    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    typedef hidden::CheckArray<Derived, structure_> Checker;
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** allocator of the memory  */
    Allocator allocator_;
    /** default constructor. */
    ICArray() : Base(), allocator_() {}
    /** constructor with specified sizes.
     *  @param sizeRows,sizeCols size of the rows and columns
     **/
    ICArray( int sizeRows, int sizeCols)
           : Base(), allocator_(sizeRows, sizeCols)
    {}
    /** constructor with specified sizes and value.
     *  @param sizeRows,sizeCols size of the rows and columns
     *  @param value the value to set
     **/
    ICArray( int sizeRows, int sizeCols, Type const& value)
           : Base(), allocator_(sizeRows, sizeCols, value)
    {}
    /** copy or wrapper constructor.
     *  @param T size of the rows
     *  @param ref is this owning its own data ?
     **/
    ICArray( Derived const& T, bool ref = false)
           : Base(), allocator_(T.allocator_, ref)
    {}
    /** wrapper constructor for 0 based C-Array.
     *  @param q pointer on the array
     *  @param sizeRows,sizeCols size of the rows and columns
     **/
    ICArray( Type* const& q, int sizeRows, int sizeCols)
           : Base(), allocator_(q, sizeRows, sizeCols)
    {}
    /** constructor by reference, ref_=1.
     *  @param allocator the allocator to wrap
     *  @param I,J range of the rows and columns to wrap
     **/
    template<class OtherAllocator>
    inline ICArray( ITContainer2D<OtherAllocator> const& allocator, Range const& I, Range const& J)
                  : Base(), allocator_(allocator.asDerived(), I, J)
    {}
    /** constructor by reference, ref_=1.
     *  @param allocator with the data
     **/
    template< class OtherAllocator>
    inline ICArray( ITContainer2D<OtherAllocator> const& allocator)
                  : Base(), allocator_(allocator.asDerived(), true)
    {}
    /**  destructor */
    ~ICArray() {}

  public:
    /** @return the Horizontal range */
    inline ColRange const&colsImpl() const { return allocator_.cols();};
    /** @return the Vertical range */
    inline RowRange const& rowsImpl() const { return allocator_.rows();}

    /** @return @c true if the container is empty, @c false otherwise */
    bool empty() const { return allocator_.empty();}
    /** @return @c true if *this is reference container, @c false otherwise */
    bool isRef() const { return allocator_.isRef();}

    /** Get a constant reference on the main allocator. */
    inline Allocator const& allocator() const { return allocator_;}
    /** Get the constant main pointer. */
    inline Type const* p_data() const { return allocator_.p_data();}
    /** Get the writable main pointer. */
    inline Type* p_data() { return allocator_.p_data();}

    /** implement the const element accessor */
    inline Type& elt2Impl( int i, int j) { return allocator_.elt(i, j);}
    /** implement the writable element accessor */
    inline Type const& elt2Impl( int i, int j) const { return allocator_.elt(i, j);}

    /** implement the const element accessor for vector/point/diagonal arrays*/
    inline Type& elt1Impl( int j) { return allocator_.elt(j);}
    /** implement the writable element accessor for vector/point/diagonal arrays*/
    inline Type const& elt1Impl( int j) const { return allocator_.elt(j);}

    /** implement the const element accessor for number arrays*/
    inline Type& elt0Impl() { return allocator_.elt();}
    /** implement the writable element accessor for number arrays*/
    inline Type const& elt0Impl() const { return allocator_.elt();}

    /** implement the row operator using a reference on the row of the allocator */
    inline Row rowImpl(int i) const { return  Row( allocator_.row(i));}
    /** implement the row operator using a reference on the row of the allocator */
    inline SubRow rowImpl(int i, Range const& J) const { return SubRow( allocator_.row( i, J));}

    /** implement the col operator using a reference on the column of the allocator */
    inline Col colImpl(int j) const { return  Col( allocator_.col(j));}
    /** implement the col operator using a reference on the column of the allocator */
    inline SubCol colImpl(Range const& I, int j) const { return SubCol( allocator_.col( I, j));}

    /** implement the row operator using a reference on the rows of the allocator */
    inline SubArray rowImpl(Range const& I) const { return SubArray( allocator_.sub(I, this->cols()));}
    /** implement the col operator using a reference on the columns of the allocator */
    inline SubArray colImpl(Range const& J) const { return SubArray( allocator_.sub( this->rows(), J));}
    /** implement the sub operator for 2D arrays using a reference on the column of the allocator */
    inline SubArray subImpl(Range const& I, Range const& J) const { return SubArray(allocator_.sub(I, J));}

    /** implement the sub operator for 1D arrays using a reference on the raw/column of the allocator */
    inline SubVector subImpl( Range const& J) const { return SubVector( allocator_.sub(J));}

    /** swap two elements: only for vectors an points. */
    void swap(int i, int  j) { std::swap(this->elt(i), this->elt(j)); }
    /** @param pos1, pos2 positions of the columns to swap */
    void swapCols(int pos1, int pos2)
    {
      if (this->beginCols() > pos1)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapCols,pos1, pos2,beginCols() >pos1);}
      if (this->lastIdxCols() < pos1)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapCols,pos1, pos2,lastIdxCols() <pos1);}
      if (this->beginCols() > pos2)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapCols,pos1, pos2,beginCols() >pos2);}
      if (this->lastIdxCols() < pos2)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapCols,pos1, pos2,lastIdxCols() <pos2);}
      // swap allocator
      allocator_.swapCols(pos1, pos2);
    }
    /** @param pos1, pos2 positions of the rows to swap */
    void swapRows(int pos1, int pos2)
    {
      if (this->beginRows() > pos1)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapRows,pos1, pos2,beginRows() >pos1);}
      if (this->lastIdxRows() < pos1)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapRows,pos1, pos2,lastIdxRows() <pos1);}
      if (this->beginRows() > pos2)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapRows,pos1, pos2,beginRows() >pos2);}
      if (this->lastIdxRows() < pos2)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapRows,pos1, pos2,lastIdxRows() <pos2);}
      // swap allocator
      allocator_.swapRows(pos1, pos2);
    }
    /** exchange this with T.
     *  @param T the container to exchange with this
     **/
    void exchange(Derived& T) { allocator_.exchange(T.allocator_);}
    /** move T to this.
     *  @param T the array to move
     **/
    void move(Derived const& T) { allocator_.move(T.allocator_);}
    /** shift the Array.
     *  @param beginRows,beginCols  first indexes of the rows and columns
     **/
    Derived& shift(int beginRows, int beginCols)
    {
      STK_STATIC_ASSERT((structure_ == (int)Arrays::array2D_)
                      ||(structure_ == (int)Arrays::lower_triangular_)
                      ||(structure_ == (int)Arrays::upper_triangular_)
                      ,YOU_CANNOT_USED_THIS_METHOD_WITH_THIS_KIND_OF_ARRAY);
      if (!Checker::shift(this->asDerived(), beginRows, beginCols)) return this->asDerived();
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(ICArray::shift,beginRows,beginCols,cannot operate on reference);}
      allocator_.shift(beginRows, beginCols);
      return this->asDerived();
    }
    /** shift the Array.
     *  @param firstIdx first index of the vector/point/diagonal/square array.
     *  @note if this method is used with arrays, upper triangular and lower triangular
     *  arrays, both indexes will be shifted.
     **/
    Derived& shift(int firstIdx)
    {
      if (!Checker::shift(this->asDerived(), firstIdx)) return this->asDerived();
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(ICArray::shift,firstIdx,cannot operate on reference);}
      allocator_.shift(firstIdx);
      return this->asDerived();
    }
    /** resize the Array.
     *  @param I, J range of the rows and columns
     **/
    Derived& resize(Range const& I, Range const& J)
    {
      STK_STATIC_ASSERT((structure_ == (int)Arrays::array2D_)
                      ||(structure_ == (int)Arrays::lower_triangular_)
                      ||(structure_ == (int)Arrays::upper_triangular_)
                      ,YOU_CANNOT_USED_THIS_METHOD_WITH_THIS_KIND_OF_ARRAY);
      if (!Checker::resize(this->asDerived(), I, J)) return this->asDerived();
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(ICArray::resize,I,J,cannot operate on reference);}
      allocator_.resize(I.size(), J.size()).shift(I.begin(), J.begin());
      return this->asDerived();
    }
    /** Resize the vector/point/diagonal/square array.
     *  @param I Range of the vector
     *  @note if this method is used with arrays, upper triangular and lower triangular
     *  arrays, both ranges will be resized.
     **/
    Derived& resize(Range const& I)
    {
      if (!Checker::resize(this->asDerived(), I)) return this->asDerived();
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(ICArray::resize,I,cannot operate on reference);}
      allocator_.resize(I.size()).shift(I.begin());
      return this->asDerived();
    }
};

} // namespace STK

#endif /* STK_DENSEARRAYBASE_H */
