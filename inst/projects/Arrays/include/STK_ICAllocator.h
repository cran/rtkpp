/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff

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
 * created on: 10 août 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ICAllocator.h
 *  @brief In this file we define the ICAllocator interface class.
 **/

#ifndef STK_ICALLOCATOR_H
#define STK_ICALLOCATOR_H

#include "STK_ITContainer.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Interface base class for 2D containers.
 *
 *  The ICAllocatorBase class is the base class for all two-dimensional containers.
 *  A two-dimensional container is defined by an horizontal range of index
 *  for the columns and a vertical range of index for the rows.
 **/
template<int SizeRows_, int SizeCols_>
class ICAllocatorBase
{
  public:
    /** Type of the Range for the rows */
    typedef TRange<SizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<SizeCols_> ColRange;

  protected:
    /** Default constructor. cols_ = 1:0 and rows_ = 1:0. */
    inline ICAllocatorBase() : rows_(), cols_() {}
    /** Constructor with specified ranges
     *  @param I the vertical range
     *  @param J the horizontal range
     **/
    inline ICAllocatorBase( RowRange const& I, ColRange const& J) : rows_(I), cols_(J) {}
    /** Copy constructor
     *  @param T the container to copy
     **/
    inline ICAllocatorBase( ICAllocatorBase const& T) : rows_(T.rows_), cols_(T.cols_) {}
    /** destructor. **/
    inline ~ICAllocatorBase() {}

    /** @return the range of the columns */
    inline ColRange const& colsImpl() const { return cols_;}
    /** @return the index of the first column */
    inline int beginColsImpl() const { return cols_.begin();}
    /** @return the ending index of the columns */
    inline int endColsImpl() const { return cols_.end();}
    /** @return the number of column */
    inline int sizeColsImpl() const { return cols_.size();}

//    /** @return the range of the rows */
//    inline Range rows() const { return Range(beginRowsImpl(), sizeRowsImpl());}
    /** @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the index of the first row */
    inline int beginRowsImpl() const { return rows_.begin();}
    /** @return the ending index of rows */
    inline int endRowsImpl() const { return rows_.end();}
    /** @return the number of rows */
    inline int sizeRowsImpl() const { return rows_.size();}

    /** @return the index of the last column */
    inline int lastIdxCols() const { return cols_.lastIdx();}
    /** @return the index of the last row */
    inline int lastIdxRows() const { return rows_.lastIdx();}

    /** @return @c true if the container is empty, @c false otherwise */
    inline bool empty() const { return (cols_.empty() || rows_.empty());}

    /** Set the first index of the rows and columns.
     *  @param rbeg the first index of the rows
     *  @param cbeg the first index of the columns
     **/
    inline void shift( int rbeg, int cbeg) { rows_.shift(rbeg); cols_.shift(cbeg);}
    /** Set the ranges of the container.
     *  @param I the vertical range
     *  @param J the horizontal range
     **/
    inline void setRanges(RowRange const& I = RowRange(), ColRange const& J = ColRange())
    { rows_ = I; cols_ =J;}
    /** Set the range of the number of rows.
     *  @param I the range of the rows number
     **/
    inline void setRows( RowRange const& I = RowRange()) { rows_ = I;}
    /** Set the first index of the rows.
     *  @param beg the first index of the rows
     **/
    inline void shiftBeginRows( int beg) { rows_.shift(beg);}
    /** Increment the range of the number of rows.
     *  @param inc the increment to apply
     **/
    inline void incRangeRows( int inc) { rows_.inc(inc);}
    /** Increment the first index of the number of rows.
     *  @param inc the increment to apply
     **/
    inline void incFirstIdxRows( int inc) { rows_.incFirst(inc);}
    /** Decrement the first index of the number of rows.
     *  @param dec the decrement to apply
     **/
    inline void decFirstIdxRows( int dec) { rows_.decFirst(dec);}
    /** Increment the end of the number of rows.
     *  @param inc the increment to apply
     **/
    inline void incLastIdxRows( int inc) { rows_.incLast(inc);}
    /** Decrement the end of the number of rows.
     *  @param dec the decrement to apply
     **/
    inline void decLastIdxRows( int dec) { rows_.decLast(dec);}

    /** Set the range of the columns.
     * @param J the range of the cols number
     **/
    inline void setCols( ColRange const& J = ColRange()) { cols_ = J;}
    /** Shift the first index of the columns to beg.
     *  @param beg the new first index
     **/
    inline void shiftBeginCols( int beg) { cols_.shift(beg);}
    /** Increment the range of the number of columns.
     *  @param inc the increment to apply
     **/
    inline void incRangeCols( int inc) { cols_.inc(inc);}
    /** increment the first index of the number of columns.
     *  @param inc the increment to apply
     **/
    inline void incbeginCols( int inc) { cols_.incFirst(inc);}
    /** Decrement the first index of the columns.
     *  @param dec the decrement to apply
     **/
    inline void decbeginCols( int dec) { cols_.decFirst(dec);}
    /** Increment the last index of the columns.
     *  @param inc the increment to apply
     **/
    inline void incLastIdxCols( int inc)  { cols_.incLast(inc);}
    /** Decrement the last index of the columns.
     *  @param dec the decrement to apply
     **/
    inline void decLastIdxCols( int dec) { cols_.decLast(dec);}
    /** exchange this container with T
     * @param T the container to exchange with this
     **/
     inline void exchange(ICAllocatorBase& T)
     {
       std::swap(T.rows_, this->rows_ );
       std::swap(T.cols_, this->cols_ );
     }
  private:
    /** Vertical range : Range of the indexes for the rows. */
    RowRange rows_;
    /** Horizontal range : Range of the indexes for the columns. */
    ColRange cols_;
};

/** @ingroup Arrays
 *
 * @brief Interface class for homogeneous 2D containers which cannot be
 * expression or part of an expression (like allocators).
 *
 * Use the curious recursive template paradigm : the template
 * parameter @c Derived is the name of the class that
 * implements the interface ICAllocator.
 * For example
 * @code
 * template<class Type>
 * class Derived : public ICAllocator< Derived<Type> >
 * {...}
 * @endcode
 *
 * @sa CAllocator
 **/
template < class Derived, int SizeRows_ = hidden::Traits<Derived>::sizeRows_
                        , int SizeCols_ = hidden::Traits<Derived>::sizeCols_>
class ICAllocator : protected ICAllocatorBase<SizeRows_, SizeCols_>
                  , public ITContainer<Derived, hidden::Traits<Derived>::structure_>
{
  public:
    /** Type of the Range for the rows */
    typedef TRange<SizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<SizeCols_> ColRange;

  protected:
    /** Type of the Base container */
    typedef ICAllocatorBase<hidden::Traits<Derived>::sizeRows_, hidden::Traits<Derived>::sizeCols_ > Base2D;
    /** Type of the Base container */
    typedef ITContainer<Derived, hidden::Traits<Derived>::structure_ > Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    /** Default constructor.*/
    inline ICAllocator() : Base2D(), Base() {}
    /** constructor with specified Range.
     *  @param I,J range of the rows and columns
     **/
    inline ICAllocator( Range const& I, Range const& J) : Base2D(I, J), Base() {}
    /** Copy constructor.
     *  @param T the container to copy
     **/
    inline ICAllocator( ICAllocator const& T) : Base2D(T), Base() {}
    /** destructor. */
    inline ~ICAllocator() {}

  public:
    /**@return the Horizontal range */
    inline ColRange const& colsImpl() const { return Base2D::colsImpl();}
    /** @return the index of the first column */
    inline int beginColsImpl() const { return Base2D::beginColsImpl();}
    /**  @return the ending index of the columns */
    inline int endColsImpl() const { return Base2D::endColsImpl();}
    /** @return the number of columns */
    inline int sizeColsImpl() const { return Base2D::sizeColsImpl();}

    /** @return the Vertical range */
    inline RowRange const& rowsImpl() const { return Base2D::rowsImpl();}
    /** @return the index of the first row */
    inline int beginRowsImpl() const { return Base2D::beginRowsImpl();}
    /** @return the ending index of the rows */
    inline int endRowsImpl() const { return Base2D::endRowsImpl();}
    /** @return the Vertical size (the number of rows) */
    inline int sizeRowsImpl() const { return Base2D::sizeRowsImpl();}
    /** @return the index of the first row */
    inline int beginRows() const { return Base2D::beginRowsImpl();}
    /** @return the ending index of the rows */
    inline int endRows() const { return Base2D::endRowsImpl();}
    /** @return the Vertical size (the number of rows) */
    inline int sizeRows() const { return Base2D::sizeRowsImpl();}

    /**  @return the index of the last column */
    inline int lastIdxCols() const { return Base2D::lastIdxCols();}
    /** @return the index of the last row */
    inline int lastIdxRows() const { return Base2D::lastIdxRows();}

    /** @return @c true if the container is empty, @c false otherwise */
    inline bool empty() const { return Base2D::empty();}
    /** @return the element (i,j) of the 2D container.
     *  @param i, j index of row and of the column
     **/
    inline Type& elt(int i, int j)
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > i)
      { STKOUT_OF_RANGE_2ARG(ICAllocator::elt, i, j, beginRows() > i);}
      if (this->endRows() <= i)
      { STKOUT_OF_RANGE_2ARG(ICAllocator::elt, i, j, endRows() <= i);}
      if (this->beginCols() > j)
      { STKOUT_OF_RANGE_2ARG(ICAllocator::elt, i, j, beginCols() > j);}
      if (this->endCols() <= j)
      { STKOUT_OF_RANGE_2ARG(ICAllocator::elt, i, j, endCols() <= j);}
#endif
      return this->asDerived().elt2Impl(i,j);
    }
    /** @return a constant reference on element (i,j) of the 2D container
     *  @param i, j indexes of the row and of the column
     **/
    inline Type const& elt(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > i)
      { STKOUT_OF_RANGE_2ARG(ICAllocator::elt, i, j, beginRows() > i);}
      if (this->endRows() <= i)
      { STKOUT_OF_RANGE_2ARG(ICAllocator::elt, i, j, endRows() <= i);}
      if (this->beginCols() > j)
      { STKOUT_OF_RANGE_2ARG(ICAllocator::elt, i, j, beginCols() > j);}
      if (this->endCols() <= j)
      { STKOUT_OF_RANGE_2ARG(ICAllocator::elt, i, j, endCols() <= j);}
#endif
      return this->asDerived().elt2Impl(i,j);
    }
    /** @return a reference on the ith element
     *  @param i index of the ith element
     **/
    inline Type& elt(int i)
    {
      STK_STATICASSERT_ONE_DIMENSION_ONLY(Derived)
      return this->asDerived().elt1Impl(i);
    }
    /** @return the constant ith element
     *  @param i index of the ith element
     **/
    inline Type const& elt(int i) const
    {
      STK_STATICASSERT_ONE_DIMENSION_ONLY(Derived)
      return this->asDerived().elt1Impl(i);
    }
    /** @return a reference on the number */
    inline Type& elt()
    {
      STK_STATICASSERT_ZERO_DIMENSION_ONLY(Derived)
      return this->asDerived().elt0Impl();
    }
    /** @return a constant reference on the number */
    inline Type const& elt() const
    {
      STK_STATICASSERT_ZERO_DIMENSION_ONLY(Derived)
      return this->asDerived().elt0Impl();
    }
};


} // namespace STK

#endif /* STK_IALLOCATOR_H */
