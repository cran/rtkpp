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

/** @file STK_IArrayBase.h
 *  @brief In this file we define the ICAllocator interface class.
 **/

#ifndef STK_IARRAYBASE_H
#define STK_IARRAYBASE_H

#include "STKernel/include/STK_Range.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Interface base class for 2D containers.
 *
 *  The IArrayBase class is the base class for all two-dimensional containers.
 *  A two-dimensional container is defined by an horizontal range of index
 *  for the columns and a vertical range of index for the rows.
 **/
template<int SizeRows_, int SizeCols_>
class IArrayBase
{
  public:
    /** Type of the Range for the rows */
    typedef TRange<SizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<SizeCols_> ColRange;

  protected:
    /** Default constructor. cols_ = 1:0 and rows_ = 1:0. */
    inline IArrayBase() : rows_(), cols_() {}
    /** Constructor with specified ranges
     *  @param I the vertical range
     *  @param J the horizontal range
     **/
    inline IArrayBase( RowRange const& I, ColRange const& J) : rows_(I), cols_(J) {}
    /** Copy constructor
     *  @param T the container to copy
     **/
    inline IArrayBase( IArrayBase const& T) : rows_(T.rows_), cols_(T.cols_) {}
    /** destructor. **/
    inline ~IArrayBase() {}

    /** @return the range of the columns */
    inline ColRange const& colsImpl() const { return cols_;}
    /** @return the index of the first column */
    inline int beginColsImpl() const { return cols_.begin();}
    /** @return the ending index of the columns */
    inline int endColsImpl() const { return cols_.end();}
    /** @return the number of column */
    inline int sizeColsImpl() const { return cols_.size();}

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
     inline void exchange(IArrayBase& T)
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


} // namespace STK

#endif /* STK_IARRAYBASE_H */
