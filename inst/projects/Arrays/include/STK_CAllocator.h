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

/** @file STK_CAllocator.h
 *  @brief In this file we define the CAllocator templated class.
 **/

#ifndef STK_CALLOCATOR_H
#define STK_CALLOCATOR_H

#include "STK_ICAllocator.h"
#include "STK_AllocatorBase.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Allocator for dense Array classes.
 *  The data are stored in two dimensions.
 *  It can be the columns or the rows allocator of any dense container.
 */
template<typename Type, Arrays::Structure Structure_, int SizeRows_, int SizeCols_, bool Orient_>
class CAllocator;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CAllocator.
 */
template< typename Scalar, Arrays::Structure Structure_, int SizeRows_, int SizeCols_, bool Orient_>
struct Traits< CAllocator<Scalar, Structure_, SizeRows_, SizeCols_, Orient_> >
{
  private:
    class Void { };

  public:
    typedef Scalar Type;
    typedef CAllocator<Scalar, Arrays::number_, 1, 1, Orient_> Number;


    typedef CAllocator<Scalar, Arrays::point_, 1, SizeCols_, Orient_> Row;
    typedef CAllocator<Scalar, Arrays::vector_, SizeRows_, 1, Orient_> Col;

    typedef CAllocator<Scalar, Arrays::point_, 1, UnknownSize, Orient_> SubRow;
    typedef CAllocator<Scalar, Arrays::vector_, UnknownSize, 1, Orient_> SubCol;

    /** If one of the Size is 1, we have a Vector (a column) or a Point (a row)
     *  (What to do if both are = 1 : Scalar or array (1,1) ?).
     **/
    typedef typename If< (SizeRows_ == 1)||(SizeCols_ == 1)   // one row or one column
                       , typename If<(SizeCols_ == 1)
                                    , typename If<SizeRows_==1, Number, SubCol>::Result
                                    , SubRow>::Result
                       , Void
                       >::Result SubVector;
    // use this as default. FIXME: Not optimal in case we just get a SubArray
    // with unmodified rows or cols size.
    typedef CAllocator<Type, Structure_, UnknownSize, UnknownSize, Orient_> SubArray;

    enum
    {
      structure_ = Structure_,
      orient_    = Orient_,
      sizeRows_  = SizeRows_,
      sizeCols_  = SizeCols_,
      storage_   = Arrays::dense_ // always dense
    };
};

} // namespace hidden

// forward declaration
template<class Derived, int SizeRows_ = hidden::Traits<Derived>::sizeRows_, int SizeCols_ = hidden::Traits<Derived>::sizeCols_>
class CAllocatorBase;
template<class Derived, int SizeRows_, int SizeCols_, bool Orient_> class OrientedCAllocator;
template<class Derived, int SizeRows_, int SizeCols_, bool Orient_> class StructuredCAllocator;


/** @ingroup Arrays
 *  @brief  Interface Base class for the CAllocator classes.
 *
 *  The allocator stores the data in a single vector by row or by column.
 *  In order to obtain the element (i,j) of the Array, we have to apply,
 *  either the formula @c i*idx_+j or @c j*idx_+i. The class deriving from this
 *  Interface cannot be an expression, thus it derive from the ICAllocator
 *  interface. Concrete implementation depend of the structure of the allocator
 *  and of the Orientation.
 *
 *  The pseudo-virtual function to implement have the following description
 *  @code
 *    // for all matrices classes except diagonal and square arrays
 *    void shift2Impl( int firstRow, int firstCol)
 *    Derived& resize2Impl(int sizeRows, int sizeCols)
 *    // for all vector classes, diagonal and square arrays
 *    void shift1Impl(int beg)
 *    Derived& resize1Impl(int size)
 *    SubVector sub1Impl(Range const& I)
 *  @endcode
 *
 *  @sa OrientedCAllocator, StructuredCAllocator, CAllocator
 **/
template<class Derived, int SizeRows_, int SizeCols_>
class CAllocatorBase: public ICAllocator<Derived, SizeRows_, SizeCols_>
{
  protected:
    /** Base class */
    typedef ICAllocator<Derived, SizeRows_, SizeCols_> Base;
    /** Default constructor */
    inline CAllocatorBase( Range const& I, Range const& J): Base(I, J) {}
    /** copy constructor */
    inline CAllocatorBase( CAllocatorBase const& A): Base(A) {}
    /** destructor */
    inline ~CAllocatorBase() {}

  public:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Row Row;
    typedef typename hidden::Traits<Derived>::Col Col;
    typedef typename hidden::Traits<Derived>::SubRow SubRow;
    typedef typename hidden::Traits<Derived>::SubCol SubCol;
    typedef typename hidden::Traits<Derived>::SubArray SubArray;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };

    /** Access to the ith row of the Allocator.
     *  @param i index of the row
     *  @return a reference on the ith row
     **/
    inline Row row(int i) const
    { return Row(this->asDerived(), Range(i,1), this->cols());}
    /** Access to the row (i,J) of the Allocator.
     *  @param i,J index of the row and range of the columns
     *  @return a reference on the ith row
     **/
    inline SubRow row(int i, Range const& J) const
    { return SubRow(this->asDerived(), Range(i,1), J);}
    /** Access to the jth column of the Allocator.
     *  @param j index of the column
     *  @return a reference on the jth column
     **/
    inline Col col(int j) const
    { return Col(this->asDerived(), this->rows(), Range(j,1));}
    /** Access to the column (I,j) of the Allocator.
     *  @param I,j range of the rows and index of the column
     *  @return a reference on the jth column
     **/
    inline SubCol col(Range const& I, int j) const
    { return SubCol(this->asDerived(), I, Range(j,1));}
    /** Access to the sub-part (I,J) of the Allocator.
     *  @param I,J range of the rows and columns
     *  @return a reference on a sub-part of the Allocator
     **/
    inline SubArray sub(Range const& I, Range const& J) const
    { return SubArray(this->asDerived(), I, J);}
    /** Access to a sub-vector. For 1D allocators only.
     *  @param I range of the rows
     *  @return a reference on a sub-part of the Allocaor
     **/
    inline SubVector sub(Range const& I) const
    {
      STK_STATICASSERT_ONE_DIMENSION_ONLY(Derived)
      return this->asDerived().sub1Impl(I);
    }
    /** shift the first indexes of the allocator.
     *  @param firstRow, firstCol indexes of the first row and first column
     **/
    inline void shift( int firstRow, int firstCol)
    { this->asDerived().shift2Impl(firstRow, firstCol);}
    /** resize the allocator
     *  @param sizeRows, sizeCols size of the rows and columns
     **/
   inline Derived& resize(int sizeRows, int sizeCols)
    {
      this->asDerived().resize2Impl(sizeRows, sizeCols);
      return this->asDerived();
    }
    /** shift the first indexes of the vector or point.
     *  @param beg the index of the first row or column
     **/
    inline void shift(int beg) { this->asDerived().shift1Impl(beg);}
    /** Resize the vector or the point
     *  @param size the size to set to the vector
     **/
    inline Derived& resize(int size)
    {
      this->asDerived().resize1Impl(size);
      return this->asDerived();
    }
    /** @param pos1, pos2 position of the first and second columns to swap
     **/
    void swapCols(int pos1, int pos2)
    {
      for (int i=this->beginRows(); i< this->endRows(); ++i)
      { std::swap(this->asDerived().elt2Impl(i, pos1),this->asDerived().elt2Impl(i, pos2));}
    }
    /** @param pos1, pos2 position of the first and second rows to swap
     **/
    void swapRows(int pos1, int pos2)
    {
      for (int j=this->beginCols(); j< this->endCols(); ++j)
      { std::swap(this->asDerived().elt2Impl(pos1, j),this->asDerived().elt2Impl(pos2, j));}
    }
};

/**  @ingroup Arrays
 *   @brief Specialization for column-oriented Allocators.*/
template<class Derived, int SizeRows_, int SizeCols_>
class OrientedCAllocator<Derived, SizeRows_, SizeCols_, Arrays::by_col_>
    : public AllocatorBase<typename hidden::Traits<Derived>::Type>
    , public CAllocatorBase<Derived, SizeRows_, SizeCols_>
{
  protected:
    typedef CAllocatorBase<Derived, SizeRows_, SizeCols_> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef AllocatorBase<Type> Allocator;
    /** default constructor */
    inline OrientedCAllocator( Range const& I, Range const& J)
                             : Allocator(prod(I, J)), Base(I, J), idx_(I.size())
    {}
    /** copy constructor */
    inline OrientedCAllocator( OrientedCAllocator const& A, bool ref)
                             :  Allocator(A, ref), Base(A), idx_(A.idx()) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherizeCols_>
    inline OrientedCAllocator( OrientedCAllocator<OtherDerived, OtherSizeRows_,OtherizeCols_, Arrays::by_col_> const& A
                             , Range const& I, Range const& J)
                             : Allocator(A, true), Base(I, J), idx_(A.idx()) {}
    /** wrapper constructor for 0 based C-Array*/
    inline OrientedCAllocator( Type* const& q, int nbRow, int nbCol)
                             : Allocator(q, nbRow * nbCol, true)
                             , Base(Range(0,nbRow), Range(0,nbCol)), idx_(nbRow)
    {}
    /** destructor */
    inline ~OrientedCAllocator() {}

  public:
    /** @return the index of the allocator*/
    inline int idx() const { return idx_;}
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i, j indexes of the element
     **/
    inline Type const& elt2Impl(int i, int j) const { return this->data(j*idx_ + i);}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i, j indexes of the element
     **/
    inline Type& elt2Impl(int i, int j) { return this->data(j*idx_ + i);}
    /** set a value to this allocator.
     *  @param v the value to set
     **/
    void setValue(Type const& v)
    {
      for (int j= this->beginCols(); j < this->endCols(); ++j)
        for (int i = this->beginRows(); i < this->endRows(); ++i)
        { this->elt(i, j) = v;}
    }
  protected:
    /** index of the data set */
    int idx_;
    /** set index of the data. */
    inline void setIdx( int idx) { idx_ = idx;}
    /** exchange T with this.
     *  @param T the container to move
     **/
    inline void exchange(OrientedCAllocator &T)
    {
      Allocator::exchange(T);
      Base::exchange(T);
      std::swap(idx_, T.idx_);
    }
    /** @brief Compute the range of the 1D Allocator when we want to
     *  allocate a 2D array with  range I for the rows and range J for the
     *  columns.
     *  @param I,J the range of the rows and columns
     *  @return The range of the 1D allocator
     **/
    static inline Range prod(Range const& I, Range const& J)
    { return Range(I.size()*J.begin()+I.begin(), I.size()*J.size()); }
    /** return the increment to apply to a zero based pointer corresponding to
     *  the actual first row and first column indexes. */
    inline int shiftInc(int firstRow, int firstCol)
    { return idx_*firstCol+firstRow; }
    /** set the index corresponding to the actual size of the allocator. */
    inline void setSizedIdx() {idx_ = this->asDerived().sizeRows();}
};

/**  @ingroup Arrays
 *   @brief Specialization for row-oriented Allocators.*/
template<class Derived, int SizeRows_, int SizeCols_>
class OrientedCAllocator<Derived, SizeRows_, SizeCols_, Arrays::by_row_>
    : public AllocatorBase<typename hidden::Traits<Derived>::Type>
    , public CAllocatorBase<Derived, SizeRows_, SizeCols_>
{
  protected:
    typedef CAllocatorBase<Derived, SizeRows_, SizeCols_> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef AllocatorBase<Type> Allocator;
    /** constructor with specified ranges */
    inline OrientedCAllocator( Range const& I, Range const& J)
                             : Allocator(prod(I, J)), Base(I, J), idx_(J.size())
    {}
    /** copy constructor */
    inline OrientedCAllocator( OrientedCAllocator const& A, bool ref)
                             : Allocator(A, ref), Base(A), idx_(A.idx()) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherizeCols_>
    inline OrientedCAllocator( OrientedCAllocator<OtherDerived, OtherSizeRows_,OtherizeCols_, Arrays::by_row_> const& A
                             , Range const& I, Range const& J)
                             : Allocator(A, true), Base(I, J), idx_(A.idx())
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline OrientedCAllocator( Type* const& q, int nbRow, int nbCol)
                             : Allocator(q, nbRow * nbCol, true)
                             , Base(Range(0,nbRow), Range(0,nbCol)), idx_(nbCol)
    {}
    inline ~OrientedCAllocator() {}

  public:
    /** @return the index of the allocator*/
    inline int idx() const { return idx_;}
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i,j indexes of the element
     **/
    inline Type const& elt2Impl(int i, int j) const { return this->data(i*idx_ + j);}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i,j indexes of the element
     **/
    inline Type& elt2Impl(int i, int j) { return this->data(i*idx_ + j);}
    /** set a value to this container.
     *  @param v the value to set
     **/
    void setValue(Type const& v)
    {
      for (int i = this->beginRows(); i < this->endRows(); ++i)
        for (int j= this->beginCols(); j < this->endCols(); ++j)
        { this->elt(i, j) = v;}
    }
  protected:
    /** index of the data set */
    int idx_;
    /** set index of the data. */
    inline void setIdx( int idx) { idx_ = idx;}
    /** exchange T with this.
     *  @param T the container to move
     **/
    inline void exchange(OrientedCAllocator &T)
    {
      Allocator::exchange(T);
      Base::exchange(T);
      std::swap(idx_, T.idx_);
    }
    /** @brief Compute the range of the 1D Allocator when we want to
     *  allocate a 2D array with I indexes in the first dimension and J indexes
     *  in the second dimension.
     *  @param I,J range of the rows and columns
     *  @return The range of the 1D allocator
     **/
    static inline Range prod(Range const& I, Range const& J)
    { return Range(J.size()*I.begin()+J.begin(), I.size()*J.size());}
    /** return the increment corresponding to the actual first row an column. */
    inline int shiftInc(int firstRow, int firstCol)
    { return idx_*firstRow+firstCol; }
    /** set the index corresponding to the actual size of the allocator. */
    inline void setSizedIdx() { idx_ = this->asDerived().sizeCols();}
};

/** @ingroup Arrays
 *  @brief  Base class for the general by_col_ structured case.
 **/
template<class Derived, int SizeRows_, int SizeCols_>
class StructuredCAllocator<Derived, SizeRows_, SizeCols_, Arrays::by_col_>
    : public OrientedCAllocator<Derived, SizeRows_, SizeCols_, Arrays::by_col_>
{
  protected:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, SizeRows_, SizeCols_, Arrays::by_col_> Base;
    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J) {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_col_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int nbRow, int nbCol)
                               : Base(q, nbRow, nbCol)
    {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T) { return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T) { Base::exchange(T);}
  public:
    /** shift the first indexes of the allocator (for square matrices).
     *  @param firstIdx the index of the first row and column
     **/
    void shift1Impl(int firstIdx)
    { this->asDerived().shift2Impl(firstIdx, firstIdx);}
};

/** @ingroup Arrays
 *  @brief  Base class for the general by_row_ structured case.
 **/
template<class Derived, int SizeRows_, int SizeCols_>
class StructuredCAllocator<Derived, SizeRows_, SizeCols_, Arrays::by_row_>
    : public OrientedCAllocator<Derived, SizeRows_, SizeCols_, Arrays::by_row_>
{
  protected:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, SizeRows_, SizeCols_, Arrays::by_row_> Base;
    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J) {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_row_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int nbRow, int nbCol)
                               : Base(q, nbRow, nbCol)
    {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T) { return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T) { Base::exchange(T);}
  public:
    /** shift the first indexes of the allocator (for square matrices).
     *  @param firstIdx the index of the first row and column
     **/
    void shift1Impl(int firstIdx)
    { this->asDerived().shift2Impl(firstIdx, firstIdx);}
};

/** @ingroup Arrays
 *  @brief specialization for the point_ case.
 **/
template<class Derived, int SizeCols_>
class StructuredCAllocator<Derived, 1, SizeCols_, Arrays::by_col_>
    : public OrientedCAllocator<Derived, 1, SizeCols_, Arrays::by_col_>
{
  protected:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, 1, SizeCols_, Arrays::by_col_> Base;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;
    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), row_(I.begin()) {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref), row_(A.row_) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_col_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J), row_(I.begin())
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int , int nbCol)
                               : Base(q, 1, nbCol), row_(0)
    {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T)
    { row_ = T.row_; return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T)
    { Base::exchange(T); std::swap(row_, T.row_);}
  public:
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param j index of the column
     **/
    inline Type const& elt1Impl( int j) const { return this->data(j*Base::idx_ + row_);}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param j index of the columns
     **/
    inline Type& elt1Impl( int j) { return this->data(j*Base::idx_ + row_);}
    /** shift the first indexes of the allocator.
     *  @param first the index of the first column and first row */
    inline void shift1Impl(int first)
    { row_ = first; this->asDerived().shift2Impl(first, first);}
    /** resize the allocator.
     *  @param sizeCols the size of the point
     **/
    void resize1Impl(int sizeCols)
    { this->asDerived().resize2Impl(1, sizeCols); row_ = this->beginRows();}
    /** @return a sub-vector in the specified range of the Allocator.
     *  @param J range of the sub-vector
     **/
    inline SubVector sub1Impl( Range const& J) const { return Base::row(row_, J);}

  private:
    /** row of the point (needed when this is a reference) */
    int row_;
};

/** @ingroup Arrays
 *  @brief specialization for the point_ case.
 **/
template<class Derived, int SizeCols_>
class StructuredCAllocator<Derived, 1, SizeCols_, Arrays::by_row_>
    : public OrientedCAllocator<Derived, 1, SizeCols_, Arrays::by_row_>
{
  protected:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, 1, SizeCols_, Arrays::by_row_> Base;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;
    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), row_(I.begin())
                               , p_start_(this->p_data() + row_*J.size())
    {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref), row_(A.row_)
                               , p_start_(this->p_data() + row_*A.idx())
    {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_row_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J), row_(I.begin())
                               , p_start_(this->p_data() + row_*A.idx())
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int , int nbCol)
                               : Base(q, 1, nbCol), row_(0)
                               , p_start_(this->p_data())
    {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T)
    { row_ = T.row_; p_start_ = T.p_start_; return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(row_, T.row_);
      std::swap(p_start_, T.p_start_);
    }
  public:
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param j index of the column
     **/
    inline Type const& elt1Impl( int j) const { return p_start_[j];}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param j index of the columns
     **/
    inline Type& elt1Impl( int j) { return p_start_[j];}
    /** shift the first indexes of the allocator.
     *  @param first the index of the first column and first row */
    inline void shift1Impl(int first)
    {
      this->asDerived().shift2Impl(first, first);
      row_ = first;
      p_start_ = this->p_data() + row_*Base::idx_;
    }
    /** resize the allocator.
     *  @param sizeCols the size of the point
     **/
    void resize1Impl(int sizeCols)
    { this->asDerived().resize2Impl(1, sizeCols);
      row_ = this->beginRows();
      p_start_ = this->p_data() + row_*Base::idx_;
    }
    /** @return a sub-vector in the specified range of the Allocator.
     *  @param J range of the sub-vector
     **/
    inline SubVector sub1Impl( Range const& J) const { return Base::row(row_, J);}

  private:
    /** row of the point (needed when this is a reference) */
    int row_;
    /** starting ptr for 1D arrays */
    Type* p_start_;
};

/** @ingroup Arrays
 *  @brief specialization for the vector_ case.
 **/
template<class Derived, int SizeRows_>
class StructuredCAllocator<Derived, SizeRows_, 1, Arrays::by_col_>
    : public OrientedCAllocator<Derived, SizeRows_, 1, Arrays::by_col_>
{
  protected:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, SizeRows_, 1, Arrays::by_col_> Base;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;

    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), col_(J.begin())
                               , p_start_(this->p_data() + col_*I.size())
    {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref), col_(A.col_)
                               , p_start_(A.p_data() + col_*A.idx())
    {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_col_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J)
                               , col_(J.begin())
                               , p_start_(A.p_data() + col_*A.idx())
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int nbRow, int)
                               : Base(q, nbRow, 1), col_(0)
                               , p_start_(this->p_data())
                               {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T)
    { col_ = T.col_; p_start_ = T.p_start_; return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(col_, T.col_);
      std::swap(p_start_, T.p_start_);
    }

  public:
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i index of the row
     **/
    inline Type const& elt1Impl( int i) const { return p_start_[i];}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i index of the row
     **/
    inline Type& elt1Impl( int i) { return p_start_[i];}
    /** shift the first indexes of the allocator.
     *  @param first the index of the first row and first column
     **/
    inline void shift1Impl(int first)
    {
      this->asDerived().shift2Impl(first, first);
      col_ = first;
      p_start_ = this->p_data() + col_*Base::idx_;
    }
    /** resize the allocator.
     *  @param sizeRow the size of the vector
     **/
    void resize1Impl(int sizeRow)
    {
      this->asDerived().resize2Impl(sizeRow, 1);
      col_ = this->beginCols();
      p_start_ = this->p_data() + col_*Base::idx_;
    }
    /** @return a sub-vector in the specified range of the Allocator.
     *  @param I range of the sub-vector
     **/
    inline SubVector sub1Impl( Range const& I) const { return Base::col(I, col_);}

  private:
    int col_;
    /** starting ptr for 1D arrays */
    Type* p_start_;
};

/** @ingroup Arrays
 *  @brief specialization for the vector_ case.
 **/
template<class Derived, int SizeRows_>
class StructuredCAllocator<Derived, SizeRows_, 1, Arrays::by_row_>
    : public OrientedCAllocator<Derived, SizeRows_, 1, Arrays::by_row_>
{
  protected:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, SizeRows_, 1, Arrays::by_row_> Base;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;
    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), col_(J.begin())
    {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref), col_(A.col_) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_row_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J), col_(J.begin())
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int nbRow, int)
                               : Base(q, nbRow, 1), col_(0)
                               {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T)
    { col_ = T.col_; return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(col_, T.col_);
    }
  public:
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i index of the row
     **/
    inline Type const& elt1Impl( int i) const { return this->data(i*Base::idx_ + col_);}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i index of the row
     **/
    inline Type& elt1Impl( int i) { return this->data(i*Base::idx_ + col_);}
    /** shift the first indexes of the allocator.
     *  @param firstCol the index of the first column
     **/
    inline void shift1Impl(int firstCol)
    { this->asDerived().shift2Impl(firstCol, firstCol);
      col_ = firstCol;
    }
    /** resize the allocator.
     *  @param sizeRow the size of the vector
     **/
    void resize1Impl(int sizeRow)
    { this->asDerived().resize2Impl(sizeRow, 1); col_ = this->beginCols();}
    /** @return a sub-vector in the specified range of the Allocator.
     *  @param I range of the sub-vector
     **/
    inline SubVector sub1Impl( Range const& I) const { return Base::col(I, col_);}
  private:
    int col_;
};

/** @ingroup Arrays
 *  @brief specialization for the number_ case.
 **/
template<class Derived>
class StructuredCAllocator<Derived, 1, 1, Arrays::by_col_>
    : public OrientedCAllocator<Derived, 1, 1, Arrays::by_col_>
{
  protected:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, 1, 1, Arrays::by_col_> Base;
    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), row_(I.begin()), col_(J.begin())
                               , start_(col_*I.size() + row_)
    {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref)
                               , row_(A.row_), col_(A.col_)
                               , start_(col_*A.idx() + row_) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_col_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J)
                               , row_(I.begin()), col_(J.begin())
                               , start_(row_*A.idx() + col_)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int , int)
                               : Base(q, 1, 1), row_(0), col_(0)
                               , start_(0)
                               {}
    inline ~StructuredCAllocator() {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T)
    { row_ = T.row_; col_ = T.col_; start_ = T.start_; return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(row_, T.row_);
      std::swap(col_, T.col_);
      std::swap(start_, T.start_);
    }
  public:
    /** @return a constant reference on the element of the Allocator. */
    inline Type const& elt0Impl() const { return this->data(start_);}
    /** @return a reference on the element of the Allocator. */
    inline Type& elt0Impl() { return this->data(start_);}
    /** @return a constant reference on the element of the Allocator. */
    inline Type const& elt1Impl(int) const { return this->data(start_);}
    /** @return a reference on the element of the Allocator. */
    inline Type& elt1Impl(int) { return this->data(start_);}

    /** shift the first indexes of the allocator.
     *  @param firstIdx the index of the first row and column
     **/
    inline void shift1Impl(int firstIdx)
    {
      this->asDerived().shift2Impl(firstIdx, firstIdx);
      row_ = firstIdx; col_ = firstIdx; start_ = col_*Base::idx_ + row_;
    }

  private:
    int row_;
    int col_;
    /** starting idx for number_ arrays */
    int start_;
};

/** @ingroup Arrays
 *  @brief specialization for the number_ case.
 **/
template<class Derived>
class StructuredCAllocator<Derived, 1, 1, Arrays::by_row_>
    : public OrientedCAllocator<Derived, 1, 1, Arrays::by_row_>
{
  protected:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, 1, 1, Arrays::by_row_> Base;
    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), row_(I.begin()), col_(J.begin())
                               , start_(row_*J.size() + col_)
    {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref)
                               , row_(A.row_), col_(A.col_)
                               , start_(row_*A.idx() + col_) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_row_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J)
                               , row_(I.begin()), col_(J.begin())
                               , start_(row_*A.idx() + col_)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int , int)
                               : Base(q, 1, 1), row_(0), col_(0)
                               , start_(0)
                               {}
    /** destructor */
    inline ~StructuredCAllocator() {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T)
    { row_ = T.row_; col_ = T.col_; start_ = T.start_; return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(row_, T.row_);
      std::swap(col_, T.col_);
      std::swap(start_, T.start_);
    }
  public:
    /** @return a constant reference on the element of the Allocator. */
    inline Type const& elt0Impl() const
    { return this->data(start_);}
    /** @return a reference on the element of the Allocator. */
    inline Type& elt0Impl() { return this->data(start_);}
    /** @return a constant reference on the element of the Allocator. */
    inline Type const& elt1Impl(int) const { return this->data(start_);}
    /** @return a reference on the element of the Allocator. */
    inline Type& elt1Impl(int) { return this->data(start_);}
    /** shift the first indexes of the allocator.
     *  @param firstIdx the index of the first row and column
     **/
    inline void shift1Impl(int firstIdx)
    {
      this->asDerived().shift2Impl(firstIdx, firstIdx);
      row_ = firstIdx; col_ = firstIdx; start_ = col_*Base::idx_ + row_;
    }
  private:
    int row_;
    int col_;
    /** starting idx for number_ arrays */
    int start_;
};

/** @ingroup Arrays
 *  @brief Allocator for the dense CArray classes.
 *  The size of the Allocator is known in both dimension
 */
template<typename Type_, Arrays::Structure Structure_, int SizeRows_, int SizeCols_, bool Orient_>
class CAllocator
      : public StructuredCAllocator<CAllocator<Type_, Structure_, SizeRows_, SizeCols_, Orient_>, SizeRows_, SizeCols_, Orient_  >
{
  public:
    typedef Type_ Type;
    typedef AllocatorBase<Type> Allocator;
    typedef StructuredCAllocator<CAllocator, SizeRows_, SizeCols_, Orient_  > Base;
    inline CAllocator(): Base(SizeRows_, SizeCols_)
    {
      STK_STATICASSERT_POINT_SIZEROWS_MISMATCH(  !((SizeRows_ != 1) && ((Structure_ == Arrays::point_)  || (Structure_ == Arrays::number_))) );
      STK_STATICASSERT_VECTOR_SIZECOLS_MISMATCH( !((SizeCols_ != 1) && ((Structure_ == Arrays::vector_) || (Structure_ == Arrays::number_))) );
      STK_STATICASSERT_SCALAR_SIZE_MISMATCH( !(((SizeRows_ != 1) || ((SizeCols_ != 1))) && (Structure_ == Arrays::number_)) );
    }
    inline CAllocator( int, int): Base(SizeRows_, SizeCols_)
    {
      STK_STATICASSERT_POINT_SIZEROWS_MISMATCH(  !((SizeRows_ != 1) && ((Structure_ == Arrays::point_)  || Structure_ == Arrays::number_)) );
      STK_STATICASSERT_VECTOR_SIZECOLS_MISMATCH( !((SizeCols_ != 1) && ((Structure_ == Arrays::vector_) || Structure_ == Arrays::number_)) );
      STK_STATICASSERT_SCALAR_SIZE_MISMATCH( !(((SizeRows_ != 1) || ((SizeCols_ != 1))) && (Structure_ == Arrays::number_)) );
    }
    inline CAllocator( int, int, Type const& v): Base(SizeRows_, SizeCols_)
    {
      STK_STATICASSERT_POINT_SIZEROWS_MISMATCH(  !((SizeRows_ != 1) && ((Structure_ == Arrays::point_)  || (Structure_ == Arrays::number_))) );
      STK_STATICASSERT_VECTOR_SIZECOLS_MISMATCH( !((SizeCols_ != 1) && ((Structure_ == Arrays::vector_) || Structure_ == Arrays::number_)) );
      STK_STATICASSERT_SCALAR_SIZE_MISMATCH( !(((SizeRows_ != 1) || ((SizeCols_ != 1))) && (Structure_ == Arrays::number_ )) );
      this->setValue(v);
    }
    inline CAllocator( CAllocator const& A, bool ref = true): Base(A, ref)
    { if (!ref) { Allocator::copy(A);} }
    template<Arrays::Structure OtherStruct_, int OtherSizeRows_, int OtherSizeCols_>
    inline CAllocator( CAllocator<Type, OtherStruct_, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                     , Range const& I, Range const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline CAllocator( Type* const& q, int , int ): Base(q, SizeRows_, SizeCols_) {}
    ~CAllocator() {}
    inline void exchange(CAllocator &T) { Base::exchange(T);}
    inline CAllocator& move(CAllocator const& T)
    {
      if (this == &T) return *this;
      Allocator::move(T);
      Base::move(T);
      Base::setRanges(T.rows(), T.cols());
      Base::setIdx(T.idx());
      return *this;
    }
    void shift2Impl(int firstRow, int firstCol)
    {
      if ((firstRow == this->beginRows())&&(firstCol == this->beginCols())) return;
      // check for reference
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(CAllocator::shift2Impl, firstRow, firstCol, cannot operate on reference);}
      // set new ranges and translate main pointer
      ICAllocatorBase<SizeRows_, SizeCols_>::shift(firstRow, firstCol);
      this->shiftData(this->shiftInc(firstRow, firstCol));
    }
    inline CAllocator& resize2Impl( int, int) { return *this;}
    inline void realloc(int, int) {}
};

/** @brief Specialized Allocator for the dense Arrays classes.
 *  The sizes of the columns and of the rows are unknown. The Orientation is
 *  either by rows or by column.
 */
template<typename Type_, Arrays::Structure Structure_, bool Orient_>
class CAllocator<Type_, Structure_, UnknownSize, UnknownSize, Orient_>
     : public StructuredCAllocator<CAllocator<Type_, Structure_, UnknownSize, UnknownSize, Orient_>, UnknownSize, UnknownSize, Orient_ >
{
  public:
    typedef Type_ Type;
    typedef AllocatorBase<Type> Allocator;
    typedef StructuredCAllocator<CAllocator, UnknownSize, UnknownSize, Orient_ > Base;
    /** Default constructor */
    inline CAllocator(): Base(0, 0)  {}
    /** Constructor with specified size.
     *  @param sizeRows, sizeCols size of the rows and columns
     **/
    inline CAllocator( int sizeRows, int sizeCols): Base(sizeRows, sizeCols)
    {}
    /** Constructor with specified size and specified value.
     *  @param sizeRows, sizeCols size of the rows and columns
     *  @param v the initial value
     **/
    inline CAllocator( int sizeRows, int sizeCols, Type const& v)
                     : Base(sizeRows, sizeCols)
    { this->setValue(v);}
    /** Copy or wrapper constructor.
     *  @param A : the array to copy
     *  @param ref : is this a wrapper of A ?
     **/
    inline CAllocator( CAllocator const& A, bool ref = true): Base(A, ref)
    { if (!ref) { Allocator::copy(A);}}
    /** Wrapper constructor. This become a reference on (some part of) the Allocator A.
     *  @param A original allocator
     *  @param I,J range of the rows and columns to wrap.
     **/
    template<Arrays::Structure OtherStruct_, int OtherSizeRows_, int OtherSizeCols_>
    inline CAllocator( CAllocator<Type, OtherStruct_, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                     , Range const& I, Range const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline CAllocator( Type* const& q, int nbRow, int nbCol)
                     : Base(q, nbRow, nbCol)
    {}
    /** Destructor */
    inline ~CAllocator() {}
    /** exchange this with T.
     *  @param T the allocator to exchange
     **/
    inline void exchange(CAllocator &T) { Base::exchange(T);}
    /** move T to this.
     *  @param T the container to move
     **/
    inline CAllocator& move(CAllocator const& T)
    {
      Allocator::move(T);
      Base::move(T);
      ICAllocatorBase<UnknownSize, UnknownSize>::setRanges(T.rows(), T.cols());
      this->setIdx(T.idx());
      return *this;
    }
    void shift2Impl(int firstRow, int firstCol)
    {
      if ((firstRow == this->beginRows())&&(firstCol == this->beginCols())) return;
      // check for reference
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG( CAllocator::shift2Impl, firstRow, firstCol, cannot operate on reference);}
      // set new ranges and  translate main pointer
      ICAllocatorBase<UnknownSize, UnknownSize>::shift(firstRow, firstCol);
      this->shiftData(this->shiftInc(firstRow, firstCol));
    }
    CAllocator& resize2Impl( int sizeRows, int sizeCols)
    {
     // check size
     if ((sizeRows <= 0)||(sizeCols<=0))
     {
       // free any allocated memory if this is not a reference
       this->free();
       // set Range values and null pointer
       this->setPtrData(0, this->prod(sizeRows, sizeCols), false);
       this->setRanges(sizeRows, sizeCols);
       this->setSizedIdx();
       return *this;
     }
     // allocate
     Allocator::malloc(this->prod(sizeRows, sizeCols));
     this->setRanges(sizeRows, sizeCols);
     this->setSizedIdx();
     return *this;
    }
   /** @brief function for memory reallocation.
    *  The function assume you want to increase or reduce the size without
    *  modification of the bases ranges.
    *  @param sizeRows, sizeCols size of the rows and columns
    **/
   void realloc(int sizeRows, int sizeCols)
   {
     if ((sizeRows == this->sizeRows())&&(sizeCols == this->sizeCols())) return;
     // create a copy the original data set
     CAllocator copy;
     this->exchange(copy);
     try
     {
       // create new container
       resize2Impl(sizeRows, sizeCols);
       shift2Impl(copy.beginRows(), copy.beginCols());
       // copy data
       const int endRow = std::min(copy.endRows(), this->endRows());
       const int endCol = std::min(copy.endCols(), this->endCols());
       for (int j= this->beginCols(); j < endCol; ++j)
         for (int i = this->beginRows(); i< endRow; ++i)
         { this->elt(i, j) = copy.elt(i, j);}
     }
     catch (std::bad_alloc & error)  // if an alloc error occur
     {
       this->exchange(copy); // restore the original container
       STKRUNTIME_ERROR_2ARG(CAllocator::realloc, sizeRows, sizeCols, memory allocation failed);
     }
   }
};

/** @brief Specialized Allocator for the dense Arrays classes.
 *  The number of rows is known, the number of columns unknown
 **/
template<typename Type_, Arrays::Structure Structure_, int SizeRows_, bool Orient_>
class CAllocator<Type_, Structure_, SizeRows_, UnknownSize, Orient_>
      : public StructuredCAllocator<CAllocator<Type_, Structure_, SizeRows_, UnknownSize, Orient_>, SizeRows_, UnknownSize, Orient_  >
{
  enum
  {
    isWith1Row_ = ((Structure_ == Arrays::point_)  || (Structure_ == Arrays::number_))
  };
  public:
    typedef Type_ Type;
    typedef StructuredCAllocator<CAllocator, SizeRows_, UnknownSize, Orient_  > Base;
    typedef AllocatorBase<Type> Allocator;
    inline CAllocator(): Base(SizeRows_, 0)
    { STK_STATICASSERT_POINT_SIZEROWS_MISMATCH( !((SizeRows_ != 1) && isWith1Row_) );}
    inline CAllocator( int, int sizeCols): Base(SizeRows_, sizeCols)
    { STK_STATICASSERT_POINT_SIZEROWS_MISMATCH( !((SizeRows_ != 1) && isWith1Row_) );}
    inline CAllocator( int, int sizeCols, Type const& v)
                     : Base(SizeRows_, sizeCols)
    { this->setValue(v);
      STK_STATICASSERT_POINT_SIZEROWS_MISMATCH( !((SizeRows_ != 1) && isWith1Row_) );
    }
    inline CAllocator( CAllocator const& A, bool ref = true)
                     : Base(A, ref)
    { if (!ref) { Allocator::copy(A);}
      STK_STATICASSERT_POINT_SIZEROWS_MISMATCH( !((SizeRows_ != 1) && isWith1Row_) );
    }
    template<Arrays::Structure OtherStruct_, int OtherSizeRows_, int OtherSizeCols_>
    inline CAllocator( CAllocator<Type, OtherStruct_, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                     , Range const& I, Range const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline CAllocator( Type* const& q, int , int nbCol)
                     : Base(q, SizeRows_, nbCol)
    {}
    inline ~CAllocator() {}
    inline void exchange(CAllocator &T) { Base::exchange(T);}
    inline CAllocator& move(CAllocator const& T)
    {
      if (this == &T) return *this;
      Allocator::move(T);
      Base::move(T);
      Base::setRanges(T.rows(), T.cols());
      this->setIdx(T.idx());
      return *this;
    }
    void shift2Impl(int firstRow, int firstCol)
    {
      if ((firstRow == this->beginRows())&&(firstCol == this->beginCols())) return;
      // check for reference
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(CAllocator::shift2Impl, firstRow, firstCol, cannot operate on reference);}
      // set new ranges and  translate main pointer
      ICAllocatorBase<SizeRows_, UnknownSize>::shift(firstRow, firstCol);
      this->shiftData(this->shiftInc(firstRow, firstCol));
    }
    CAllocator& resize2Impl( int, int sizeCols)
    {
      // check size
      if (sizeCols<=0)
      {
        // free any allocated memory if it is not a reference
        this->free();
        // set Range values and null pointer
        this->setPtrData(0, this->prod(SizeRows_, sizeCols), false);
        this->setRanges(SizeRows_, sizeCols);
        this->setSizedIdx();
        return *this;
      }
      // allocate
      Allocator::malloc(this->prod(SizeRows_, sizeCols));
      this->setRanges(SizeRows_, sizeCols);
      this->setSizedIdx();
      return *this;
    }
    void realloc(int, int sizeCols)
    {     // create a copy the original data set
      CAllocator copy;
      this->exchange(copy);
      try
      {
        // create new container
        resize2Impl(SizeRows_, sizeCols);
        shift2Impl(copy.beginRows(), copy.beginCols());
        // copy data
        const int endCol = std::min(copy.endCols(), this->endCols());
        for (int j= this->beginCols(); j < endCol; ++j)
          for (int i = this->beginRows(); i < this->endRows(); ++i)
        { this->elt(i, j) = copy.elt(i, j);}

      }
      catch (std::bad_alloc const& error)  // if an alloc error occur
      {
        this->exchange(copy); // restore the original container
        STKRUNTIME_ERROR_2ARG(CAllocator::realloc, SizeRows_, sizeCols, memory allocation failed);
      }
    }
};

/** @brief Specialized Allocator for the dense Arrays classes.
 *  The sizes of the columns is known, the number of rows is unknown
 */
template<typename Type_, Arrays::Structure Structure_, bool Orient_, int SizeCols_>
class CAllocator<Type_, Structure_, UnknownSize, SizeCols_, Orient_>
      : public StructuredCAllocator<CAllocator<Type_, Structure_, UnknownSize, SizeCols_, Orient_>, UnknownSize, SizeCols_, Orient_  >
{
  public:
    typedef Type_ Type;
    typedef AllocatorBase<Type> Allocator;
    typedef StructuredCAllocator<CAllocator, UnknownSize, SizeCols_, Orient_  > Base;
    inline CAllocator() : Base(0, SizeCols_)
    { STK_STATICASSERT_POINT_SIZEROWS_MISMATCH(  !((SizeCols_ != 1) && ((Structure_ == Arrays::vector_)  || (Structure_ == Arrays::number_))) );}
    inline CAllocator( int sizeRows, int): Base(sizeRows, SizeCols_)
    { STK_STATICASSERT_POINT_SIZEROWS_MISMATCH(  !((SizeCols_ != 1) && ((Structure_ == Arrays::vector_)  || (Structure_ == Arrays::number_))) );}
    inline CAllocator( int sizeRows, int, Type const& v)
                     : Base(sizeRows, SizeCols_)
    {
      STK_STATICASSERT_POINT_SIZEROWS_MISMATCH(  !((SizeCols_ != 1) && ((Structure_ == Arrays::vector_)  || (Structure_ == Arrays::number_))) );
      this->setValue(v);
    }
    inline CAllocator( CAllocator const& A, bool ref = true)
                     : Base(A, ref)
    { if (!ref) { Allocator::copy(A);} }
    /** wrap other allocator */
    template<Arrays::Structure OtherStruct_, int OtherSizeRows_, int OtherSizeCols_>
    inline CAllocator( CAllocator<Type, OtherStruct_, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                     , Range const& I, Range const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline CAllocator( Type* const& q, int nbRow, int )
                     : Base(q, nbRow, SizeCols_)
    {}
    ~CAllocator() {}

    inline void exchange(CAllocator &T) { Base::exchange(T);}
    inline CAllocator& move(CAllocator const& T)
    {
      if (this == &T) return *this;
      Allocator::move(T);
      Base::move(T);
      Base::setRanges(T.rows(), T.cols());
      this->setIdx(T.idx());
      return *this;
    }
    void shift2Impl(int firstRow, int firstCol)
    {
      // check if there is something to do
      if ((firstRow == this->beginRows())&&(firstCol == this->beginCols())) return;
      // check for reference
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(CAllocator::shift2Impl, firstRow, firstCol, cannot operate on reference.);}
      // set new ranges and  translate main pointer
      ICAllocatorBase<UnknownSize, SizeCols_>::shift(firstRow, firstCol);
      this->shiftData(this->shiftInc(firstRow, firstCol));
    }
    CAllocator& resize2Impl( int sizeRows, int)
    {
      // check size
      if (sizeRows <= 0)
      {
        // free any allocated memory if it is not a reference
        this->free();
        // set Range values and null pointer
        this->setPtrData(0, this->prod(sizeRows, SizeCols_), false);
        this->setRanges(sizeRows, SizeCols_);
        this->setSizedIdx();
        return *this;
     }
     // allocate
     Allocator::malloc(this->prod(sizeRows,  SizeCols_));
     this->setRanges(sizeRows,  SizeCols_);
     this->setSizedIdx();
     return *this;
   }
   void realloc(int sizeRows, int)
   {
     // create a copy the original data set
     CAllocator copy;
     this->exchange(copy);
     try
     {
       // create new container
       resize2Impl(sizeRows, SizeCols_);
       shift2Impl(copy.beginRows(), copy.beginCols());
       // copy data
       const int endRow = std::min(copy.endRows(), this->endRows());
       for (int j= this->beginCols(); j < this->endCols(); ++j)
         for (int i = this->beginRows(); i< endRow; ++i)
         { this->elt(i, j) = copy.elt(i, j);}
     }
     catch (std::bad_alloc & error)  // if an alloc error occur
     {
       this->exchange(copy); // restore the original container
       STKRUNTIME_ERROR_2ARG(CAllocator::realloc, sizeRows, SizeCols_, memory allocation failed);
     }
   }
};

/** @ingroup Arrays
 *  ostream for CArray.
 *  @param s the output stream
 *  @param V the CArray to write
 **/
template <typename Type, Arrays::Structure Structure_, int SizeRows_, int SizeCols_, bool Orient_>
ostream& operator<<(ostream& s, const CAllocator<Type, Structure_, SizeRows_, SizeCols_, Orient_>& V)
{ return out2D(s,V);}


} // namespace STK

#endif /* STK_CALLOCATOR_H */
