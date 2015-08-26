/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015 Serge Iovleff

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
 *  @brief In this file we define the CAllocator class.
 **/

#ifndef STK_CALLOCATOR_H
#define STK_CALLOCATOR_H

#include "STK_ITContainer2D.h"
#include "STK_AllocatorBase.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Allocator for dense Array classes.
 *  The data are stored in two dimensions.
 *  It can be the columns or the rows allocator of any dense container.
 */
template<typename Type, int SizeRows_, int SizeCols_, bool Orient_>
class CAllocator;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CAllocator.
 */
template< typename Type_, int SizeRows_, int SizeCols_, bool Orient_>
struct Traits< CAllocator<Type_, SizeRows_, SizeCols_, Orient_> >
{
  private:
    class Void { };

  public:
    enum
    {
      orient_    = Orient_,
      sizeRows_  = SizeRows_,
      sizeCols_  = SizeCols_,
      sizeProd_  = ProductSizeRowsBySizeCols<SizeRows_, SizeCols_>::prod_,
      storage_   = Arrays::dense_ // always dense
    };

    typedef CAllocator<Type_, 1, 1, Orient_> Number;

    typedef CAllocator<Type_, 1, SizeCols_, Orient_> Row;
    typedef CAllocator<Type_, SizeRows_, 1, Orient_> Col;

    typedef CAllocator<Type_, 1, UnknownSize, Orient_> SubRow;
    typedef CAllocator<Type_, UnknownSize, 1, Orient_> SubCol;

    /** If one of the Size is 1, we have a Vector (a column) or a Point (a row)
     *  (What to do if both are = 1 : Type_ or array (1,1) ?).
     **/
    typedef typename If< (SizeRows_ == 1)||(SizeCols_ == 1)   // one row or one column
                       , typename If<(SizeCols_ == 1)
                                    , typename If<SizeRows_==1, Number, SubCol>::Result
                                    , SubRow>::Result
                       , Void
                       >::Result SubVector;
    /** If one of the Size is 1, we have a Vector (a column) or a Point (a row)
     *  (What to do if both are = 1 : Type_ or array (1,1) ?).
     **/
    typedef typename If< (SizeRows_ >= SqrtUnknownSize)||(SizeCols_ >= SqrtUnknownSize)   // one row or one column
                       , AllocatorBase<Type_, UnknownSize>
                       , AllocatorBase<Type_, sizeProd_>
                       >::Result Allocator;
    typedef Type_  Type;
    typedef typename RemoveConst<Type_>::Type const& ReturnType;

    // use this as default. FIXME: Not optimal in case we just get a SubArray
    // with unmodified rows or cols size.
    typedef CAllocator<Type, UnknownSize, UnknownSize, Orient_> SubArray;
};

} // namespace hidden

// forward declaration
template<class Derived, bool Orient_> class OrientedCAllocator;
template<class Derived, int SizeRows_, int SizeCols_> class StructuredCAllocator;

/**  @ingroup Arrays
 *   @brief Specialization for column-oriented Allocators.*/
template<class Derived>
class OrientedCAllocator<Derived, Arrays::by_col_>: public ITContainer2D<Derived>
{
  public:
    typedef ITContainer2D<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    enum
    {
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    typedef typename hidden::Traits<Derived>::Allocator Allocator;
    typedef IContainer2D<sizeRows_, sizeCols_> LowBase;

  protected:
    /** default constructor */
    OrientedCAllocator( Range const& I, Range const& J)
                             : Base(I, J), ldx_(I.size()), allocator_(prod(I, J))
    {}
    /** copy constructor */
    OrientedCAllocator( OrientedCAllocator const& A, bool ref)
                             : Base(A), ldx_(A.ldx()), allocator_(A.allocator_, ref)
    { if (!ref) allocator_.copy(A.allocator_);}
    /** Reference constructor */
    template<class OtherDerived>
    inline OrientedCAllocator( OrientedCAllocator<OtherDerived, Arrays::by_col_> const& A
                             , Range const& I, Range const& J)
                             : Base(I, J), ldx_(A.ldx()), allocator_(A.allocator(), true)  {}
    /** wrapper constructor for 0 based C-Array*/
    OrientedCAllocator( Type* const& q, int nbRow, int nbCol)
                             : Base(Range(0,nbRow), Range(0,nbCol)), ldx_(nbRow)
                             , allocator_(q, nbRow * nbCol, true)
    {}
    /** destructor */
    ~OrientedCAllocator() {}

  public:
    inline Type* p_data() { return allocator_.p_data();}
    inline Type const* p_data() const { return allocator_.p_data();}
    bool isRef() const { return allocator_.isRef();}
    inline Allocator& allocator() { return allocator_;}
    inline Allocator const& allocator() const { return allocator_;}

    /** @return the index of the allocator*/
    int ldx() const { return ldx_;}
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i, j indexes of the element
     **/
    inline Type const& elt2Impl(int i, int j) const { return p_data()[j*ldx_ + i];}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i, j indexes of the element
     **/
    inline Type& elt2Impl(int i, int j) { return p_data()[j*ldx_ + i];}
    /** set a value to this allocator.
     *  @param v the value to set
     **/
    void setValue(Type const& v)
    {
      for (int j= this->beginCols(); j < this->endCols(); ++j)
        for (int i = this->beginRows(); i < this->endRows(); ++i)
        { this->elt(i, j) = v;}
    }
    /** move T to this without copying data*/
    Derived& move(OrientedCAllocator const& T)
    {
      if (this->asPtrDerived() != T.asPtrDerived())
      {
        allocator_.move(T.allocator_);
        this->asDerived().moveImpl(T.asDerived()); // move any other data fields
        this->setRanges(T.rows(), T.cols());
        ldx_ = T.ldx_;
      }
      return this->asDerived();
    }
  protected:
    /** do nothing by default */
    void moveImpl(OrientedCAllocator const& T) {}
    /** index of the data set */
    int ldx_;
    /** set index of the data. */
    void setLdx( int idx) { ldx_ = idx;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    void exchange(OrientedCAllocator &T)
    {
      allocator_.exchange(T.allocator_);
      Base::exchange(T);
      std::swap(ldx_, T.ldx_);
    }
    /** @brief Compute the range of the 1D Allocator when we want to
     *  allocate a 2D array with  range I for the rows and range J for the
     *  columns.
     *  @param I,J the range of the rows and columns
     *  @return The range of the 1D allocator
     **/
    static Range prod(Range const& I, Range const& J)
    { return Range(I.size()*J.begin()+I.begin(), I.size()*J.size()); }
    /** return the increment to apply to a zero based pointer corresponding to
     *  the actual first row and first column indexes. */
    int shiftInc(int beginRows, int beginCols)
    { return ldx_*beginCols+beginRows; }
    /** set the index corresponding to the actual size of the allocator. */
    void setSizedIdx() {ldx_ = this->sizeRows();}
    /** manager of the memory */
    Allocator allocator_;
};

/**  @ingroup Arrays
 *   @brief Specialization for row-oriented Allocators.*/
template<class Derived>
class OrientedCAllocator<Derived, Arrays::by_row_>: public ITContainer2D<Derived>
{
  public:
    typedef ITContainer2D<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    enum
    {
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    typedef typename hidden::Traits<Derived>::Allocator Allocator;
    typedef IContainer2D<sizeRows_, sizeCols_> LowBase;

  protected:
    /** constructor with specified ranges */
    OrientedCAllocator( Range const& I, Range const& J)
                             : Base(I, J), ldx_(J.size()), allocator_(prod(I, J))
    {}
    /** copy constructor */
    OrientedCAllocator( OrientedCAllocator const& A, bool ref)
                             : Base(A), ldx_(A.ldx_), allocator_(A.allocator_, ref)
    { if (!ref) allocator_.copy(A.allocator_);}
    /** Reference constructor */
    template<class OtherDerived>
    inline OrientedCAllocator( OrientedCAllocator<OtherDerived, Arrays::by_row_> const& A
                             , Range const& I, Range const& J)
                             : Base(I, J), ldx_(A.ldx()), allocator_(A.allocator(), true)
    {}
    /** wrapper constructor for 0 based C-Array*/
    OrientedCAllocator( Type* const& q, int nbRow, int nbCol)
                             : Base(Range(0,nbRow), Range(0,nbCol)), ldx_(nbCol)
                             , allocator_(q, nbRow * nbCol, true)
    {}
    /** destructor */
    ~OrientedCAllocator() {}

  public:
    inline Type* p_data() { return allocator_.p_data();}
    inline Type const* p_data() const { return allocator_.p_data();}
    bool isRef() const { return allocator_.isRef();}
    inline Allocator& allocator() { return allocator_;}
    inline Allocator const& allocator() const { return allocator_;}

    /** @return the index of the allocator*/
    int ldx() const { return ldx_;}
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i,j indexes of the element
     **/
    inline Type const& elt2Impl(int i, int j) const  { return p_data()[i*ldx_ + j];}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i,j indexes of the element
     **/
    inline Type& elt2Impl(int i, int j) { return p_data()[i*ldx_ + j];}
    /** set a value to this container.
     *  @param v the value to set
     **/
    void setValue(Type const& v)
    {
      for (int i = this->beginRows(); i < this->endRows(); ++i)
        for (int j= this->beginCols(); j < this->endCols(); ++j)
        { this->elt(i, j) = v;}
    }
    /** move T to this without copying data*/
    Derived& move(OrientedCAllocator const& T)
    {
      if (this->asPtrDerived() != T.asPtrDerived())
      {
        allocator_.move(T.allocator_);
        this->asDerived().moveImpl(T.asDerived()); // move any other data fields
        this->setRanges(T.rows(), T.cols());
        ldx_ = T.ldx_;
      }
      return this->asDerived();
    }
  protected:
    /** index of the data set */
    int ldx_;
    /** set index of the data. */
    void setLdx( int idx) { ldx_ = idx;}
    /** exchange T with this.
     *  @param T the container to move
     **/
    void exchange(OrientedCAllocator &T)
    {
      allocator_.exchange(T.allocator_);
      Base::exchange(T);
      std::swap(ldx_, T.ldx_);
    }
    /** @brief Compute the range of the 1D Allocator when we want to
     *  allocate a 2D array with I indexes in the first dimension and J indexes
     *  in the second dimension.
     *  @param I,J range of the rows and columns
     *  @return The range of the 1D allocator
     **/
    static Range prod(Range const& I, Range const& J)
    { return Range(J.size()*I.begin()+J.begin(), I.size()*J.size());}
    /** return the increment corresponding to the actual first row an column. */
    int shiftInc(int beginRows, int beginCols)
    { return ldx_*beginRows+beginCols; }
    /** set the index corresponding to the actual size of the allocator. */
    void setSizedIdx() { ldx_ = this->sizeCols();}
    /** manager of the memory */
    Allocator allocator_;
};

/** @ingroup Arrays
 *  @brief  Base class for the general by_col_ structured case.
 **/
template<class Derived, int SizeRows_, int SizeCols_>
class StructuredCAllocator: public OrientedCAllocator<Derived, hidden::Traits<Derived>::orient_>
{
  public:
    enum
    {
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    typedef OrientedCAllocator<Derived, orient_> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Allocator Allocator;
    typedef IContainer2D<sizeRows_, sizeCols_> LowBase;

  protected:
    /** Default constructor */
    StructuredCAllocator( Range const& I, Range const& J): Base(I, J) {}
    /** copy constructor */
    StructuredCAllocator( StructuredCAllocator const& A, bool ref): Base(A, ref) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    StructuredCAllocator( Type* const& q, int nbRow, int nbCol)
                               : Base(q, nbRow, nbCol)
    {}

  public:
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    void exchange(StructuredCAllocator &T) { Base::exchange(T);}
    /** shift the first indexes of the allocator
     *  @param beginRows, beginCols indexes to set
     **/
    void shift2Impl(int beginRows, int beginCols)
    {
      if ((beginRows == this->beginRows())&&(beginCols == this->beginCols())) return;
      LowBase::shift(beginRows, beginCols);
      this->allocator_.shiftData(this->shiftInc(beginRows, beginCols));
    }
    /** shift the first indexes of the allocator (for square and diagonal matrices).
     *  @param begin the index of the first row and column
     **/
    void shift1Impl(int begin) { shift2Impl(begin, begin);}
    /** resize the allocator (for square and diagonal matrices).
     *  @param size the size of the point
     **/
    Derived& resize1Impl(int size)
    { return this->asDerived().resize2Impl(size, size);}
};


/** @ingroup Arrays
 *  @brief specialization for the point_ case.
 **/
template<class Derived, int SizeCols_>
class StructuredCAllocator<Derived, 1, SizeCols_>
    : public OrientedCAllocator<Derived, hidden::Traits<Derived>::orient_>
{
  public:
    enum
    {
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    typedef OrientedCAllocator<Derived, orient_> Base;
    typedef IContainer2D<sizeRows_, sizeCols_> LowBase;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Allocator Allocator;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;

  protected:
    /** Default constructor */
    StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), row_(I.begin())
    {}
    /** copy constructor */
    StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref), row_(A.row_)
    {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J), row_(I.begin())
    {}
    /** wrapper constructor for 0 based C-Array*/
    StructuredCAllocator( Type* const& q, int , int nbCol)
                               : Base(q, 1, nbCol), row_(0)
    {}
  public:
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(row_, T.row_);
    }
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param j index of the column
     **/
    inline Type const& elt1Impl( int j) const { return this->elt2Impl(row_,j);}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param j index of the columns
     **/
    inline Type& elt1Impl( int j) { return this->elt2Impl(row_,j);}
    /** shift the first indexes of the allocator
     *  @param beginRows, beginCols indexes to set
     **/
    void shift2Impl(int beginRows, int beginCols)
    {
      if ((beginRows == this->beginRows())&&(beginCols == this->beginCols())) return;
      LowBase::shift(beginRows, beginCols);
      this->allocator_.shiftData(this->shiftInc(beginRows, beginCols));
      row_ = beginRows;
    }
    /** shift the first column index of the allocator.
     *  @param beginCols the index of the first column */
    void shift1Impl(int beginCols)
    { shift2Impl(beginCols, beginCols);
      row_ = beginCols;
    }
    /** resize the allocator.
     *  @param sizeCols the size of the point
     **/
    Derived& resize1Impl(int sizeCols)
    { this->asDerived().resize2Impl(1, sizeCols);
      row_ = this->beginRows();
      return this->asDerived();
    }
    /** @return a sub-vector in the specified range of the Allocator.
     *  @param J range of the sub-vector
     **/
    SubVector sub1Impl( Range const& J) const { return Base::row(row_, J);}
    /** move T to this.
     *  @param T the container to move
     **/
    void moveImpl(StructuredCAllocator const& T)
    { row_ = T.row_;}

  private:
    /** row of the point (needed when this is a reference) */
    int row_;
};

/** @ingroup Arrays
 *  @brief specialization for the vector_ case.
 **/
template<class Derived, int SizeRows_>
class StructuredCAllocator<Derived, SizeRows_, 1>
    : public OrientedCAllocator<Derived, hidden::Traits<Derived>::orient_>
{
  public:
    enum
    {
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    typedef OrientedCAllocator<Derived, orient_> Base;
    typedef IContainer2D<sizeRows_, sizeCols_> LowBase;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Allocator Allocator;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;

  protected:
    /** Default constructor */
    StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), col_(J.begin())
    {}
    /** copy constructor */
    StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref), col_(A.col_)
    {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J)
                               , col_(J.begin())
    {}
    /** wrapper constructor for 0 based C-Array*/
    StructuredCAllocator( Type* const& q, int nbRow, int)
                               : Base(q, nbRow, 1), col_(0)
                               {}
  public:
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(col_, T.col_);
    }
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i index of the row
     **/
    inline Type const& elt1Impl( int i) const { return this->elt2Impl(i, col_);}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i index of the row
     **/
    inline Type& elt1Impl( int i) { return this->elt2Impl(i, col_);}
    /** shift the first indexes of the allocator
     *  @param beginRows, beginCols indexes to set
     **/
    void shift2Impl(int beginRows, int beginCols)
    {
      if ((beginRows == this->beginRows())&&(beginCols == this->beginCols())) return;
      LowBase::shift(beginRows, beginCols);
      this->allocator_.shiftData(this->shiftInc(beginRows, beginCols));
      col_ = beginCols;
    }
    /** shift the first row index of the allocator.
     *  @param beginRows the index of the first row
     **/
    void shift1Impl(int beginRows)
    { shift2Impl(beginRows, beginRows);
      col_ = beginRows;
    }
    /** resize the allocator.
     *  @param sizeRow the size of the vector
     **/
    Derived& resize1Impl(int sizeRow)
    {
      this->asDerived().resize2Impl(sizeRow, 1);
      col_ = this->beginCols();
      return this->asDerived();
    }
    /** @return a sub-vector in the specified range of the Allocator.
     *  @param I range of the sub-vector
     **/
    SubVector sub1Impl( Range const& I) const { return Base::col(I, col_);}
    /** move T to this.
     *  @param T the container to move
     **/
    void moveImpl(StructuredCAllocator const& T)
    { col_ = T.col_;}

  private:
    int col_;
};


/** @ingroup Arrays
 *  @brief specialization for the number_ case.
 **/
template<class Derived>
class StructuredCAllocator<Derived, 1, 1>
    : public OrientedCAllocator<Derived, hidden::Traits<Derived>::orient_>
{
  public:
    enum
    {
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    typedef OrientedCAllocator<Derived, orient_> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Allocator Allocator;
    typedef IContainer2D<sizeRows_, sizeCols_> LowBase;

  protected:
    /** Default constructor */
    StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), row_(I.begin()), col_(J.begin())
    {}
    /** copy constructor */
    StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref)
                               , row_(A.row_), col_(A.col_)
    {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J)
                               , row_(I.begin()), col_(J.begin())
    {}
    /** wrapper constructor for 0 based C-Array*/
    StructuredCAllocator( Type* const& q, int , int)
                               : Base(q, 1, 1), row_(0), col_(0)
                               {}
    ~StructuredCAllocator() {}
    /** move T to this.
     *  @param T the container to move
     **/
    StructuredCAllocator& move(StructuredCAllocator const& T)
    { row_ = T.row_; col_ = T.col_; return *this;}

  public:
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(row_, T.row_);
      std::swap(col_, T.col_);
    }
    /** @return a constant reference on the element of the Allocator. */
    inline Type const& elt0Impl() const { return this->elt2Impl(row_, col_);}
    /** @return a reference on the element of the Allocator. */
    inline Type& elt0Impl() { return this->elt2Impl(row_, col_);}
    /** @return a constant reference on the element of the Allocator. */
    inline Type const& elt1Impl(int) const { return this->elt2Impl(row_, col_);}
    /** @return a reference on the element of the Allocator. */
    inline Type& elt1Impl(int) { return this->elt2Impl(row_, col_);}

    /** shift the first indexes of the allocator
     *  @param beginRows, beginCols indexes to set
     **/
    void shift2Impl(int beginRows, int beginCols)
    {
      if ((beginRows == this->beginRows())&&(beginCols == this->beginCols())) return;
      LowBase::shift(beginRows, beginCols);
      this->allocator_.shiftData(this->shiftInc(beginRows, beginCols));
      row_ = beginRows;
      col_ = beginCols;
   }
    /** shift the first indexes of the allocator.
     *  @param begin the index of the first row and column
     **/
    void shift1Impl(int begin)
    {
      this->shift2Impl(begin, begin);
      row_ = begin; col_ = begin;
    }
    /** move T to this.
     *  @param T the container to move
     **/
    void moveImpl(StructuredCAllocator const& T)
    { row_ = T.row_; col_ = T.col_;}

  private:
    int row_;
    int col_;
};

/** @ingroup Arrays
 *  @brief Allocator for the dense CArray classes.
 *  The size of the Allocator is known in both dimension
 */
template<typename Type_, int SizeRows_, int SizeCols_, bool Orient_>
class CAllocator: public StructuredCAllocator<CAllocator<Type_, SizeRows_, SizeCols_, Orient_>, SizeRows_, SizeCols_>
{
  public:
    typedef Type_ Type;
    typedef AllocatorBase<Type, SizeRows_* SizeCols_> Allocator;
    typedef StructuredCAllocator<CAllocator, SizeRows_, SizeCols_> Base;
    typedef IContainer2D<SizeRows_, SizeCols_> LowBase;
    using Base::allocator_;
    /** default constructor */
    CAllocator(): Base(SizeRows_, SizeCols_) {}
    CAllocator( int, int): Base(SizeRows_, SizeCols_) {}
    CAllocator( int, int, Type const& v): Base(SizeRows_, SizeCols_) { this->setValue(v);}
    CAllocator( CAllocator const& A, bool ref = true): Base(A, ref)
    { if (!ref) { allocator_.copy(A.allocator_);} }
    template< int OtherSizeRows_, int OtherSizeCols_>
    inline CAllocator( CAllocator<Type, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                     , Range const& I, Range const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    CAllocator( Type* const& q, int , int ): Base(q, SizeRows_, SizeCols_) {}
    ~CAllocator() {}
    CAllocator& move(CAllocator const& T)
    {
      if (this == &T) return *this;
      allocator_.move(T);
      Base::move(T);
      this->setRanges(T.rows(), T.cols());
      this->setLdx(T.ldx());
      return *this;
    }
    CAllocator& resize2Impl( int, int) { return *this;}
    void realloc(int, int) {}
};

/** @brief Specialized Allocator for the dense Arrays classes.
 *  The sizes of the columns and of the rows are unknown. The Orientation is
 *  either by rows or by column.
 */
template<typename Type_, bool Orient_>
class CAllocator<Type_, UnknownSize, UnknownSize, Orient_>
     : public StructuredCAllocator<CAllocator<Type_, UnknownSize, UnknownSize, Orient_>, UnknownSize, UnknownSize>
{
  public:
    typedef Type_ Type;
    typedef AllocatorBase<Type, UnknownSize> Allocator;
    typedef StructuredCAllocator<CAllocator, UnknownSize, UnknownSize> Base;
    typedef IContainer2D<UnknownSize, UnknownSize> LowBase;
    using Base::allocator_;

    /** Default constructor */
    CAllocator(): Base(0, 0)  {}
    /** Constructor with specified size.
     *  @param sizeRows, sizeCols size of the rows and columns
     **/
    CAllocator( int sizeRows, int sizeCols): Base(sizeRows, sizeCols)
    {}
    /** Constructor with specified size and specified value.
     *  @param sizeRows, sizeCols size of the rows and columns
     *  @param v the initial value
     **/
    CAllocator( int sizeRows, int sizeCols, Type const& v)
                     : Base(sizeRows, sizeCols)
    { this->setValue(v);}
    /** Copy or wrapper constructor.
     *  @param A : the array to copy
     *  @param ref : is this a wrapper of A ?
     **/
    CAllocator( CAllocator const& A, bool ref = true): Base(A, ref)
    { if (!ref) { allocator_.copy(A.allocator_);}}
    /** Wrapper constructor. This become a reference on (some part of) the Allocator A.
     *  @param A original allocator
     *  @param I,J range of the rows and columns to wrap.
     **/
    template< int OtherSizeRows_, int OtherSizeCols_>
    inline CAllocator( CAllocator<Type, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                     , Range const& I, Range const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    CAllocator( Type* const& q, int nbRow, int nbCol)
                     : Base(q, nbRow, nbCol)
    {}
    /** Destructor */
    ~CAllocator() {}
    CAllocator& resize2Impl( int sizeRows, int sizeCols)
    {
      if (this->sizeRows() == sizeRows && this->sizeCols() == sizeCols) return *this;
     // check size
     if ((sizeRows <= 0)||(sizeCols<=0))
     {
       // free any allocated memory if this is not a reference
       allocator_.free();
       this->setRanges(sizeRows, sizeCols);
       this->setSizedIdx();
       return *this;
     }
     // allocate
     allocator_.malloc(this->prod(sizeRows, sizeCols));
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
       this->shift2Impl(copy.beginRows(), copy.beginCols());
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
template<typename Type_, int SizeRows_, bool Orient_>
class CAllocator<Type_, SizeRows_, UnknownSize, Orient_>
      : public StructuredCAllocator<CAllocator<Type_, SizeRows_, UnknownSize, Orient_>, SizeRows_, UnknownSize>
{
  public:
    typedef Type_ Type;
    typedef AllocatorBase<Type, UnknownSize> Allocator;
    typedef StructuredCAllocator<CAllocator, SizeRows_, UnknownSize> Base;
    typedef IContainer2D<SizeRows_, UnknownSize> LowBase;
    using Base::allocator_;

    CAllocator(): Base(SizeRows_, 0) {}
    CAllocator( int, int sizeCols): Base(SizeRows_, sizeCols) {}
    CAllocator( int, int sizeCols, Type const& v): Base(SizeRows_, sizeCols)
    { this->setValue(v);}
    CAllocator( CAllocator const& A, bool ref = true)
                     : Base(A, ref)
    { if (!ref) { allocator_.copy(A.allocator_);}}
    template< int OtherSizeRows_, int OtherSizeCols_>
    inline CAllocator( CAllocator<Type, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                     , Range const& I, Range const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    CAllocator( Type* const& q, int , int nbCol)
                     : Base(q, SizeRows_, nbCol)
    {}
    ~CAllocator() {}
    CAllocator& resize2Impl( int, int sizeCols)
    {
      if (this->sizeCols() == sizeCols) return *this;
      // check size
      if (sizeCols<=0)
      {
        // free any allocated memory if it is not a reference
        allocator_.free();
        this->setRanges(SizeRows_, sizeCols);
        this->setSizedIdx();
        return *this;
      }
      // allocate
      allocator_.malloc(this->prod(SizeRows_, sizeCols));
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
        this->shift2Impl(copy.beginRows(), copy.beginCols());
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
template<typename Type_, bool Orient_, int SizeCols_>
class CAllocator<Type_, UnknownSize, SizeCols_, Orient_>
      : public StructuredCAllocator<CAllocator<Type_, UnknownSize, SizeCols_, Orient_>, UnknownSize, SizeCols_>
{
  public:
    typedef Type_ Type;
    typedef AllocatorBase<Type, UnknownSize> Allocator;
    typedef StructuredCAllocator<CAllocator, UnknownSize, SizeCols_> Base;
    typedef IContainer2D<UnknownSize, SizeCols_> LowBase;
    using Base::allocator_;

    CAllocator(): Base(0, SizeCols_) {}
    CAllocator( int sizeRows, int): Base(sizeRows, SizeCols_) {}
    CAllocator( int sizeRows, int, Type const& v)
                     : Base(sizeRows, SizeCols_)
    { this->setValue(v);}
    CAllocator( CAllocator const& A, bool ref = true)
                     : Base(A, ref)
    { if (!ref) { allocator_.copy(A.allocator_);} }

    /** wrap other allocator */
    template< int OtherSizeRows_, int OtherSizeCols_>
    inline CAllocator( CAllocator<Type, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                     , Range const& I, Range const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    CAllocator( Type* const& q, int nbRow, int )
                     : Base(q, nbRow, SizeCols_)
    {}
    ~CAllocator() {}
    CAllocator& resize2Impl( int sizeRows, int)
    {
      if (this->sizeRows() == sizeRows) return *this;
      // check size
      if (sizeRows <= 0)
      {
        // free any allocated memory if it is not a reference
        allocator_.free();
        this->setRanges(sizeRows, SizeCols_);
        this->setSizedIdx();
        return *this;
     }
     // allocate
     allocator_.malloc(this->prod(sizeRows,  SizeCols_));
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
       this->shift2Impl(copy.beginRows(), copy.beginCols());
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

} // namespace STK

#endif /* STK_CALLOCATOR_H */
