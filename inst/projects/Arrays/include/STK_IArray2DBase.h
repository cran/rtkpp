/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2007  Serge Iovleff

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
 * Purpose:  Define the Interface for the Array classes.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_IArray2DBase.h
 *  @brief Interface base class for the Array2D classes, this is an internal
 *  header file, included by other containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array2D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_IARRAY2DBASE_H
#define STK_IARRAY2DBASE_H

#include "STK_IArrayBase.h"
#include "STK_ArrayBase.h"

#include "STK_ExprBaseProduct.h"
#include "STK_ExprBaseDot.h"
#include "STK_ExprBaseVisitor.h"
#include "STK_ArrayBaseApplier.h"
#include "STK_ArrayBaseAssign.h"
#include "STK_ArrayBaseInitializer.h"

#include "STK_Array1D.h"

namespace STK
{
// forward declaration
template < class PTRCOL, class Derived
         , int SizeRows_ = hidden::Traits<Derived>::sizeRows_
         , int SizeCols_ = hidden::Traits<Derived>::sizeCols_>
class IArray2DBase;
/** @ingroup Arrays
 *  @brief Templated interface base class for two-dimensional arrays.
 *
 * A IArray2DBase is an interface class for two-dimensional Arrays
 * stored in columns and having flexible dimensions. It is possible
 * to add, remove easily columns and rows to the Derived class.
 *
 * Each column has a Range stored in the array @c rangeCols_ and a
 * capacity stored in the array @c availableRows_. It should be worth
 * noting that we should have
 * @code
 *   (rangeCols_[j].size() <= availableRows_[j]) == true;
 *   (rangeCols_[j].isIn(this->rows()) == true;
 * @endcode
 *
 * @tparam PTRCOL is the type of the ptr of column in a two-dimensional
 * array: for exemple @c TYPE*, @c Array1D<TYPE>*, @c DBACCESS*....
 * @tparam Derived is the name of the class that implements @c IArray2DBase.
 **/
template < class PTRCOL, class Derived, int SizeRows_, int SizeCols_>
class IArray2DBase :  protected IArrayBase<SizeRows_, SizeCols_>, public ArrayBase<Derived>
{
  template <class OTHERPTRCOL, class OtherDerived, int OtherSizeRows, int OtherSizeCols>
  friend    class IArray2DBase;

  public:
    /** Type of the Range for the rows */
    typedef TRange<SizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<SizeCols_> ColRange;
    /** Type for the IArrayBase base Class. */
    typedef IArrayBase<SizeRows_, SizeCols_ > Base2D;
    /** Type for the Base Class. */
    typedef AllocatorBase<PTRCOL> Allocator;
    /** type of the Base Container Class. */
    typedef ArrayBase<Derived> Base;

    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Row Row;
    typedef typename hidden::Traits<Derived>::Col Col;
    typedef typename hidden::Traits<Derived>::SubRow SubRow;
    typedef typename hidden::Traits<Derived>::SubCol SubCol;
    typedef typename hidden::Traits<Derived>::SubArray SubArray;
    // for 1D container
    typedef typename hidden::Traits<Derived>::SubVector SubVector;


    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };

  protected:
    /** Default constructor */
    IArray2DBase(): Base2D(), Base(), allocator_()
                  , availableRows_(), rangeCols_()
                  , availableCols_(0), capacityByCols_(0)
    { mallocHo(this->cols());}
    /** constructor with specified ranges
     *  @param I range of the Rows
     *  @param J range of the columns
     **/
    IArray2DBase( Range const& I, Range const& J)
                : Base2D(I, J), Base(), allocator_()
                , availableRows_(), rangeCols_()
                , availableCols_(0), capacityByCols_(0)
    { mallocHo(this->cols());}
    /** Copy constructor If we want to wrap T, the main ptr will be wrapped
     *  in AllocatorBase class. If we want to copy  T, Allocator is
     *  initialized to default values.
     *  @note bug correction, we have to use a copy of T.rangeCols_. in case
     *  we are using the code
     *  @code
     *  Array2DVector<TYPE> Dref(D.sub(J), true)
     *  @endcode
     *  we get
     *  @param T the container to copy
     *  @param ref true if we wrap T
     **/
    IArray2DBase( IArray2DBase const& T, bool ref =false)
                : Base2D(T), Base(), allocator_(T.allocator_, ref)
                , availableRows_(T.availableRows_, ref)
                , rangeCols_(T.rangeCols_) // we have to copy it again, in case T is a temporary
                , availableCols_(T.availableCols_), capacityByCols_(T.capacityByCols_)
    { if (!ref) mallocHo(this->cols());}
    /** constructor by reference, ref_=1.
     *  @param T the container to copy
     *  @param I,J ranges of the rows and columns to wrap
     **/
    template<class OtherDerived>
    IArray2DBase( IArray2DBase<PTRCOL, OtherDerived> const& T, Range const& I, Range const& J)
                : Base2D(I, J), Base(), allocator_(T.allocator(), true)
                , availableRows_(T.availableRows(), J) // just a reference
                , rangeCols_(T.rangeCols())        // we have to create it again
                , availableCols_(J.size()), capacityByCols_(I.size())
    {
      for (int j=J.begin(); j<J.end(); j++)
      { rangeCols_[j] = Range::inf(I, T.rangeCols()[j]);}
    }
    /** Wrapper constructor We get a reference of the data.
     *  @param q pointer on data
     *  @param I,J range of the rows and columns to wrap
     **/
    IArray2DBase( PTRCOL* q, Range const& I, Range const& J)
                : Base2D(I, J), Base(), allocator_(q, J, true)
                , availableRows_(J, I.size()), rangeCols_(J, I)
                , availableCols_(I.size()), capacityByCols_(J.size())
    {}
    /** destructor. Allocated horizontal memory (the array with the pointers
     *  on the columns) is liberated by the Allocator.
     **/
    ~IArray2DBase() {}

  public:
    /**@return the Horizontal range */
    inline ColRange const& colsImpl() const { return Base2D::colsImpl();}
    /** @return the index of the first column */
    inline int beginColsImpl() const { return Base2D::beginColsImpl();}
    /**  @return the ending index of columns */
    inline int endColsImpl() const { return Base2D::endColsImpl();}
    /** @return the number of columns */
    inline int sizeColsImpl() const { return Base2D ::sizeColsImpl();}

    /** @return the Vertical range */
    inline RowRange const& rowsImpl() const { return Base2D::rowsImpl();}
    /** @return the index of the first row */
    inline int beginRowsImpl() const { return Base2D::beginRowsImpl();}
    /** @return the ending index of the rows */
    inline int endRowsImpl() const { return Base2D::endRowsImpl();}
    /** @return the number of rows */
    inline int sizeRowsImpl() const { return Base2D::sizeRowsImpl();}

    /**  @return the index of the last column */
    inline int lastIdxCols() const { return Base2D::lastIdxCols();}
    /** @return the index of the last row */
    inline int lastIdxRows() const { return Base2D::lastIdxRows();}

    /**  @return @c true if the container is empty, @c false otherwise */
    inline bool empty() const { return Base2D::empty();}

    /** access to an element.
     *  @note Take care that @c PTRCOL can be accessed using @c operator[]
     *  @param i,j indexes of the row and of the column
     *  @return a reference on the (i,j) element
     **/
    inline Type& elt2Impl( int i, int j) { return this->data(j)[i];}
    /** constant access to an element.
     *  @note Take care that @c PTRCOL can be accessed using @c operator[]
     *  @param i,j indexes of the row and of the column
     *  @return a constant reference on the (i,j) element
     **/
    inline Type const& elt2Impl( int i, int j) const { return this->data(j)[i];}
    /** access to a part of a column.
     *  @param j index of the column
     *  @return a reference in the range I of the column j of this
     **/
    inline Col colImpl( int j) const
    { return Col( this->asDerived(), this->rangeRowsInCol(j), j);}
    /** access to a part of a column.
     *  @param I range of the rows
     *  @param j index of the col
     *  @return a reference in the range I of the column j of this
     **/
    inline SubCol colImpl(Range const& I, int j) const
    { return SubCol( this->asDerived(), Range::inf(I, this->rangeRowsInCol(j)), j);}
    /** access to many columns.
     *  @param J range of the index of the cols
     *  @return a 2D array containing the Container in the Horizontal range @c J
     **/
    inline SubArray colImpl(Range const& J) const
    { return SubArray( this->asDerived(), this->rows(), J);}
    /** access to a part of a row.
     *  @param i index of the row
     *  @return a reference of the row i.
     **/
    inline Row rowImpl( int i) const
    { return Row( this->asDerived(), this->rangeColsInRow(i), i);}
    /** access to a part of a row.
     *  @param i index of the row
     *  @param J range of the columns
     *  @return a reference of the row i.
     **/
    inline SubRow rowImpl(int i, Range const& J) const
    { return SubRow( this->asDerived(), Range::inf(J, this->rangeColsInRow(i)), i);}
    /** access to many rows.
     *  @param I range of the index of the rows
     *  @return a 2D array containing the Container in the vertical range @c I
     **/
    inline SubArray rowImpl(Range const& I) const
    { return SubArray(this->asDerived(), I, this->cols());}
    /** @return  many elements.
     *  @param J Range of the elements
     **/
    inline SubVector subImpl(Range const& J) const
    { return SubVector(this->asDerived(), J);}
    /** access to a sub-array.
     *  @param I,J range of the rows and of the columns
     **/
    inline SubArray subImpl(Range const& I, Range const& J) const
    { return SubArray(this->asDerived(), I, J);}

    /** @return a constant pointer on the j-th column of the container
     *  @param j the index of the column
     **/
    inline PTRCOL const& data(int j) const { return allocator_.data(j);}
    /** @return a constant pointer on the main pointer of the container */
    inline PTRCOL* const& p_data() const { return allocator_.p_data();}
    /**  @return @c true if the container is empty, @c false otherwise */
    inline bool isRef() const { return allocator_.isRef();}
    /** @return a reference on the element (i,j) of the 2D container.
     *  @param i,j indexes of the row and of the column
     **/
    inline Type& operator()(int i, int j) { return this->elt(i,j);}
    /** @return a constant reference on the element (i,j) of the 2D container.
     *  @param i,j indexes of the row and of the column
     **/
    inline Type const operator()(int i, int j) const { return this->elt(i,j);}
    /** @param I range of the index of the rows
     *  @param j index of the column
     *  @return a Vertical container containing the column @c j of this
     *  in the range @c I
     **/
    inline SubCol operator()(Range const& I, int j) const
    { return this->asDerived().col(I, j);}
    /** @param i index of the row
     *  @param J range of the columns
     *  @return an Horizontal container containing the row @c i of this
     *  in the range @c J
     **/
    inline SubRow operator()(int i, Range const& J) const
    { return this->asDerived().row(i, J);}
    /** @param I,J range of the rows and of the columns
     *  @return a 2D container containing this in the range @c I, @c J
     **/
    inline SubArray operator()(Range const& I, Range const& J) const
    { return this->asDerived().sub(I, J);}
    /** @return the column j.
     *  @param j index of the column
     **/
    inline SubCol atCol(int j) const
    {
      if (this->beginCols() > j)
      { STKOUT_OF_RANGE_1ARG(IArray2DBase::atCol, j, beginCols() > j);}
      if (this->endCols() <= j)
      { STKOUT_OF_RANGE_1ARG(IArray2DBase::atCol, j, endCols() <= j);}
      return this->asDerived().col(j);
    }
    /** @return the row i.
     *  @param i the index of the row
     **/
    inline Row atRow(int i) const
    {
      if (this->beginRows() > i)
      { STKOUT_OF_RANGE_1ARG(IArray2DBase::atRow, i, beginRows() > i);}
      if (this->endRows() <= i)
      { STKOUT_OF_RANGE_1ARG(IArray2DBase::at, i, lastIdxRows() < i);}
      return this->asDerived().row(i);
    }
    /** @return the allocator. */
    inline Allocator const& allocator() const { return allocator_;}
    /** @return the maximum possible number of columns without reallocation. */
    inline int availableCols() const { return availableCols_;}
    /** @return the maximum possible number of rows without reallocation for all Cols. */
    inline const Array1D<int>& availableRows() const { return availableRows_;}
    /** @return the capacity to used in a column. */
    inline int capacityByCols() const { return capacityByCols_;}
    /** @return the capacity of the column @c col.
     *  @param col index of the column we want the capacity
     **/
    inline int capacityCol(int col) const { return availableRows_[col];}
    /** @return the range of each columns. */
    inline Array1D<Range> const& rangeCols() const { return rangeCols_;}
    /** @return the range of a column.
     *  @param col index of the column we want the range
     **/
    inline Range const rangeCol(int col) const { return rangeCols_[col];}
    /** internal method for reserving memory for the columns.
     *  @param sizeCols the size to reserve.
     **/
    void reserveCols(int sizeCols)
    {
      if (availableCols_ >= sizeCols) return;
      Range J(this->beginCols(), sizeCols);
      // try to allocate memory
      try
      {
        // re-allocate memory for the columns
        allocator_.realloc(J);
        // initialize availableRows_
        availableRows_.resize(J);
        // initialize this->rangeCols_
        rangeCols_.resize(J);
      }
      catch (runtime_error & error)   // if an error occur
      {
        // set default capacity (0)
        setAvailableCols();
        // set default range
        Base2D::setCols();
        // clear this->availableRows_
        this->availableRows_.clear();
        // clear this->rangeCols_
        this->rangeCols_.clear();
        // throw the error
        throw error;
      }
      // set new capacity if no error occur
      setAvailableCols(sizeCols);
    }
    /** New beginning index for the columns of the object.
     *  @param cbeg the index of the first column to set
     **/
    void shiftBeginCols(int cbeg)
    {
      // if there is something to do
      if ((cbeg - this->beginCols()) != 0)
      {
        // is this structure just a pointer?
        if (this->isRef())
        { STKRUNTIME_ERROR_1ARG(IArray2DBase::shiftBeginCols,cbeg,cannot operate on references);}
        // translate data
        allocator_.shiftData(cbeg);
        // tranlate availableRows_
        availableRows_.shift(cbeg);
        // translate rangeCols_
        rangeCols_.shift(cbeg);
        // adjust dimensions
        Base2D::shiftBeginCols(cbeg);
      }
    }
    /** Swapping two columns.
     *  @param pos1, pos2 positions of the columns to swap
     **/
    void swapCols(int pos1, int pos2)
    {
      if (this->beginCols() > pos1)
      { STKOUT_OF_RANGE_2ARG(IArray2D::swapCols,pos1, pos2,beginCols() >pos1);}
      if (this->endCols() <= pos1)
      { STKOUT_OF_RANGE_2ARG(IArray2D::swapCols,pos1, pos2,endCols() <= pos1);}
      if (this->beginCols() > pos2)
      { STKOUT_OF_RANGE_2ARG(IArray2D::swapCols,pos1, pos2,beginCols() >pos2);}
      if (this->endCols() <= pos2)
      { STKOUT_OF_RANGE_2ARG(IArray2D::swapCols,pos1, pos2,endCols() <=pos2);}
      // swap allocator part
      allocator_.swap(pos1, pos2);
      // swap availableRows_
      std::swap(availableRows_[pos1], availableRows_[pos2]);
      // swap rangeCols_
      std::swap(rangeCols_[pos1], rangeCols_[pos2]);
    }
    /** exchange this container with T.
     *  @param T the container to exchange with this
     **/
    void exchange(IArray2DBase &T)
    {
      // swap AllocatorBase part
      allocator_.exchange(T.allocator_);
      // swap IArrayBase part
      Base2D::exchange(T);
      // swap this part
      std::swap(availableCols_, T.availableCols_);
      std::swap(capacityByCols_, T.capacityByCols_);
      availableRows_.exchange(T.availableRows_);
      rangeCols_.exchange(T.rangeCols_);
    }
    /** swap two elements: only for vectors and points
     * @param i,j indexes of the elemet to swap
     **/
    inline void swap(int i, int  j) { std::swap(this->elt(i), this->elt(j)); }
    /** Append the container @c other to @c this without copying the data
     *  explicitly. The column of @c other are appended to this and
     *  @c other will become a reference container. Observe that the @c const
     *  keyword is not respected in this method: but it is useful to
     *  define this method even for constant objects. The data in itself are not
     *  altered, the Array2D become a reference on its own data.
     *  @param other the container to merge with this
     **/
    template<class Other>
    void merge(IArray2DBase<PTRCOL, Other> const& other)
    {
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_NO_ARG(IArray2DBase::merge(other),*this is a reference.);}
      // is other just a pointer?
      if (other.isRef())
      { STKRUNTIME_ERROR_NO_ARG(IArray2DBase::merge(other),other is a reference.);}
      // if there is no columns, we can safely modify the vertical range
      if (this->sizeCols() <= 0) this->setRows(other.rows());
      // Are ranges corrects ?
      if (this->rows() != other.rows())
      { STKRUNTIME_ERROR_NO_ARG(IArray2DBase::merge(other),this->rows() != other.rows());}
      // break const reference
      IArray2DBase<PTRCOL, Other>& Tref = const_cast<IArray2DBase<PTRCOL, Other>&>(other);
      // compute horizontal range of the container after insertion
      Range cols(this->cols());
      // compute first index of the first column added
      const int first = cols.end();
      // reallocate memory for the columns
      cols.incLast(Tref.sizeCols());
      reallocCols(cols);
      Base2D::setCols(cols);
      // align other range
      Tref.shiftBeginCols(first); // easiest like that
      // copy data from other
      for (int j=first; j< cols.end(); j++) { transferColumn(Tref, j, j);}
      // delete and set view on the data
      Tref.allocator().free();
      Tref.allocator().setPtrData(allocator_.p_data(), allocator_.rangeData(), true);
    }
    /** Append the vector @c other to @c this without copying the data
     *  explicitly. @c other is appended to this and
     *  @c other will become a reference container. The data in itself are not
     *  altered, the Array1D become a reference on its own data.
     *  @param other the container to merge with this
     **/
    template<class Other>
    void merge(IArray1D<Other> const& other)
    {
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_NO_ARG(IArray2DBase::merge(IArray1D),*this is a reference.);}
      // is other just a pointer?
      if (other.isRef())
      { STKRUNTIME_ERROR_NO_ARG(IArray2DBase::merge(IArray1D),other is a reference.);}
      // if there is no columns, we can safely modify the vertical range
      if (this->sizeCols() <= 0) this->setRows(other.range());
      // Are ranges corrects ?
      if (this->rows() != other.range())
      { STKRUNTIME_ERROR_NO_ARG(IArray2DBase::merge(IArray1D),this->rows() != other.range());}
      // compute horizontal range of the container after insertion
      Range cols(this->cols());
      // reallocate memory for the columns
      cols.incLast(1);
      this->reallocCols(cols);
      this->setCols(cols);
      // set column
      data(cols.lastIdx()) = other.p_data();
      availableRows_[cols.lastIdx()] = other.sizeData();
      rangeCols_[cols.lastIdx()] = other.range();
      // set other as reference
      other.setRef(true);
    }
  protected:
    /** allocator of the column data set */
    Allocator allocator_;
    /** capacity of the columns of the container (for each column: number of
     *  available Rows without reallocation in this column)
    **/
    Array1D<int> availableRows_;
    /** range of the index of the columns of the container. **/
    Array1D<Range> rangeCols_;
    /** set the maximum possible number of columns without reallocation.
     *  @param capacity the maximum number of columns
     **/
    inline void setAvailableCols(int capacity = 0) { availableCols_ = capacity;}
    /** set the default capacity of a column.
     *  @param capacity the capacity to used
     **/
    inline void setCapacityByCols(int capacity) { capacityByCols_ = capacity;}
    /** @return the allocator. */
    inline Allocator& allocator() { return allocator_;}
    /** @return a pointer on the j-th column of the container
     *  @param j the index of the column
     **/
    inline PTRCOL& data(int j) { return allocator_.data(j);}
    /** move T to this
     *  @param T the container to move
     **/
    inline void move(Derived const& T)
    {
      allocator_.move(T.allocator_);
      // move this part
      availableRows_.move(T.availableRows_);
      rangeCols_.move(T.rangeCols_);
      availableCols_ = T.availableCols_;
      // Set IArrayBase part
      this->setCols(T.cols());
      this->setRows(T.rows());
    }
    /** copy the column pos2 of the container T to the column
     *  pos1 of this. One of the container (either this or T but not both)
     *  have to be a reference otherwise, user will experiment a memory leak.
     *
     *  @param T the container with the column to transfer
     *  @param pos1 index of the column to initialize
     *  @param pos2 the column in the container T to transfer in this
     **/
    template<class Other>
    void copyColumn( IArray2DBase<PTRCOL, Other> const& T, int pos1, int pos2)
    {
      // copy column pos2 of T in pos1 of this
      data(pos1) = T.data(pos2);
      // set availableRows_
      availableRows_[pos1] = T.availableRows()[pos2];
      // set rangeCols_
      rangeCols_[pos1] = T.rangeCols()[pos2];
    }
    /** Transfer the column pos2 of the container T to the column
     *  pos1 of this. Set the column pos2 in T to a default value.
     *  The column pos1 should not exists or should be deleted previously
     *  otherwise user will experiment a memory leak.
     *
     *  @param T the container with the column to transfer
     *  @param pos1 index of the column to initialize
     *  @param pos2 the column in the container T to transfer in this
     **/
    template<class Other>
    void transferColumn( IArray2DBase<PTRCOL, Other>& T, int pos1, int pos2)
    {
      // copy column pos2 of T in pos1 of this
      data(pos1) = T.data(pos2);
      // set availableRows_
      availableRows_[pos1] = T.availableRows_[pos2];
      // set rangeCols_
      rangeCols_[pos1] = T.rangeCols_[pos2];
      // set column of T to default
      T.setDefaultCol(pos2);
    }
    /** Method for memory allocation and initialization of the horizontal
     *  range of the container.
     *  The vertical range is not set in this method. If an
     *  error occur, we set the cols_ of the container to default.
     *  @param J horizontal range
     **/
    void mallocHo(Range const& J)
    {
      // compute the size necessary (can be 0)
      int size = Arrays::evalSizeCapacity(J.size());
      // try to allocate memory
      try
      {
        // resize availableRows_ and rangeCols
        availableRows_.resize(J);
        rangeCols_.resize(J);
        // allocate memory for the columns
        allocator_.malloc(Range(J.begin(), size));
        setAvailableCols(size);
      }
      catch (runtime_error & error)   // if an error occur
      {
        // set default capacity (0)
        setAvailableCols();
        // set default range
        Base2D::setCols();
        // clear this->availableRows_
        availableRows_.clear();
        // clear this->rangeCols_
        rangeCols_.clear();
        // throw the error
        throw error;
      }
    }
    /** Method for memory reallocation and initialization of the horizontal
     *  range of the container.
     *  The vertical range is not set in this method. If an
     *  error occur, we set the cols_ of the container to default.
     *  @param J horizontal range
     **/
    void reallocCols(Range const& J)
    {
      // compute the size necessary (can be 0)
      int size = Arrays::evalSizeCapacity(J.size());
      // try to allocate memory
      try
      {
        // allocate memory for the columns
        allocator_.realloc(Range(J.begin(), size));
        // initialize this->availableRows_
        availableRows_.resize(J);
        // initialize this->rangeCols_
        rangeCols_.resize(J);
      }
      catch (runtime_error & error)   // if an error occur
      {
        // set default capacity (0)
        setAvailableCols();
        // set default range
        Base2D::setCols();
        // clear this->availableRows_
        availableRows_.clear();
        // clear this->rangeCols_
        rangeCols_.clear();
        // throw the error
        throw error;
      }
      // set new capacity if no error occur
      setAvailableCols(size);
    }
    /** Horizontal Memory deallocation.
     *  This method clear all allocated memory. The range of the columns
     *  is set to (firstCol_:firstCol_-1). The range of the Rows remain
     *  unmodified. If there is allocated memory for the columns, it
     *  should be liberated prior to this method.
     **/
    void freeRows()
    {
      // Nothing to do for reference
      if (this->isRef()) return;
      // free memory allocated in AllocatorBase
      allocator_.free();
      // set capacity size to default
      setAvailableCols(0);
      // set range of the columns to default
      Base2D::setCols(Range(this->beginCols(), 0));
      // clear arrays
      availableRows_.clear();
      rangeCols_.clear();
    }
  private:
    /** Horizontal capacity of the container (number of available
     *  columns without reallocation)
     **/
    int availableCols_;
    /** default capacity of a column */
    int capacityByCols_;
    /** set the default parameters and dimension to a column of the container.
     *  @param col the position of the column to initialize to a default value.
     *  @note if data is allocated, it will be lost
     **/
    inline void setDefaultCol(int col)
    {
      // set column of T to default
      data(col) = 0;
      // set availableRows_
      availableRows_[col] = 0;
      // set rangeCols_
      rangeCols_[col] = Range();
    }
};

} // namespace STK

#endif
// STK_ITARRAY2DBASE_H
