/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2011  Serge Iovleff

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
 * created on: 26 nov. 2011
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ITContainer.h
 *  @brief This is an internal header file, included by other
 *  Containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array2D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_ITCONTAINER2D_H
#define STK_ITCONTAINER2D_H

#include "Sdk/include/STK_Macros.h"
#include "Sdk/include/STK_IRecursiveTemplate.h"
#include "Sdk/include/STK_StaticAssert.h"

#include "STK_Traits.h"
#include "STK_Arrays_Util.h"

namespace STK
{

namespace hidden
{
/** @ingroup hidden
 * classes inheriting NoAssignOperator should not generate a default operator=.
 **/
class NoAssignOperator
{ private:
    NoAssignOperator& operator=(const NoAssignOperator&);
};

} // hidden namespace

/** @ingroup Arrays
 * @brief Interface base class for 2D containers.
 *
 * Use the curious recursive template paradigm : the template
 * parameter @c Derived is the name of the class that
 * implements the interface ITContainer.
 * For example
 * @code
 * template<class Type>
 * class Derived : public ITContainer< Derived<Type> >
 * {...}
 * @endcode
 *
 * The ITContainerBase class is the base class for all two-dimensional containers.
 * A two-dimensional container is defined by an horizontal range of index
 * for the columns and a vertical range of index for the rows.
 *
 * The pseudo virtual function defined in this interface and to implement
 * in the derived classes have the following definitions for the dimensions:
 * @code
 * ColRange const& colsImpl() const;
 * int beginColsImpl() const;
 * int endColsImpl() const;
 * int sizeColsImpl() const;
 * RowRange const& rowsImpl() const;
 * int beginRowsImpl() const;
 * int endRowsImpl() const;
 * int sizeRowsImpl() const;
 * @endcode
 * and for the accessors the following definitions:
 * @code
 *   Type const elt2Impl( int i, int j) const;
 *   Type& elt2Impl( int i, int j);
 * @endcode
 * for the 2D arrays (@ref Arrays::array2D_, @ref Arrays::square_, @ref Arrays::upper_triangular_
 * and @ref Arrays::lower_triangular_),
 * @code
 *   Type const elt1Impl( int i) const;
 *   Type& elt1Impl( int i);
 * @endcode
 * for the diagonal arrays (@ref Arrays::diagonal_) and the 1D arrays
 * (@ref Arrays::vector_, @ref Arrays::point_)
 * and
 * @code
 *   Type const elt0Impl() const;
 *   Type& elt0Impl();
 * @endcode
 * for the @ref Arrays::number_ arrays.
 **/
template <class Derived>
class ITContainerBase: public IRecursiveTemplate<Derived>, hidden::NoAssignOperator
{
  public:
    typedef IRecursiveTemplate<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** Default constructor */
    inline ITContainerBase() : Base() {}
    /** destructor. **/
    inline ~ITContainerBase() {}

  public:
    /** @return the range of the columns */
    inline ColRange const& cols() const { return this->asDerived().colsImpl();};
    /**  @return the index of the first column */
    inline int beginCols() const { return this->asDerived().beginColsImpl();}
    /**  @return the ending index of the columns */
    inline int endCols() const { return this->asDerived().endColsImpl();}
    /** @return the Horizontal size (the number of column) */
    inline int sizeCols() const { return this->asDerived().sizeColsImpl();}

    /** @return the range of the rows */
    inline RowRange const& rows() const { return this->asDerived().rowsImpl();}
    /** @return the index of the first row*/
    inline int beginRows() const { return this->asDerived().beginRowsImpl();}
    /** @return the ending index of the rows*/
    inline int endRows() const { return this->asDerived().endRowsImpl();}
    /** @return the Vertical size (the number of rows) */
    inline int sizeRows() const { return this->asDerived().sizeRowsImpl();}

    // for backward compatibility
    /** @return the index of the first column */
    inline int firstIdxCols() const { return beginCols();}
    /** @return the index of the first row */
    inline int firstIdxRows() const { return beginRows();}
    /** @return the index of the last column */
    inline int lastIdxCols() const { return endCols()-1;}
    /** @return the index of the last row */
    inline int lastIdxRows() const { return endRows()-1;}

    /** @return @c true if the container is empty, @c false otherwise */
    inline bool empty() const { return (sizeCols()<=0 || sizeRows()<=0);}

    /** @return the size of the container (the number of rows by the number of columns */
    inline int sizeArray() const { return sizeRows()*sizeCols();}
    /** @return a constant reference on element (i,j) of the 2D container
     *  @param i,j indexes of the row and of the column
     **/
    inline Type const elt(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > i)
      { STKOUT_OF_RANGE_2ARG(ITContainerBase::elt, i, j, beginRows() > i);}
      if (this->endRows() <= i)
      { STKOUT_OF_RANGE_2ARG(ITContainerBase::elt, i, j, endRows() <= i);}
      if (this->beginCols() > j)
      { STKOUT_OF_RANGE_2ARG(ITContainerBase::elt, i, j, beginCols() > j);}
      if (this->endCols() <= j)
      { STKOUT_OF_RANGE_2ARG(ITContainerBase::elt, i, j, endCols() <= j);}
#endif
      return this->asDerived().elt2Impl(i,j);
    }
    /** @return a constant reference on the element (i,j) of the 2D container.
     *  @param i,j indexes of the row and column
     **/
    inline Type const operator()(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > i)
      { STKOUT_OF_RANGE_2ARG(ITContainerBase::elt, i, j, beginRows() > i);}
      if (this->endRows() <= i)
      { STKOUT_OF_RANGE_2ARG(ITContainerBase::elt, i, j, endRows() <= i);}
      if (this->beginCols() > j)
      { STKOUT_OF_RANGE_2ARG(ITContainerBase::elt, i, j, beginCols() > j);}
      if (this->endCols() <= j)
      { STKOUT_OF_RANGE_2ARG(ITContainerBase::elt, i, j, endCols() <= j);}
#endif
      return this->asDerived().elt2Impl(i,j);
    }
    /** @return safely the constant element (i, j).
     *  @param i,j indexes of the row and column
     **/
    Type const at(int i, int j) const
    {
      if (this->beginRows() > i)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, beginRows() > i);}
      if (this->lastIdxRows() < i)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, lastIdxRows() < i);}
      if (this->beginCols() > j)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, beginCols() > j);}
      if (this->lastIdxCols() < j)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, lastIdxCols() < j);}
      return this->asDerived().elt2Impl(i, j);
    }
    /** @return the constant ith element
     *  @param i index of the ith element
     **/
    inline Type const elt(int i) const { return this->asDerived().elt1Impl(i);}
    /** @return the ith element
     *  @param i index of the ith element
     **/
    inline Type const operator[](int i) const { return this->asDerived().elt1Impl(i);}
    /** @return safely the constant ith element
     *  @param i index of the element
     **/
    Type const at(int i) const
    {
      if (this->asDerived().begin() > i)
      { STKOUT_OF_RANGE_1ARG(ITContainer::at, i, begin() > i);}
      if (this->asDerived().end() <= i)
      { STKOUT_OF_RANGE_1ARG(ITContainer::at, i, end() <= i);}
      return this->asDerived().elt1Impl(i);
    }
    /** @return a constant reference on the number */
    inline Type const elt() const { return this->asDerived().elt0Impl();}
    /** @return a constant reference on the number */
    inline Type const operator()() const { return this->asDerived().elt0Impl();}
};

/** @ingroup Arrays
 *  @brief Specialized interface class for homogeneous arrays that can be
 *  either a 2D arrays and 1D arrays.
 **/
template < class Derived, int Structure_ = hidden::Traits<Derived>::structure_ >
class ITContainer;

/** @ingroup Arrays
 *  Specialization for array2D_ */
template <class Derived>
class ITContainer<Derived, Arrays::array2D_> : public ITContainerBase<Derived>
{
  public:
    typedef ITContainerBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };

    /** compute the range of the stored elements in a column. */
    inline Range rangeRowsInCol(int) const { return this->rows();}
    /** compute the range of the stored elements in a row. */
    inline Range rangeColsInRow(int) const { return this->cols();}

  protected:
    /** default constructor. */
    inline ITContainer() : Base() {}
    /** destructor. */
    inline ~ITContainer() {}
};

/** @ingroup Arrays Specialization for square_ */
template <class Derived>
class ITContainer<Derived, Arrays::square_> : public ITContainerBase<Derived>
{
  public:
    typedef ITContainerBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** default constructor. */
    inline ITContainer() : Base() {}
    /** destructor. */
    inline ~ITContainer() {}

  public:
    /** @return the range of the square container. */
    inline RowRange const& range() const { return this->asDerived().rows();}
    /** @return the first index of the square container. */
    inline int begin() const { return this->asDerived().beginRowsImpl();}
    /** @return the ending index of the square container. */
    inline int end() const { return this->asDerived().endRowsImpl();}
    /** @return the size of the rows and columns of the container. */
    inline int size() const { return this->asDerived().sizeRowsImpl();}

    // for backward compatibility
    /** @return the first index of the square container */
    inline int firstIdx() const { return this->beginRows();}
    /** @return the last index of the square container. */
    inline int lastIdx() const { return this->endRows()-1;}

    /** compute the range of the stored elements in a column. */
    inline Range rangeRowsInCol(int) const { return this->asDerived().rows();}
    /** compute the range of the stored elements in a row. */
    inline Range rangeColsInRow(int) const { return this->asDerived().rows();}

    Type const trace() const
    {
      Type sum = 0.0;
      for (int k = begin(); k< end(); k++) sum += this->elt(k, k);
      return sum;
    }
};

/** @ingroup Arrays
 *  Specialization for lower_triangular_ */
template <class Derived>
class ITContainer<Derived, Arrays::lower_triangular_> : public ITContainerBase<Derived>
{
  public:
    typedef ITContainerBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;

  protected:
    /** default constructor. */
    inline ITContainer() : Base() {}
    /** destructor. */
    inline ~ITContainer() {}

  public:
    /** @return safely the constant element (i, j).
     *  @param i,j indexes of the row and column
     **/
    Type const at(int i, int j) const
    {
      if (this->beginRows() > i)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, beginRows() > i);}
      if (this->lastIdxRows() < i)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, lastIdxRows() < i);}
      if (this->beginCols() > j)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, beginCols() > j);}
      if (this->lastIdxCols() < j)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, lastIdxCols() < j);}
      return (j>i) ? Type() : this->asDerived().elt2Impl(i, j);
    }
    /** @return the Range of the effectively stored elements in the column @c icol.
     *  @param icol the number of the column to compute the range
     **/
    inline Range rangeRowsInCol( int const& icol) const
    { return Range(icol, this->lastIdxRows(), 0);}
    /** compute the range of the effectively stored elements in the row @c irow.
     *  @param irow the index of the row
     **/
    inline Range rangeColsInRow( int const& irow) const
    { return Range(this->beginCols(), std::min(irow, this->lastIdxCols()), 0);}
};

/** @ingroup Arrays
 *  Specialization for upper_triangular_ */
template <class Derived>
class ITContainer<Derived, Arrays::upper_triangular_> : public ITContainerBase<Derived>
{
  public:
    typedef ITContainerBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;

  protected:
    /** default constructor. */
    inline ITContainer() : Base() {}
    /** destructor. */
    inline ~ITContainer() {}

  public:
    /** @return safely the constant element (i, j).
     *  @param i, j indexes of the row and of the column
     **/
    Type const at(int i, int j) const
    {
      if (this->beginRows() > i)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, beginRows() > i);}
      if (this->lastIdxRows() < i)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, lastIdxRows() < i);}
      if (this->beginCols() > j)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, beginCols() > j);}
      if (this->lastIdxCols() < j)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, lastIdxCols() < j);}
      return (i>j) ? Type() : this->asDerived().elt2Impl(i, j);
    }
    /** @return the range of the effectively stored elements in the row @c irow.
     *  @param irow the index of the row we want to compute the range
     **/
    inline Range rangeColsInRow( int irow) const
    { return Range(irow, this->lastIdxCols(), 0);}
    /** compute the range of the effectively stored elements in the col
     *  icol.
     *  @param icol the number of the column to compute the range
     **/
    inline Range rangeRowsInCol( int icol) const
    { return Range(this->beginRows(), std::min(icol, this->lastIdxRows()), 0);}
};

/** @ingroup Arrays
 *  Specialization for diagonal_ */
template <class Derived>
class ITContainer<Derived, Arrays::diagonal_> : public ITContainerBase<Derived>
{
  public:
    typedef ITContainerBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** default constructor. */
    inline ITContainer() : Base() {}
    /** destructor. */
    inline ~ITContainer() {}

  public:
    /** @return the range of the diagonal container. */
    inline RowRange const& range() const { return this->asDerived().rows();}
    /** @return the first index of the diagonal container. */
    inline int begin() const { return this->asDerived().beginRowsImpl();}
    /** @return the ending index of the diagonal container. */
    inline int end() const { return this->asDerived().endRowsImpl();}
    /** @return the size of the rows and columns of the container. */
    inline int size() const { return this->asDerived().sizeRowsImpl();}

    // for backward compatibility
    /** @return the first index of the diagonal container */
    inline int firstIdx() const { return this->beginRows();}
    /** @return the last index of the diagonal container. */
    inline int lastIdx() const { return this->endRows()-1;}

    /** @return the first diagonal element */
    inline Type const front() const { return this->asDerived().elt(begin());}
    /** @return the last diagonal element */
    inline Type const back() const { return this->asDerived().elt(lastIdx());}
    /** @return safely the constant element (i, j).
     *  @param i index of the row
     *  @param j index of the col
     **/
    Type const at(int i, int j) const
    {
      if (this->begin() > i)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, begin() > i);}
      if (this->lastIdx() < i)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, lastIdx() < i);}
      if (this->begin() > j)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, begin() > j);}
      if (this->lastIdx() < j)
      { STKOUT_OF_RANGE_2ARG(ITContainer::at, i, j, lastIdx() < j);}
      return (i==j) ? this->asDerived().elt2Impl(i, i) : Type();
    }
    /** @return safely the constant ith diagonal element.
     *  @param i index of the diagonal element
     **/
    inline Type const at(int i) const { return at(i, i);}
    /** @return the Range of the column pos. */
    inline Range rangeRowsInCol(int pos) const { return Range(pos,1);}
    /** @return the Range of the row pos. */
    inline Range rangeColsInRow(int pos) const { return Range(pos,1);}
};

/** @ingroup Arrays
 *  Specialization for vector_ */
template <class Derived>
class ITContainer<Derived, Arrays::vector_> : public ITContainerBase<Derived>
{
  public:
    typedef ITContainerBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** default constructor. */
    inline ITContainer() : Base() {}
    /** destructor. */
    inline ~ITContainer() {}

  public:
    /** @return the range of the container */
    inline RowRange const& range() const  { return this->asDerived().rows();}
    /** @return the index of the first element */
    inline int begin() const { return this->asDerived().beginRowsImpl();}
    /**  @return the ending index of the elements */
    inline int end() const  { return this->asDerived().endRowsImpl();}
    /**  @return the size of the vector */
    inline int size() const  { return this->asDerived().sizeRowsImpl();}

    // for backward compatibility
    /** @return the first index of the vector */
    inline int firstIdx() const { return this->beginRows();}
    /** @return the last index of the vector */
    inline int lastIdx() const { return this->endRows()-1;}

    /** @return the first element */
    inline Type const front() const { return this->asDerived().elt(begin());}
    /** @return the last element */
    inline Type const back() const { return this->asDerived().elt(lastIdx());}

    /** @return the index of the column of the vector */
    inline int colIdx() const { return this->beginCols();}

    /** @return the range of the effectively stored elements in the column. */
    inline Range rangeRowsInCol(int) const { return this->rows();}
    /** @return the range of the effectively stored elements in the row. */
    inline Range rangeColsInRow(int) const { return this->cols();}
};

/** @ingroup Arrays
 *  Specialization for point_ */
template <class Derived>
class ITContainer<Derived, Arrays::point_> : public ITContainerBase<Derived>
{
  public:
    typedef ITContainerBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** default constructor. */
    inline ITContainer() : Base() {}
    /** destructor. */
    inline ~ITContainer() {}

  public:
    /** @return the range of the container */
    inline ColRange const& range() const  { return this->asDerived().cols();}
    /** @return the index of the first element */
    inline int begin() const { return this->asDerived().beginColsImpl();}
    /**  @return the ending index of the elements */
    inline int end() const  { return this->asDerived().endColsImpl();}
    /**  @return the size of the container */
    inline int size() const  { return this->asDerived().sizeColsImpl();}

    /** @return the index of the row of the point */
    inline int rowIdx() const { return this->beginRows();}

    // for backward compatibility
    /** @return the first index of the point */
    inline int firstIdx() const { return this->beginCols();}
    /** @return the last index of the vector */
    inline int lastIdx() const { return this->endCols()-1;}

    /** @return the first element */
    inline Type const front() const { return this->asDerived().elt(begin());}
    /** @return the last element */
    inline Type const back() const { return this->asDerived().elt(lastIdx());}

    /** @return the range of the effectively stored elements in the column. */
    inline Range rangeRowsInCol(int) const { return this->rows();}
    /** @return the range of the effectively stored elements in the row. */
    inline Range rangeColsInRow(int) const { return this->cols();}
};

/** @ingroup Arrays
 *  Specialization for number_ */
template <class Derived>
class ITContainer<Derived, Arrays::number_> : public ITContainerBase<Derived>
{
  public:
    typedef ITContainerBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;

  protected:
    /** default constructor. */
    inline ITContainer() : Base() {}
    /** destructor. */
    inline ~ITContainer() {}

  public:
    /** @return the range of the effectively stored elements in the column. */
    inline Range rangeRowsInCol(int) const { return this->rows();}
    /** @return the range of the effectively stored elements in the row. */
    inline Range rangeColsInRow(int) const { return this->cols();}
    /** Conversion to scalar */
    inline operator Type const() const {return this->asDerived().elt0Impl();}
};

} // namespace STK

#endif
// STK_ITCONTAINER_H
