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
 * Project:  stkpp::Arrays
 * created on: 21 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_SlicingOperators.h
 *  @brief In this file we implement the RowOperator, ColOperator and SubOperator classes.
 **/

#ifndef STK_SLICINGOPERATORS_H
#define STK_SLICINGOPERATORS_H

namespace STK
{

// forward declaration
template< typename Derived> class SlicingOperatorBase;
template< typename Array> class RowOperator;
template< typename Array> class ColOperator;
// only for vectors, points and the like
template< typename Array, int structure_ = hidden::Traits<Array>::structure_>
class SubOperator;

namespace hidden
{
/** @ingroup hidden
  * Helper Trait class for column operator.
  **/
template<int Structure_> struct ColumnTraits
{ enum { structure_ = Arrays::vector_ }; };
/** specialization for point_ */
template<> struct ColumnTraits<Arrays::point_>
{ enum { structure_ = Arrays::number_}; };
/** specialization for number_ */
template<> struct ColumnTraits<Arrays::number_>
{ enum { structure_ = Arrays::number_}; };

/** @ingroup hidden
  * Helper Trait class for row operator.
  **/
template<int Structure_> struct RowTraits
{ enum { structure_ = Arrays::point_ }; };
/** specialization for vector_ */
template<> struct RowTraits<Arrays::vector_>
{ enum { structure_ = Arrays::number_}; };
/** specialization for number_ */
template<> struct RowTraits<Arrays::number_>
{ enum { structure_ = Arrays::number_}; };

} // end namespace hidden

// forward declaration
template< typename Lhs> class RowOperatorBase;
template< typename Lhs> class ColOperatorBase;


namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the row operator
 */
template<typename Lhs>
struct Traits< RowOperator <Lhs> >
{
  typedef typename Lhs::Type Type;
  enum
  {
    structure_ = RowTraits<Lhs::structure_>::structure_,
    orient_    = Lhs::orient_,
    sizeRows_  = 1,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
};

} // namespace hidden
/** @ingroup Arrays
  * @class RowOperator
  *
  * \brief Generic expression when the row of an expression is accessed
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * row operator.
  *
  * This class represents an expression where a row operator is applied to
  * an expression. It is the return type of the row operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name RowOperator type explicitly.
  */
template< typename Lhs>
class RowOperator: public RowOperatorBase< Lhs>, public TRef<1>
{
  public:
    typedef RowOperatorBase< Lhs> Base;
    typedef typename hidden::Traits< RowOperator<Lhs> >::Type Type;
    enum
    {
        structure_ = hidden::Traits< RowOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< RowOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< RowOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< RowOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< RowOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline RowOperator( Lhs const& lhs, int i)
                      : Base(i), lhs_(lhs), rows_(i,1), cols_(lhs_.rangeColsInRow(i)) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the first index of the rows */
    inline int const beginRowsImpl() const { return rows_.begin();}
    /** @return the ending index of the rows */
    inline int const endRowsImpl() const { return rows_.end();}
    /** @return the number of rows */
    inline int const sizeRowsImpl() const { return 1;}

    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return cols_;}
    /** @return the first index of the columns */
    inline int const beginColsImpl() const { return cols_.begin();}
    /** @return the ending index of the columns */
    inline int const endColsImpl() const { return cols_.end();}
    /** @return the number of columns */
    inline int const sizeColsImpl() const { return cols_.size();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }

  protected:
    Lhs const& lhs_;
    RowRange rows_;
    ColRange cols_;
};

/** @ingroup Arrays
  * @brief implement the access to the elements in the general case.
  **/
template< typename Lhs>
class RowOperatorBase: public ExprBase< RowOperator< Lhs> >
{
  public:
    typedef typename hidden::Traits< RowOperator<Lhs> >::Type Type;
    typedef ExprBase< RowOperator< Lhs> > Base;
    /** constructor. */
    inline RowOperatorBase(int i) : Base(), i_(i) {}
    /** @return the element (i,j)
     *  @param i, j index of the row and column
     **/
    inline Type const elt2Impl(int i, int j) const
    {
#ifdef STK_DEBUG
      if (i != i_)
      { STKRUNTIME_ERROR_2ARG(RowOperatorBase::elt2Impl,i,j,row index is not valid);}
#endif
      return (this->asDerived().lhs().elt(i_, j));
    }
    /** @return the element jth element
     *  @param j index of the jth element
     **/
    inline Type const elt1Impl(int j) const
    { return (this->asDerived().lhs().elt(i_, j));}
    /** accesses to the element */
    inline Type const elt0Impl() const
    { return (this->asDerived().lhs().elt());}
  protected:
    int i_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the column operator
 */
template<typename Lhs>
struct Traits< ColOperator <Lhs> >
{
  typedef typename Lhs::Type Type;
  enum
  {
    structure_ = ColumnTraits<Lhs::structure_>::structure_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = 1,
    storage_   = Lhs::storage_
  };
};

} // namespace hidden

/** @ingroup Arrays
  * @class ColOperator
  *
  * @brief Generic expression when the column of an expression is accessed
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * column operator.
  *
  * This class represents an expression where a column operator is applied to
  * an expression. It is the return type of the column operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name ColOperator type explicitly.
  */
template< typename Lhs>
class ColOperator: public ColOperatorBase< Lhs>, public TRef<1>
{
  public:
    typedef ColOperatorBase< Lhs> Base;
    typedef typename hidden::Traits< ColOperator<Lhs> >::Type Type;
    enum
    {
        structure_ = hidden::Traits< ColOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< ColOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< ColOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< ColOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< ColOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline ColOperator( Lhs const& lhs, int j)
                      : Base(j), lhs_(lhs), rows_(lhs_.rangeRowsInCol(j)), cols_(j,1) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the first index of the rows */
    inline int const beginRowsImpl() const { return rows_.begin();}
    /** @return the ending index of the rows */
    inline int const endRowsImpl() const { return rows_.end();}
    /** @return the number of rows */
    inline int const sizeRowsImpl() const { return rows_.size();}

    /** @return the range of the columns */
    inline ColRange const& colsImpl() const { return cols_;}
    /** @return the first index of the columns */
    inline int const beginColsImpl() const { return cols_.begin();}
    /** @return the ending index of the columns */
    inline int const endColsImpl() const { return cols_.end();}
    /** @return the number of columns */
    inline int const sizeColsImpl() const { return 1;}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }

  protected:
    Lhs const& lhs_;
    RowRange rows_;
    ColRange cols_;
};

/** @ingroup Arrays
  * @brief implement the access to the elements in the general case.
  **/
template< typename Lhs>
class ColOperatorBase : public ExprBase< ColOperator< Lhs> >
{
  public:
    typedef typename hidden::Traits< ColOperator<Lhs> >::Type Type;
    typedef ExprBase< ColOperator< Lhs> > Base;
    /** constructor. */
    inline ColOperatorBase(int j) : Base(), j_(j) {}
    /** @return the element (i,j)
     *  @param i, j indexes of the row and column
     **/
    inline Type const elt2Impl(int i, int j) const
    {
#ifdef STK_DEBUG
      if (j != j_)
      { STKRUNTIME_ERROR_2ARG(ColOperatorBase::elt2Impl,i,j,column index is not valid);}
#endif
      return (this->asDerived().lhs().elt(i, j_));
    }
    /** @return the element ith element
     *  @param i index of the ith element
     **/
    inline Type const elt1Impl(int i) const
    { return (this->asDerived().lhs().elt(i, j_));}
    /** accesses to the element */
    inline Type const elt0Impl() const
    { return (this->asDerived().lhs().elt());}
  protected:
    int j_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the sub operator
 */
template<typename Lhs>
struct Traits< SubOperator <Lhs> >
{
  typedef typename Lhs::Type Type;
  enum
  {
    structure_ = Lhs::structure_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
};

} // namespace hidden

/** @ingroup Arrays
  * @class SubOperator
  *
  * @brief Generic expression when the sub-part of an expression is accessed
  * (specialization for vectors)
  *
  * @tparam Lhs the type of the column vector expression to which we are
  * applying the sub operator.
  *
  * This class represents an expression where a sub operator is applied to
  * a column-vector expression. It is the return type of the sub operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name SubOperator type explicitly.
  */
template< typename Lhs>
class SubOperator<Lhs, Arrays::vector_ > : public ExprBase< SubOperator< Lhs> >, public TRef<1>
{
  public:
    typedef typename hidden::Traits< SubOperator<Lhs> >::Type Type;
    typedef ExprBase< SubOperator< Lhs> > Base;
    enum
    {
      structure_ = hidden::Traits< SubOperator<Lhs> >::structure_,
      orient_    = hidden::Traits< SubOperator<Lhs> >::orient_,
      sizeRows_  = hidden::Traits< SubOperator<Lhs> >::sizeRows_,
      sizeCols_  = hidden::Traits< SubOperator<Lhs> >::sizeCols_,
      storage_   = hidden::Traits< SubOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<UnknownSize> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<1> ColRange;

    /** Constructor */
    inline SubOperator( Lhs const& lhs, RowRange const& I) : Base(), lhs_(lhs), rows_(I) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the first row index of the sub expression */
    inline int const beginRowsImpl() const { return rows_.begin();}
    /** @return the ending row index of the sub expression */
    inline int const endRowsImpl() const { return rows_.end();}
    /** @return the number of rows of the sub expression */
    inline int const sizeRowsImpl() const { return rows_.size();}

    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs_.cols();}
    /** @return the first column index of the sub expression */
    inline int const beginColsImpl() const { return lhs_.begin();}
    /** @return the ending column index of the sub expression */
    inline int const endColsImpl() const { return lhs_.end();}
    /** @return the number of columns of the sub expression */
    inline int const sizeColsImpl() const { return 1;}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the element ith element
     *  @param i index of the ith element
     **/
    inline Type const elt1Impl(int i) const
    { return (this->asDerived().lhs().elt1Impl(i));}
    /** accesses to the element */
    inline Type const elt0Impl() const
    { return (this->asDerived().lhs().elt());}

  protected:
    Lhs const& lhs_;
    RowRange rows_;
};

/** @ingroup Arrays
  * @class SubOperator
  *
  * @brief Generic expression when the sub-part of an expression is accessed
  * (specialization for points)
  *
  * @tparam Lhs the type of the row-vector expression to which we are applying
  * the sub operator.
  *
  * This class represents an expression where a sub operator is applied to
  * an row-vector expression. It is the return type of the sub operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name SubOperator type explicitly.
  */
template< typename Lhs>
class SubOperator<Lhs, Arrays::point_ > : public ExprBase< SubOperator< Lhs> >, public TRef<1>
{
  public:
    typedef typename hidden::Traits< SubOperator<Lhs> >::Type Type;
    typedef ExprBase< SubOperator< Lhs> > Base;
    enum
    {
      structure_ = hidden::Traits< SubOperator<Lhs> >::structure_,
      orient_    = hidden::Traits< SubOperator<Lhs> >::orient_,
      sizeRows_  = hidden::Traits< SubOperator<Lhs> >::sizeRows_,
      sizeCols_  = hidden::Traits< SubOperator<Lhs> >::sizeCols_,
      storage_   = hidden::Traits< SubOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<1> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<UnknownSize> ColRange;
    /** Constructor */
    inline SubOperator( Lhs const& lhs, ColRange const& J) : Base(), lhs_(lhs), cols_(J) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs_.rows();}
    /** @return the first row index of the sub expression */
    inline int const beginRowsImpl() const { return lhs_.begin();}
    /** @return the ending row index of the sub expression */
    inline int const endRowsImpl() const { return lhs_.end();}
    /** @return the number of rows */
    inline int const sizeRowsImpl() const { return 1;}

    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return cols_;}
    /** @return the first column index of the sub expression */
    inline int const beginColsImpl() const { return cols_.begin();}
    /** @return the ending column index of the sub expression */
    inline int const endColsImpl() const { return cols_.end();}
    /** @return the number of columns */
    inline int const sizeColsImpl() const { return cols_.size();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the element ith element
     *  @param i index of the ith element
     **/
    inline Type const elt1Impl(int i) const
    { return (this->asDerived().lhs().elt1Impl(i));}
    /** accesses to the element */
    inline Type const elt0Impl() const
    { return (this->asDerived().lhs().elt());}

  protected:
    Lhs const& lhs_;
    ColRange cols_;
};

} // namespace STK

#endif /* STK_SLICINGOPERATORS_H */
