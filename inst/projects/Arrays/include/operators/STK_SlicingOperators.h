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

/** @ingroup hidden
 *  @brief Traits class for the row operator
 */
template<typename Lhs>
struct Traits< RowOperator <Lhs> >
{
  enum
  {
    structure_ = RowTraits<Lhs::structure_>::structure_,
    orient_    = Lhs::orient_,
    sizeRows_  = 1,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ReturnType ReturnType;
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
class RowOperator: public ExprBase< RowOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< RowOperator< Lhs> > Base;
    typedef typename hidden::Traits< RowOperator< Lhs> >::Type Type;
    typedef typename hidden::Traits< RowOperator< Lhs> >::ReturnType ReturnType;
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
                      : Base(), i_(i), lhs_(lhs), rows_(i,1), cols_(lhs_.rangeColsInRow(i)) {}
    /** Copy constructor */
    inline RowOperator( RowOperator const& row, bool ref)
                      : Base(row), i_(row.i_), lhs_(row.lhs_), rows_(row.rows_), cols_(row.cols_) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return cols_;}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the element (i,j)
     *  @param i, j index of the row and column
     **/
    inline ReturnType elt2Impl(int i, int j) const
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
    inline ReturnType elt1Impl(int j) const { return (this->asDerived().lhs().elt(i_, j));}
    /** accesses to the element */
    inline ReturnType elt0Impl() const { return (this->asDerived().lhs().elt());}

  protected:
    int i_;
    Lhs const& lhs_;
    RowRange rows_;
    ColRange cols_;
};


namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the column operator
 */
template<typename Lhs>
struct Traits< ColOperator <Lhs> >
{
  enum
  {
    structure_ = ColumnTraits<Lhs::structure_>::structure_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = 1,
    storage_   = Lhs::storage_
  };
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ReturnType ReturnType;
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
class ColOperator: public ExprBase< ColOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< ColOperator< Lhs> > Base;
    typedef typename hidden::Traits< ColOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< ColOperator<Lhs> >::ReturnType ReturnType;
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
                      : Base(), j_(j), lhs_(lhs), rows_(lhs_.rangeRowsInCol(j)), cols_(j,1) {}
    /** Copy constructor */
    inline ColOperator( ColOperator const& col, bool ref)
                      : Base(col), j_(col.j_), lhs_(col.lhs_), rows_(col.rows_), cols_(col.cols_) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the columns range */
    inline ColRange const& colsImpl() const { return cols_;}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the ith element
     *  @param i,j indexes of the element
     **/
    inline ReturnType elt2Impl(int i, int j) const
    {
#ifdef STK_DEBUG
      if (j != j_)
      { STKRUNTIME_ERROR_2ARG(ColOperatorBase::elt2Impl,i,j,column index is not valid);}
#endif
      return (lhs_.elt(i, j_));
    }
    /** @return the element ith element
     *  @param i index of the ith element
     **/
    inline ReturnType elt1Impl(int i) const { return (lhs_.elt(i, j_));}
    /** access to the element */
    inline ReturnType elt0Impl() const { return (lhs_.elt(j_));}

  protected:
    int j_;
    Lhs const& lhs_;
    RowRange rows_;
    ColRange cols_;
};


namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the sub operator
 */
template<typename Lhs>
struct Traits< SubOperator <Lhs> >
{
  enum
  {
    structure_ = Lhs::structure_,
    orient_    = Lhs::orient_,
    sizeRows_  = (orient_ == Arrays::point_) ? 1 : UnknownSize,
    sizeCols_  = (orient_ == Arrays::vector_) ? 1 : UnknownSize,
    storage_   = Lhs::storage_
  };
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ReturnType ReturnType;
};

} // namespace hidden

/** @ingroup Arrays
  * @class SubOperator
  *
  * @brief Generic expression when the sub-part of an expression is accessed
  * (specialization for vectors)
  *
  * @tparam Lhs the type of the array or expression to which we are
  * applying the sub operator.
  *
  * This class represents an expression where a sub operator is applied to
  * an array expression. It is the return type of the sub operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name SubOperator type explicitly.
  */
template< typename Lhs, int Orient_>
class SubOperator: public ExprBase< SubOperator<Lhs, Orient_> >, public TRef<1>
{
  public:
    typedef ExprBase< SubOperator<Lhs, Orient_> > Base;
    typedef typename hidden::Traits< SubOperator< Lhs, Orient_> >::Type Type;
    typedef typename hidden::Traits< SubOperator< Lhs, Orient_> >::ReturnType ReturnType;
    enum
    {
      structure_ = hidden::Traits< SubOperator<Lhs, Arrays::vector_> >::structure_,
      orient_    = hidden::Traits< SubOperator<Lhs, Arrays::vector_> >::orient_,
      sizeRows_  = hidden::Traits< SubOperator<Lhs, Arrays::vector_> >::sizeRows_,
      sizeCols_  = hidden::Traits< SubOperator<Lhs, Arrays::vector_> >::sizeCols_,
      storage_   = hidden::Traits< SubOperator<Lhs, Arrays::vector_> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<UnknownSize> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<UnknownSize> ColRange;

    /** Constructor */
    inline SubOperator( Lhs const& lhs, Range const& I, Range const& J)
                      : Base(), lhs_(lhs), rows_(I), cols_(J) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the range of the Columns */
    inline ColRange const&colsImpl() const { return cols_;}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the element ith element
     *  @param i index of the ith element
     **/
    inline ReturnType elt1Impl(int i) const { return (lhs_.elt1Impl(i));}
    /** accesses to the element */
    inline ReturnType elt0Impl() const { return (lhs_.elt());}

  protected:
    Lhs const& lhs_;
    RowRange rows_;
    ColRange cols_;
};

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
class SubOperator<Lhs, Arrays::vector_ >: public ExprBase< SubOperator<Lhs, Arrays::vector_> >
                                        , public TRef<1>
{
  public:
    typedef ExprBase< SubOperator<Lhs, Arrays::vector_> > Base;
    typedef typename hidden::Traits< SubOperator< Lhs, Arrays::vector_> >::Type Type;
    typedef typename hidden::Traits< SubOperator< Lhs, Arrays::vector_> >::ReturnType ReturnType;
    enum
    {
      structure_ = hidden::Traits< SubOperator<Lhs, Arrays::vector_> >::structure_,
      orient_    = hidden::Traits< SubOperator<Lhs, Arrays::vector_> >::orient_,
      sizeRows_  = hidden::Traits< SubOperator<Lhs, Arrays::vector_> >::sizeRows_,
      sizeCols_  = hidden::Traits< SubOperator<Lhs, Arrays::vector_> >::sizeCols_,
      storage_   = hidden::Traits< SubOperator<Lhs, Arrays::vector_> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<UnknownSize> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<1> ColRange;

    /** Constructor */
    inline SubOperator( Lhs const& lhs, Range const& I) : Base(), lhs_(lhs), rows_(I) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs_.cols();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the element ith element
     *  @param i index of the ith element
     **/
    inline ReturnType elt1Impl(int i) const { return (lhs_.elt1Impl(i));}
    /** accesses to the element */
    inline ReturnType elt0Impl() const { return (lhs_.elt());}

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
class SubOperator<Lhs, Arrays::point_ > : public ExprBase< SubOperator<Lhs, Arrays::point_> >
                                        , public TRef<1>
{
  public:
    typedef ExprBase< SubOperator<Lhs, Arrays::point_> > Base;
    typedef typename hidden::Traits< SubOperator<Lhs, Arrays::point_> >::Type Type;
    typedef typename hidden::Traits< SubOperator<Lhs, Arrays::point_> >::ReturnType ReturnType;
    enum
    {
      structure_ = hidden::Traits< SubOperator<Lhs, Arrays::point_> >::structure_,
      orient_    = hidden::Traits< SubOperator<Lhs, Arrays::point_> >::orient_,
      sizeRows_  = hidden::Traits< SubOperator<Lhs, Arrays::point_> >::sizeRows_,
      sizeCols_  = hidden::Traits< SubOperator<Lhs, Arrays::point_> >::sizeCols_,
      storage_   = hidden::Traits< SubOperator<Lhs, Arrays::point_> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<1> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<UnknownSize> ColRange;
    /** Constructor */
    inline SubOperator( Lhs const& lhs, Range const& J) : Base(), lhs_(lhs), cols_(J) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs_.rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return cols_;}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the element ith element
     *  @param i index of the ith element
     **/
    inline ReturnType elt1Impl(int i) const
    { return (this->asDerived().lhs().elt1Impl(i));}
    /** accesses to the element */
    inline ReturnType elt0Impl() const
    { return (this->asDerived().lhs().elt());}

  protected:
    Lhs const& lhs_;
    ColRange cols_;
};

} // namespace STK

#endif /* STK_SLICINGOPERATORS_H */
