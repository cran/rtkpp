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
 * created on: 17 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_DiagonalOperator.h
 *  @brief In this file we implement the DiagonalOperator and DiagonalizeOperator classes.
 **/

#ifndef STK_DIAGOPERATOR_H
#define STK_DIAGOPERATOR_H


#include "Sdk/include/STK_StaticAssert.h"
#include "STK_SlicingOperators.h"

namespace STK
{

// forward declaration
template< typename Array> class DiagonalizeOperator;

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the diagonalize operator
 */
template<typename Lhs>
struct Traits< DiagonalizeOperator <Lhs> >
{
  enum
  {
    structure_ = Arrays::diagonal_,
    orient_    = Lhs::orient_,
    sizeRows_  = ( (Lhs::sizeRows_ != UnknownSize) && (Lhs::structure_!= (int)Arrays::point_) )
                 ?  Lhs::sizeRows_  : UnknownSize,
    sizeCols_  = ( (Lhs::sizeCols_ != UnknownSize) && (Lhs::structure_!= (int)Arrays::vector_) )
                 ?  Lhs::sizeCols_  : UnknownSize,
    // this is safe as we can use diagonalize operator only on 1D container
    size_      = (sizeRows_ != UnknownSize) ? sizeRows_ : sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<DiagonalizeOperator < Lhs> > Row;
  typedef ColOperator<DiagonalizeOperator < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ReturnType ReturnType;
};
} // end namespace hidden

/** @ingroup Arrays
 *  @class DiagonalizeOperator
  *
  * @brief Generic expression when a vector expression is "diagonalized".
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * diagonalize operator.
  *
  * This class represents an expression where a diagonalize operator is applied to
  * a vector expression. It is the return type of the diagonalize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name DiagonalizeOperator type explicitly.
  */
template< typename Lhs>
class DiagonalizeOperator: public ExprBase< DiagonalizeOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< DiagonalizeOperator< Lhs> > Base;
    typedef typename hidden::Traits< DiagonalizeOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< DiagonalizeOperator<Lhs> >::ReturnType ReturnType;

    typedef typename hidden::Traits< DiagonalizeOperator<Lhs> >::Row Row;
    typedef typename hidden::Traits< DiagonalizeOperator<Lhs> >::Col Col;
    enum
    {
        structure_ = hidden::Traits< DiagonalizeOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< DiagonalizeOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< DiagonalizeOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< DiagonalizeOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< DiagonalizeOperator<Lhs> >::storage_,
        size_      = hidden::Traits< DiagonalizeOperator<Lhs> >::size_
    };
    /** Type of the Range for the rows */
    typedef TRange<size_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<size_> ColRange;
    /** Constructor */
    inline DiagonalizeOperator( Lhs const& lhs)
                       : Base(), lhs_(lhs)
                       , rows_(lhs_.range())
                       , cols_(lhs_.range())
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Lhs);
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the range of the Columns */
    inline ColRange const&colsImpl() const { return cols_;}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }

    /** @return the element (i,j) of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ReturnType elt2Impl(int i, int j) const { return (lhs_.elt(i, j));}
    /** @return the element ith element of the expression
     *  @param i index of the ith element
     **/
    inline ReturnType elt1Impl(int i) const { return (lhs_.elt(i));}
    /** accesses to the element of the expression */
    inline ReturnType elt0Impl() const { return (lhs_.elt());}

  protected:
    Lhs const& lhs_;
    RowRange rows_;
    ColRange cols_;
};


// forward declaration
template< typename Array> class DiagonalOperator;

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the diag operator
 */
template<typename Lhs>
struct Traits< DiagonalOperator <Lhs> >
{
  enum
  {
    structure_ = Arrays::diagonal_,
    orient_    = Lhs::orient_,
    sizeRows_  = ((Lhs::sizeRows_ < Lhs::sizeCols_)) ?  Lhs::sizeRows_ : Lhs::sizeCols_,
    sizeCols_  = sizeRows_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<DiagonalOperator < Lhs> > Row;
  typedef ColOperator<DiagonalOperator < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ReturnType ReturnType;
};

} // end namespace hidden

/** @ingroup Arrays
 *  @class DiagonalOperator
  *
  * @brief Generic expression when we get the .
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * diagonal operator.
  *
  * This class represents an expression where a diagonal operator is applied to
  * an expression. It is the return type of the diagonal operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name DiagonalOperator type explicitly.
  */
template< typename Lhs>
class DiagonalOperator: public ExprBase< DiagonalOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< DiagonalOperator< Lhs> > Base;
    typedef typename hidden::Traits< DiagonalOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< DiagonalOperator<Lhs> >::ReturnType ReturnType;

    typedef typename hidden::Traits< DiagonalOperator<Lhs> >::Row Row;
    typedef typename hidden::Traits< DiagonalOperator<Lhs> >::Col Col;
    enum
    {
        structure_ = hidden::Traits< DiagonalOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< DiagonalOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< DiagonalOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< DiagonalOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< DiagonalOperator<Lhs> >::storage_,
        // get size of the
        size_      = (sizeCols_ < sizeRows_) ? sizeCols_ : sizeRows_
    };
    /** Type of the diagonal Range */
    typedef TRange<size_> DiagRange;

    /** Constructor */
    inline DiagonalOperator( Lhs const& lhs)
                           : Base(), lhs_(lhs)
                           , range_( lhs_.beginRows(), (size_ != UnknownSize) ? size_ : lhs_.sizeRows())
    {
      if (lhs.rows()!=lhs.cols())
        STKRUNTIME_ERROR_NO_ARG(DiagonalOperatorBase,lhs.rows()!=lhs.cols());
    }
    /**  @return the range of the rows */
    inline DiagRange const& rowsImpl() const { return range_;}
    /** @return the range of the Columns */
    inline DiagRange const& colsImpl() const { return range_;}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }

    /** @return the element (i,j) of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ReturnType elt2Impl(int i, int j) const { return (this->asDerived().lhs().elt(i, j));}
    /** @return the element ith element of the transposed expression
     *  @param i index of the ith element
     **/
    inline ReturnType elt1Impl(int i) const { return (this->asDerived().lhs().elt(i,i));}
    /** accesses to the element of the transposed expression */
    inline ReturnType elt0Impl() const { return (this->asDerived().lhs().elt());}

  protected:
    Lhs const& lhs_;
    DiagRange range_;
};

} // namespace STK

#endif /* STK_DIAGOPERATOR_H */
