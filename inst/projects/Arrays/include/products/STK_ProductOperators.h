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
 * created on: 20 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ProductOperators.h
 *  @brief In this file we implement the ProductOperator class.
 **/

#ifndef STK_PRODUCTOPERATOR_H
#define STK_PRODUCTOPERATOR_H


#define EGAL(arg1, arg2) ((arg1::structure_ == int(Arrays::arg2)))

#include "../STK_CAllocator.h"
#include "STK_ProductImpl.h"

namespace STK
{

// forward declarations (needed because there is no CAllocator with this structure)
template<typename> class Array2DUpperTriangular;
template<typename> class Array2DLowerTriangular;

namespace hidden
{

/** @ingroup hidden
 *  @brief Traits class to get the correct returned Structure, Type, allocator,...
 *  of operator*. This Traits class is used by the functors classes operating
 *  on the Array2D  classes.
 *  @note the impossible cases are tracked in ArrayByArrayProduct class.
 **/
template<typename Lhs, typename Rhs, int LStructure_, int RStructure_>
struct ProductTraits;

//------------------------------------------------------------
// general lhs case. result is general except if rhs is vector
/** specialization for lhs is array2D */
template<typename Lhs, typename Rhs, int RStructure_>
struct ProductTraits<Lhs, Rhs, Arrays::array2D_, RStructure_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_ = Traits<Lhs>::sizeRows_,
    sizeCols_ = Traits<Rhs>::sizeCols_,
    // if there is more block on the left side, we use bp --> result is by_col_
    orient_   =  ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : (Traits<Lhs>::sizeRows_)/blockSize < (Traits<Rhs>::sizeCols_)/blockSize
                   ? int(Arrays::by_row_) : int(Arrays::by_col_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
              ? int(Arrays::dense_) : int(Arrays::sparse_)
  };

  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ReturnType;
  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

//----------------------------------
// square lhs case. result is general except if rhs square
template<typename Lhs, typename Rhs, int RStructure_>
struct ProductTraits<Lhs, Rhs, Arrays::square_, RStructure_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_ = Traits<Lhs>::sizeRows_,
    sizeCols_ = Traits<Rhs>::sizeCols_,
    // if there is more block on the left side, we use bp --> result is by_col_
    orient_   =  ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : (Traits<Lhs>::sizeRows_)/blockSize < (Traits<Rhs>::sizeCols_)/blockSize
                   ? int(Arrays::by_row_) : int(Arrays::by_col_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
              ? int(Arrays::dense_) : int(Arrays::sparse_)
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ReturnType;

  typedef CAllocator<Type, sizeRows_,sizeCols_, orient_> Allocator;
};

template<typename Lhs, typename Rhs>
struct ProductTraits<Lhs, Rhs, Arrays::square_, Arrays::square_>
{
  enum
  {
    structure_ = Arrays::square_,
    sizeRows_ = (int(Traits<Lhs>::sizeRows_) < int(Traits<Rhs>::sizeCols_)) ? int(Traits<Lhs>::sizeRows_) : int(Traits<Rhs>::sizeCols_),
    sizeCols_ = (int(Traits<Lhs>::sizeCols_) < int(Traits<Rhs>::sizeRows_)) ? int(Traits<Lhs>::sizeCols_) : int(Traits<Rhs>::sizeRows_),
    // if there is more block on the left side, we use bp --> result is by_col_
    orient_   =  ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : (Traits<Lhs>::sizeRows_)/blockSize < (Traits<Rhs>::sizeCols_)/blockSize
                   ? int(Arrays::by_row_) : int(Arrays::by_col_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
              ? int(Arrays::dense_) : int(Arrays::sparse_)
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ReturnType;

  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

//----------------------------------
// lower triangular case
template<typename Lhs, typename Rhs, int RStructure_>
struct ProductTraits<Lhs, Rhs, Arrays::lower_triangular_, RStructure_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_ = Traits<Lhs>::sizeRows_,
    sizeCols_ = Traits<Rhs>::sizeCols_,
    // if there is more block on the left side, we use bp --> result is by_col_
    orient_   =  ( sizeRows_ == UnknownSize || sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : sizeRows_/blockSize < sizeCols_/blockSize
                   ? int(Arrays::by_row_) : int(Arrays::by_col_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
              ? int(Arrays::dense_) : int(Arrays::sparse_)
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ReturnType;

  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

template<typename Lhs, typename Rhs>
struct ProductTraits<Lhs, Rhs, Arrays::lower_triangular_, Arrays::lower_triangular_>
{
  enum
  {
    structure_ = Arrays::lower_triangular_,
    sizeRows_ = ((int)Traits<Lhs>::sizeRows_ < (int)Traits<Rhs>::sizeCols_) ? Traits<Lhs>::sizeRows_ : Traits<Rhs>::sizeCols_,
    sizeCols_ = ((int)Traits<Lhs>::sizeCols_ < (int)Traits<Rhs>::sizeRows_) ? Traits<Lhs>::sizeCols_ : Traits<Rhs>::sizeRows_,
    // if there is more block on the left side, we use bp --> result is by_col_
    orient_   =  ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : (Traits<Lhs>::sizeRows_)/blockSize < (Traits<Rhs>::sizeCols_)/blockSize
                   ? int(Arrays::by_row_) : int(Arrays::by_col_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
              ? int(Arrays::dense_) : int(Arrays::sparse_)
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ReturnType;

  typedef Array2DLowerTriangular<Type> Allocator; // no CAllocator
};

//----------------------------------
// upper triangular case
template<typename Lhs, typename Rhs, int RStructure_>
struct ProductTraits<Lhs, Rhs, Arrays::upper_triangular_, RStructure_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_ = Traits<Lhs>::sizeRows_,
    sizeCols_ = Traits<Rhs>::sizeCols_,
    // if there is more block on the left side, we use bp --> result is by_col_
    orient_   =  ( sizeRows_ == UnknownSize || sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : (Traits<Lhs>::sizeRows_)/blockSize < (Traits<Rhs>::sizeCols_)/blockSize
                   ? int(Arrays::by_row_) : int(Arrays::by_col_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
              ? int(Arrays::dense_) : int(Arrays::sparse_)
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ReturnType;

  typedef CAllocator<Type, sizeRows_ , sizeCols_, orient_> Allocator;
};

template<typename Lhs, typename Rhs>
struct ProductTraits<Lhs, Rhs, Arrays::upper_triangular_, Arrays::upper_triangular_>
{
  enum
  {
    structure_ = Arrays::upper_triangular_,
    sizeRows_ = ((int)Traits<Lhs>::sizeRows_ < (int)Traits<Rhs>::sizeCols_) ? Traits<Lhs>::sizeRows_ : Traits<Rhs>::sizeCols_,
    sizeCols_ = ((int)Traits<Lhs>::sizeCols_ < (int)Traits<Rhs>::sizeRows_) ? Traits<Lhs>::sizeCols_ : Traits<Rhs>::sizeRows_,
    // if there is more block on the left side, we use bp --> result is by_col_
    orient_   =  ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : (Traits<Lhs>::sizeRows_)/blockSize < (Traits<Rhs>::sizeCols_)/blockSize
                   ? int(Arrays::by_row_) : int(Arrays::by_col_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
              ? int(Arrays::dense_) : int(Arrays::sparse_)
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ReturnType;

  typedef Array2DUpperTriangular<Type> Allocator;  // no CAllocator
};

} // namespace hidden

/** @ingroup Arrays
  * @class PointByArrayProduct
  * @brief Generic expression where a product operator is
  * applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side
  * @tparam Rhs the type of the right-hand side
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions. The left hand Side is a point, the Right Hand Side
  * is any compatible exprssion or matrix.
  **/
template<typename Lhs, typename Rhs> class PointByArrayProduct;

/** @ingroup Arrays
  * @class ArrayByVectorProduct
  * @brief Generic expression where a product operator is
  * applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side
  * @tparam Rhs the type of the right-hand side
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions. The left hand Side can be of any kind, the Right Hand Side
  * is a vector or vectorial expression.
  **/
template<typename Lhs, typename Rhs> class ArrayByVectorProduct;

/** @ingroup Arrays
  * @class ArrayByDiagonalProduct
  * @brief Generic expression where a product operator is
  * applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side array expression
  * @tparam Rhs the type of the right-hand side diagonal expression
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions. The left hand Side can be of any kind, the Left Hand Side
  * is a diagonal matrix or expression.
  **/
template<typename Lhs, typename Rhs> class ArrayByDiagonalProduct;

/** @ingroup Arrays
  * @class DiagonalByArrayProduct
  * @brief Generic expression where a product operator is
  * applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side, a diagonal expression
  * @tparam Rhs the type of the right-hand side, a vector expression
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions. The left hand Side is a diagonal expression,
  * the Right Hand Side is a vector expression.
  **/
template<typename Lhs, typename Rhs> class DiagonalByArrayProduct;

/** @ingroup Arrays
  * @class VectorByPointProduct
  * @brief Generic expression where a product operator is
  * applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side, a vector expression
  * @tparam Rhs the type of the right-hand side, a row-vector expression
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions. The left hand Side is a vector expression,
  * the Right Hand Side is a row-vector expression.
  **/
template<typename Lhs, typename Rhs> class VectorByPointProduct;

/** @ingroup Arrays
  * @class ArrayByArrayProduct
  * @brief Generic expression where a product operator is
  * applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side
  * @tparam Rhs the type of the right-hand side
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions.
  * It is the return type of product operators, by which we mean only those
  * product operators where both the left-hand side and the right-hand side
  * are expressions. For example, the return type of matrix1*matrix2 is a
  * ArrayByArrayProduct.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name ArrayByArrayProduct types explicitly.
  **/
template<typename Lhs, typename Rhs> class ArrayByArrayProduct;

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for array2d by diagonal product
 */
template< typename Lhs, typename Rhs>
struct Traits< ArrayByDiagonalProduct < Lhs, Rhs> >
{
  enum
  {
    structure_ = Lhs::structure_,
    sizeRows_  = Lhs::sizeRows_ ,
    sizeCols_  = Rhs::sizeCols_,
    orient_    = Lhs::orient_,
    storage_   = Lhs::storage_
  };
  typedef typename Lhs::Type Type;
  typedef typename RemoveConst<Type>::Type ReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for the DiagonalByArrayProduct class
 */
template< typename Lhs, typename Rhs>
struct Traits< DiagonalByArrayProduct < Lhs, Rhs> >
{
  enum
  {
    structure_ = Rhs::structure_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Rhs::sizeCols_,
    orient_    = Rhs::orient_,
    storage_   = Rhs::storage_
  };
  typedef typename Promote<typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type ReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for the PointByArrayProduct
 */
template< typename Lhs, typename Rhs>
struct Traits< PointByArrayProduct < Lhs, Rhs> >
{
  enum
   {
     structure_ = Arrays::point_,
     sizeRows_  = 1,
     sizeCols_  = Traits<Rhs>::sizeCols_,
     orient_    = Arrays::by_row_,
     storage_   = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                  ? int(Arrays::dense_) : int(Arrays::sparse_)
   };

  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ReturnType;

  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

/** @ingroup hidden
 *  @brief Traits class for the ArrayByVectorProduct class
 */
template< typename Lhs, typename Rhs>
struct Traits< ArrayByVectorProduct < Lhs, Rhs> >
{
  enum
   {
     structure_ = Arrays::vector_,
     orient_    = Arrays::by_col_,
     sizeRows_  = Traits<Lhs>::sizeRows_,
     sizeCols_  = 1,
     storage_   = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                  ? int(Arrays::dense_) : int(Arrays::sparse_)
   };

  typedef ProductTraits<Lhs, Rhs, Lhs::structure_, Rhs::structure_> Base;
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ReturnType;

  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};
/** @ingroup hidden
 *  @brief Traits class for the DiagonalByArrayProduct class
 */
template< typename Lhs, typename Rhs>
struct Traits< VectorByPointProduct < Lhs, Rhs> >
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_  = Traits<Lhs>::sizeRows_,
    sizeCols_  = Traits<Rhs>::sizeCols_,
    orient_    = Arrays::by_col_,
    storage_   = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                 ? int(Arrays::dense_) : int(Arrays::sparse_)
  };
  typedef typename Promote<typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef Type ReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for the ArrayByArrayProduct class
 */
template< typename Lhs, typename Rhs>
struct Traits< ArrayByArrayProduct<Lhs, Rhs> >
{
  typedef ProductTraits<Lhs, Rhs, Traits<Lhs>::structure_, Traits<Rhs>::structure_> Base;

  enum
  {
    structure_ = Base::structure_,
    sizeRows_ = Traits<Lhs>::sizeRows_,
    sizeCols_ = Traits<Rhs>::sizeCols_,
    // if there is more block on the left side, we use bp --> result is by_col_
    orient_   =  ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : (Traits<Lhs>::sizeRows_)/blockSize < (Traits<Rhs>::sizeCols_)/blockSize
                   ? int(Arrays::by_row_) : int(Arrays::by_col_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                  ? int(Arrays::dense_) : int(Arrays::sparse_)
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ReturnType;
  typedef typename Base::Allocator Allocator;
};

} // end namespace hidden

template<typename Lhs, typename Rhs>
class ArrayByDiagonalProduct : public ExprBase< ArrayByDiagonalProduct<Lhs, Rhs> >
                             , public TRef<1>
{
  public:
    typedef ExprBase< ArrayByDiagonalProduct<Lhs, Rhs> > Base;
    typedef typename hidden::Traits<ArrayByDiagonalProduct>::Type Type;
    typedef typename hidden::Traits<ArrayByDiagonalProduct>::ReturnType ReturnType;

    enum
    {
      structure_ = hidden::Traits<ArrayByDiagonalProduct>::structure_,
      orient_    = hidden::Traits<ArrayByDiagonalProduct>::orient_,
      sizeRows_  = hidden::Traits<ArrayByDiagonalProduct>::sizeRows_,
      sizeCols_  = hidden::Traits<ArrayByDiagonalProduct>::sizeCols_,
      storage_   = hidden::Traits<ArrayByDiagonalProduct>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    inline ArrayByDiagonalProduct( const Lhs& lhs, const Rhs& rhs)
                                  : Base(), lhs_(lhs), rhs_(rhs)
    {
      if (lhs.cols() != rhs.rows())
      { STKRUNTIME_ERROR_NO_ARG(ArrayByDiagonalProduct, sizes mismatch);}
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs_.rows();}
    /** @return the first index of the rows */
    inline int beginRowsImpl() const { return(lhs_.beginRows());}
    /** @return the ending index of the rows */
    inline int endRowsImpl() const { return(lhs_.endRows());}
    /** @return the number of rows */
    inline int sizeRowsImpl() const { return lhs_.sizeRows();}

    /** @return the range of the columns */
    inline ColRange const& colsImpl() const { return rhs_.cols();}
    /** @return the first index of the columns */
    inline int beginColsImpl() const { return(rhs_.beginCols());}
    /** @return the ending index of the columns */
    inline int endColsImpl() const { return(rhs_.endCols());}
    /** @return the number of columns */
    inline int sizeColsImpl() const { return rhs_.sizeCols();}

    inline ReturnType elt2Impl(int i, int j) const { return lhs_.elt(i,j)*rhs_.elt(j);}
    /** access to the ith element */
    inline ReturnType elt1Impl(int i) const { return lhs_.elt(i)*rhs_.elt(i);}
    /** access to the element */
    inline ReturnType elt0Impl() const { return lhs_.elt()*rhs_.elt();}


    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;
};

template<typename Lhs, typename Rhs>
class DiagonalByArrayProduct : public ExprBase< DiagonalByArrayProduct<Lhs, Rhs> >
                              , public TRef<1>
{
  public:
    typedef ExprBase< DiagonalByArrayProduct<Lhs, Rhs> > Base;
    typedef typename hidden::Traits<DiagonalByArrayProduct>::Type Type;
    typedef typename hidden::Traits<DiagonalByArrayProduct>::ReturnType ReturnType;

    enum
    {
      structure_ = hidden::Traits<DiagonalByArrayProduct>::structure_,
      orient_    = hidden::Traits<DiagonalByArrayProduct>::orient_,
      sizeRows_  = hidden::Traits<DiagonalByArrayProduct>::sizeRows_,
      sizeCols_  = hidden::Traits<DiagonalByArrayProduct>::sizeCols_,
      storage_   = hidden::Traits<DiagonalByArrayProduct>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    inline DiagonalByArrayProduct( const Lhs& lhs, const Rhs& rhs)
                                 : Base(), lhs_(lhs), rhs_(rhs)
    {
      if (lhs.cols() != rhs.rows())
      { STKRUNTIME_ERROR_NO_ARG(DiagonalByArrayProduct, sizes mismatch);}
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs_.rows();}
    /** @return the first index of the rows */
    inline int beginRowsImpl() const { return(lhs_.beginRows());}
    /** @return the ending index of the rows */
    inline int endRowsImpl() const { return(lhs_.endRows());}
    /** @return the number of rows */
    inline int sizeRowsImpl() const { return lhs_.sizeRows();}

    /** @return the range of the columns */
    inline ColRange const& colsImpl() const { return rhs_.cols();}
    /** @return the first index of the columns */
    inline int beginColsImpl() const { return(rhs_.beginCols());}
    /** @return the ending index of the columns */
    inline int endColsImpl() const { return(rhs_.endCols());}
    /** @return the number of columns */
    inline int sizeColsImpl() const { return rhs_.sizeCols();}

    /** access to the element (i,j) */
    inline ReturnType elt2Impl(int i, int j) const { return lhs_.elt(i)*rhs_.elt(i,j);}
    /** access to the ith element */
    inline ReturnType elt1Impl(int i) const { return lhs_.elt(i)*rhs_.elt(i);}
    /** access to the element */
    inline ReturnType elt0Impl() const { return lhs_.elt()*rhs_.elt();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;
};

template<typename Lhs, typename Rhs>
class PointByArrayProduct : public ExprBase< PointByArrayProduct<Lhs, Rhs> >
                           , public TRef<1>
{
  public:
    typedef ExprBase< PointByArrayProduct<Lhs, Rhs> > Base;
    typedef typename hidden::Traits<PointByArrayProduct>::Type Type;
    typedef typename hidden::Traits<PointByArrayProduct>::ReturnType ReturnType;
    typedef typename hidden::Traits<PointByArrayProduct>::Allocator Allocator;

    enum
    {
      structure_ = hidden::Traits<PointByArrayProduct>::structure_,
      orient_    = hidden::Traits<PointByArrayProduct>::orient_,
      sizeRows_  = hidden::Traits<PointByArrayProduct>::sizeRows_,
      sizeCols_  = hidden::Traits<PointByArrayProduct>::sizeCols_,
      storage_   = hidden::Traits<PointByArrayProduct>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    inline PointByArrayProduct( const Lhs& lhs, const Rhs& rhs)
                              : Base(), lhs_(lhs), rhs_(rhs)
                              , result_(1, rhs.sizeCols(), Type(0))
    {
      if (lhs.range() != rhs.rows())
      { STKRUNTIME_ERROR_2ARG(PointByArrayProduct, lhs.range(), rhs.rows(), sizes mismatch);}
      result_.shift(rhs_.beginCols());
      hidden::ProductImpl<Lhs, Rhs, Allocator>::run(lhs, rhs, result_);
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return result_.rows();}
    /** @return the first index of the rows */
    inline int beginRowsImpl() const { return result_.beginRows();}
    /** @return the ending index of the rows */
    inline int endRowsImpl() const { return result_.endRows();}
    /** @return the size of the vector */
    inline int sizeRowsImpl() const { return 1;}

    /** @return the range of the columns */
    inline ColRange const& colsImpl() const { return result_.cols();}
    /** @return the first index of the columns */
    inline int beginColsImpl() const { return result_.beginCols();}
    /** @return the ending index of the columns */
    inline int endColsImpl() const { return result_.endCols();}
    /** @return the fixed size type if available to enable compile time optimizations */
    inline int sizeColsImpl() const { return result_.sizeCols();}

    /** @return the element (i,j) */
    inline ReturnType elt2Impl(int i, int j) const { return result_.elt(i, j);}
    /** @return the ith element */
    inline ReturnType elt1Impl(int i) const { return result_.elt(i);}
    /** @return the element */
    inline ReturnType elt0Impl() const { return result_.elt();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& result() const { return result_; }

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;

  private:
    Allocator result_;
};



template<typename Lhs, typename Rhs>
class ArrayByVectorProduct : public ExprBase< ArrayByVectorProduct<Lhs, Rhs> >
                           , public TRef<1>
{
  public:
    typedef ExprBase< ArrayByVectorProduct<Lhs, Rhs> > Base;
    typedef typename hidden::Traits<ArrayByVectorProduct>::Type Type;
    typedef typename hidden::Traits<ArrayByVectorProduct>::ReturnType ReturnType;
    typedef typename hidden::Traits<ArrayByVectorProduct>::Allocator Allocator;
    enum
    {
      structure_ = hidden::Traits<ArrayByVectorProduct>::structure_,
      orient_    = hidden::Traits<ArrayByVectorProduct>::orient_,
      sizeRows_  = hidden::Traits<ArrayByVectorProduct>::sizeRows_,
      sizeCols_  = hidden::Traits<ArrayByVectorProduct>::sizeCols_,
      storage_   = hidden::Traits<ArrayByVectorProduct>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    inline ArrayByVectorProduct( const Lhs& lhs, const Rhs& rhs)
                              : Base(), lhs_(lhs), rhs_(rhs)
                              , result_(lhs.sizeRows(), 1, Type(0))
    {
      if (lhs.cols() != rhs.range())
      { STKRUNTIME_ERROR_NO_ARG(ArrayByVectorProduct, sizes mismatch);}
      result_.shift(lhs_.beginRows());
      hidden::ProductImpl<Lhs, Rhs, Allocator>::run(lhs, rhs, result_);
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return result_.rows();}
    /** @return the first index of the rows */
    inline int beginRowsImpl() const { return(result_.beginRows());}
    /** @return the ending index of the rows */
    inline int endRowsImpl() const { return(result_.endRows());}
    /** @return the number of rows */
    inline int sizeRowsImpl() const { return result_.sizeRows();}

    /** @return the range of the columns */
    inline ColRange const& colsImpl() const { return result_.cols();}
    /** @return the first index of the columns */
    inline int beginColsImpl() const { return(result_.beginCols());}
    /** @return the ending index of the columns */
    inline int endColsImpl() const { return(result_.endCols());}
    /** @return the number of columns */
    inline int sizeColsImpl() const { return int(1);}

    /** @return the element (i,j) */
    inline ReturnType elt2Impl(int i, int j) const { return result_.elt(i, j);}
    /** @return the ith element */
    inline ReturnType elt1Impl(int i) const { return result_.elt(i);}
    /** @return the element */
    inline ReturnType elt0Impl() const { return result_.elt();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& result() const { return result_; }

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;

  private:
    Allocator result_;
};

template<typename Lhs, typename Rhs>
class VectorByPointProduct : public ExprBase< VectorByPointProduct<Lhs, Rhs> >
                           , public TRef<1>
{
  public:
    typedef ExprBase< VectorByPointProduct<Lhs, Rhs> > Base;
    typedef typename hidden::Traits<VectorByPointProduct>::Type Type;
    typedef typename hidden::Traits<VectorByPointProduct>::ReturnType ReturnType;

    enum
    {
      structure_ = hidden::Traits<VectorByPointProduct>::structure_,
      orient_    = hidden::Traits<VectorByPointProduct>::orient_,
      sizeRows_  = hidden::Traits<VectorByPointProduct>::sizeRows_,
      sizeCols_  = hidden::Traits<VectorByPointProduct>::sizeCols_,
      storage_   = hidden::Traits<VectorByPointProduct>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    inline VectorByPointProduct( const Lhs& lhs, const Rhs& rhs)
                               : Base(), lhs_(lhs), rhs_(rhs)
    {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs_.rows();}
    /** @return the first index of the rows */
    inline int beginRowsImpl() const { return(lhs_.beginRows());}
    /** @return the ending index of the rows */
    inline int endRowsImpl() const { return(lhs_.endRows());}
    /** @return the number of rows */
    inline int sizeRowsImpl() const { return lhs_.sizeRows();}

    /** @return the range of the columns */
    inline ColRange const& colsImpl() const { return rhs_.cols();}
    /** @return the first index of the columns */
    inline int beginColsImpl() const { return(rhs_.beginCols());}
    /** @return the ending index of the columns */
    inline int endColsImpl() const { return(rhs_.endCols());}
    /** @return the number of columns */
    inline int sizeColsImpl() const { return rhs_.sizeCols();}

    /** @return the element (i,j) */
    inline ReturnType elt2Impl(int i, int j) const { return lhs_.elt(i)*rhs_.elt(j);}
    /** @return the element */
    inline ReturnType elt0Impl() const { return lhs_.elt()*rhs_.elt();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side expression */
    inline Rhs const& rhs() const { return rhs_; }

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;
};
// forward declaration
template< typename Lhs, typename Rhs> class ArrayByArrayProductBase;

template<typename Lhs, typename Rhs>
class ArrayByArrayProduct : public ArrayByArrayProductBase< Lhs, Rhs >, public TRef<1>
{
  public:
    typedef ArrayByArrayProductBase< Lhs, Rhs> Base;
    typedef typename hidden::Traits<ArrayByArrayProduct>::Type Type;
    typedef typename hidden::Traits<ArrayByArrayProduct>::ReturnType ReturnType;

    typedef typename hidden::Traits<ArrayByArrayProduct < Lhs, Rhs> >::Allocator Allocator;

    enum
    {
      // All the valid cases for ArrayByArray operator
     isValid_ =(  (EGAL(Lhs,array2D_) && !EGAL(Rhs,point_)  && !EGAL(Rhs,number_) && !EGAL(Rhs,diagonal_) )
               || (EGAL(Lhs,square_)  && !EGAL(Rhs,point_)  && !EGAL(Rhs,number_) && !EGAL(Rhs,diagonal_) )
               || (EGAL(Lhs,lower_triangular_) && !EGAL(Rhs,point_)  && !EGAL(Rhs,number_) && !EGAL(Rhs,diagonal_) )
               || (EGAL(Lhs,upper_triangular_) && !EGAL(Rhs,point_)  && !EGAL(Rhs,number_) && !EGAL(Rhs,diagonal_) )
               || (EGAL(Lhs,vector_) && EGAL(Rhs,point_) )
               || (EGAL(Lhs,point_)  && !EGAL(Rhs,point_) && !EGAL(Rhs,vector_) && !EGAL(Rhs,number_) && !EGAL(Rhs,diagonal_) )
               ),

      structure_ = hidden::Traits<ArrayByArrayProduct>::structure_,
      orient_    = hidden::Traits<ArrayByArrayProduct>::orient_,
      sizeRows_  = hidden::Traits<ArrayByArrayProduct>::sizeRows_,
      sizeCols_  = hidden::Traits<ArrayByArrayProduct>::sizeCols_,
      storage_   = hidden::Traits<ArrayByArrayProduct>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    inline ArrayByArrayProduct( const Lhs& lhs, const Rhs& rhs)
                              : Base(), lhs_(lhs), rhs_(rhs)
                              , result_(lhs.sizeRows(), rhs.sizeCols(), Type(0))
    {
      STK_STATICASSERT_PRODUCT_OPERATOR_MISMATCH( isValid_ );
      if (lhs.cols() != rhs.rows())
      { STKRUNTIME_ERROR_NO_ARG(ArrayByArrayProduct,sizes mismatch for 2D array);}
      result_.shift(lhs_.beginRows(), rhs_.beginCols());
      (orient_) ? hidden::ProductImpl<Lhs, Rhs, Allocator>::runpb(lhs, rhs, result_)
                : hidden::ProductImpl<Lhs, Rhs, Allocator>::runbp(lhs, rhs, result_);
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs_.rows();}
    /** @return the first index of the rows */
    inline int beginRowsImpl() const { return(lhs_.beginRows());}
    /** @return the ending index of the rows */
    inline int endRowsImpl() const { return(lhs_.endRows());}
    /** @return the number of rows */
    inline int sizeRowsImpl() const { return(lhs_.sizeRows());}

    /** @return the range of the columns */
    inline ColRange const& colsImpl() const { return rhs_.cols();}
    /** @return the first index of the columns */
    inline int beginColsImpl() const { return(rhs_.beginCols());}
    /** @return the ending index of the columns */
    inline int endColsImpl() const { return(rhs_.endCols());}
    /** @return the number of columns */
    inline int sizeColsImpl() const { return rhs_.sizeCols();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }
    /** @return the result */
    inline Allocator const& result() const { return result_; }

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;

  private:
    Allocator result_;
};

/** @ingroup Arrays
  * @brief implement the access to the elements in the (2D) general case.
  **/
template< typename Lhs, typename Rhs>
class ArrayByArrayProductBase : public ExprBase< ArrayByArrayProduct<Lhs, Rhs> >
{
  public:
    typedef ArrayByArrayProduct<Lhs, Rhs> Derived;
    typedef ExprBase< Derived> Base;

    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::ReturnType ReturnType;

    /** constructor. */
    inline ArrayByArrayProductBase() : Base() {}
    /** access to the element (i,j) */
    inline Type const& elt2Impl(int i, int j) const
    { return this->asDerived().result().elt(i,j);}
    /** access to the element i */
    inline Type const& elt1Impl(int i) const
    { return this->asDerived().result().elt(i);}
    /** access to the element */
    inline Type const& elt0Impl() const
    { return this->asDerived().result().elt();}
};

}  // namespace STK

#undef EGAL

#endif /* STK_PRODUCTOPERATOR_H */
