/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  rtkpp
 * created on: 25 juil. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_RMatrix.h
 *  @brief In this file we define the wrapper of the the Rcpp matrices.
 **/


#ifndef STK_RMATRIX_H
#define STK_RMATRIX_H

#include "../../Arrays/include/STK_ExprBaseVisitor.h"
#include "../../Arrays/include/STK_ExprBaseDot.h"
#include "../../Arrays/include/STK_ExprBaseProduct.h"

#include "../../Arrays/include/STK_ArrayBaseApplier.h"
#include "../../Arrays/include/STK_ArrayBaseAssign.h"
#include "../../Arrays/include/STK_ArrayBaseInitializer.h"


namespace STK
{

// forward declaration
template <typename Type_> class RMatrix;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for the  STK::RMatrix class.
 **/
template<typename Type_>
struct Traits< RMatrix<Type_> >
{
  private:
    class Void {};
  public:
    typedef RowOperator< RMatrix<Type_> > Row;
    typedef RowOperator< RMatrix<Type_> > SubRow;
    typedef ColOperator< RMatrix<Type_> > Col;
    typedef ColOperator< RMatrix<Type_> > SubCol;
    typedef Void SubVector;
    typedef Void SubArray;

    typedef Type_ Type;
    typedef Type const& ReturnType;
    enum
    {
      structure_ = Arrays::array2D_,
      orient_    = Arrays::by_col_,
      sizeRows_  = UnknownSize,
      sizeCols_  = UnknownSize,
      storage_   = Arrays::dense_
    };
};

} // namespace hidden

template <typename Type_>
class RMatrix : public ArrayBase< RMatrix<Type_> >
{
  public:
    typedef typename hidden::Traits<RMatrix<Type_> >::Row Row;
    typedef typename hidden::Traits<RMatrix<Type_> >::Col Col;
    typedef typename hidden::Traits<RMatrix<Type_> >::SubRow SubRow;
    typedef typename hidden::Traits<RMatrix<Type_> >::SubCol SubCol;

    typedef typename hidden::Traits<RMatrix<Type_> >::Type Type;
    typedef typename hidden::Traits<RMatrix<Type_> >::ReturnType ReturnType;
    enum
    {
      structure_ = hidden::Traits<RMatrix<Type_> >::structure_,
      orient_    = hidden::Traits<RMatrix<Type_> >::orient_,
      sizeRows_  = hidden::Traits<RMatrix<Type_> >::sizeRows_,
      sizeCols_  = hidden::Traits<RMatrix<Type_> >::sizeCols_,
      storage_   = hidden::Traits<RMatrix<Type_> >::storage_,

      Rtype_ = hidden::RcppTraits<Type_>::Rtype_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Default Constructor. */
    inline RMatrix() : matrix_(),rows_(), cols_() {}
    /** Constructor with SEXP. */
    inline RMatrix( SEXP robj)
                     : matrix_(robj),rows_(0, matrix_.rows()), cols_(0, matrix_.cols())
    {}
    /** Constructor */
    inline RMatrix( Rcpp::Matrix<Rtype_> matrix)
                     : matrix_(matrix),rows_(0, matrix_.rows()), cols_(0, matrix_.cols())
    {}
    /** set Matrix .
     *  @param matrix the Rcpp matrix to wrap
     *  @note cannot be passed as const& due to a bug from the Rcpop side
     * */
    inline void setMatrix( Rcpp::Matrix<Rtype_> matrix)
    { matrix_ = matrix;
      rows_ = RowRange(0, matrix_.rows());
      cols_ = RowRange(0, matrix_.cols());
    }
    /** cast operator */
    inline operator Rcpp::Matrix<Rtype_>() const { return matrix_;}

    /** @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the index of the first row */
    inline int beginRowsImpl() const { return 0;}
    /** @return the ending index of the rows */
    inline int endRowsImpl() const { return matrix_.rows();}
    /** @return the number of rows */
    inline int sizeRowsImpl() const { return matrix_.rows();}

    /**@return the range of the columns */
    inline ColRange const& colsImpl() const { return cols_;}
    /** @return the index of the first column */
    inline int beginColsImpl() const { return 0;}
    /**  @return the ending index of the columns */
    inline int endColsImpl() const { return matrix_.cols();}
    /** @return the number of columns */
    inline int sizeColsImpl() const { return matrix_.cols();}

    /** @return the j-th column of this. */
    inline Col colImpl(int j) const { return Col(this->asDerived(), j);}
    /** @return the i-th row of this. */
    inline Row rowImpl(int i) const { return Row(this->asDerived(), i);}
   /** @return a constant reference on element (i,j) of the 2D container
     *  @param i, j indexes of the row and of the column
     **/
    inline Type const& elt2Impl(int i, int j) const { return static_cast<Type const&>(matrix_(i,j));}
    /** @return a reference on the element (i,j) of the 2D container.
     *  @param i, j indexes of the row and of the column
     **/
    inline Type& elt2Impl(int i, int j) { return (matrix_(i,j));}
    /** overwrite the RMatrix with mat using Rcpp::operator=.
     *  @param mat the RMatrix to copy
     **/
    inline RMatrix& operator=( RMatrix const& mat)
    {
      matrix_ = mat.matrix_;  rows_ = mat.rows_; cols_ = mat.cols_;
      return *this;
    }
    /** overwrite the RMatrix with mat, using Rcpp::operator=.
     *  @param mat the matrix to copy
     **/
    template<class OtherType>
    inline RMatrix& operator=( RMatrix<OtherType> const& mat)
    {
      matrix_ = mat.matrix_;  rows_ = mat.rows_; cols_ = mat.cols_;
      return *this;
    }
    /** overwrite the RMatrix with mat, using Rcpp::operator=.
     *  @param mat the matrix to copy
     **/
    template<int OtherRtype>
    inline RMatrix& operator=( Rcpp::Matrix<OtherRtype> const& mat)
    {
      matrix_ = mat; rows_ = mat.rows_; cols_ = mat.cols_;
      return *this;
    }
  private:
    Rcpp::Matrix<Rtype_> matrix_;
    RowRange rows_;
    ColRange cols_;
};

} // namespace STK


#endif /* STK_RMATRIX_H */
