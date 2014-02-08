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

/** @file RcppMatrix.h
 *  @brief In this file we define the wrapper of the the Rcpp matrices.
 **/


#ifndef RCPPMATRIX_H
#define RCPPMATRIX_H

#include "../../Arrays/include/STK_ExprBaseVisitor.h"
#include "../../Arrays/include/STK_ExprBaseDot.h"
#include "../../Arrays/include/STK_ExprBaseProduct.h"

#include "../../Arrays/include/STK_ArrayBaseApplier.h"
#include "../../Arrays/include/STK_ArrayBaseAssign.h"
#include "../../Arrays/include/STK_ArrayBaseInitializer.h"


namespace STK
{

// forward declaration
template <typename Type_> class RcppMatrix;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for the  STK::IntegerMatrix class.
 **/
template<typename Type_>
struct Traits< RcppMatrix<Type_> >
{
  private:
    class Void {};
  public:
    typedef Type_ Type;
    typedef RowOperator< RcppMatrix<Type_> > Row;
    typedef RowOperator< RcppMatrix<Type_> > SubRow;
    typedef ColOperator< RcppMatrix<Type_> > Col;
    typedef ColOperator< RcppMatrix<Type_> > SubCol;
    typedef Void SubVector;
    typedef Void SubArray;

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
class RcppMatrix : public ArrayBase< RcppMatrix<Type_> >
{
  public:
    typedef typename hidden::Traits<RcppMatrix >::Type Type;
    typedef typename hidden::Traits<RcppMatrix >::Row Row;
    typedef typename hidden::Traits<RcppMatrix >::Col Col;
    typedef typename hidden::Traits<RcppMatrix >::SubRow SubRow;
    typedef typename hidden::Traits<RcppMatrix >::SubCol SubCol;
    enum
    {
      Rtype_ = hidden::RcppTraits<Type_>::Rtype_
    };
    /** Default Constructor. */
    inline RcppMatrix() : matrix_() {}
    /** Constructor with SEXP. */
    inline RcppMatrix(SEXP robj) : matrix_(robj) {}
    /** Constructor */
    inline RcppMatrix( Rcpp::Matrix<Rtype_> matrix) : matrix_(matrix) {}
    /** set Matrix .
     *  @param matrix the Rcpp matrix to wrap
     *  @note cannot be passed as const& due to a bug from the Rcpop side
     * */
    inline void setMatrix( Rcpp::Matrix<Rtype_> matrix) { matrix_ = matrix;}
    /** cast operator */
    inline operator Rcpp::Matrix<Rtype_>() const { return matrix_;}

    /** @return the Vertical range */
    inline Range rows() const { return Range(0, matrix_.rows());}
    /** @return the index of the first row */
    inline int beginRowsImpl() const { return 0;}
    /** @return the ending index of the rows */
    inline int endRowsImpl() const { return matrix_.rows();}
    /** @return the number of rows */
    inline int sizeRowsImpl() const { return matrix_.rows();}

    /**@return the Horizontal range */
    inline Range cols() const { return Range(0, matrix_.cols());}
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
    inline Type const elt2Impl(int i, int j) const { return (matrix_(i,j));}
    /** @return a reference on the element (i,j) of the 2D container.
     *  @param i, j indexes of the row and of the column
     **/
    inline Type& elt2Impl(int i, int j) { return (matrix_(i,j));}

  private:
    Rcpp::Matrix<Rtype_> matrix_;
};

} // namespace STK


#endif /* RCPPMATRIX_H */
