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

/** @file RcppVector.h
 *  @brief In this file we define the wrapper of the the Rcpp vectors.
 **/


#ifndef RCPPVECTOR_H
#define RCPPVECTOR_H

#include "../projects/Arrays/include/STK_ExprBaseVisitor.h"
#include "../projects/Arrays/include/STK_ExprBaseDot.h"
#include "../projects/Arrays/include/STK_ExprBaseProduct.h"

#include "../projects/Arrays/include/STK_ArrayBaseApplier.h"
#include "../projects/Arrays/include/STK_ArrayBaseAssign.h"
#include "../projects/Arrays/include/STK_ArrayBaseInitializer.h"

namespace STK
{

// forward declaration
template <typename Type > class RcppVector;

namespace hidden
{

/** @ingroup hidden
 *  @brief Specialization of the Traits class for the  STK::IntegerVector class.
 **/
template<typename Type_>
struct Traits< RcppVector<Type_> >
{
  private:
    class Void {};
  public:
    typedef Type_ Type;
    typedef RowOperator< RcppVector<Type_> > Row;
    typedef ColOperator< RcppVector<Type_> > Col;
    typedef RowOperator< RcppMatrix<Type_> > SubRow;
    typedef ColOperator< RcppMatrix<Type_> > SubCol;
    typedef SubOperator< RcppMatrix<Type_> > SubVector;
    typedef Void SubArray;

    enum
    {
      structure_ = Arrays::vector_,
      orient_    = Arrays::by_col_,
      sizeRows_  = UnknownSize,
      sizeCols_  = UnknownSize,
      storage_   = Arrays::dense_
    };
};

} // namespace hidden

template <typename Type_>
class RcppVector : public ArrayBase< RcppVector<Type_> >
{
  public:
    typedef typename hidden::Traits<RcppVector >::Type Type;
    typedef typename hidden::Traits<RcppVector >::Row Row;
    typedef typename hidden::Traits<RcppVector >::Col Col;
    typedef typename hidden::Traits<RcppVector >::SubRow SubRow;
    typedef typename hidden::Traits<RcppVector >::SubCol SubCol;
    typedef typename hidden::Traits<RcppVector >::SubVector SubVector;

    enum
    {
      Rtype_ = hidden::RcppTraits<Type_>::Rtype_
    };
    /** Constructor */
    inline RcppVector(Rcpp::Vector<Rtype_> vector) : vector_(vector) {}
    /** Constructor with SEXP. */
    inline RcppVector(SEXP robj) : vector_(robj) {}
    /** set Vector .
     *  @param vector the Rcpp matrix to wrap
     *  @note cannot be passed as const& due to a bug from the (old versions of) Rcpp side
     * */
    inline void setVector( Rcpp::Vector<Rtype_> vector) { vector_ = vector;}

    /** @return the Vertical range */
    inline Range rows() const { return Range(0, vector_.size());}
    /** @return the index of the first row */
    inline int beginRowsImpl() const { return 0;}
    /** @return the ending index of the rows */
    inline int endRowsImpl() const { return vector_.size();}
    /** @return the number of rows */
    inline int sizeRowsImpl() const { return vector_.size();}

    /**@return the Horizontal range */
    inline Range cols() const { return Range(0,1);}
    /** @return the index of the first column */
    inline int beginColsImpl() const { return 0;}
    /**  @return the ending index of the columns */
    inline int endColsImpl() const { return 1;}
    /** @return the number of columns */
    inline int sizeColsImpl() const { return 1;}

    /** @return the j-th column of this. */
    inline Col colImpl(int j) const { return Col(this->asDerived(), j);}
    /** @return the i-th row of this. */
    inline Row rowImpl(int i) const { return Row(this->asDerived(), i);}
    /** @return the i-th row of this. */
    inline SubVector sub(Range I) const { return SubVector(this->asDerived(), I);}
    /** @return the ith element of the operator
     *  @param i index of the ith element
     **/
    inline Type const elt1Impl(int i) const { return (vector_[i]);}
    /** overwrite the RcppMatrix with T, using Rcpp operator=.
     *  @param T the container to copy
     **/
    inline RcppVector& operator=( RcppVector const& T)
    { vector_ = T.vector_; return *this;}
    /** overwrite the RcppMatrix with T, using Rcpp operator=.
     *  @param T the container to transfer
     **/
    template<class OtherType>
    inline RcppVector& operator=( RcppVector<OtherType> const& T)
    { vector_ = T.vector_; return *this;}
    /** overwrite the RcppMatrix with T, using Rcpp operator=.
     *  @param vector the container to copy
     **/
    template<int OtherRtype>
    inline RcppVector& operator=( Rcpp::Vector<OtherRtype> const& vector)
    { vector_ = vector; return *this;}
  private:
    Rcpp::Vector<Rtype_> vector_;
};

} // namespace STK


#endif /* RCPPVECTOR_H */
