/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2015  Serge Iovleff, University Lille 1, Inria

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

/** @file STK_RVector.h
 *  @brief In this file we define the wrapper of the the Rcpp vectors.
 **/


#ifndef STK_RVECTOR_H
#define STK_RVECTOR_H

#include "../projects/Arrays/include/STK_ExprBaseVisitor.h"
#include "../projects/Arrays/include/STK_ExprBaseDot.h"
#include "../projects/Arrays/include/STK_ExprBaseProduct.h"

#include "../projects/Arrays/include/STK_ArrayBaseApplier.h"
#include "../projects/Arrays/include/STK_ArrayBaseAssign.h"
#include "../projects/Arrays/include/STK_ArrayBaseInitializer.h"

namespace STK
{

// forward declaration
template <typename Type > class RVector;

namespace hidden
{

/** @ingroup hidden
 *  @brief Specialization of the Traits class for the  STK::IntegerVector class.
 **/
template<typename Type_>
struct Traits< RVector<Type_> >
{
  private:
    class Void {};
  public:
    typedef RowOperator< RVector<Type_> > Row;
    typedef ColOperator< RVector<Type_> > Col;
    typedef RowOperator< RMatrix<Type_> > SubRow;
    typedef ColOperator< RMatrix<Type_> > SubCol;
    typedef SubOperator< RMatrix<Type_> > SubVector;
    typedef Void SubArray;

    typedef Type_ Type;
    typedef Type const& ReturnType;

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
class RVector : public ArrayBase< RVector<Type_> >, public TRef<1>
{
  public:
    typedef typename hidden::Traits<RVector<Type_> >::Row Row;
    typedef typename hidden::Traits<RVector<Type_> >::Col Col;
    typedef typename hidden::Traits<RVector<Type_> >::SubRow SubRow;
    typedef typename hidden::Traits<RVector<Type_> >::SubCol SubCol;
    typedef typename hidden::Traits<RVector<Type_> >::SubVector SubVector;

    typedef typename hidden::Traits<RVector<Type_> >::Type Type;
    typedef typename hidden::Traits<RVector<Type_> >::ReturnType ReturnType;
    enum
    {
      structure_ = hidden::Traits<RVector<Type_> >::structure_,
      orient_    = hidden::Traits<RVector<Type_> >::orient_,
      sizeRows_  = hidden::Traits<RVector<Type_> >::sizeRows_,
      sizeCols_  = hidden::Traits<RVector<Type_> >::sizeCols_,
      storage_   = hidden::Traits<RVector<Type_> >::storage_,

      Rtype_ = hidden::RcppTraits<Type_>::Rtype_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Default Constructor. */
    inline RVector(): vector_(),rows_(), cols_(0,1) {}
    /** Constructor */
    inline RVector( Rcpp::Vector<Rtype_> vector)
                  : vector_(vector), rows_(), cols_(0,1) {}
    /** Constructor with SEXP. */
    inline RVector( SEXP robj)
                  : vector_(robj), rows_(0, vector_.size()), cols_(0,1) {}
    /** Constructor with SEXP. */
    inline RVector( RVector robj, bool ref)
                  : vector_(robj), rows_(0, robj.size()), cols_(0,1) {}

    /** @return the underlying Rcpp matrix */
    inline Rcpp::Vector<Rtype_> const& vector() const { return vector_;}
    /** set Vector .
     *  @param vector the Rcpp matrix to wrap
     *  @note cannot be passed as const& due to a bug from the (old versions of) Rcpp side
     * */
    inline void setVector( Rcpp::Vector<Rtype_> vector)
    { vector_ = vector; rows_ = RowRange(0, vector_.size());}

    /** @return the Vertical range */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the index of the first row */
    inline int beginRowsImpl() const { return 0;}
    /** @return the ending index of the rows */
    inline int endRowsImpl() const { return vector_.size();}
    /** @return the number of rows */
    inline int sizeRowsImpl() const { return vector_.size();}

    /**@return the Horizontal range */
    inline ColRange const& colsImpl() const { return cols_;}
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
    /** @return a constant reference on ith element
     *  @param i index of the ith element
     **/
    inline Type const& elt1Impl(int i) const { return static_cast<Type const&>(vector_[i]);}
    /** @return the ith element of the operator
     *  @param i index of the ith element
     **/
    inline Type& elt1Impl(int i) { return static_cast<Type&>(vector_[i]);}
    /** overwrite the RVector with vec using Rcpp::operator=.
     *  @param vec the vector to copy
     **/
    inline RVector& operator=( RVector const& vec)
    {
      vector_ = vec.vector_RVector; rows_ = vec.rows_;
      return *this;
    }
    /** overwrite the RVector with vec using Rcpp::operator=.
     *  @param vec the vector to transfer
     **/
    template<class OtherType>
    inline RVector& operator=( RVector<OtherType> const& vec)
    {
      vector_ = vec.vector_; rows_ = vec.rows_;
      return *this;
    }
    /** overwrite the RVector with vec using Rcpp::operator=.
     *  @param vec the vector to copy
     **/
    template<int OtherRtype>
    inline RVector& operator=( Rcpp::Vector<OtherRtype> const& vec)
    {
      vector_ = vec;  rows_ = RowRange(0, vec.size());
      return *this;
    }
  private:
    Rcpp::Vector<Rtype_> vector_;
    RowRange rows_;
    const ColRange cols_;
};

} // namespace STK


#endif /* STK_RVECTOR_H */
