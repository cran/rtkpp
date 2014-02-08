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
 * Project:  stkpp::Algebra
 * Purpose:  Define The Interface SymEigen Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IQr.h
 *  @brief In this file we define the IQr class (for a general matrix).
 **/
 
#ifndef STK_IQR_H
#define STK_IQR_H

#include "STKernel/include/STK_Real.h"
#include "Sdk/include/STK_IRunner.h"

#include "Arrays/include/STK_CArray.h"
#include "Arrays/include/STK_Array2D.h"
#include "Arrays/include/STK_Array2DVector.h"
#include "Arrays/include/STK_Array2DUpperTriangular.h"

#include "../include/STK_Givens.h"

namespace STK
{
/** @ingroup Algebra
 *  @brief The class IQr is an interface class for the method
 *  computing the QR Decomposition of a Matrix.
 * 
 *  The QR decomposition of a matrix require
 *  - Input:  A matrix of size M-by-N
 *  - Output:
 *     -# Q Matrix of size M-by-N.
 *     -# R Upper Triangular matrix of dimension min(M,N)-by-N
 *     -# \f$ A = QR \f$
 **/
class IQr : public IRunnerBase
{
  protected:
    typedef IRunnerBase Base;
    /** @brief Constructor
     *  @param data reference on a matrix expression
     *  @param ref do we use data as reference for Q_ ?
     */
    inline IQr( Array2D<Real> const& data, bool ref = false)
              : Base(), Q_(data, ref), compq_(false)
    {
      if (data.beginRows() != data.beginCols())
        STKRUNTIME_ERROR_NO_ARG(IQR::IQR,Wrong data set: beginRows row must be equal to beginCols);
    }
    /** @brief Constructor
     *  @param data reference on a matrix expression
     */
    template<class Derived>
    IQr( ExprBase<Derived> const& data) : Base(), compq_(false)
    {
      if (data.beginRows() != data.beginCols())
      { STKRUNTIME_ERROR_NO_ARG(IQR::IQR,Wrong data set: beginRows row must be equal to beginCols);}
      Q_ = data.asDerived();
    }
    /** Copy constructor.
     *  @param decomp the decomposition  to copy
     **/
    IQr( IQr const& decomp);

  public :
    /** virtual destructor */
    inline virtual ~IQr() {}
    /** Operator = : overwrite the this with decomp. */
    IQr& operator=(IQr const& decomp);
    /** Is Q computed ?
     *  @return @c true if Q_ is computed, @c false otherwise
     */
    inline bool isCompQ() const { return compq_;}
    /** give the matrix Q of the QR decomposition.
     * @return the matrix Q of the QR decomposition
     **/
    inline Matrix const& Q() const  { return Q_;}
    /** give the matrix R of the QR decomposition.
     * @return the matrix R of the QR decomposition
     **/
    inline MatrixUpperTriangular const& R() const { return R_;}
    /** Compute Q (to use after run). After the run process, Q_ store
     *  the householder vector in its column. Call compQ, if you want to
     *  obtain Q in its true form.
     *  Without effect if (compq_ == true)
     **/
    void compQ();
    /** Delete the n last columns and update the QR decomposition.
     *  @param n number of column to delete
     **/
    void popBackCols(int n =1);
    /** Delete the column pos and update the QR decomposition.
     *  @param pos the position of the column to delete
     **/
    void eraseCol(int pos);
    /** Add a column with value T and update th QR decomposition.
     *  @param T the column to add
     **/
    template<class ColVector>
    void pushBackCol(ColVector const& T);
    /** Add a column with value T at the position pos and update the QR
     *  decomposition.
     *  @param T the column to insert
     *  @param pos the position of the column to insert
     **/
    template<class ColVector>
    void insertCol(ColVector const& T, int pos);

    /* TODO : Delete the ith row and update the QR decomposition :
     *  default is the last row.
     **/
    //Qr& popBackRows();
    //Qr& eraseRows(int i);

    /* TODO : Add a row with value T and update th QR decomposition :
     *  default is the last row.
     **/
    //Qr& pushBackRows(const Array2DPoint<double> &T);
    //Qr& insertRows(const Array2DPoint<double> &T, int i);

    /** overloading of setData.
     * @param data the data set to set.
     **/
    template<class Derived>
    void setData( ExprBase<Derived> const& data)
    { Q_ = data.asDerived(); R_.clear(); compq_ = false;}

  protected :
    /** Q Matrix of the QR decomposition */
    Matrix Q_;
    /** R Matrix of th QR decomposition */
    MatrixUpperTriangular R_;
    /// is Q computed ?
    bool compq_;
};

/* Adding the last column and update the QR decomposition. */
template<class ColVector>
void IQr::pushBackCol(ColVector const& T)
{
  STK_STATICASSERT(ColVector::structure_==(int)Arrays::vector_||ColVector::structure_==(int)Arrays::point_,YOU_HAVE_TO_USE_A_VECTOR_OR_POINT_IN_THIS_METHOD)
  // check conditions
  if (T.range() != Q_.rows())
  { STKRUNTIME_ERROR_NO_ARG(Qr::pushBackCol,T.range() != Q_.rows());}
  // if Q_ is not computed yet
  if (!compq_) compQ();
  // Adding a column to R
  int lastColR = R_.endCols();
  // Create an auxiliary container
  Vector Rncolr = Q_.transpose() * T; // Rncolr of size Q_.cols()
  // update Q_
  for (int iter = Q_.lastIdxCols()-1, iter1 = Q_.lastIdxCols(); iter>=lastColR; iter--, iter1--)
  {
    Real sinus, cosinus;
    // compute the Givens rotation
    Rncolr[iter] = compGivens( Rncolr[iter], Rncolr[iter1], cosinus, sinus);
    // apply Givens rotation if necessary
    if (sinus)
    { rightGivens(Q_, iter, iter1, cosinus, sinus);}
  }
  // update R_
  R_.pushBackCols();
  R_.col(lastColR).copy(Rncolr.sub(R_.rangeRowsInCol(lastColR)));
}


/* Add a column with value T at the position pos and update the QR
 *  decomposition.
 *  @param T the column to insert
 *  @param pos the position of the column to insert
 **/
template<class ColVector>
void IQr::insertCol(ColVector const& T, int pos)
{
  STK_STATICASSERT(ColVector::structure_==(int)Arrays::vector_||ColVector::structure_==(int)Arrays::point_,YOU_HAVE_TO_USE_A_VECTOR_OR_POINT_IN_THIS_METHOD)
  if (pos < R_.beginCols())
  { STKOUT_OF_RANGE_1ARG(Qr::insertCol,pos,pos<R_.beginCols());}
  if (R_.lastIdxCols() < pos)
  { STKOUT_OF_RANGE_1ARG(Qr::insertCol,pos,pos<R_.lastIdxCols()<pos);}
  if (T.range() != Q_.rows())
  { STKRUNTIME_ERROR_1ARG(Qr::insertCol,pos,T.range() != Q_.rows());}
  // if Q_ is not computed yet
  if (!compq_) compQ();
  // Adding a column to R
  R_.insertCols(pos);
  // update the range of the remaining cols of R_
  R_.update( _R(pos+1, std::min(R_.lastIdxRows(), R_.lastIdxCols())) );
  for (int i=pos+1; i< std::min(R_.endRows(), R_.endCols()); ++i) R_(i,i) = 0.0;

  Vector Rpos =  Q_.transpose() * T;
  // Zeroed the unwanted elements
  for (int iter= Q_.lastIdxCols(), iterm1= Q_.lastIdxCols()-1; iter>pos; iterm1--, iter--)
  {
    Real sinus, cosinus;
    // compute the Givens rotation
    Rpos[iterm1]  = compGivens(Rpos[iterm1], Rpos[iter], cosinus, sinus);
    // apply Givens rotation if necessary
    if (sinus)
    {
      // create a reference on the sub-Matrix
      MatrixUpperTriangular Rsub(R_.col(_R(iter, R_.lastIdxCols())), true);
      // Update the next rows (iter:ncolr_) of R_
      leftGivens( Rsub, iterm1, iter, cosinus, sinus);
      // Update the cols of Q_
      rightGivens(Q_, iterm1, iter, cosinus, sinus);
    }
  }
  // update R_
  R_.col(pos) = Rpos.sub(R_.rangeRowsInCol(pos));
  R_.update(pos);
}

} // namespace STK

#endif //STK_IQR_H
