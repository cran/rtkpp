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
 * Project:  Algebra
 * Purpose:  Implement the IIQr Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_IQr.cpp
 *  @brief In this file we implement the Interface IIQr.
 **/
 
#include "../include/STK_IQr.h"

#include "Arrays/include/STK_Array2DPoint.h"

#include "../include/STK_Householder.h"
#include "../include/STK_Givens.h"

#ifdef STK_ALGEBRA_DEBUG
#include "../../Arrays/include/STK_Display.h"

template< class Container2D >
void print(Container2D const& A, STK::String const& name)
{
  stk_cout << "print: " << name << _T("\n";);
  stk_cout << name << _T(".isRef() =")        << A.isRef()  << _T("\n");
  stk_cout << name << _T(".capacityHo() =")   << A.capacityHo()  << _T("\n");
  stk_cout << name << _T(".cols() =")      << A.cols()  << _T("\n");
  stk_cout << name << _T(".rows() =")      << A.rows()  << _T("\n\n");
  stk_cout << name << _T(".rangeCols().isRef() =")  << A.rangeCols().isRef() << _T("\n");
  stk_cout << name << _T(".rangeCols() =\n")  << A.rangeCols() << _T("\n");
  stk_cout << name << _T(".capacityCols().isRef() =") << A.capacityCols().isRef()  << _T("\n");
  stk_cout << name << _T(".capacityCols() =\n") << A.capacityCols()  << _T("\n");
}
#endif

namespace STK
{

/* Copy constructor.
 *  @param eigen the EigenValue to copy
 **/
IQr::IQr( IQr const& decomp): Q_(decomp.Q_), R_(decomp.R_), compq_(decomp.compq_)
{}

/* Operator = : overwrite the this with decomp. */
IQr& IQr::operator=(const IQr& decomp)
{
  Q_ = decomp.Q_;
  R_ = decomp.R_;
  compq_ = decomp.compq_;

  return *this;
}

/* Computation of Q. */
void IQr::compQ()
{
#ifdef STK_ALGEBRA_VERBOSE
  stk_cout << _T("Entering IQr::compQ()") << _T("\n");
#endif
  // if Q_ is computed yet
  if (compq_) return;
  // number of non zero cols of Q_  
  int ncol  = std::min(Q_.sizeRows(), Q_.sizeCols()), lastCol;
  // add or remove the column
  if (ncol < Q_.sizeCols())
  {
    Q_.popBackCols(Q_.sizeCols() - ncol);
    lastCol = Q_.lastIdxCols();
  }
  else
  {
    lastCol = Q_.lastIdxCols();
    if (ncol < Q_.sizeRows())
    {
      Q_.pushBackCols(Q_.sizeRows() -ncol);
      // Initialize added columns
      Q_.col( _R( lastCol+1, Q_.lastIdxCols()) ).setValue(0);
      for (int i=lastCol+1; i< Q_.endCols(); ++i) { Q_(i, i) = 1.0;}
    }
  }
  // compute other columns
  for (int i=lastCol, i1= lastCol +1; i>=Q_.beginCols(); i--, i1--)
  {
    //stk_cout << _T("i= ") << i << ", Q_ =\n"<< Q_ << "\n";
    // get current householder vector
    Vector u(Q_, _R(i, Q_.lastIdxRows()), i);
    // Apply Householder vector to the right of the matrix
    leftHouseholder( Q_( _R(i, Q_.lastIdxRows()), _R(i1, Q_.lastIdxCols())), u);
    // Get beta
    Real beta = Q_(i,i);
    // update the column i
    Q_( _R(Q_.beginRows(),i-1) , i) = 0.0;     //     0:(i-1)
    Q_(i,i) += 1.0;                            //     i:i
    Q_( _R(i1, Q_.lastIdxRows()), i ) *= beta; // (i+1):M
    // update the column i
  }
  // Q_ is now computed
  compq_ = true;
#ifdef  STK_ALGEBRA_VERBOSE
  stk_cout << _T("Terminating IQr::compQ().") << _T("\n");
#endif
}


/* Delete the jth column and update the QR decomposition : default
 * is the last col
 **/    
void IQr::popBackCols(int n) { R_.popBackCols(n);}

void IQr::eraseCol(int pos)
{
  if (pos < R_.beginCols())
  { STKOUT_OF_RANGE_1ARG(Qr::eraseCol,pos,pos<R_.beginCols());}
  if (R_.lastIdxCols() < pos)
  { STKOUT_OF_RANGE_1ARG(Qr::eraseCol,pos,pos<R_.lastIdxCols()<pos);}
  // if Q_ is not computed yet
  if (!compq_) compQ();
  // compute the number of iteration for updating to zeroed
  int niter = std::min(R_.lastIdxCols(), R_.lastIdxRows());//R_.beginCols()-1+std::min(R_.sizeRows(), R_.sizeCols());
  // Zeroed the remaining elements (z)
  for (int iter = pos+1; iter<= niter; iter++)
  {
    Real sinus, cosinus;
    // compute the Givens rotation
    R_(iter-1, iter) = compGivens( R_(iter-1, iter), R_(iter, iter), cosinus, sinus);
    R_(iter, iter)   = 0.0;
    // if necessary update R_ and Q_
    if (sinus)
    {
      // create a reference on the sub-Matrix
      MatrixUpperTriangular Rsub(R_.col( _R(iter+1, R_.lastIdxCols()) ), true);
      // Update the next rows (iter1:ncolr_) of R_
      leftGivens(Rsub, iter-1, iter, cosinus, sinus);
      // Update the cols of Q_
      rightGivens(Q_, iter-1, iter, cosinus, sinus);
    }
  }
  // erase the column pos
  R_.eraseCols(pos);
  // update the range of the remaining cols of the container
  R_.update(Range(pos, std::min(R_.lastIdxRows(), R_.lastIdxCols()), 0));
}

} // namespace STK

