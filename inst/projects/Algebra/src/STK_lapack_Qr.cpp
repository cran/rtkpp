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
 * Purpose:  Implement the Qr Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Qr.cpp
 *  @brief In this file we implement the Qr Class (QR decomposition).
 **/
 
#include "../include/STK_lapack_Qr.h"

#include "Arrays/include/STK_Array2DPoint.h"
#include "Arrays/include/STK_CArray.h"

namespace STK
{

namespace lapack
{
/* Constructor */
Qr::Qr( Matrix const& data, bool ref):  IQr(data, ref) {}
/* copy constructor */
Qr::Qr( Qr const& decomp): IQr(decomp) {}

/* Computing the QR decomposition of the matrix Q_. */
bool Qr::run()
{
  if (Q_.empty())     // if the container is empty
  {
    compq_ = true;
    return true;
  }
  // start qr iterations
  int lwork =-1, m = Q_.sizeRows(), n= Q_.sizeCols();
  int info, size = std::min(m, n);
  Real iwork;
  Real *p_work, *p_tau;
  p_tau  = new Real[size];
  CArrayXX a = Q_;
  a.shift(0, 0);
  // get size for work
  info = geqrf(m, n, a.p_data(), m, p_tau, &iwork, lwork);
  // check any error
  if (info!=0)
  {
    delete[] p_tau;
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(SymEigen::run,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(SymEigen::run,-info,error parameter);
    return false;
  }
  // create
  lwork = (int)iwork;
  p_work = new Real[lwork];
  info = geqrf(m, n, a.p_data(), m, p_tau, p_work, lwork);
  // check any error
  if (info!=0)
  {
    delete[] p_tau;
    delete[] p_work;
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(SymEigen::run,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(SymEigen::run,-info,error parameter);
    return false;
  }
  // get results
  a.shift(Q_.beginRows(), Q_.beginCols());
  R_.resize(Q_.rows(), Q_.cols());
  for (int i=Q_.beginRows(); i< Q_.endRows(); ++i)
  {
    int jmax = std::min(i, Q_.endCols());
    for (int j=Q_.beginCols(); j< jmax; ++j) { Q_(i,j) = a(i,j);}
    for (int j=i; j< Q_.endCols(); ++j)      { R_(i,j) = a(i,j);}
  }
  for (int i=0, j=Q_.beginCols(); i< size; ++i, ++j) { Q_(j,j) = -p_tau[i];}
  // release
  delete[] p_work;
  delete[] p_tau;

  return true;
}

/* Destructor */
Qr::~Qr() {}

/* Operator = : overwrite the Qr with S. */
Qr& Qr::operator=(Qr const& decomp)
{
  IQr::operator =(decomp);
  return *this;
}

int Qr::geqrf(int m, int n, double* a, int lda, double* tau, double *work, int lwork)
{
  int info = 1;

#ifdef STKUSELAPACK
#ifdef STKREALAREFLOAT
  sgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
#else
  dgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
#endif

#endif

  return info;
}

} // namespace lapack

} // namespace STK

