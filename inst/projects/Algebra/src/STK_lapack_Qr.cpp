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
 * Project:  stkpp::Algebra
 * created on: 20 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_lapack_Qr.h
 *  @brief In this file we define the enclosing class of the dgeqrf lapack routine.
 **/

#include "../include/STK_lapack_Qr.h"

namespace STK
{

namespace lapack
{

/* private method for computing the Qr decomposition using a CArrayXX array */
bool Qr::computeQr(CArrayXX& a, CVectorX& tau)
{
  // start qr iterations
  int lWork =-1, m = a.sizeRows(), n= a.sizeCols();
  int info;
  Real iwork;
  Real *p_work;

  // get size for work
  info = geqrf(m, n, a.p_data(), m, tau.p_data(), &iwork, lWork);
  // check any error
  if (info!=0)
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(lapack::Qr::computeQr,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(lapack::Qr::computeQr get,-info,error parameter);
    return false;
  }
  // create
  lWork = (int)iwork;
  p_work = new Real[lWork];
  info = geqrf(m, n, a.p_data(), m, tau.p_data(), p_work, lWork);
  // release working set
  delete[] p_work;
  // check any error
  if (info!=0)
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(lapack::Qr::computeQr,internal error);
     return false;
    }
    this->msg_error_= STKERROR_1ARG(lapack::Qr::computeQr get,-info,error parameter);
    return false;
  }
  return true;
}
/* Computing the QR decomposition of the matrix Q_. */
bool Qr::runImpl()
{
  int begin = Q_.beginRows(); // same first index checked at construction of QR
  int end   = std::min(Q_.endRows(), Q_.endCols());
  int size  = std::min(Q_.sizeRows(), Q_.sizeCols());
  CVectorX tau(size);
  CArrayXX a = Q_;
  a.shift(0, 0);
  tau.shift(0);
  // compute results
  if (!computeQr(a, tau)) { return false;};
  // get results
  a.shift(begin, begin);
  tau.shift(begin);
  R_.resize(Q_.rows(), Q_.cols());
  for (int i=Q_.beginRows(); i< Q_.endRows(); ++i)
  {
    int jmax = std::min(i, Q_.endCols());
    for (int j=Q_.beginCols(); j< jmax; ++j) { Q_(i,j) = a(i,j);}
    for (int j=i; j< Q_.endCols(); ++j)      { R_(i,j) = a(i,j);}
  }
  for (int j=Q_.beginCols(); j<end; ++j) { Q_(j,j) = -tau[j];}
  return true;
}

inline int Qr::geqrf(int m, int n, Real* a, int lda, Real* tau, Real *work, int lWork)
{
  int info = 1;

#ifdef STKUSELAPACK
#ifdef STKREALAREFLOAT
  sgeqrf_(&m, &n, a, &lda, tau, work, &lWork, &info);
#else
  dgeqrf_(&m, &n, a, &lda, tau, work, &lWork, &info);
#endif
#endif

  return info;
}

} // namespace lapack

} // namespace STK

