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

/** @file STK_lapack_Svd.h
 *  @brief In this file we implement the enclosing class of the dgeqrf lapack routine.
 **/

#include "../include/STK_lapack_Svd.h"

namespace STK
{
namespace lapack
{
/** @brief Run svd decomposition */
bool Svd::runImpl()
{
  if (jobz_ == 'A' || jobz_ == 'S')
  {
    msg_error_ = _T("In lapack::Svd::runImpl, the options 'A' and 'S' are not available");
    return false;
  }
  int beginRow = U_.beginRows(), beginCol = U_.beginCols();
  bool result = true;
  // compute results
  if ( (jobz_ == 'O' && U_.sizeRows() >= U_.sizeCols()) || jobz_ == 'N')
  {
    if (!computeSvd(U_, U_, D_, V_)) { result = false;}
  }
  else
  { // jobz_ == 'O' and m<n, V_ will contain u and U_ will contain vt
    if (!computeSvd(U_, V_, D_, V_)) { return false;};
    U_.exchange(V_); // U_ is (m,m) and V_ is (m,n)
  }
  U_.shift(beginRow, beginCol);
  D_.shift(beginCol);
  V_.shift(beginCol); // u*s.diagonalize()*vt work
  return result;
}


/* Computing the Svd decomposition of the matrix U_. */
bool Svd::computeSvd(CArrayXX& a, CArrayXX& u, CVectorX& s, CArrayXX& vt)
{
  int m = a.sizeRows(), n = a.sizeCols(), nbSv = std::min(m,n);
  a.shift(0,0);
  // Workspace and status variables:
  double workSize;
  double *work = &workSize;
  int* iwork = new int[8*nbSv];
  int lwork = -1;
  int info;
  //
  // Call dgesdd_ with lwork = -1 to query optimal workspace size:
  info = gesdd( jobz_, m, n
               , a.p_data(), m, s.p_data(), u.p_data(), m, vt.p_data(), n
               , work, lwork, iwork);
  // check any error
  if (info!=0)
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(lapack::Svd::computeSvd,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(lapack::Svd::computeSvd get,-info,error parameter);
    return false;
  }
  // optimal storage dimension is returned in work[0]
  lwork = workSize;
  work = new Real[lwork];
  // resize u
  if ( !((jobz_ == 'O' && m >= n) || (jobz_ == 'N')) )
  {
    int ucol = (jobz_ == 'A' || jobz_ == 'O') ? m : nbSv;
    u.resize(m, ucol);
  }
  if ( !((jobz_ == 'O' && m < n)  || (jobz_ == 'N')) )
  {
    int ldvt = (jobz_ == 'A' || jobz_ == 'O') ? n : nbSv;
    vt.resize(ldvt,n);
  }
  s.resize(nbSv).shift(0);
  u.shift(0,0);
  vt.shift(0,0);

  // Call dgesdd_ to do the actual computation:
  info = gesdd( jobz_, m, n
              , a.p_data(), a.sizeRows(), s.p_data()
              , u.p_data(), u.sizeRows()
              , vt.p_data(), vt.sizeRows()
              , work, lwork, iwork);
  // clean
  delete[] work;
  delete[] iwork;
  // check any error
  if (info!=0)
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(lapack::Svd::computeSvd,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(lapack::Svd::computeSvd get,-info,error parameter);
    return false;
  }
  return true;
}
//
inline int Svd::gesdd( char jobz, int m, int n, Real *a, int lda
                     , Real *s, Real *u, int ldu, Real *vt, int ldvt
                     , Real *work, int lWork, int *iWork
         )
{
  int info = 1;

#ifdef STKUSELAPACK
#ifdef STKREALAREFLOAT
  sgesdd_(&m, &n, a, &lda, tau, work, &lWork, &info);
#else
  dgesdd_( &jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lWork, iWork, &info);
#endif
#endif

  return info;
}

} // namespace lapack

} // namespace STK

