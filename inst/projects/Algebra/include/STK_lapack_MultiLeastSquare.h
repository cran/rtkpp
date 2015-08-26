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
 * created on: 1 avr. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_lapack_MultiLeastSquare.h
 *  @brief In this file we define the class MultiLeastSQquare using lapack.
 **/


#ifndef STK_LAPACK_MULTILEASTSQUARE_H
#define STK_LAPACK_MULTILEASTSQUARE_H

#include <STKernel/include/STK_Real.h>
#include <Arrays/include/STK_CArray.h>
#include <Arrays/include/STK_CArrayVector.h>

#include "STK_ILeastSquare.h"

#ifdef STKUSELAPACK

extern "C"
{
/** lapack routine computing the optimal block size */
int ilaenv_(int *, char *, char *, int *, int *,int *, int *);

#ifdef STKREALAREFLOAT
/** LAPACK routine in float to compute the least square solution */
extern
int sgelsd_( int *m, int *n, int *nrhs
           , float *a, int *lda, float *b, int *ldb, float *s, float *rcond
           , int *rank, float *work, int *lWork
           , int *iwork, int *info);
#else
/** LAPACK routine in double to compute the least square solution */
extern int dgelsd_( int *m, int *n, int *nrhs
                  , double *a, int *lda, double *b, int *ldb, double *s, double *rcond
                  , int *rank, double *work, int *lWork
                  , int *iwork, int *info);
#endif
}

#endif // STKUSELAPACK


namespace STK
{

namespace lapack
{
// forward declaration
template<class ArrayB, class ArrayA> class MultiLeastSquare;
}

namespace hidden
{
/** @ingroup hidden
 *  Specialization for the lapack::MultiLeastSquare class.
 **/
template<class ArrayB_, class ArrayA_>
struct AlgebraTraits< lapack::MultiLeastSquare<ArrayB_, ArrayA_> >
{
  typedef ArrayA_ ArrayA;
  typedef ArrayB_ ArrayB;
};

} // namespace hidden


namespace lapack
{
/** @ingroup Algebra
 *  @brief The class MultiLeastSQquare solve the least square problem when
 *  the response @e b is multidimensional.
 *  computes the minimum-norm solution to a real linear least squares
 *  problem: minimize 2-norm(| b - A*x |) using the singular value
 *  decomposition (SVD) of A.
 *  A is an M-by-N matrix which may be rank-deficient.
 **/
template<class ArrayB, class ArrayA>
class MultiLeastSquare: public ILeastSquare< MultiLeastSquare<ArrayB, ArrayA> >
{
  public:
    typedef ILeastSquare< MultiLeastSquare<ArrayB, ArrayA> > Base;
    using Base::b_;
    using Base::a_;
    using Base::x_;
    using Base::rank_;
    /** @brief constructor
     *  @param b,a the left hand side and the right hand side of the least square problem.
     *  @param isBref,isAref are the left hand side and the right hand side references ?
     */
    MultiLeastSquare( ArrayB const& b, ArrayA const& a, bool isBref=false, bool isAref=false)
                    : Base(b, a, isBref, isAref), rcond_(-1) {};
    /** @brief templated constructor
     *  @param b,a the left hand side and the right hand side of the least square
     *  problem.
     */
    template<class ArrayB_, class ArrayA_>
    MultiLeastSquare( ExprBase<ArrayB_> const& b, ExprBase<ArrayA_> const& a)
                    : Base(b, a), rcond_(-1) {}
    /** Destructor */
    virtual ~MultiLeastSquare() {};
    /** @return the condition number */
    inline Real rcond() const { return rcond_;}
    /** return the array with the singular values of A */
    inline CVectorX const& s() const { return s_;}
    /** @param rcond the condition number. If rcond<0, the machine precision is
     *  used.*/
    inline void setRcond(Real rcond) { rcond_ =rcond;}
    /** @brief solve the multi-linear least square problem.
     *  Launch gelsd LAPACK routine to perform the  decomposition.
     *  @return @c true if no error occur, @c false otherwise
     */
    bool runImpl();
    /** @brief solve the weighted least square problem.
     *  Launch gelsd LAPACK routine to perform the  decomposition.
     *  @return @c true if no error occur, @c false otherwise
     */
    template<class Weights>
    bool runImpl(Weights const& weights);
    /** @brief wrapper of the LAPACK GELSD routine.
     *  GELSD computes the minimum-norm solution to a real linear least squares
     *  problem: minimize 2-norm(| b - A*x |) using the singular value
     *  decomposition (SVD) of A. A is an M-by-N matrix which may be
     *  rank-deficient.
     *
     *  @param[in] m  The number of rows of A. <em> m >= 0 </em>.
     *  @param[in] n  The number of columns of A. <em> n>= 0 </em>.
     *  @param[in] nrhs The number of right hand sides, i.e., the number of
     *  columns of the matrices B and X. <em>nrhs >= 0</em>.
     *  @param[in] a On entry, the M-by-N matrix A. On exit, A has been destroyed.
     *  @param[in] lda The leading dimension of the array A.  <em>lda >= max(1,m)</em>.
     *  @param[in,out] b On entry, the M-by-NRHS right hand side matrix B.
     *  @verbatim
     *    On exit, B is overwritten by the N-by-NRHS solution matrix X.
     *    If m >= n and RANK = n, the residual sum-of-squares for the solution in
     *    the i-th column is given by the sum of squares of elements n+1:m in that
     *    column.
     *  @endverbatim
     *  @param[in] ldb The leading dimension of the array B.
     *  <em> ldb >= max(1,max(m,n)) </em>.
     *
     *  @param[out] s The singular values of A in decreasing order.
     *  @verbatim
     *    The condition number of A in the 2-norm = s(1)/s(min(m,n)).
     *  @endverbatim
     *
     *  @param[out] rcond used to determine the effective rank of A.
     *  @verbatim
     *    Singular values s(i) <= RCOND*s(1) are treated as zero.
     *    If RCOND < 0, machine precision is used instead.
     *  @endverbatim
     *
     *  @param[out] rank The effective rank of A, i.e., the number of singular values
     *  which are greater than RCOND*S(1).
     *
     *  @param[out] work On exit, if info = 0, work(1) returns the optimal lWork.
     *
     *  @param[in] lWork The dimension of the array @c work. @c lWork must be at least 1.
     *  @verbatim
     *    The exact minimum amount of workspace needed depends on M, N and NRHS.
     *    As long as lWork is at least 12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,
     *    if M is greater than or equal to N or 12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,
     *    if M is less than N, the code will execute correctly.
     *    SMLSIZ is returned by ILAENV and is equal to the maximum size of
     *    the subproblems at the bottom of the computation tree (usually about 25),
     *    and NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
     *    For good performance, lWork should generally be larger.
     *    If lWork = -1, then a workspace query is assumed; the routine only
     *    calculates the optimal size of the work array, returns this value as the
     *    first entry of the work array, and no error message related to lWork is
     *    issued by XERBLA.
     *  @endverbatim
     *
     *  @param iwork  array of integer of dimension (MAX(1,LiWork))
     *  @verbatim
     *    LiWork >= 3 * MINMN * NLVL + 11 * MINMN, where MINMN = MIN( M,N ).
     *    On exit, if info = 0, iWork(1) returns the minimum LiWork
     *  @endverbatim
     *
     *  @return info
     *  @verbatim
     *    = 0:  successful exit
     *    < 0:  if info = -i, the i-th argument had an illegal value.
     *    > 0:  the algorithm for computing the SVD failed to converge;
     *    if info = i, i off-diagonal elements of an intermediate bidiagonal
     *    form did not converge to zero.
     *  @endverbatim
     *
     *  @verbatim
     *    Further Details
     *    ===============
     *      Based on contributions by
     *      Ming Gu and Ren-Cang Li, Computer Science Division, University of
     *      California at Berkeley, USA
     *      Osni Marques, LBNL/NERSC, USA
     * @endverbatim
     **/
    static int gelsd( int m, int n, int nrhs
                    , Real *a, int lda
                    , Real *b, int ldb
                    , Real *s
                    , Real *rcond, int* rank
                    , Real *work, int lWork, int* iwork);

  protected:
    /** condition number used for determining the effective rank of A */
    Real rcond_;
    /** Array of the singular values */
    CVectorX s_;

  private:
    /** private method for computing the LS solution */
    bool computeLS(CArrayXX& b, CArrayXX& a);
};

/* @brief Run LS solution
 *  Launch gelsd LAPACK routine to compute the solution.
 *  @return @c true if no error occur, @c false otherwise
 */
template<class ArrayB, class ArrayA>
bool MultiLeastSquare<ArrayB, ArrayA>::runImpl()
{
  Range brows = (b_.sizeRows()< a_.sizeCols()) ? a_.cols() : b_.rows();
  // local arrays, b is resized if necessary
  CArrayXX a(a_), b(brows, b_.cols());
  b.sub(b_.rows(), b_.cols()) = b_;
  // start
  return computeLS(b,a);
}
/* @brief Run LS solution
 *  Launch gelsd LAPACK routine to compute the solution.
 *  @return @c true if no error occur, @c false otherwise
 */
template<class ArrayB, class ArrayA>
template<class Weights>
bool MultiLeastSquare<ArrayB, ArrayA>::runImpl(Weights const& w)
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
  Range brows = (b_.sizeRows()< a_.sizeCols()) ? a_.cols() : b_.rows();
  // local arrays, b is resized if necessary
  CArrayXX a(w.sqrt().diagonalize() * a_), b(brows, b_.cols());
  b.sub(b_.rows(), b_.cols()) = w.sqrt().diagonalize() * b_;
  // start
  return computeLS(b,a);
}

/* @brief Run LS solution
 *  Launch gelsd LAPACK routine to compute the solution.
 *  @return @c true if no error occur, @c false otherwise
 */
template<>
inline bool MultiLeastSquare<CArrayXX, CArrayXX>::runImpl()
{ return computeLS(b_,a_);}


/* private method for computing the LS solution */
template<class ArrayB, class ArrayA>
bool MultiLeastSquare<ArrayB, ArrayA>::computeLS(CArrayXX& b, CArrayXX& a)
{
  Range arows = a.rows(), acols = a.cols(), brows = b.rows(), bcols = b.cols();
  int m = arows.size(), n= acols.size(), nrhs = bcols.size();
  // resize if necessary
  if (brows.size()<n)
  {
    if (!b.isRef())
    {
      CArrayXX tmp(b);
      b.resize(acols,bcols);
      b.sub(brows, bcols) = tmp;
    }
    else
    { STKRUNTIME_ERROR_NO_ARG(MultiLeastSquare::computeLS,b has not enough rows);}
  }
  int lda = a.sizeRows(), ldb = b.sizeRows();
  // shift to 0
  a.shift(0,0);
  b.shift(0,0);
  s_.resize(m); s_ = 0.;
  Range srange = s_.range();
  s_.shift(0);
  // get sizes
  Real wkopt;
  int lWork = -1, iwkopt;
  int info = gelsd( m, n, nrhs
                  , a.p_data(), lda
                  , b.p_data(), ldb
                  , s_.p_data()
                  , &rcond_, &rank_
                  , &wkopt, lWork, &iwkopt);
  if( info != 0 )
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(lapack::MultiLeastSquare::computeLS,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(lapack::MultiLeastSquare::computeLS get,-info,error parameter);
    return false;
  }
  // get working sizes
  lWork = (int)wkopt;
  Real* work = new Real[lWork];
  int* iwork = new int[iwkopt];
  // Solve the least square problem
  info = gelsd( m, n, nrhs
              , a.p_data(), lda
              , b.p_data(), ldb
              , s_.p_data()
              , &rcond_, &rank_
              , work, lWork, iwork);
  delete[] iwork;
  delete[] work;
  if( info != 0 )
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(lapack::MultiLeastSquare::computeLS,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(lapack::MultiLeastSquare::computeLS get,-info,error parameter);
    return false;
  }
  // shift back
  a.shift(arows.begin(), acols.begin());
  b.shift(brows.begin(), bcols.begin());
  s_.shift(srange.begin());
  x_ = b.sub(acols,bcols);
  return true;
}

template<class ArrayB, class ArrayA>
int MultiLeastSquare<ArrayB, ArrayA>::gelsd( int m, int n, int nrhs
                                           , Real * a, int lda
                                           , Real * b, int ldb
                                           , Real * s
                                           , Real *rcond, int *rank
                                           , Real *work, int lWork, int* iwork)
{
  int info = 1;
#ifdef STKUSELAPACK
#ifdef STKREALAREFLOAT
  sgelsd_( &m, &n, &nrhs, a, &lda, b, &ldb, s, rcond, rank, work, &lWork, iwork, &info);
#else
  dgelsd_( &m, &n, &nrhs, a, &lda, b, &ldb, s, rcond, rank, work, &lWork, iwork, &info);
#endif
#endif
  return info;
}

} // namespace lapack

} // namespace STK

#endif /* STK_LAPACK_MULTILEASTSQUARE_H */
