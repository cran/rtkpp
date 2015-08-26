/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

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

/** @file STK_lapack_SymEigen.h
 *  @brief In this file we define the enclosing class of the syevr lapck routine.
 **/


#ifndef STK_LAPACK_SYMEIGEN_H
#define STK_LAPACK_SYMEIGEN_H

#include "STK_ISymEigen.h"

#ifdef STKUSELAPACK

extern "C"
{
#ifdef STKREALAREFLOAT
/** LAPACK routine in float to compute the eigenvalues */
extern void ssyevr_( char *, char *, char *, int *, float *, int *, float *,
                     float *, int *, int *, float *, int *, float *, float *, int *,
                    int *, float *, int *, int *, int *, int *);
#else
/** LAPACK routine in double to compute the eigenvalues */
extern void dsyevr_( char *, char *, char *, int *, double *, int *, double *,
                     double *, int *, int *, double *, int *, double *, double *, int *,
                    int *, double *, int *, int *, int *, int *);
#endif
}

#endif // STKUSELAPACK

namespace STK
{

// forward declaration
namespace lapack
{
template<class SquareArray> class SymEigen;
}

namespace hidden
{
/** @ingroup hidden
 *  Specialization for the Qr class.
 **
 **/
template<class SquareArray_>
struct AlgebraTraits< lapack::SymEigen<SquareArray_> >
{
  typedef SquareArray_ SquareArray;
};

} // namespace hidden


namespace lapack
{
/** @ingroup Algebra
 *  {
 *    @class SymEigen
 *    @brief  SymEigen computes the eigenvalues and optionally the
 *    eigenvectors of a symmetric real matrix using the syevr Lapack routine.
 */
template<class SquareArray>
class SymEigen : public ISymEigen<SymEigen<SquareArray> >
{
  public:
    typedef ISymEigen<SymEigen> Base;
    using Base::eigenValues_;
    using Base::eigenVectors_;
    using Base::range_;
    using Base::norm_;
    using Base::rank_;
    using Base::det_;
    /** @brief Constructor
     *  @param data reference on a symmetric square matrix
     *  @param ref @c true if we overwrite the data set, @c false otherwise
     *  @note data can be a reference if and only if it is a CSquareX
     */
    SymEigen( SquareArray const& data, bool ref =false);
    /** @brief Constructor
     *  @param data reference on a symmetric square expression
     */
    template<class Derived>
    SymEigen( ArrayBase<Derived> const& data)
            : Base(data)
            , jobz_('V'), RANGE_('A'), UPLO_('U')
            , VL_(0.0), VU_(0.0), IL_(0), IU_(0)
    {}
    /** @brief copy constructor
     *  @param eigen the SymEigen to copy
     */
    inline SymEigen( SymEigen const& eigen)
                   : Base(eigen)
                   , jobz_(eigen.jobz_), RANGE_(eigen.RANGE_), UPLO_(eigen.UPLO_)
                   , VL_(eigen.VL_), VU_(eigen.VU_), IL_(eigen.IL_), IU_(eigen.IU_)
     {}
    /** Destructor. */
    inline virtual ~SymEigen() {}
    /** @param jobz If jobz ='N': Compute eigenvalues only; If jobz = 'V': Compute
     * eigenvalues and eigenvectors.
     **/
    inline void setJobz(char jobz) { jobz_ = jobz;}
    /** @param range range of the eigenvalues to be found.
     *  If range = 'A': all eigenvalues will be found.
     *  If range = 'V': all eigenvalues in the half-open interval  (VL_,VU_] will be found.
     *  If range = 'I': the IL_-th through IU_-th eigenvalues will be found.
     **/
    inline void setRange(char range) { RANGE_ = range;}
    /** @param uplo values to used in A.
     * If uplo = 'U':  Upper triangle of A is stored;
     * If uplo = 'L':  Lower triangle of A is stored.
     **/
    inline void setUplo(char uplo) { UPLO_ = uplo;}
    /** @param[in] vl,vu lower and upper bounds of the interval to be searched
     * for eigenvalues.
     * Not referenced if RANGE_ = 'A' or 'I'.
    **/
    inline void setVlAndVu(Real const& vl, Real const& vu) { VL_ = vl; VU_ = vu;}
    /** @param il, iu
     *  If RANGE_='I', the indices (in ascending order)
     *  of the smallest and largest eigenvalues to be returned.
     *  1 <= IL_ <= IU_ <= NL, if NL > 0; IL_ = 1 and IU_ = 0 if NL = 0. Not
     *  referenced if RANGE_ = 'A' or 'V'.
    **/
    inline void setIlAndIu(int il, int iu) { IL_ = il; IU_ = iu;}
    /** @brief clone pattern */
    inline virtual SymEigen* clone() const { return new SymEigen(*this);}
    /** @brief Run eigenvalues decomposition
     *  Launch SYEVR LAPACK routine to perform the eigenvalues decomposition.
     *  @return @c true if no error occur, @c false otherwise
     */
    bool runImpl();

    /** wrapper of the LAPACK SYEVR routine. Compute the eigenvalues of a symmetric
     *  square matrix.
     *
     *  @param[in] jobz
     * @verbatim
     *  CHARACTER*1
     *  = 'N':  Compute eigenvalues only;
     *  = 'V':  Compute eigenvalues and eigenvectors.
     * @endverbatim
     *
     *  @param[in] range
     *  @verbatim
     *  CHARACTER*1
     *  = 'A': all eigenvalues will be found.
     *  = 'V': all eigenvalues in the half-open interval  (VL_,VU_]  will be found.
     *  = 'I': the IL_-th through IU_-th eigenvalues will be found.
     * @endverbatim
     *
     *  @param[in] uplo
     * @verbatim
     * CHARACTER*1
     * = 'U':  Upper triangle of A is stored;
     * = 'L':  Lower triangle of A is stored.
     * @endverbatim
     *
     * @param[in] n The order of the matrix A.  N >= 0.
     *
     * @param[in,out] a Real array, dimension (lda, N)
     * @verbatim
     * On entry, the symmetric matrix A.  If UPLO_ = 'U', the leading
     * N-by-N upper triangular part of A contains the upper  triangular part
     * of the  matrix  A.   If  UPLO_  = 'L', the leading N-by-N lower triangular
     * part of A contains the lower triangular part of the matrix A.  On
     * exit, the lower triangle (if UPLO_='L') or the upper triangle
     * (if UPLO_='U') of A, including the diagonal, is destroyed.
     * @endverbatim
     *
     * @param[in] lda The leading dimension of the array A.  lda >= max(1,N).
     *
     * @param[in] vl,vu
     * @verbatim
     *  Real If RANGE_='V', the lower and  upper  bounds
     *  of  the  interval to be searched for eigenvalues. VL_ < VU_.  Not
     *  referenced if RANGE_ = 'A' or 'I'.
     * @endverbatim
     *
     * @param[in] il, iu
     * @verbatim
     *  Integer If RANGE_='I', the indices (in ascending order)
     *  of the smallest and largest eigenvalues to be returned.
     *  1 <= IL_ <= IU_ <= NL, if NL > 0; IL_ = 1 and IU_ = 0 if NL = 0. Not
     *  referenced if RANGE_ = 'A' or 'V'.
     * @endverbatim
     *
     * @param[in] abstol
     * @verbatim
     *  The  absolute error tolerance for the eigenvalues.  An approximate
     *  eigenvalue is accepted as converged when it is  determined
     *  to lie in an interval [a,b] of width less than or equal to
     *  ABSTOL + EPS *   max( |a|,|b| ) ,
     *  where  EPS is the machine precision.  If ABSTOL is less than or
     *  equal to zero, then  EPS*|T|  will be used in its place,  where
     *  |T|  is the 1-norm of the tridiagonal matrix obtained by reducing A
     *  to tridiagonal form.
     *  If high relative accuracy is important, set ABSTOL  to  SLAMCH(
     *  'Safe minimum' ).  Doing so will guarantee that eigenvalues are
     *  computed to high relative  accuracy  when  possible  in  future
     *  releases.   The current code does not make any guarantees about
     *  high relative accuracy, but future releases will. See J. Barlow
     *  and J. Demmel, "Computing Accurate Eigensystems of Scaled Diagonally
     *  Dominant Matrices", LAPACK Working Note #7, for  a  discussion of
     *  which matrices define their eigenvalues to high relative accuracy.
     * @endverbatim
     *  @see "Computing Small Singular  Values  of  Bidiagonal  Matrices
     *  with  Guaranteed  High Relative Accuracy," by Demmel and Kahan,
     *  LAPACK Working Note #3.
     *
     * @param[out] m
     * @verbatim
     *   The  total number of eigenvalues found.  0 <= M <= NL.  If RANGE_
     *   = 'A', M = NL, and if RANGE_ = 'I', M = IU_-IL_+1.
     * @endverbatim
     *
     * @param[out] w
     * @verbatim
     *  array, dimension (NL)
     *  The first  M  elements  contain  the  selected  eigenvalues  in
     *  ascending order.
     * @endverbatim
     *
     * @param[out] z
     * @verbatim
     *  array, dimension (LDZ, max(1,M))
     *  If  jobz_ = 'V', then if info = 0, the first M columns of Z contain
     *  the orthonormal eigenvectors of the matrix A corresponding
     *  to  the selected eigenvalues, with the i-th column of Z holding
     *  the eigenvector associated with W(i).  If jobz_ = 'N', then Z is
     *  not  referenced.   Note:  the  user  must  ensure that at least
     *  max(1,M) columns are supplied in the array Z; if RANGE_  =  'V',
     *  the exact value of M is not known in advance and an upper bound
     *  must be used.  Supplying N columns is always safe.
     * @endverbatim
     *
     * @param[in] ldz
     * @verbatim
     *  The leading dimension of the array Z.  LDZ >= 1, and if jobz_  =
     *  'V', LDZ >= max(1,N).
     * @endverbatim
     *
     * @param[out] isuppz array, dimension ( 2*max(1,M) )
     * @verbatim
     *  The  support  of the eigenvectors in Z, i.e., the indices indicating
     *  the nonzero elements  in  Z.  The  i-th  eigenvector  is
     *  nonzero only in elements ISUPPZ( 2*i-1 ) through ISUPPZ( 2*i ).
     * @endverbatim
     *
     * @param[in,out] work Real array, dimension (MAX(1,lWork))
     * @verbatim
     *   On exit, if info = 0, work(1) returns the optimal lWork.
     * @endverbatim
     *
     * @param[in] lWork The  dimension  of  the array work
     * @verbatim
     *  lWork >= max(1,26*N). For optimal efficiency, lWork >= (NB+6)*N, where
     *  NB is the max of the  blocksize for SSYTRD and SORMTR returned by ILAENV.
     *  If lWork = -1, then a workspace query is assumed; the routine only
     *  calculates  the  optimal  sizes  of  the work and iWork arrays,
     *  returns these values as the first entries of the work and iWork
     *  arrays,  and  no  error  message  related to lWork or LiWork is
     *  issued by XERBLA.
     * @endverbatim
     *
     * @param[in,out] iwork array, dimension (MAX(1,LiWork))
     * @verbatim
     *  On exit, if info = 0, iWork(1) returns the optimal lWork.
     * @endverbatim
     *
     * @param[in] liwork The dimension of the array iWork.
     * @verbatim
     *  LiWork >=  max(1,10*N). If LiWork  =  -1, then a workspace query is
     *  assumed; the routine only calculates the optimal sizes of the work and
     *  iWork arrays, returns these values as the first entries of the work and
     *  iWork arrays, and no error message related  to  lWork  or  LiWork  is
     *  issued by XERBLA.
     * @endverbatim
     *
     * @return info
     * @verbatim
     *  = 0:  successful exit
     *  < 0:  if info = -i, the i-th argument had an illegal value
     *  > 0:  Internal error
     * @endverbatim
     */
    static int syevr( char jobz, char range, char uplo, int n
                    , Real* a, int lda
                    , Real vl, Real vu, int il, int iu
                    , Real abstol, int *m, Real *w
                    , Real *z, int ldz, int *isuppz
                    , Real *work, int lWork, int *iwork, int liwork
                    );
  private:
    /** Lapack pptions */
    char jobz_, RANGE_, UPLO_;
    Real VL_, VU_;
    int IL_, IU_;
};
/** @} */



/* @brief Constructor
 *  @param data reference on a symmetric square matrix
 *  @param ref @c true if we overwrite the data set, @c false otherwise
 */
template<class SquareArray>
inline SymEigen<SquareArray>::SymEigen( SquareArray const& data, bool ref)
                                      : Base(data)
                                      , jobz_('V'), RANGE_('A'), UPLO_('U')
                                      , VL_(0.0), VU_(0.0), IL_(0), IU_(0)
{}
/* @brief Constructor
 *  @param data reference on a symmetric square matrix
 *  @param ref @c true if we overwrite the data set, @c false otherwise
 */
template<>
inline SymEigen<CSquareX>::SymEigen( CSquareX const& data, bool ref)
                                   : Base(data, ref)
                                   , jobz_('V'), RANGE_('A'), UPLO_('U')
                                   , VL_(0.0), VU_(0.0), IL_(0), IU_(0)
{}

/* @brief Run eigen decomposition
 *  Launch SYEVR LAPACK routine to perform the eigenvalues decomposition.
 *  @return @c true if no error occur, @c false otherwise
 */
template<class SquareArray>
bool SymEigen<SquareArray>::runImpl()
{
#ifdef STK_ALGEBRA_VERY_VERBOSE
  stk_cout << _T("Enter in SymEigen::run\n");
#endif
  /* copy square matrix with the original data set. */
  CSquareX data_ = eigenVectors_;
  // shift data sets
  data_.shift(0);
  eigenVectors_.shift(0);
  eigenValues_.shift(0);
  this->SupportEigenVectors_.shift(0);
  /* set default behavior */
  Real absTol = 0.0; // let Lapack chose the correct tolerance
  // get optimal size necessary for work
  Real work; // work is just one place, get the optimal size for work
  int iwork; // iwork is just on place, get the optimal size for iWork
  int lWork =-1, liwork =-1; // workspace variable

  int info = 1;
#ifdef STK_ALGEBRA_DEBUG
  stk_cout << _T("Data dimensions: ") << data_.rows() << " " << data_.cols() << "\n";
  stk_cout << _T("eigenValues_ dimensions: ") << eigenValues_.rows() << " " << eigenValues_.cols() << "\n";
  stk_cout << _T("eigenVectors_ dimensions: ") << eigenVectors_.rows() << " " << eigenVectors_.cols() << "\n";
  stk_cout << _T("Options: ") << jobz_ << " " << RANGE_ << " " << UPLO_ << "\n";
#endif
  info = syevr( jobz_, RANGE_, UPLO_
              , range_.size(), data_.p_data(), range_.size()
              , VL_, VU_, IL_, IU_
              , absTol, &rank_,  eigenValues_.p_data()
              , eigenVectors_.p_data(), range_.size(), this->SupportEigenVectors_.p_data()
              , &work, lWork, &iwork, liwork);
  // check any error
  if (info!=0)
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(SymEigen::run,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(SymEigen::run,-info,error parameter);
    return false;
  }
#ifdef STK_ALGEBRA_DEBUG
  stk_cout << _T("Size needed:") << (int)work << " " << iwork << " " << "\n";
#endif
  // get results and allocate space
  lWork = (int)work;
  liwork = iwork;
  Real* p_work = new Real[lWork];
  int* p_iwork = new int[liwork];

  // Call SYEVR with the optimal block size
  info = syevr( jobz_, RANGE_, UPLO_
              , range_.size(), data_.p_data(), range_.size()
              , VL_, VU_, IL_, IU_
              , absTol, &rank_, eigenValues_.p_data()
              , eigenVectors_.p_data(), range_.size(), this->SupportEigenVectors_.p_data()
              , p_work, lWork, p_iwork, liwork);
  // recover memory
  delete[] p_work;
  delete[] p_iwork;

  // finalize
  data_.shift(range_.begin());
  eigenVectors_.shift(range_.begin());
//  stk_cout << _T("eigenValues_  (before shift) =\n") << eigenValues_ << "\n";
  eigenValues_.shift(range_.begin());
//  stk_cout << _T("eigenValues_  (after shift) =\n") << eigenValues_ << "\n";
  this->SupportEigenVectors_.shift(range_.begin());
  this->finalizeStep();
  // return the result of the computation
  if (!info) return true;
  if (info>0)
  { this->msg_error_ = STKERROR_NO_ARG(SymEigen ::run,internal error);
    return false;
  }
  this->msg_error_= STKERROR_1ARG(SymEigen ::run,-info,error parameter);
  return false;
}


/* wrapper of the LAPACK routine to compute the eigenvalues */
template<class SquareArray>
int SymEigen<SquareArray>::syevr( char jobz, char range, char uplo
                   , int n, Real* a, int lda
                   , Real vl, Real vu, int il, int iu
                   , Real abstol, int *m, Real *w
                   , Real *z, int ldz, int *isuppz
                   , Real *work, int lWork, int *iwork, int liwork
                   )
{
  int info = 0;
#ifdef STKUSELAPACK
#ifdef STKREALAREFLOAT
  ssyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il,
          &iu, &abstol, m, w, z, &ldz, isuppz, work,
          &lWork, iwork, &liwork, &info);
#else
  dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il,
          &iu, &abstol, m, w, z, &ldz, isuppz, work,
          &lWork, iwork, &liwork, &info);
#endif
#endif
  return info;
}


} // namespace lapack

} // namespace STK

#endif /* STK_LAPACK_SYMEIGEN_H */
