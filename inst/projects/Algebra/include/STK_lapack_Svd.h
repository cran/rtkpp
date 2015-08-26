/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015  Serge Iovleff

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR a PARTICULAR PURPOSE.  See the
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
 *  @brief In this file we define the enclosing class of the dgeqrf lapack routine.
 **/


#ifndef STK_LAPACK_SVD_H
#define STK_LAPACK_SVD_H

#include <Arrays/include/STK_CArray.h>
#include <Arrays/include/STK_CArrayVector.h>
#include "STK_ISvd.h"

#ifdef STKUSELAPACK

extern "C"
{
#ifdef STKREALAREFLOAT
/* LAPACK routine in float to compute the SVD decomposition */
void sgesdd_( char *jobz, int *M, int *N, float *A, int *lda, float *S, float *U, int* ldu
            , float *vt, int *ldvt, float *work, int *lWork, int *iwork, int *info);
#else
/* LAPACK routine in double to compute the SVD decomposition */
void dgesdd_( char *jobz, int *M, int *N, double *A, int *lda, double *S, double *U, int* ldu
            , double *vt, int *ldvt, double *work, int *lWork, int *iwork, int *info);
#endif
}

#endif // STKUSELAPACK

namespace STK
{

namespace lapack
{
class Svd;
}

namespace hidden
{
/** @ingroup hidden
 *  Specialization for the Svd class.
 **/
template<>
struct AlgebraTraits< lapack::Svd >
{
  typedef CArrayXX ArrayU;
  typedef CVectorX ArrayD;
  typedef CArrayXX ArrayV;
};

} // namespace hidden


namespace lapack
{
/** @ingroup Algebra
 *  {
 *    @class Svd
 *    @brief Svd computes the SVD decomposition of a real matrix using the
 *    Lapack routine dgeqrf.
 */
class Svd : public ISvd<Svd>
{
  public:
    typedef ISvd< Svd > Base;
    typedef CArrayXX::Col ColVector;
    typedef CArrayXX::Row RowVector;
    using Base::U_;
    using Base::D_;
    using Base::V_;
    using Base::withU_;
    using Base::withV_;
    using Base::nrowU;
    using Base::ncolU;
    using Base::ncolV;
    using Base::norm_;
    using Base::rank_;

    /** Default constructor.
     *  @param A the matrix to decompose
     *  @param ref if true, U_ is a reference of A.
     *  @param withU if @c true save the left housolder transforms in @c U_.
     *  @param withV if @c true save the right housolder transforms in @c V_.
     **/
    inline Svd( CArrayXX const&  A, bool ref = false, bool withU= true, bool withV= true)
              : Base(A, ref, withU, withV)
              , jobz_( (withU|withV) ? 'O':'N') {}
    /** constructor with other kind of array/expression
     *  @param A the matrix/expression to decompose.
     *  @param withU if @c true save the left housolder transforms in @c U_.
     *  @param withV if @c true save the right housolder transforms in @c V_.
     */
    template<class OtherArray>
    inline Svd( ArrayBase<OtherArray> const& A, bool withU= true, bool withV= true)
              : Base(A, withU, withV)
              , jobz_( (withU|withV) ? 'O':'N') {}
    /** Copy constructor.
     *  @param decomp the decomposition to copy
     **/
    inline Svd( Svd const& decomp): Base(decomp), jobz_(decomp.jobz_) {}
    /** virtual destructor */
    inline virtual ~Svd() {}
    /** @return the option chosen for the svd */
    char jobz() const { return jobz_;}
    /** set the option chosen for the svd */
    void setJobz(char jobz) { jobz_ = jobz;}
    /** @brief clone pattern */
    inline virtual Svd* clone() const { return new Svd(*this);}
    /** Operator = : overwrite the Svd with decomp. */
    inline Svd& operator=(Svd const& decomp)
    {
      Base::operator=(decomp);
      jobz_ = decomp.jobz_;
      return *this;
    }
    /** @brief Run svd decomposition */
    bool runImpl();

  protected:
    /** wrapper of the LAPACK DGESDD routine. Compute the Svd decomposition
     *  of a matrix.
     *
     * @verbatim
     * DGESDD computes the singular value decomposition (SVD) of a real
     * m-by-n matrix a, optionally computing the left and right singular
     * vectors.  If singular vectors are desired, it uses a
     * divide-and-conquer algorithm.
     *
     * The SVD is written
     *
     *      a = u * SIGMA * transpose(v)
     *
     * where SIGMA is an m-by-n matrix which is zero except for its
     * min(m,n) diagonal elements, u is an m-by-m orthogonal matrix, and
     * v is an n-by-n orthogonal matrix. The diagonal elements of SIGMA
     * are the singular values of a; they are real and non-negative, and
     * are returned in descending order. The first min(m,n) columns of
     * u and v are the left and right singular vectors of a.
     *
     * Note that the routine returns vt = v**T, not v.
     *
     * @endverbatim
     *
     *  Arguments:
     *  ==========
     *
     * @param[in] jobz
     * @verbatim
     *          jobz is Char*1
     *          Specifies options for computing all or part of the matrix u:
     *          = 'A':  all m columns of u and all n rows of v**T are
     *                  returned in the arrays u and vt;
     *          = 'S':  the first min(m,n) columns of u and the first
     *                  min(m,n) rows of v**T are returned in the arrays U
     *                  and vt;
     *          = 'O':  If m >= n, the first n columns of u are overwritten
     *                  on the array a and all rows of v**T are returned in
     *                  the array vt;
     *                  otherwise, all columns of u are returned in the
     *                  array u and the first m rows of v**T are overwritten
     *                  in the array a;
     *          = 'N':  no columns of u or rows of v**T are computed.
     * @endverbatim
     *
     * @param[in] m
     * @verbatim
     *          m is Integer
     *          The number of rows of the input matrix a.  m >= 0.
     * @endverbatim
     *
     * @param[in] n
     * @verbatim
     *          n is Integer
     *          The number of columns of the input matrix a.  n >= 0.
     * @endverbatim
     *
     * @param[in,out] a
     * @verbatim
     *          a is STK::Real array, dimension (lda,n)
     *          On entry, the m-by-n matrix a.
     *          On exit,
     *          if jobz = 'O',  a is overwritten with the first n columns
     *                          of u (the left singular vectors, stored
     *                          columnwise) if m >= n;
     *                          a is overwritten with the first m rows
     *                          of v**T (the right singular vectors, stored
     *                          rowwise) otherwise.
     *          if jobz .ne. 'O', the contents of a are destroyed.
     * @endverbatim
     *
     * @param[in] lda
     * @verbatim
     *          lda is Integer
     *          The leading dimension of the array a.  lda >= max(1,m).
     * @endverbatim
     *
     * @param[out] s
     * @verbatim
     *          s is STK::Real array, dimension (min(m,n))
     *          The singular values of a, sorted so that s[i] >= s[i+1].
     * @endverbatim
     *
     * @param[out] u
     * @verbatim
     *          u is STK::Real array, dimension (ldu,ucol)
     *          ucol = m if jobz = 'A' or jobz = 'O' and m < n;
     *          ucol = min(m,n) if jobz = 'S'.
     *          If jobz = 'A' or jobz = 'O' and m < n, u contains the m-by-m
     *          orthogonal matrix u;
     *          if jobz = 'S', u contains the first min(m,n) columns of u
     *          (the left singular vectors, stored columnwise);
     *          if jobz = 'O' and m >= n, or jobz = 'N', u is not referenced.
     * @endverbatim
     *
     * @param[in] ldu
     * @verbatim
     *          ldu is Integer
     *          The leading dimension of the array U.  ldu >= 1; if
     *          jobz = 'S' or 'A' or jobz = 'O' and m < n, ldu >= m.
     * @endverbatim
     *
     * @param[out] vt
     * @verbatim
     *          vt is STK::Real array, dimension (ldvt,n)
     *          If jobz = 'A' or jobz = 'O' and m >= n, vt contains the
     *          N-by-N orthogonal matrix v**T;
     *          if jobz = 'S', vt contains the first min(m,n) rows of
     *          v**T (the right singular vectors, stored rowwise);
     *          if jobz = 'O' and m < n, or jobz = 'N', vt is not referenced.
     * @endverbatim
     *
     * @param[in] ldvt
     * @verbatim
     *          ldvt is Integer
     *          The leading dimension of the array vt.  ldvt >= 1; if
     *          jobz = 'A' or jobz = 'O' and m >= n, ldvt >= n;
     *          if jobz = 'S', ldvt >= min(m,n).
     * @endverbatim
     *
     * @param[out] work
     * @verbatim
     *          work is STK::Real array, dimension (MAX(1,lWork))
     *          On exit, if info = 0, work(1) returns the optimal lWork;
     * @endverbatim
     *
     * @param[in] lWork
     * @verbatim
     *          lWork is Integer
     *          The dimension of the array work. lWork >= 1.
     *          If jobz = 'N',
     *            lWork >= 3*min(m,n) + max(max(m,n),7*min(m,n)).
     *          If jobz = 'O',
     *            lWork >= 3*min(m,n) +
     *                     max(max(m,n),5*min(m,n)*min(m,n)+4*min(m,n)).
     *          If jobz = 'S' or 'A'
     *            lWork >= min(m,n)*(6+4*min(m,n))+max(m,n)
     *          For good performance, lWork should generally be larger.
     *          If lWork = -1 but other input arguments are legal, work(1)
     *          returns the optimal lWork.
     * @endverbatim
     *
     * @param[out] iWork
     * @verbatim
     *          iWork is Integer array, dimension (8*min(m,n))
     * @endverbatim
     *
     * @return info
     * @verbatim
     *          info is Integer
     *          = 0:  successful exit.
     *          < 0:  if info = -i, the i-th argument had an illegal value.
     *          > 0:  DBDSDC did not converge, updating process failed.
     * @endverbatim
     *
     *  Authors:
     *  ========
     *
     * @author Univ. of Tennessee
     * @author Univ. of California Berkeley
     * @author Univ. of Colorado Denver
     * @author NAG Ltd.
     *
     * Contributors:
     * ============
     *
     *     Ming Gu and Huan Ren, Computer Science Division, University of
     *     California at Berkeley, USA
     *
     **/
    static int gesdd( char jobz, int m, int n
                    , Real *a, int lda, Real *s, Real *u, int ldu, Real *vt, int ldvt
                    , Real *work, int lWork, int *iWork
                    );
  private:
     /** option */
     char jobz_;
     /** compute the svd decomposition. a contains either u, vt (if jobz_=='O')
      *  or is destroyed at the end of the oputput. */
     bool computeSvd(CArrayXX& a, CArrayXX& u, CVectorX& s, CArrayXX& v);
};


/* @} */

} // namespace lapack

} // namespace STK

#endif /* STK_LAPACK_SVD_H */
