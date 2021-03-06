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


#ifndef STK_LAPACK_QR_H
#define STK_LAPACK_QR_H

#include "STK_IQr.h"
#include <Arrays/include/STK_CArray.h>
#include <Arrays/include/STK_CArrayVector.h>
#include <Arrays/include/STK_Array2D.h>

#ifdef STKUSELAPACK

extern "C"
{
#ifdef STKREALAREFLOAT
/** LAPACK routine in float to compute the QR decomposition */
extern void sgeqrf_(int* M, int* N, float* A, int* lda, float* TAU, float* work, int* lWork, int* info );
#else
/** LAPACK routine in double to compute the QR decomposition */
extern void dgeqrf_(int* M, int* N, double* A, int* lda, double* TAU, double* work, int* lWork, int* info );
#endif
}

#endif // STKUSELAPACK

namespace STK
{

namespace lapack
{
class Qr;
}

namespace hidden
{
/** @ingroup hidden
 *  Specialization for the Qr class.
 **
 **/
template<>
struct AlgebraTraits< lapack::Qr >
{
  typedef ArrayXX Array;
};

} // namespace hidden


namespace lapack
{
/** @ingroup Algebra
 *  {
 *    @class Qr
 *    @brief Qr computes the QR decomposition of a real matrix using the
 *    Lapack routine dgeqrf.
 */
class Qr : public IQr<Qr >
{
  public:
    typedef IQr<Qr > Base;
    using Base::Q_;
    using Base::R_;
    /** Default constructor.
     *  @param data the matrix to decompose
     *  @param ref true if we overwrite A
     **/
    inline Qr( ArrayXX const&  data, bool ref = false): Base(data, ref) {}
    /** @brief Constructor
     *  @param data reference on a matrix expression
     */
    template<class Derived>
    Qr( ArrayBase<Derived> const& data): Base(data){}
    /** Copy constructor.
     *  @param decomp the decomposition  to copy
     **/
    inline Qr( Qr const& decomp): Base(decomp) {}
    /** virtual destructor */
    inline virtual ~Qr() {}
    /** @brief clone pattern */
    inline virtual Qr* clone() const { return new Qr(*this);}
    /** Operator = : overwrite the Qr with decomp. */
    inline Qr& operator=(Qr const& decomp)
    {
      Base::operator =(decomp);
      return *this;
    }
    /** @brief Run qr decomposition
     *  Launch geqrf LAPACK routine to perform the qr decomposition.
     *  @return @c true if no error occur, @c false otherwise
     */
    bool runImpl();
    /** wrapper of the LAPACK DGEQRF routine. Compute the Qr decomposition
     *  of a matrix.
     *
     * @param[in] m The number of rows of the matrix A.  M >= 0.
     *
     * @param[in] n The number of columns of the matrix A.  N >= 0.
     *
     * @param[in,out] a Real array, dimension (lda, N)
     * \verbatim
     *     On entry, the M-by-N matrix A.
     *     On exit, the elements on and above the diagonal of the array
     *     contain the min(M,N)-by-N upper trapezoidal matrix R (R is
     *     upper triangular if m >= n); the elements below the diagonal,
     *     with the array TAU, represent the orthogonal matrix Q as a
     *     product of min(m,n) elementary reflectors (see Further Details).
     * \endverbatim
     *
     * @param[in] lda The leading dimension of the array A.  lda >= max(1,M).
     *
     * @param[out] tau Real array, dimension min(M,N)
     * The scalar factors of the elementary reflectors (see Further Details).
     *
     * @param[in,out] work Real array, dimension (MAX(1,lWork))
     * \verbatim
     *   On exit, if info = 0, work(1) returns the optimal lWork.
     * \endverbatim
     *
     * @param[in] lWork The  dimension  of  the array work
     * \verbatim
     *  lWork >= max(1,N).
     *  For optimum performance lWork >= N*NB, where NB is the optimal blocksize.
     *
     *  If lWork = -1, then a workspace query is assumed; the routine
     *  only calculates the optimal size of the work array, returns
     *  this value as the first entry of the work array, and no error
     *  message related to lWork is issued by XERBLA.
     * \endverbatim
     *
     * @return info
     * \verbatim
     *  = 0:  successful exit
     *  < 0:  if info = -i, the i-th argument had an illegal value
     * \endverbatim
     *
     * @verbatim
     *  Further Details
     *  ===============
     *
     *  The matrix Q is represented as a product of elementary reflectors
     *
     *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
     *
     *  Each H(i) has the form
     *
     *     H(i) = I - tau * v * v'
     *
     *  where tau is a real scalar, and v is a real vector with
     *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
     *  and tau in TAU(i).
     * @endverbatim
     */
    static int geqrf(int m, int n, Real* a, int lda, Real* tau, Real *work, int lWork);

  private:
    /** private method for computing the Qr decomposition using a CArrayXX array */
    bool computeQr(CArrayXX& a, CVectorX& tau);
};
/** @} */


} // namespace lapack

} // namespace STK

#endif /* STK_LAPACK_QR_H */
