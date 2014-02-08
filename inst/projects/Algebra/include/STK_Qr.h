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
 * Purpose:  Define the Qr Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Qr.h
 *  @brief This file define the Qr class (QR decomposition).
 **/
 
#ifndef STK_QR_H
#define STK_QR_H

#include "STK_IQr.h"

namespace STK
{

/** @ingroup Algebra
 *  @brief The class Qr perform the QR decomposition of a Matrix.
 * 
 *  - Input:  A matrix (nrow,ncol)
 *  - Output:
 *    -# Q  matrix (nrow,ncol) with the Housholder vectors in the min(nrow, ncol) first columns.
 *    -# R  matrix (nrow,ncol) upper triangular.
 */
class Qr: public IQr
{
  public :
    /** Default constructor.
     *  @param A the matrix to decompose
     *  @param ref true if we overwrite A
     **/
    Qr( Matrix const&  A =Matrix(), bool ref = false);
    /** @brief Constructor
     *  @param data reference on a matrix expression
     */
    template<class Derived>
    Qr( ExprBase<Derived> const& data) : IQr(data) {}
    /** Copy constructor.
     *  @param decomp the decomposition  to copy
     **/
    Qr( Qr const& decomp);
    /** virtual destructor */
    virtual ~Qr();
    /** clone pattern */
    inline virtual Qr* clone() const { return new Qr(*this);}
    /** Operator = : overwrite the Qr with decomp. */
    Qr& operator=(Qr const& decomp);
    /** Compute the QR decomposition. **/
    virtual bool run();

  private:
  /** Compute the qr decoposition of the matrix Q_ */
    void qr();
};

} // namespace STK

#endif
// STK_QR_H

