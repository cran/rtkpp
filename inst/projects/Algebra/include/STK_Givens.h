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
 * Purpose:  Define Givens rotation methods for Algebra classes.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Givens.h
 *  @brief In this file we define Givens methods used by the Algebra
 *  classes.
 **/
 
#ifndef STK_GIVENS_H
#define STK_GIVENS_H

#include "STKernel/include/STK_Real.h"
#include "Arrays/include/STK_ArrayBase.h"

namespace STK
{
/** @ingroup Algebra
 *  @brief Compute Givens rotation.
 * 
 * Compute the Givens rotation
 * \f[
 * \begin{bmatrix}
 *  c & s \\
 * -s & c
 * \end{bmatrix}
 * \begin{bmatrix}
 *  a \\
 *  b
 * \end{bmatrix}
 * =
 * \begin{bmatrix}
 *  r \\
 *  0
 * \end{bmatrix}.
 * \f]
 *  in order to eliminate the coefficient z and
 * return the value r of the rotated element.
 * 
 *  @see http://en.wikipedia.org/wiki/Givens_rotation
 * 
 * @param y The coefficient to rotate (input)
 * @param z the coefficient to eliminate (input)
 * @param cosinus the cosinus of the Givens rotation (output)
 * @param sinus the sinus of the Givens rotation rotation (output)
 **/
Real compGivens( Real const& y, Real const& z, Real& cosinus, Real& sinus);

/** @ingroup Algebra
 *  @brief Apply Givens rotation.
 * 
 * Perform a right multiplication of the Container M with a
 * Givens Matrix on the col1 and col2. col1 should be less than col2.
 * The Matrix M is passed as const as we are using reference
 * on the two cols we want to rotate.
 * 
 * @see http://en.wikipedia.org/wiki/Givens_rotation
 * 
 * @param M the Container to multiply
 * @param j1,j2 the first and second columns
 * @param cosinus,sinus the cosinus and sinus of the givens rotation
 **/
template < class TContainer2D>
void rightGivens( ArrayBase<TContainer2D>& M, int j1, int j2
                , Real const& cosinus, Real const& sinus
                )
{
  // Apply givens rotation
  for (int i = M.beginRows(); i < M.endRows(); i++)
  {
    const Real aux1 = M.elt(i, j1), aux2 = M.elt(i, j2);
    M.elt(i, j1) = cosinus * aux1 + sinus * aux2;
    M.elt(i, j2) = cosinus * aux2 - sinus * aux1;
  }
}

/** @ingroup Algebra
 *  @brief left multiplication by a Givens Matrix.
 * 
 * Perform a left multiplication of the matrix M with a Givens Matrix
 * on the row1 and row2. row1 should be less than row2.
 * The Matrix M is passed as const as we are using reference
 * on the two rows we want to rotate.
 * 
 * @see http://en.wikipedia.org/wiki/Givens_rotation
 * 
 * @param M the matix to multiply
 * @param i1,i2 the first and second rows
 * @param cosinus,sinus the cosinus and sinus of the givens rotation
 **/
template < class TContainer2D>
void leftGivens( ArrayBase<TContainer2D>& M, int i1, int i2
               , Real const& cosinus, Real const& sinus
               )
{
  // apply left Givens rotation
  for (int j = M.beginCols(); j< M.endCols(); j++)
  {
    const Real aux1 = M.elt(i1, j), aux2 = M.elt(i2, j);
    M.elt(i1, j) = cosinus * aux1 + sinus * aux2;
    M.elt(i2, j) = cosinus * aux2 - sinus * aux1;
  }
}

} // namespace STK

#endif /*STK_GIVENS_H*/
