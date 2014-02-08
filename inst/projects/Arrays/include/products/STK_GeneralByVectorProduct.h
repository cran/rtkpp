/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2012  Serge Iovleff

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
 * Project:  stkpp::
 * created on: 30 déc. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_GeneralByVectorProduct.h
 *  @brief In this file we implement the General Matrix by Vector product.
 **/


#ifndef STK_GENERALBYVECTORPRODUCT_H
#define STK_GENERALBYVECTORPRODUCT_H

namespace STK
{

namespace hidden
{
/** @ingroup hidden
 *  Methods to use for C=AB with A a general matrix and B a vector.
 *  The structure bv contains only static methods and typedef and should normally
 *  not be used directly.
 **/
template<typename Lhs, typename Rhs, typename Result>
struct bv
{
  typedef typename Result::Type Type;
  typedef hidden::MultImpl<Type> Cmult;
  /** Main method for Matrices by vector multiplication implementation
   *  res have been resized and initialized to zero outside
   *  this method.
   **/
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int nbInnerLoop = lhs.sizeCols()/vectorSize; // = rhs.sizeRows()/blockSize;
    int vSize = lhs.sizeCols() - nbInnerLoop * vectorSize;
    Type* p_lhs  = new Type[vectorSize];
    Type* p_rhs  = new Type[vectorSize];
    for (int k=0, iPos = lhs.beginCols(); k< nbInnerLoop; ++k, iPos+= vectorSize)
    {
      for (int j=0; j<vectorSize; j++) p_rhs[j]  = rhs.elt(iPos+j);
      for (int iRow=lhs.beginRows(); iRow< lhs.endRows(); ++iRow)
      {
        for (int k=0; k<vectorSize; k++) p_lhs[k]  = lhs.elt(iRow, iPos+k);
        res.elt(iRow) += Cmult::vectorByVector(p_lhs, p_rhs);
      }
    }
    int iPos= lhs.beginCols()+vectorSize*nbInnerLoop;
    for (int j=0; j<vSize; j++) p_rhs[j]  = rhs.elt(iPos+j);
    for (int iRow=lhs.beginRows(); iRow< lhs.endRows(); ++iRow)
    {
      for (int j=0; j<vSize; j++) p_lhs[j]  = lhs.elt(iRow, iPos+j);
      res.elt(iRow) += Cmult::vectorByVector(p_lhs, p_rhs, vSize);
    }
    delete[] p_lhs;
    delete[] p_rhs;
  }
}; // struct bv

/** @ingroup hidden
 *  Methods to use for C=AB with A a point and B a matrix.
 *  The structure vb contains only static method and typedef and should normally
 *  not be used directly.
 **/
template<typename Lhs, typename Rhs, typename Result>
struct vb
{
  typedef typename Result::Type Type;
  typedef hidden::MultImpl<Type> Cmult;
  /** Main method for point by Matrices multiplication implementation.
   *  res have been resized and initialized to zero outside
   *  this method.
   **/
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    int nbInnerLoop = rhs.sizeRows()/vectorSize; // = rhs.sizeRows()/blockSize;
    int vSize = rhs.sizeRows() - nbInnerLoop * vectorSize;
    Type* p_lhs  = new Type[vectorSize];
    Type* p_rhs  = new Type[vectorSize];
    for (int k = 0, iPos = lhs.begin(); k<nbInnerLoop; ++k, iPos+=vectorSize)
    {
      for (int j=0; j<vectorSize; ++j) p_lhs[j]  = lhs.elt(iPos+j);
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      {
        for (int i=0; i<vectorSize; ++i) p_rhs[i]  = rhs.elt(iPos+i, j);
        res.elt(j) += Cmult::vectorByVector(p_lhs, p_rhs);
      }
    } // k loop
    int iPos = lhs.begin()+vectorSize*nbInnerLoop;
    for (int j=0; j<vSize; ++j) p_lhs[j]  = lhs.elt(iPos+j);
    for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
    {
      for (int i=0; i<vSize; ++i) p_rhs[i]  = rhs.elt(iPos+i, j);
      res.elt(j) += Cmult::vectorByVector(p_lhs, p_rhs, vSize);
    } // j loop
    delete[] p_lhs;
    delete[] p_rhs;
  }
}; // struct pb

} // namespace hidden



} // namespace STK

#endif /* STK_GENERALBYVECTOR_H */
