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
 * Project:  stkpp::Arrays
 * created on: 30 déc. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ArrayByVectorProduct.h
 *  @brief In this file we implement the General Array by Vector product.
 **/


#ifndef STK_ARRAYBYVECTORPRODUCT_H
#define STK_ARRAYBYVECTORPRODUCT_H

namespace STK
{

namespace hidden
{

/** @ingroup hidden
 *  this structure regroup all the methods using only pointers on the Type
 **/
template<typename Type>
struct MultImpl
{
  /** multiplication of two vectors */
  static Type vectorByVector(Type const* p_lhs, Type const* p_rhs)
  {
    Type sum = Type(0);
    for (int k=0; k< vectorSize; ++k) sum += p_lhs[k] * p_rhs[k];
    return(sum);
  }
  /** multiplication of two vectors */
  static Type PanelByVector(Type const* p_lhs, Type const* p_rhs)
  {
    Type sum = Type(0);
    for (int k=0; k< vectorSize; ++k) sum += p_lhs[k] * p_rhs[k];
    return(sum);
  }
};

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
    for (int j= lhs.beginCols(); j<lhs.endCols(); ++j)
    {
      for (int i= lhs.beginRows(); i< lhs.endRows(); ++i)
      { res.elt(i) += lhs(i,j) * rhs[j];}
    }
    return;
//    // compute dimensions
//    int nbRawVec = rhs.size()/panelSize;
//    int nbPanels = lhs.sizeRows()/blockSize;
//    // create panels and blocks
//    Panel<Type>* tabPanel = new Panel<Type>[nbPanels];
//    RawVec<Type> vec;
//    // start panels by RawVec
//    for (int k = 0, kPos = rhs.begin(); k<nbRawVec; ++k, kPos += panelSize)
//    {
//      vectorToRawVec(rhs, vec, kPos);
//      for (int i = 0, iRow = lhs.beginRows(); i<nbPanels; ++i, iRow += blockSize)
//      { arrayToPanel( lhs, tabPanel[i], iRow, kPos);}
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//      for (int i = 0; i<nbPanels; ++i)
//      {
//        int iRow = lhs.beginRows() + i * blockSize;
//        PanelByRawVec(tabPanel[i],vec,res, iRow);
//      }
//      // last rows
//      for (int iRow=lhs.beginRows() + nbPanels * blockSize; iRow<rhs.endRows(); ++iRow)
//      { ArrayByRawVec(lhs,vec,res, iRow, kPos);}
//    }
//    delete[] tabPanel;
//    // remaining columns
//    for (int j= lhs.beginCols() + nbRawVec * panelSize; j<lhs.endCols(); ++j)
//    {
//      for (int i= lhs.beginRows(); i< lhs.endRows(); ++i)
//      { res.elt(i) += lhs(i,j) * rhs[j];}
//    }
  }
  /** Default dimension */
  static void vectorToRawVec( Rhs const& rhs, RawVec<Type>& vec, int iRow)
  {
    for (int j=0; j<panelSize; ++j)
    { vec[j] = rhs[iRow + j];}
  }
  /** Default dimension */
  static void arrayToPanel( Lhs const& lhs, Panel<Type>& panel, int iRow, int jCol)
  {
    for (int j=0; j<panelSize; ++j)
    { panel[j] = lhs(iRow  , jCol+j);}
    for (int j=0; j<panelSize; ++j)
    { panel[ panelSize+j] = lhs(iRow+1, jCol+j);}
    for (int j=0; j<panelSize; ++j)
    { panel[2*panelSize+j] = lhs(iRow+2, jCol+j);}
    for (int j=0; j<panelSize; ++j)
    { panel[3*panelSize+j] = lhs(iRow+3, jCol+j);}
  }
  /** Default dimension */
  static void PanelByRawVec( Panel<Type> const& panel, RawVec<Type> const& vec
                           , Result& res, int iRow)
  {
    for (int j=0; j<panelSize; ++j)
    { res.elt(iRow)   += panel[            j] * vec[j];}
    for (int j=0; j<panelSize; ++j)
    { res.elt(iRow+1) += panel[  panelSize+j] * vec[j];}
    for (int j=0; j<panelSize; ++j)
    { res.elt(iRow+2) += panel[2*panelSize+j] * vec[j];}
    for (int j=0; j<panelSize; ++j)
    { res.elt(iRow+3) += panel[3*panelSize+j] * vec[j];}
  }
  /** Default dimension */
  static void ArrayByRawVec( Lhs const& lhs, RawVec<Type> const& vec
                           , Result& res, int iRow, int jCol)
  {
    for (int j=0; j<panelSize; ++j)
    { res.elt(iRow) += lhs(iRow, jCol+j) * vec[j];}
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
    for (int j= rhs.beginCols(); j<rhs.endCols(); ++j)
    {
      for (int i= rhs.beginRows(); i< rhs.endRows(); ++i)
      { res.elt(j) += lhs[i] * rhs(i,j);}
    }
//    int nbInnerLoop = rhs.sizeRows()/vectorSize; // = rhs.sizeRows()/blockSize;
//    int vSize = rhs.sizeRows() - nbInnerLoop * vectorSize;
//    Type p_lhs[vectorSize];
//    Type p_rhs[vectorSize];
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//    for (int k = 0; k<nbInnerLoop; ++k)
//    {
//      int iPos = lhs.begin() + k * vectorSize;
//      for (int j=0; j<vectorSize; ++j) p_lhs[j]  = lhs.elt(iPos+j);
//      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
//      {
//        for (int i=0; i<vectorSize; ++i) p_rhs[i]  = rhs.elt(iPos+i, j);
//        res.elt(j) += Cmult::vectorByVector(p_lhs, p_rhs);
//      }
//    } // k loop
//    int iPos = lhs.begin()+vectorSize*nbInnerLoop;
//    for (int j=0; j<vSize; ++j) p_lhs[j]  = lhs.elt(iPos+j);
//    for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
//    {
//      for (int i=0; i<vSize; ++i) p_rhs[i]  = rhs.elt(iPos+i, j);
//      for (int k=0; k<vSize; ++k) res.elt(j) += p_lhs[k] * p_rhs[k];
//    } // j loop
  }
}; // struct pb

} // namespace hidden



} // namespace STK

#endif /* STK_ARRAYBYVECTORPRODUCT_H */
