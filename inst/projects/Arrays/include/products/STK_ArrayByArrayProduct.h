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

/** @file STK_ArrayByArrayProduct.h
 *  @brief In this file we implement the General Array by Array product.
 **/


#ifndef STK_ARRAYBYARRAYPRODUCT_H
#define STK_ARRAYBYARRAYPRODUCT_H

namespace STK
{
namespace hidden
{
template<typename Lhs, typename Rhs, typename Result, bool Orient_> struct bp;

/** @ingroup hidden
 *  Methods to use for C=AB with A divided in blocks and B divided in panels.
 *  The structure bp contains only static method and typedef and should normally
 *  not be used directly.
 **/
template<typename Lhs, typename Rhs, typename Result>
struct bp<Lhs, Rhs, Result, (bool)Arrays::by_row_>
{
  typedef typename Result::Type Type;
  typedef hidden::MultImpl<Type> Cmult;
  typedef hidden::MultCoefImpl<Lhs, Rhs, Result> MultCoeff;

  /** Main method for Matrices multiplication implementation.
   *  @note res have been resized and initialized to zero outside this method.
   **/
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res)
  {
#ifdef STK_DEBUG
     stk_cout << _T("Entering bp::run()\n");
#endif
    // check small sizes
    switch (lhs.sizeRows())
    {
      case 0:
        return; break;
      case 1:
        MultCoeff::mult1RowOuterCol(lhs, rhs, res, lhs.beginRows());
        return; break;
      case 2:
        MultCoeff::mult2RowOuterCol(lhs, rhs, res, lhs.beginRows());
        return; break;
      case 3:
        MultCoeff::mult3RowOuterCol(lhs, rhs, res, lhs.beginRows());
        return; break;
      default: break;
    }
    switch (rhs.sizeCols())
    {
      case 0:
        return; break;
      case 1:
        MultCoeff::mult1ColOuterCol(lhs, rhs, res, rhs.beginCols());
        return; break;
      case 2:
        MultCoeff::mult2ColOuterCol(lhs, rhs, res, rhs.beginCols());
        return; break;
      case 3:
        MultCoeff::mult3ColOuterCol(lhs, rhs, res, rhs.beginCols());
        return; break;
      default: break;
    }
    // compute dimensions
    int nbInnerLoop = lhs.sizeCols()/blockSize; // = rhs.sizeRows()/blockSize;
    int nbBlocks    = lhs.sizeRows()/blockSize;
    int nbPanels    = rhs.sizeCols()/panelSize;
    // remaining sizes in the matrices
    int pSize = rhs.sizeCols() - panelSize*nbPanels;
    int bSize = lhs.sizeRows() - blockSize*nbBlocks;
    int tSize = lhs.sizeCols() - blockSize*nbInnerLoop;
              // = rhs.sizeRows() -  rhs.sizeRows()/blockSize
    int iLastRow = lhs.beginRows() + nbBlocks * blockSize;
    int jLastCol = rhs.beginCols() + nbPanels * panelSize;
    int kLastPos = lhs.beginCols() + blockSize * nbInnerLoop;
    if (nbInnerLoop)
    {
      // create panels and blocks
      Panel<Type>* tabPanel = new Panel<Type>[nbPanels+1];
      Block<Type>* tabBlock = new Block<Type>[nbBlocks+1];
      // start blocks by panel
      for (int k = 0; k<nbInnerLoop; ++k)
      {
        int kPos = lhs.beginCols() + k * blockSize;
        for (int i = 0; i<nbBlocks; ++i)
        {
          int iRow = lhs.beginRows() + i * blockSize;
          arrayToBlock( lhs, tabBlock[i], iRow, kPos);
        }
        arrayToBlock( lhs, tabBlock[nbBlocks], iLastRow, kPos, bSize);
        for (int j = 0; j<nbPanels; ++j)
        {
          int jCol = rhs.beginCols() + j*panelSize;
          arrayToPanel( rhs, tabPanel[j], kPos, jCol);
        }
        arrayToPanel( rhs, tabPanel[nbPanels], kPos, jLastCol, pSize);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i<nbBlocks; ++i)
        {
          int iRow = lhs.beginRows() + i * blockSize;
          for (int j = 0; j<nbPanels; ++j)
          {
            int jCol = rhs.beginCols() + j * panelSize;
            blockByPanel( tabBlock[i], tabPanel[j], res, iRow, jCol);
          }
        }
        for (int i = 0; i<nbBlocks; ++i)
        {
          int iRow = lhs.beginRows() + i * blockSize;
          blockByPanel( tabBlock[i], tabPanel[nbPanels], res, iRow, jLastCol, pSize);
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int j = 0; j<nbPanels; ++j)
        {
          int jCol = rhs.beginCols() + j * panelSize;
          blockByPanel( tabBlock[nbBlocks], tabPanel[j], res, iLastRow, jCol, panelSize, bSize);
        }
        blockByPanel( tabBlock[nbBlocks], tabPanel[nbPanels], res, iLastRow, jLastCol, pSize, bSize);
      } // InnerLoop
      delete[] tabBlock;
      delete[] tabPanel;
    } // if IneerLoop
    // treat the remaining rows, columns
    switch (tSize)
    {
      case 1:
        MultCoeff::mult1Col(lhs, rhs, res, kLastPos, kLastPos);
        break;
      case 2:
        MultCoeff::mult2Col(lhs, rhs, res, kLastPos, kLastPos);
        break;
      case 3:
        MultCoeff::mult3Col(lhs, rhs, res, kLastPos, kLastPos);
        break;
      default:
        break;
    }
  }
  /** Default size block */
  static void arrayToBlock( Lhs const& lhs, Block<Type>& block, int iRow, int jCol)
  {
     block[0]  = lhs.elt(iRow, jCol);
     block[1]  = lhs.elt(iRow, jCol+1);
     block[2]  = lhs.elt(iRow, jCol+2);
     block[3]  = lhs.elt(iRow, jCol+3);
     block[4]  = lhs.elt(iRow+1, jCol);
     block[5]  = lhs.elt(iRow+1, jCol+1);
     block[6]  = lhs.elt(iRow+1, jCol+2);
     block[7]  = lhs.elt(iRow+1, jCol+3);
     block[8]  = lhs.elt(iRow+2, jCol);
     block[9]  = lhs.elt(iRow+2, jCol+1);
     block[10] = lhs.elt(iRow+2, jCol+2);
     block[11] = lhs.elt(iRow+2, jCol+3);
     block[12] = lhs.elt(iRow+3, jCol);
     block[13] = lhs.elt(iRow+3, jCol+1);
     block[14] = lhs.elt(iRow+3, jCol+2);
     block[15] = lhs.elt(iRow+3, jCol+3);
  }
  /** with block row-size given */
  static void arrayToBlock( Lhs const& lhs, Block<Type>& block, int iRow, int jCol, int bSize)
  {
    for (int i=0; i<bSize; ++i)
    {
      block[i*blockSize]    = lhs.elt(iRow+i, jCol);
      block[i*blockSize+1]  = lhs.elt(iRow+i, jCol+1);
      block[i*blockSize+2]  = lhs.elt(iRow+i, jCol+2);
      block[i*blockSize+3]  = lhs.elt(iRow+i, jCol+3);
    }
  }
  /** Default dimension */
  static void arrayToPanel( Rhs const& rhs, Panel<Type>& panel, int iRow, int jCol)
  {
    for (int j=0; j<panelSize; ++j)
    {
      panel[j*blockSize]   = rhs.elt(iRow,   jCol+j);
      panel[j*blockSize+1] = rhs.elt(iRow+1, jCol+j);
      panel[j*blockSize+2] = rhs.elt(iRow+2, jCol+j);
      panel[j*blockSize+3] = rhs.elt(iRow+3, jCol+j);
    }
  }
  /** with panel size given */
  static void arrayToPanel( Rhs const& rhs, Panel<Type>& panel, int iRow, int jCol, int pSize)
  {
    for (int j=0; j<pSize; ++j)
    {
      panel[j*blockSize]   = rhs.elt(iRow,   jCol+j);
      panel[j*blockSize+1] = rhs.elt(iRow+1, jCol+j);
      panel[j*blockSize+2] = rhs.elt(iRow+2, jCol+j);
      panel[j*blockSize+3] = rhs.elt(iRow+3, jCol+j);
    }
  }
  /** Default dimension */
  static void blockByPanel( Block<Type> const& block, Panel<Type> const& panel
                          , Result& res, int iRow, int jCol)
  {
    for (int j=0; j<panelSize; ++j)
    {
      res.elt(iRow  ,jCol+j) += panel[j*blockSize]    * block[0]
                              + panel[j*blockSize+ 1] * block[1]
                              + panel[j*blockSize+ 2] * block[2]
                              + panel[j*blockSize+ 3] * block[3];
      res.elt(iRow+1,jCol+j) += panel[j*blockSize]    * block[4]
                              + panel[j*blockSize+ 1] * block[5]
                              + panel[j*blockSize+ 2] * block[6]
                              + panel[j*blockSize+ 3] * block[7];
      res.elt(iRow+2,jCol+j) += panel[j*blockSize]    * block[8]
                              + panel[j*blockSize+ 1] * block[9]
                              + panel[j*blockSize+ 2] * block[10]
                              + panel[j*blockSize+ 3] * block[11];
      res.elt(iRow+3,jCol+j) += panel[j*blockSize]    * block[12]
                              + panel[j*blockSize+ 1] * block[13]
                              + panel[j*blockSize+ 2] * block[14]
                              + panel[j*blockSize+ 3] * block[15];
    }
  }
  /** with panel size given */
  static void blockByPanel( Block<Type> const& block, Panel<Type> const& panel
                          , Result& res, int iRow, int jCol, int pSize)
  {
    for (int j=0; j<pSize; ++j)
    {
      res.elt(iRow  ,jCol+j) += panel[j*blockSize]    * block[0]
                              + panel[j*blockSize+ 1] * block[1]
                              + panel[j*blockSize+ 2] * block[2]
                              + panel[j*blockSize+ 3] * block[3];
      res.elt(iRow+1,jCol+j) += panel[j*blockSize]    * block[4]
                              + panel[j*blockSize+ 1] * block[5]
                              + panel[j*blockSize+ 2] * block[6]
                              + panel[j*blockSize+ 3] * block[7];
      res.elt(iRow+2,jCol+j) += panel[j*blockSize]    * block[8]
                              + panel[j*blockSize+ 1] * block[9]
                              + panel[j*blockSize+ 2] * block[10]
                              + panel[j*blockSize+ 3] * block[11];
      res.elt(iRow+3,jCol+j) += panel[j*blockSize]    * block[12]
                              + panel[j*blockSize+ 1] * block[13]
                              + panel[j*blockSize+ 2] * block[14]
                              + panel[j*blockSize+ 3] * block[15];
    }
  }
  /** with panel size given */
  static void blockByPanel( Block<Type> const& block, Panel<Type> const& panel
                          , Result& res, int iRow, int jCol, int pSize, int bSize)
  {
    for (int j=0; j<pSize; ++j)
      for (int i=0; i<bSize; ++i)
      { res.elt(iRow+i,jCol+j) += panel[j*blockSize]   * block[i*blockSize]
                                + panel[j*blockSize+1] * block[i*blockSize+1]
                                + panel[j*blockSize+2] * block[i*blockSize+2]
                                + panel[j*blockSize+3] * block[i*blockSize+3];}
  }
}; // struct bp

/** @ingroup hidden
 *  Methods to use for C=AB with A divided in panels and B divided in blocks.
 * The structure pb contains only static method and typedef and should normally
 * not be used directly.
 **/
template<typename Lhs, typename Rhs, typename Result>
struct bp<Lhs, Rhs, Result, (bool)Arrays::by_col_>
{
  typedef typename Result::Type Type;
  typedef hidden::MultImpl<Type> Cmult;
  typedef hidden::MultCoefImpl<Lhs, Rhs, Result> MultCoeff;
  /** Main method for Matrices multiplication implementation.
   *  @note res have been resized and initialized to zero outside this method.
   **/
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res)
  {
    // check small sizes
    switch (lhs.sizeRows())
    {
      case 0:
        return; break;
      case 1:
        MultCoeff::mult1RowOuterRow(lhs, rhs, res, lhs.beginRows());
        return; break;
      case 2:
        MultCoeff::mult2RowOuterRow(lhs, rhs, res, lhs.beginRows());
        return; break;
      case 3:
        MultCoeff::mult3RowOuterRow(lhs, rhs, res, lhs.beginRows());
        return; break;
      default: break;
    }
    switch (rhs.sizeCols())
    {
      case 0:
        return; break;
      case 1:
        MultCoeff::mult1ColOuterRow(lhs, rhs, res, rhs.beginCols());
        return; break;
      case 2:
        MultCoeff::mult2ColOuterRow(lhs, rhs, res, rhs.beginCols());
        return; break;
      case 3:
        MultCoeff::mult3ColOuterRow(lhs, rhs, res, rhs.beginCols());
        return; break;
      default: break;
    }
    // compute dimensions
    int nbInnerLoop = lhs.sizeCols()/blockSize; // = rhs.sizeRows()/blockSize;
    int nbBlocks = rhs.sizeCols()/blockSize;
    int nbPanels = lhs.sizeRows()/panelSize;

    // remaining sizes in the matrices
    int pSize = lhs.sizeRows() - panelSize*nbPanels;
    int bSize = rhs.sizeCols() - blockSize*nbBlocks;
    int tSize = lhs.sizeCols() - blockSize*nbInnerLoop;
    //
    int lastCol = rhs.beginCols() + blockSize * nbBlocks;
    int lastRow = lhs.beginRows() + panelSize * nbPanels;
    int lastPos = rhs.beginRows() + blockSize * nbInnerLoop;

    if (nbInnerLoop)
    {
      // create panels
      Panel<Type>* tabPanel = new Panel<Type>[nbPanels+1];
      Block<Type>* tabBlock = new Block<Type>[nbBlocks+1];
      // start blocks by panel
      for (int k = 0; k<nbInnerLoop; ++k)
      {
        int kPos = rhs.beginRows() + k * blockSize;
        // get panels
        for (int i = 0, iRow= lhs.beginRows(); i<nbPanels; ++i, iRow+= panelSize)
        { arrayToPanel( lhs, tabPanel[i], iRow, kPos);}
        arrayToPanel( lhs, tabPanel[nbPanels], lastRow, kPos, pSize);
        // get blocks
        for (int j = 0, jCol = rhs.beginCols(); j<nbBlocks; ++j, jCol+=blockSize)
        { arrayToBlock( rhs, tabBlock[j], kPos, jCol);}
        arrayToBlock( rhs, tabBlock[nbBlocks], kPos, lastCol, bSize);
        // perform the products blocks * panel
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int j = 0; j<nbBlocks; ++j)
        {
          int jCol = rhs.beginCols() + j * blockSize;
          for (int i = 0; i<nbPanels; ++i)
          {
            int iRow = lhs.beginRows() + i * panelSize;
            panelByBlock( tabPanel[i], tabBlock[j], res, iRow, jCol);
          }
        }
        for (int j = 0; j<nbBlocks; ++j)
        {
          int jCol = rhs.beginCols() + j * blockSize;
          panelByBlock( tabPanel[nbPanels], tabBlock[j], res, lastRow, jCol, pSize);
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i<nbPanels; ++i)
        {
          int iRow= lhs.beginRows() + i * panelSize;
          panelByBlock( tabPanel[i],  tabBlock[nbBlocks], res, iRow, lastCol, panelSize, bSize);
        }
        panelByBlock( tabPanel[nbPanels],  tabBlock[nbBlocks], res, lastRow, lastCol, pSize, bSize);
      } // k loop
      delete[] tabPanel;
      delete[] tabBlock;
    }
    // treat the remaining rows, columns
    switch (tSize)
    {
      case 1:
        MultCoeff::mult1Col(lhs, rhs, res, lastPos, lastPos);
        break;
      case 2:
        MultCoeff::mult2Col(lhs, rhs, res, lastPos, lastPos);
        break;
      case 3:
        MultCoeff::mult3Col(lhs, rhs, res, lastPos, lastPos);
        break;
      default:
        break;
    }
  }
  /** default dimensions */
  static void arrayToBlock( Rhs const& rhs, Block<Type>& block, int iRow, int jCol)
  {
    block[0]  = rhs.elt(iRow  , jCol);
    block[1]  = rhs.elt(iRow+1, jCol);
    block[2]  = rhs.elt(iRow+2, jCol);
    block[3]  = rhs.elt(iRow+3, jCol);
    block[4]  = rhs.elt(iRow  , jCol+1);
    block[5]  = rhs.elt(iRow+1, jCol+1);
    block[6]  = rhs.elt(iRow+2, jCol+1);
    block[7]  = rhs.elt(iRow+3, jCol+1);
    block[8]  = rhs.elt(iRow  , jCol+2);
    block[9]  = rhs.elt(iRow+1, jCol+2);
    block[10] = rhs.elt(iRow+2, jCol+2);
    block[11] = rhs.elt(iRow+3, jCol+2);
    block[12] = rhs.elt(iRow  , jCol+3);
    block[13] = rhs.elt(iRow+1, jCol+3);
    block[14] = rhs.elt(iRow+2, jCol+3);
    block[15] = rhs.elt(iRow+3, jCol+3);
  }
  /** with block size given */
  static void arrayToBlock( Rhs const& rhs, Block<Type>& block, int iRow, int jCol, int bSize)
  {
    for (int j=0; j<bSize; ++j)
    {
      block[j*blockSize]    = rhs.elt(iRow,   jCol+j);
      block[j*blockSize+1]  = rhs.elt(iRow+1, jCol+j);
      block[j*blockSize+2]  = rhs.elt(iRow+2, jCol+j);
      block[j*blockSize+3]  = rhs.elt(iRow+3, jCol+j);
    }
  }
  /** default dimensions */
  static void arrayToPanel( Lhs const& lhs, Panel<Type>& panel, int iRow, int kPos)
  {
    for (int i=0; i<panelSize; ++i)
    {
      panel[i*blockSize]   = lhs.elt(iRow+i,kPos);
      panel[i*blockSize+1] = lhs.elt(iRow+i,kPos+1);
      panel[i*blockSize+2] = lhs.elt(iRow+i,kPos+2);
      panel[i*blockSize+3] = lhs.elt(iRow+i,kPos+3);
    }
  }
  /** with panel size dimension given */
  static void arrayToPanel( Lhs const& lhs, Panel<Type>& panel, int iRow, int kPos, int pSize)
  {
    for (int i=0; i<pSize; ++i)
    {
      panel[i*blockSize]   = lhs.elt(iRow+i,kPos);
      panel[i*blockSize+1] = lhs.elt(iRow+i,kPos+1);
      panel[i*blockSize+2] = lhs.elt(iRow+i,kPos+2);
      panel[i*blockSize+3] = lhs.elt(iRow+i,kPos+3);
    }
  }
  /** Default dimension */
  static void panelByBlock( Panel<Type> const& panel, Block<Type> const& block
                          , Result& res, int iRow, int jCol)
  {
    for (int i=0; i<panelSize; ++i)
    {
      res.elt(iRow+i,jCol)   += panel[i*blockSize]    * block[0]
                              + panel[i*blockSize+ 1] * block[1]
                              + panel[i*blockSize+ 2] * block[2]
                              + panel[i*blockSize+ 3] * block[3];
      res.elt(iRow+i,jCol+1) += panel[i*blockSize]    * block[4]
                              + panel[i*blockSize+ 1] * block[5]
                              + panel[i*blockSize+ 2] * block[6]
                              + panel[i*blockSize+ 3] * block[7];
      res.elt(iRow+i,jCol+2) += panel[i*blockSize]    * block[8]
                              + panel[i*blockSize+ 1] * block[9]
                              + panel[i*blockSize+ 2] * block[10]
                              + panel[i*blockSize+ 3] * block[11];
      res.elt(iRow+i,jCol+3) += panel[i*blockSize]    * block[12]
                              + panel[i*blockSize+ 1] * block[13]
                              + panel[i*blockSize+ 2] * block[14]
                              + panel[i*blockSize+ 3] * block[15];
    }
  }
  static void panelByBlock( Panel<Type> const& panel, Block<Type> const& block
                          , Result& res, int iRow, int jCol, int pSize)
  {
    for (int i=0; i<pSize; ++i)
    {
      res.elt(iRow+i,jCol)   += panel[i*blockSize]    * block[0]
                              + panel[i*blockSize+ 1] * block[1]
                              + panel[i*blockSize+ 2] * block[2]
                              + panel[i*blockSize+ 3] * block[3];
      res.elt(iRow+i,jCol+1) += panel[i*blockSize]    * block[4]
                              + panel[i*blockSize+ 1] * block[5]
                              + panel[i*blockSize+ 2] * block[6]
                              + panel[i*blockSize+ 3] * block[7];
      res.elt(iRow+i,jCol+2) += panel[i*blockSize]    * block[8]
                              + panel[i*blockSize+ 1] * block[9]
                              + panel[i*blockSize+ 2] * block[10]
                              + panel[i*blockSize+ 3] * block[11];
      res.elt(iRow+i,jCol+3) += panel[i*blockSize]    * block[12]
                              + panel[i*blockSize+ 1] * block[13]
                              + panel[i*blockSize+ 2] * block[14]
                              + panel[i*blockSize+ 3] * block[15];
    }
  }
  /** with panel size dimension given */
  static void panelByBlock( Panel<Type> const& panel, Block<Type> const&  block
                          , Result& res, int iRow, int jCol, int pSize, int bSize)
  {
    for (int i=0; i<pSize; ++i)
      for (int j=0; j<bSize; ++j)
        res.elt(iRow+i,jCol+j) += panel[i*blockSize]   * block[j*blockSize]
                                + panel[i*blockSize+1] * block[j*blockSize+1]
                                + panel[i*blockSize+2] * block[j*blockSize+2]
                                + panel[i*blockSize+3] * block[j*blockSize+3];
  }
}; // struct pb

} // namespace hidden

} // namespace STK

#endif /* STK_ARRAYBYARRAYPRODUCT_H */
