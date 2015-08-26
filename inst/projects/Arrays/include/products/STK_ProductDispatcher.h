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
 * created on: 25 déc. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ProductDispatcher.h
 *  @brief In this file we select the product method to use.
 **/

#ifndef STK_PRODUCTDISPATCHER_H
#define STK_PRODUCTDISPATCHER_H

#include "STK_ProductRaw.h"
#include "STK_ArrayByVectorProduct.h"
#include "STK_ArrayByArrayProduct.h"

namespace STK
{

namespace hidden
{
/** In the general case, e.g. when some of the structures are triangular for
 *  examples, we use the usual matrix multiplication formula.
 **/
template < class Lhs, class Rhs, class Result
         , int lhsStructure_ = hidden::Traits<Lhs>::structure_
         , int RhsStructure_ = hidden::Traits<Rhs>::structure_ >
struct ProductDispatcher
{
  enum
  {
    // structure_ = Traits<Result>::structure_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  typedef MultCoefImpl<Lhs, Rhs, Result> MultCoeff;

  static void runbp(Lhs const& lhs, Rhs const& rhs, Result& res )
  {
    for (int j=rhs.beginCols(); j< rhs.endCols(); j++)
    {
      Integer const end = res.rangeRowsInCol(j).end();
      for (int i=res.rangeRowsInCol(j).begin(); i< end; i++)
      { MultCoeff::dot(lhs, rhs, res, i, j);}
    }
  }
  static void runpb(Lhs const& lhs, Rhs const& rhs, Result& res )
  {
    for (int i=lhs.beginRows(); i< lhs.endRows(); i++)
    {
      Integer const end = res.rangeColsInRow(i).end();
      for (int j=res.rangeColsInRow(i).begin(); j< end; j++)
      { MultCoeff::dot(lhs, rhs, res, i, j);}
    }
  }
};

/** Specialization for the array2d by array2D case. */
template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::array2D_, Arrays::array2D_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void runbp(Lhs const& lhs, Rhs const& rhs, Result& res )
  { bp<Lhs,Rhs,Result,orient_>::run(lhs, rhs, res);}
  static void runpb(Lhs const& lhs, Rhs const& rhs, Result& res )
  { bp<Lhs,Rhs,Result,orient_>::run(lhs, rhs, res);}
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::array2D_, Arrays::square_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void runbp(Lhs const& lhs, Rhs const& rhs, Result& res )
  { bp<Lhs,Rhs,Result,orient_>::run(lhs, rhs, res);}
  static void runpb(Lhs const& lhs, Rhs const& rhs, Result& res )
  { bp<Lhs,Rhs,Result,orient_>::run(lhs, rhs, res);}
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::square_, Arrays::square_>
{
  enum
  {
    structure_ = Arrays::square_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void runbp(Lhs const& lhs, Rhs const& rhs, Result& res )
  {
    bp<Lhs,Rhs,Result,orient_>::run(lhs, rhs, res);
  }
  static void runpb(Lhs const& lhs, Rhs const& rhs, Result& res )
  { bp<Lhs,Rhs,Result,orient_>::run(lhs, rhs, res);}
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::square_, Arrays::array2D_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void runbp(Lhs const& lhs, Rhs const& rhs, Result& res )
  { bp<Lhs,Rhs,Result, orient_>::run(lhs, rhs, res);}
  static void runpb(Lhs const& lhs, Rhs const& rhs, Result& res )
  { bp<Lhs,Rhs,Result,orient_>::run(lhs, rhs, res);}
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::array2D_, Arrays::vector_>
{
  enum
  {
    structure_ = Arrays::vector_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  { bv<Lhs,Rhs,Result>::run(lhs, rhs, res);}
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::square_, Arrays::vector_>
{
  enum
  {
    structure_ = Arrays::vector_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  { bv<Lhs,Rhs,Result>::run(lhs, rhs, res);}
};

template <class Lhs, class Rhs, class Result, int lhsStructure_>
struct ProductDispatcher<Lhs, Rhs, Result, lhsStructure_, Arrays::vector_>
{
  enum
  {
    structure_ = Arrays::vector_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  typedef MultCoefImpl<Lhs, Rhs, Result> MultCoeff;
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  {
    for (int i=lhs.beginRows(); i< lhs.endRows(); i++)
    { MultCoeff::dot(lhs, rhs, res, i);}
  }
};


template <class Lhs, class Rhs, class Result, int RhsStructure_>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::point_, RhsStructure_>
{
  enum
  {
    structure_ = Arrays::point_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  typedef MultCoefImpl<Lhs, Rhs, Result> MultCoeff;
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  {
    for (int j=rhs.beginCols(); j< rhs.endCols(); j++)
    { MultCoeff::dot(lhs, rhs, res, j);}
  }
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::point_, Arrays::array2D_>
{
  enum
  {
    structure_ = Arrays::point_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  { vb<Lhs,Rhs,Result>::run(lhs, rhs, res);}
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::point_, Arrays::square_>
{
  enum
  {
    structure_ = Arrays::point_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  { vb<Lhs,Rhs,Result>::run(lhs, rhs, res);}
};

} // namespace hidden

} // namespace STK

#endif /* STK_PRODUCTDISPATCHER_H */
