/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2011  Serge Iovleff

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
 * created on: 25 nov. 2011
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_CArrayPoint.h
 *  @brief In this file we implement the final class CArrayPoint.
 **/

#ifndef STK_CARRAY2DPOINT_H
#define STK_CARRAY2DPOINT_H

#include "STK_ICArray.h"

namespace STK
{
// forward declarations
template< typename Type, int SizeCols_=UnknownSize, bool Orient_ = Arrays::by_col_>
class CArrayPoint;
template< typename Type, int SizeRows_, int SizeCols_, bool Orient_> class CArray;
template< typename Type, int Size_, bool Orient_> class CArraySquare;
template< typename Type, int SizeRows_, bool Orient_> class CArrayVector;
template< typename Type, bool Orient_> class CArrayNumber;

// useful typedef
typedef CArrayPoint<Real, UnknownSize, Arrays::by_col_>   CPointX;
typedef CArrayPoint<Real, 2, Arrays::by_col_>             CPoint2;
typedef CArrayPoint<Real, 3, Arrays::by_col_>             CPoint3;
typedef CArrayPoint<double, UnknownSize, Arrays::by_col_> CPointXd;
typedef CArrayPoint<double, 2, Arrays::by_col_>           CPoint2d;
typedef CArrayPoint<double, 3, Arrays::by_col_>           CPoint3d;
typedef CArrayPoint<int, UnknownSize, Arrays::by_col_>    CPointXi;
typedef CArrayPoint<int, 2, Arrays::by_col_>              CPoint2i;
typedef CArrayPoint<int, 3, Arrays::by_col_>              CPoint3i;

typedef CArrayPoint<Real, UnknownSize, Arrays::by_row_>   CPointByRowX;
typedef CArrayPoint<Real, 2, Arrays::by_row_>             CPointByRow2;
typedef CArrayPoint<Real, 3, Arrays::by_row_>             CPointByRow3;
typedef CArrayPoint<double, UnknownSize, Arrays::by_row_> CPointByRowXd;
typedef CArrayPoint<double, 2, Arrays::by_row_>           CPointByRow2d;
typedef CArrayPoint<double, 3, Arrays::by_row_>           CPointByRow3d;
typedef CArrayPoint<int, UnknownSize, Arrays::by_row_>    CPointByRowXi;
typedef CArrayPoint<int, 2, Arrays::by_row_>              CPointByRow2i;
typedef CArrayPoint<int, 3, Arrays::by_row_>              CPointByRow3i;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CArray class.
 */
template<typename Type_, int SizeCols_, bool Orient_>
struct Traits< CArrayPoint<Type_, SizeCols_, Orient_> >
{
    typedef CArrayNumber<Type_, Orient_> Number;
//    typedef CArrayNumber<Type_, Orient_> Col;
//    typedef CArrayNumber<Type_, Orient_> SubCol;
    typedef CArrayVector<Type_, 1, Orient_> Col;
    typedef CArrayVector<Type_, 1, Orient_> SubCol;

    typedef CArrayPoint<Type_, SizeCols_, Orient_> Row;
    typedef CArrayPoint<Type_, UnknownSize, Orient_> SubRow;

    /* Type or array (1,1) ? */
    typedef typename If<(SizeCols_ == 1), Number, SubRow>::Result SubVector;
    typedef typename If<(SizeCols_ == 1), Number, SubRow>::Result SubArray;
    // The CAllocator have to have the same structure than the CArray
    typedef CAllocator<Type_, 1, SizeCols_, Orient_> Allocator;

    typedef Type_                Type;
    typedef typename RemoveConst<Type_>::Type const& ReturnType;
    enum
    {
      structure_ = Arrays::point_,
      orient_    = Orient_,
      sizeRows_  = 1,
      sizeCols_  = SizeCols_,
      storage_   = Arrays::dense_
    };
};

} // namespace hidden


/** @ingroup Arrays
 * @brief declaration of  the point case.
 */
template <typename Type_, int SizeCols_, bool Orient_>
class CArrayPoint : public ICArray < CArrayPoint<Type_, SizeCols_, Orient_> >
{
  public:
    typedef ICArray < CArrayPoint<Type_, SizeCols_, Orient_> > Base;
    typedef ArrayBase < CArrayPoint<Type_, SizeCols_, Orient_> > LowBase;

    typedef typename hidden::Traits<CArrayPoint<Type_, SizeCols_, Orient_> >::Type Type;
    typedef typename hidden::Traits<CArrayPoint<Type_, SizeCols_, Orient_> >::ReturnType ReturnType;

    enum
    {
      structure_ = Arrays::point_,
      orient_    = Orient_,
      sizeRows_  = 1,
      sizeCols_  = SizeCols_,
      storage_   = Arrays::dense_
    };

    /** Default constructor. */
    inline CArrayPoint(): Base() {}
    /** constructor with specified dimension.
     *  @param sizeCols range of the columns
     **/
    inline CArrayPoint( int sizeCols): Base(1, sizeCols) {}
    /** constructor with rbeg, rend, cbeg and cend specified,
     *  initialization with a constant.
     *  @param sizeCols range of the columns
     *  @param v initial value of the container
     **/
    inline CArrayPoint( int sizeCols, Type const& v): Base(1, sizeCols, v) {}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    inline CArrayPoint( CArrayPoint const& T, bool ref=false): Base(T, ref) {}
    /** wrapper constructor for 0 based C-Array.
     *  @param q pointer on the array
     *  @param nbCol number of columns
     **/
    inline CArrayPoint( Type* const& q, int nbCol): Base(q, 1, nbCol) {}
    /** constructor by reference.
     *  @param allocator the allocator to wrap
     **/
    template<class OtherAllocator>
    inline CArrayPoint( CAllocatorBase<OtherAllocator> const& allocator)
                      : Base(allocator.asDerived())
    {}
    /** Copy constructor using an expression.
     *  @param T the container to wrap
     **/
    template<class OtherDerived>
    inline CArrayPoint( ExprBase<OtherDerived> const& T): Base(1,T.size()) { LowBase::operator=(T);}
    /** destructor. */
    inline ~CArrayPoint() {}
    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    inline CArrayPoint& operator=(Type const& v) { return LowBase::setValue(v);}
    /** operator = : overwrite the CArrayPoint with the Right hand side T.
     *  @param T the container to copy
     **/
    template<class Rhs>
    inline CArrayPoint& operator=( ExprBase<Rhs> const& T) { return LowBase::assign(T);}
    /** operator = : overwrite the CArray with the Right hand side rhs.
     *  @param rhs the container to copy
     **/
    inline CArrayPoint& operator=(CArrayPoint rhs) { return LowBase::assign(rhs);}
};

} // namespace STK


#endif /* STK_CARRAY2DPOINT_H */
