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
 * Purpose:  Define the Interface 1D templated Container class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_ITContainer1D.h
 *  @brief in this file we define the interface classes IContainer1D and ITContainer1D.
 **/

#ifndef STK_ITCONTAINER1D_H
#define STK_ITCONTAINER1D_H


#include <Sdk/include/STK_IRecursiveTemplate.h>
#include <Sdk/include/STK_Macros.h>

#include "STK_Traits.h"
#include "STK_Arrays_Util.h"

namespace STK
{
/** @ingroup Arrays
 *  @brief Interface base class for 1D containers.
 *
 * The IContainer1D class is an interface base class for all
 * one dimensional containers.
 **/
template<int Size_>
class IContainer1D
{
  public:
    typedef TRange<Size_> Range1D;

    /** Default constructor */
    IContainer1D() : range_() {}
    /** constructor with specified range
     *  @param I the Range of the container
     **/
    IContainer1D(Range const& I) : range_(I) {}
    /** copy constructor
     *  @param V the container to copy
     **/
    IContainer1D(IContainer1D const& V) : range_(V.range_) {}
    /** destructor. */
    ~IContainer1D() {}

    /**  @return the range of the container */
    Range1D const& range() const  { return range_;}
    /** @return the index of the first element */
    inline int begin() const { return range_.begin();}
    /**  @return the ending index of the elements */
    inline int end() const { return range_.end();}
    /**  @return the size of the container */
    inline int size() const  { return range_.size();}
    /**  @return the index of the last element */
    inline int lastIdx() const  { return range_.lastIdx();}

    /** Is there some data ?
     *  @return @c true if the container is empty, @c false otherwise
     **/
    bool empty() const { return range_.empty();}

  protected:
    /** exchange this container with T
     * @param T the container to exchange with T
     **/
     void exchange(IContainer1D& T) { std::swap(T.range_, range_ );}
    /** Set the beginning of the range
     *  @param beg the first index of the range
     **/
    void shift( int const& beg) { range_.shift(beg);}
    /** Set range of the rows of the container.
     *  @param I the range to set (default empty)
     **/
    void setRange(Range const& I = Range1D()) { range_ = I;}
    /** increment the range of the container (can be negative).
     *  @param inc increment to apply to the range
     **/
    void incRange(int const& inc) { range_.inc(inc);}
    /** increment the beginning of the container (can be negative).
     *  @param inc the increment to apply to the beginning of the range
     **/
    void incFirst(int const& inc) { range_.incFirst(inc);}
    /** decrement the beginning of the container.
     *  @param inc the decrement to apply to the beginning of the range
     **/
    void decFirst(int const& inc) { range_.decFirst(inc);}
    /** increment the end of the container (can be negative).
     *  @param inc the increment to apply to the end of the range
     **/
    void incLast(int const& inc =1) { range_.incLast(inc);}
    /** decrement the end of the container.
     *  @param inc the decrement to apply to the end of the range
     **/
    void decLast(int const& inc =1) { range_.decLast(inc);}

  private:
    /** range of the array. */
    Range1D range_;
};

/** @ingroup Arrays
 *  @brief Interface base class for homogeneous 1D containers.
 *
 * The ITContainer1D class is the templated base class for all
 * homogeneous one-dimensional containers containing element of type @c Type
 * where Type is note necessarily a scalar. It assumes that the derived class
 * cannot be part of an expression and is not constant, so that it can be
 * #- shifted,
 * #- resized
 * #- accessed in modification.
 *
 * Implement the curious recursive template paradigm : the template
 * parameter @c Derived is the name of the class that
 * implements @c ITContainer1D. For example
 * <code>
 * template<class Type>
 * class Derived : public ITContainer1D< Derived<Type> >
 * {...}
 * </code>
 *
 * The pseudo virtual function defined in this interface have the
 * following definition:
 * @code
 *   Type& elt1Impl(int const& pos);
 *   Type const& elt1Impl(int const& pos) const;
 *   void shiftImpl(int const& beg);
 *   Derived& resizeImpl(Range const& I);
 * @endcode
 * All these methods have to be implemented in the Derived class.
 *
 * @sa List1D, Array1D
 *
 * @note The constant getter @c elt1Impl(pos) have to return a reference as we
 * are using derived classes for storing any kind of data.
 **/
template <class Derived>
class ITContainer1D : public IContainer1D< hidden::Traits<Derived>::size_ >
                    , public IRecursiveTemplate<Derived>
{
  public:
    typedef IContainer1D< hidden::Traits<Derived>::size_ > Base;
    using Base::begin;
    using Base::end;
    using Base::lastIdx;

  protected:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;

    /** Default constructor */
    ITContainer1D(): Base() {}
    /** constructor with a specified range.
     *  @param I : the range of the container
     **/
    ITContainer1D( Range const& I): Base(I) {}
    /** Copy constructor
     *  @param T the container to copy
     **/
    ITContainer1D( ITContainer1D const& T): Base(T) {}
    /** destructor. */
    ~ITContainer1D() {}

  public:
    /** @return the ith element for vector_, point_ and diagonal_ containers
     *  @param i index of the ith element
     **/
    inline Type& elt(int i)
    {
#ifdef STK_BOUNDS_CHECK
      if (this->asDerived().begin() > i)
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::elt, i, begin() > i);}
      if (this->asDerived().end() <= i)
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::elt, i, end() <= i);}
#endif
      return this->asDerived().elt1Impl(i);
    }
    /** @return a constant reference on the ith element for vector_, point_ and diagonal_ containers
     *  @param i index of the ith element
     **/
    inline Type const& elt(int i) const
    {
#ifdef STK_BOUNDS_CHECK
      if (this->asDerived().begin() > i)
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::elt, i, begin() > i);}
      if (this->asDerived().end() <= i)
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::elt, i, end() <= i);}
#endif
      return this->asDerived().elt1Impl(i);
    }
    /** @return the element ith element
     *  @param i index of the ith element
     **/
    inline Type& operator[](int i) { return elt(i);}
    /** @return a constant reference on the ith  element
     *  @param i index of the ith element
     **/
    inline Type const& operator[](int i) const { return elt(i);}
    /** @return safely the jth element
     *  @param i index of the element
     **/
    inline Type& at(int i)
    {
      if (this->asDerived().begin() > i)
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::at, i, begin() > i);}
      if (this->asDerived().end() <= i)
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::at, i, lastIdx() <= i);}
      return elt(i);
    }
    /** @return safely the constant jth element
     *  @param i index of the element
     **/
    Type const& at(int i) const
    {
      if (this->asDerived().begin() > i)
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::at, i, begin() > i);}
      if (this->asDerived().lastIdx() < i)
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::at, i, lastIdx() < i);}
      return elt(i);
    }
    /** Access to many elements.
     *  @param I the range of the elements
     *  @return a reference container with the elements of this in the range I
     **/
    inline SubVector sub(Range const& I) const
    {
#ifdef STK_BOUNDS_CHECK
      if ((I.begin()<begin()))
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::sub,I,I.begin()<begin());}
      if ((I.end()>end()))
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::sub,I,I.end()>end());}
#endif
      return SubVector(this->asDerived(),I);
    }

    /** @return a reference on the first element. */
    inline Type& front() { return elt(begin());}
    /** @return a constant reference on the first element */
    inline Type const& front() const { return elt(begin());}
    /** @return a reference on the last element */
    inline Type& back() { return elt(lastIdx());}
    /** @return a constant reference on the last element */
    inline Type const& back() const { return elt(lastIdx());}

    /**  @param beg the index of the first column to set */
    void shift(int beg) { this->asDerived().shiftImpl(beg);}
    /** @return the resized container.
     *  @param I the range of the container
     **/
    Derived& resize(Range const& I) { return this->asDerived().resizeImpl(I);}
};

} // namespace STK

#endif // STK_ITCONTAINER1D_H
