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
 * Purpose:  Define the Base Interface for the Array classes.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_AllocatorBase.h
  * @brief In this file we define the class AllocatorBase
 **/

#ifndef STK_ALLOCATORBASE_H
#define STK_ALLOCATORBASE_H

#include <STKernel/include/STK_Range.h>
#include <Sdk/include/STK_Macros.h>
#include "STK_IContainerRef.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Templated base class for all Allocator classes.
 *
 *  The AllocatorBase class is the base class of all Allocators string data
 *  in memory. This class manages the main pointer on the data.
 *  It derives from the IContainerRef class as an array stored in
 *  memory can always be wrapped in some way or be a wrapper of
 *  data stored in memory.
 * 
 *  This class can also be used as a concrete class.
 *  @tparam Type can be any type of data that can be stored in memory.
 **/
template<typename Type, int Size>
class AllocatorBase: public IContainerRef
{
  public:
    typedef TRange<Size> AllocRange;

    /** Default constructor. */
    AllocatorBase() : IContainerRef(false)
                           , p_data_(0)
                           , rangeData_()
    {}
    /** constructor with specified Range
     *  @param I range of the data
     **/
    AllocatorBase( Range const& I) : IContainerRef(false)
                                   , p_data_(0)
                                   , rangeData_()
    { malloc(I);}
    /** Copy constructor.We don't know, how the user classes want to copy the
     *  data if this is not a reference.
     *  @param T : the array to copy
     *  @param ref : is this a wrapper of T ?
     **/
    AllocatorBase( AllocatorBase const& T, bool ref = false)
                        : IContainerRef(ref)
                        , p_data_(ref ? T.p_data_: 0)
                        , rangeData_(T.rangeData_)
    {/* derived class have to copy the data if ref==false */}
    template<int OtherSize>
    inline AllocatorBase( AllocatorBase<Type, OtherSize> const& T, bool ref = false)
                        : IContainerRef(ref)
                        , p_data_(ref ? T.p_data(): 0)
                        , rangeData_(T.rangeData())
    {/* derived class have to copy the data if ref==false */}
    /** @brief Wrapper or copy constructor. It is set as protected as we don't
     *  know, how the end-user want to copy the data.
     *  @param q ptr to the data to wrap
     *  @param I range of the data wrapped
     *  @param ref is this a wrapper ?
     **/
    AllocatorBase( Type* const& q, Range const& I, bool ref)
                        : IContainerRef(ref)
                        , p_data_(ref ? q : 0)
                        , rangeData_(I)
    { /* derived class have to copy the data if ref==false */}
    /** Wrapper or copy constructor : second form. This constructor assume the
     *  data as a C-like array. Thus the first index is 0 and the last
     *  index is size-1.
     *  @param q ptr to the data to wrap
     *  @param size of the data to wrap
     *  @param ref is this a wrapper ?
     **/
    AllocatorBase( Type* const& q, int size, bool ref)
                        : IContainerRef(ref)
                        , p_data_(ref ? q : 0)
                        , rangeData_(Range(0,size-1))
    { /* derived class have to copy the data if ref==false */}

    /** destructor. */
    ~AllocatorBase() { free(); }

    /** @return the range of the data*/
    AllocRange const& rangeData() const { return rangeData_;}
    /** @return the first index of the data. */
    inline int firstData() const { return rangeData_.begin();}
    /**@return the ending index of the data */
    inline int endData() const { return rangeData_.end();}
   /**@return the last index of the data */
    inline int lastData() const { return rangeData_.lastIdx();}
    /** @return the size of the data */
    inline int sizeData() const { return rangeData_.size();}
    /** @return a pointer on the constant data set*/
    inline Type* const& p_data() const { return p_data_;}
    /** @return a pointer on the data set */
    inline Type* p_data() { return p_data_;}
    /** Get the const element number pos.
     *  @param pos the position of the element we get 
     **/
    inline Type const& data( int pos) const
    {
#ifdef STK_BOUNDS_CHECK
      if (pos < firstData())
      { STKOUT_OF_RANGE_1ARG(AllocatorBase::data,pos,AllocatorBase::firstData() > pos);}
      if (pos > lastData())
      { STKOUT_OF_RANGE_1ARG(AllocatorBase::data,pos,AllocatorBase::lastData() < pos);}
#endif
      return p_data_[pos];
    }
    /** Get the element number pos.
     *  @param pos the position of the element we get 
     **/
    inline Type& data(int pos)
    {
#ifdef STK_BOUNDS_CHECK
      if (pos < firstData())
      { STKOUT_OF_RANGE_1ARG(AllocatorBase::data,pos,AllocatorBase::firstData() > pos);}
      if (pos > lastData())
      { STKOUT_OF_RANGE_1ARG(AllocatorBase::data,pos,AllocatorBase::lastData() < pos);}
#endif
      return p_data_[pos];
    }
    /** @brief main ptr memory allocation.
     *  @param I range of the data allocated
     **/
    void malloc( Range const& I);
    /** @brief function for main ptr memory reallocation.
     *
     *  If the size requested is greater than the allocated size,
     *  the Type stored are saved and copied using the operator=. the Type
     *  class have to provide this operator.
     *
     *  If the size requested is lesser than the allocated size, only
     *  the first elements fitting in the container are copied.
     *  @param I range of the data to reserve
     **/
    void realloc( Range const& I);
    /** function for main ptr memory deallocation. */
    void free();
    /** exchange this with T.
     *  @param T the container to exchange with T
     **/
    AllocatorBase& exchange(AllocatorBase &T)
    {
      std::swap(p_data_, T.p_data_);
      std::swap(rangeData_, T.rangeData_);
      IContainerRef::exchange(T);
      return *this;
    }
    /** swap two elements of the Allocator.
     *  @param pos1, pos2 the positions of the first and second element
     **/
    void swap(int pos1, int pos2)
    { std::swap(p_data_[pos1], p_data_[pos2]);}
    /** @brief copy the Allocator T by value.
     *  The memory is free and the Allocator T is physically copied in this.
     *  @param T the allocator to copy by value
     *  @return a copy of this
     **/
    AllocatorBase& copy( AllocatorBase const& T)
    {
      // allocate memory if necessary
      malloc(T.rangeData_);
      // copy values
      for (int pos= T.firstData(); pos < T.endData(); ++pos)
      { p_data_[pos] = T.p_data_[pos];}
      return *this;
    }
    /** @brief move the Allocator T to this.
     *  The memory is free and T become a reference of this. This method allow
     *  to steal the data of T without physical copy.
     *
     *  @return this object.
     *
     *  @note the data member ref_ is mutable so that T can be passed as a
     *  constant reference.
     *
     *  @param T the allocator to copy as reference
     **/
    AllocatorBase& move( AllocatorBase const& T)
    {
      if (this == &T) return *this;
      free();
      setPtrData(T.p_data_, T.rangeData_, T.isRef());
      T.setRef(true);
      return *this;
    }
    /** shift the first index of the data to first.
     *  @param first the index of the first data to set
     **/
    void shiftData(int const& first)
    {
      // check if there is something to do
      if (first == firstData()) return;
      // check for reference
      if (this->isRef())
        STKRUNTIME_ERROR_1ARG(AllocatorBase::shiftData, first, cannot operate on reference.);
      // compute increment
      int inc = first - firstData();
      // translate data
      decPtrData(inc);
    }
    /** @brief Set the address and the range of the data.
     * This method is to be used when the memory have been allocated outside.
     * Allocated memory is not freed.
     *
     *  @param p_data the address to set
     *  @param rangeData the range of the data
     *  @param ref is p_data_ a wrapper ?
     **/
    void setPtrData( Type* p_data,  Range const& rangeData, bool ref)
    { p_data_ = p_data; rangeData_ = rangeData; this->setRef(ref);}

  private:
    /** Set the address of the data : this method is not destined
     *  to the end-user.
     *  @param p_data the address to set
     **/
    void setPtrData( Type* p_data = 0) { p_data_ = p_data;}
    /** Set the index of the first data : this method is not destined
     *  to the end-user.
     *  @param rangeData the range of the data to set
     **/
    void setRangeData( Range const& rangeData = Range())
    { rangeData_ = rangeData;}
    /** Set array members to default values. */
    void setDefault()
    {
      setPtrData();
      setRangeData();
    }
    /** Increment the address of the data.
     *  @param inc the increment to apply
     **/
    void incPtrData( int const& inc)
    {
      if (p_data_) { p_data_ += inc;}
      rangeData_.dec(inc);
    }
    /** Decrement the address of the data.
     *  @param dec the increment to apply
     **/
    void decPtrData( int const& dec)
    {
      if (p_data_) { p_data_ -= dec;}
      rangeData_.inc(dec);
    }
  private:
    /** Main pointer on the data. */
    Type* p_data_;
    /** Range of the data */
    AllocRange rangeData_;
};

template<typename Type, int Size>
void AllocatorBase<Type,Size>::malloc( Range const& I)
{
  {
    if ((this->rangeData() == I)&&(p_data_)&&(!this->isRef())) return;
    // free any existing data
    free();
    // check size
    if (I.size() <= 0)
    {
      setPtrData(0, I, false);
      this->setRef(false);
      return;
    }
    // allocate memory
    try
    {
      setPtrData(new Type[I.size()], Range(0, I.size()), false);
      decPtrData(I.begin());
    }
    catch (std::bad_alloc const& error)
    {
      setDefault();
      STKRUNTIME_ERROR_1ARG(AllocatorBase::malloc, I, memory allocation failed);
    }
    this->setRef(false);
  }
}

template<typename Type, int Size>
void AllocatorBase<Type,Size>::realloc( Range const& I)
{
  if ((this->rangeData() == I)&&(p_data_)&&(!this->isRef())) return;
  // check size
  if (I.size() <= 0)
  {
    free();
    setPtrData(0, I, false);
    return;
  }
  try
  {
    // allocate memory and apply increment
    Type* p  = new Type[I.size()];
     p -= I.begin();
    // no error: copy data
    const int begin = std::max(rangeData_.begin(), I.begin())
            , end  = std::min(rangeData_.end(), I.end());
    for (int i = begin; i<end; ++i) { p[i] = p_data_[i];}
    // liberate old memory
    free();
    // set default values
    setPtrData(p, I, false);
  }
  catch (std::bad_alloc const& error)
  { STKRUNTIME_ERROR_1ARG(AllocatorBase::realloc, I, memory allocation failed);}
}
/** function for main ptr memory deallocation. */
template<typename Type, int Size>
void AllocatorBase<Type,Size>::free()
{
  // nothing to do for reference
  if (this->isRef()) return;
  // if there is elts
  if (p_data_)
  {
    incPtrData(firstData());  // translate
    delete [] p_data_;        // erase
    setDefault();             // set default values
  }
}

} // namespace STK

#endif /* STK_ALLOCATORBASE_H */
