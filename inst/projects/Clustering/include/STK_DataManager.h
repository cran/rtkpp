/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013 Serge Iovleff

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 15 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_DataManager.h
 *  @brief In this file we define the data manager class associated to the MixtureBridge class.
 **/

#ifndef STK_DATAMANAGER_H
#define STK_DATAMANAGER_H

#include "STK_IDataManager.h"

namespace STK
{

namespace hidden
{
/** @ingroup hidden
 *  Implementation of the safeValue method. Default implementation.
 */
template<int Id, class Data>
struct SafeValueImpl
{
  // type of the data
  typedef typename Data::Type Type;

  /** @return a safe value for the jth variable
   *  @param  m_dataij the matrix of the data set
   *  @param j index of the column with the safe value needed */
  static Type run(Data const& m_dataij, int j)
  { return m_dataij.col(j).safe().mean();}
};
/** @ingroup hidden
 *  Implementation of the safeValue method. Specialization for Gamma_ajk_bjk_ models
 */
template<class Data>
struct SafeValueImpl<Clust::Gamma_, Data>
{
  // type of the data
  typedef typename Data::Type Type;

  /** @return a safe value for the jth variable
   *  @param  m_dataij the matrix of the data set
   *  @param j index of the column with the safe value needed */
  static Type run(Data const& m_dataij, int j)
  { return m_dataij.col(j).safe(1).mean();}
};
/** @ingroup hidden
 *  Implementation of the safeValue method. Speciualization for Categorical_pjk_ models
 */
template<class Data>
struct SafeValueImpl<Clust::Categorical_, Data>
{
  // type of the data
  typedef typename Data::Type Type;

  /** @return a safe value for the jth variable
   *  @param  m_dataij the matrix of the data set
   *  @param j index of the column with the safe value needed */
  static Type run(Data const& m_dataij, int j)
  {
     int lmin = m_dataij.col(j).safe().minElt(), lmax = m_dataij.col(j).safe().maxElt();
     Array2DVector<int> count(Range(lmin, lmax, 0), 0);
     for (int i= m_dataij.beginRows(); i < m_dataij.endRows(); ++i)
     {
       if (Arithmetic<int>::isFinite(m_dataij(i,j)))
         count[m_dataij(i,j)]++;
     }
     int l; count.maxElt(l);
     return l;
  }
};

} // namespace hidden

/** @ingroup Clustering
 *  @brief bridge a data set with a MixtureBridge.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface class STK::IMixtureModel.
 */
template<Clust::MixtureClass Id, class Data>
class DataManager : public IDataManager
{
  public:
    typedef std::vector<std::pair<int,int> >::const_iterator ConstIterator;
    // type of the data
    typedef typename Data::Type Type;

    /** default constructor. */
    inline DataManager(std::string const& idData) : IDataManager(idData), m_dataij_() {}
    /** copy constructor (Warning: will copy the data set)
     *  @param manager the DataManager to copy
     **/
    DataManager( DataManager const& manager)
               : IDataManager(manager), m_dataij_(manager.m_dataij_) {}
    /** getter. @return a constant reference on the data set */
    Data const& m_dataij() const { return m_dataij_;}
    /** data set (public) */
    Data m_dataij_;
    /** utility function for lookup the data set and find missing values
     *  coordinates. */
   virtual void findMissing()
   {
#ifdef STK_MIXTURE_VERBOSE
     stk_cout << _T("Entering findMissing()\n");
#endif
     for (int j=m_dataij_.beginCols(); j< m_dataij_.endCols(); ++j)
     {
       for (int i=m_dataij_.beginRows(); i< m_dataij_.endRows(); ++i)
       {
         if (Arithmetic<Type>::isNA(m_dataij_(i,j)))
         { v_missing_.push_back(std::pair<int,int>(i,j));}
       }
     }
#ifdef STK_MIXTURE_VERBOSE
     stk_cout << _T("findMissing() terminated, nbMiss= ") << v_missing_.size() << _T("\n");
#endif
   }
   /** utility function for lookup the data set and remove the missing values.*/
   virtual void removeMissing()
   {
     Type value = Type();
     int j, old_j = Arithmetic<int>::NA();
     for(ConstIterator it = v_missing_.begin(); it!= v_missing_.end(); ++it)
     {
        j = it->second; // get column
        if (j != old_j)
        {
          old_j = j;
          value = safeValue(j);
        }
        m_dataij_(it->first, j) = value;
      }
   }
   /** @return a safe value for the jth variable
    *  @param j index of the column with the safe value needed
    **/
   Type safeValue( int j) const
   { return hidden::SafeValueImpl<Id, Data>::run(m_dataij_, j);}
   /** get the missing values of a data set.
    *  @param data the array to return with the missing values
    **/
   template<typename Type_>
   void getMissingValues( std::vector< std::pair<std::pair<int,int>, Type_ > >& data) const
   {
     data.resize(v_missing_.size());
     for(size_t i = 0; i< v_missing_.size(); ++i)
     {
       data[i].first  = v_missing_[i];
       data[i].second = m_dataij_(v_missing_[i].first, v_missing_[i].second);
     }
   }
};

} // namespace STK

#endif /* STK_DATAMANAGER_H */
