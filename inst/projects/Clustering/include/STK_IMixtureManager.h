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
 * Project:  stkpp::Clustering
 * created on: 15 mars 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IMixtureManager.h
 *  @brief In this file we define the Interface IMixtureManager class.
 **/


#ifndef STK_IMIXTUREMANAGER_H
#define STK_IMIXTUREMANAGER_H

#include <DManager/include/STK_DataHandlerBase.h>
#include <Arrays/include/STK_Array2D.h>

#include "STK_IMixture.h"
#include "STK_MixtureData.h"

namespace STK
{
/** @ingroup Clustering
 *  @brief Interface base class for mixture managers.
 *
 *  A mixture manager is a factory class for injection dependency in the
 *  STK++ derived class of the IMixtureComposer. It handles all the creation
 *  and initialization stuff needed by the mixture models. It allows to get the
 *  parameters and imputed missing values from specific mixtures. It allows also
 *  to set the parameters to a specific mixture.
 *
 *  @tparam DataHandler any concrete class from the interface IDataHandler
 */
template<class DataHandler>
class IMixtureManager
{
  public:
    typedef std::vector<IMixtureData*>::const_iterator ConstDataIterator;
    typedef std::vector<IMixtureData*>::iterator DataIterator;

    /** Default constructor, need an instance of a DataHandler.  */
    IMixtureManager(DataHandler const* const p_handler);
    /** destructor */
    virtual ~IMixtureManager();
    /** @return constant pointer on the data handler */
    DataHandler const* const p_handler() const { return p_handler_;}

    /** Utility function allowing to find the idModel from the idData
     *  @param idData the id name of the data we want the idModel
     *  @return the idModel
     **/
    Clust::Mixture getIdModel( String const& idData) const;
    /** Utility function allowing to find the idModel name from the idData
     *  @param idData the id name of the data we want the idModel
     *  @return the idModel name
     **/
    String getIdModelName( String const& idData) const;
    /** @brief create a mixture and initialize it.
     *  This method get the idModelName from the DataHandler and then delegate
     *  the concrete creation to derived class using the pure virtual method
     *   @c createMixtureImpl.
     *  @param idData name of the model
     *  @param nbCluster number of cluster of the model
     *  @return 0 if the idData is not find, the result of
     *  @c createMixture( idModelName, idData, nbCluster) otherwise.
     **/
    IMixture* createMixture(String const& idData, int nbCluster);
    /** @brief register a data manager to the IMixtureManager.
     *  For each mixture created and registered, a data manager is created
     *  and registered so that it will be deleted when the mixture itself is
     *  deleted.
     *  @param p_data a pointer on the data manager
     **/
    void registerMixtureData(IMixtureData* p_data);
    /** release a data manager from v_data_.
     *  @param idData name of the data set to release
     **/
    void releaseMixtureData(String const& idData);
    // templated methods
    /** get the missing values of a data set.
     *  @param idData Id name of the data set attached to the mixture
     *  @param data the array to return with the missing values
     **/
    template<typename Type>
    void getMissingValues( String const& idData, std::vector< std::pair< std::pair<int,int>, Type > >& data) const;
    /** get the wrapper for any kind of data set using its Id
     *  @param idData Id name of the data set attached to the mixture
     *  @return a constant reference on the array with the data set
     **/
    template<typename Type>
    typename hidden::DataHandlerTraits<DataHandler, Type>::Data const& getData( String const& idData) const;
    // pure virtual methods
    /** get the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param data the array to return with the parameters
     **/
    virtual void getParameters(IMixture* p_mixture, ArrayXX& data) const =0;
    /** set the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param data the array with the parameters to set
     **/
    virtual void setParameters(IMixture* p_mixture, ArrayXX const& data) const =0;

  protected:
    /** Utility lookup function allowing to find a MixtureData from its idData
     *  @param idData the id name of the mixture we want to get
     *  @return a pointer on the MixtureData
     **/
    IMixtureData* getMixtureData( String const& idData) const;

  private:
    /** create a concrete mixture and initialize it.
     *  @param idModelName, idData strings with the Id name of the model and of the data
     *  @param nbCluster number of cluster of the model
     **/
    virtual IMixture* createMixtureImpl(String const& idModelName, String const& idData, int nbCluster) =0;
    /** A pointer on the concrete instance of the data handler */
    DataHandler const* const p_handler_;
    /** vector of pointers to the data components */
    std::vector<IMixtureData*> v_data_;
};

/* Utility function allowing to find the idModel from the idData
 *  @param idData the id name of the data we want the idModel
 *  @return the idModel
 **/
template<class DataHandler>
Clust::Mixture IMixtureManager<DataHandler>::getIdModel( String const& idData) const
{
  std::string idModelName;
  if (!p_handler()->getIdModelName( idData, idModelName))
  {
#ifdef STK_MIXTURE_VERY_VERBOSE
    stk_cout << _T("In IMixtureManager::getIdModel, fail to get idData = ") << idData << _T("\n");
#endif
    return Clust::unknown_mixture_;
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("In IMixtureManager::getIdModel, success to get idData = ") << idData << _T("\n");
  stk_cout << _T("In IMixtureManager::getIdModel, idModelName = ") << idModelName << _T("\n");
#endif
  return Clust::stringToMixture(idModelName);
}

/* Default constructor, need an instance of a DataHandler.  */
template<class DataHandler>
IMixtureManager<DataHandler>::IMixtureManager( DataHandler const* const p_handler)
                                             : p_handler_(p_handler)
{}
/* destructor */
template<class DataHandler>
IMixtureManager<DataHandler>::~IMixtureManager()
{
  for (DataIterator it = v_data_.begin() ; it != v_data_.end(); ++it)
  { delete (*it);}
}
/* Utility function allowing to find the idModel name from the idData
 *  @param idData the id name of the data we want the idModel
 *  @return the idModel name
 **/
template<class DataHandler>
String IMixtureManager<DataHandler>::getIdModelName( String const& idData) const
{
  std::string idModelName;
  if (!p_handler_->getIdModelName( idData, idModelName))
  {
#ifdef STK_MIXTURE_VERY_VERBOSE
    stk_cout << _T("In IMixtureManager::getIdModelName, fail to get idData = ") << idData << _T("\n");
#endif
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("In IMixtureManager::getIdModeName, success to get idData = ") << idData << _T("\n");
  stk_cout << _T("In IMixtureManager::getIdModel, idModelName = ") << idModelName << _T("\n");
#endif
  return idModelName;
}
/* create a mixture and initialize it.
 *  @param idData name of the model
 *  @param nbCluster number of cluster of the model
 *  @return 0 if the idData is not find, the result of @c createMixture( idModel, idData, nbCluster)
 *  otherwise.
 **/
template<class DataHandler>
IMixture* IMixtureManager<DataHandler>::createMixture(String const& idData, int nbCluster)
{
  std::string idModelName;
  if (!p_handler_->getIdModelName( idData, idModelName)) { return 0;};
  return createMixtureImpl( idModelName, idData, nbCluster);
}
/* @brief register a data manager to the IMixtureManager.
 *  For each mixture created and registered, a data manager is created
 *  and registered so that it will be deleted when the mixture itself is
 *  deleted.
 *  @param p_data a pointer on the data manager
 **/
template<class DataHandler>
void IMixtureManager<DataHandler>::registerMixtureData(IMixtureData* p_data)
{ v_data_.push_back(p_data);}
/* release a data set from v_data_.
 *  @param idData name of the data set to release
 **/
template<class DataHandler>
void IMixtureManager<DataHandler>::releaseMixtureData(String const& idData)
{
  for (DataIterator it = v_data_.begin(); it != v_data_.end(); ++it)
  { if ((*it)->idData() == idData) {delete (*it); v_data_.erase(it); break;}}
}

/* Utility lookup function allowing to find a MixtureData from its idData
 *  @param idData the id name of the mixture we want to get
 *  @return a pointer on the MixtureData
 **/
template<class DataHandler>
IMixtureData* IMixtureManager<DataHandler>::getMixtureData( String const& idData) const
{
  for (ConstDataIterator it = v_data_.begin(); it != v_data_.end(); ++it)
  {  if ((*it)->idData() == idData) return (*it);}
  return 0;
}

/** get the missing values of a data set.
 *  @param idData Id name of the data set attached to the mixture
 *  @param data the array to return with the missing values
 **/
template<class DataHandler>
template<typename Type>
void IMixtureManager<DataHandler>::getMissingValues( String const& idData, std::vector< std::pair< std::pair<int,int>, Type > >& data) const
{
  typedef typename hidden::DataHandlerTraits<DataHandler, Type>::Data DataType;
  typedef MixtureData<DataType> MixtureDataType;

  IMixtureData* p_data = getMixtureData(idData);
  // up-cast... (Yes it's bad....;)...)
  if (p_data)
  { static_cast<MixtureDataType const*>(p_data)->getMissingValues(data);}
}
/** get the wrapper for any kind of data set using its Id
 *  @param idData Id name of the data set attached to the mixture
 *  @return a constant reference on the array with the data set
 **/
template<class DataHandler>
template<typename Type>
typename hidden::DataHandlerTraits<DataHandler, Type>::Data const& IMixtureManager<DataHandler>::getData( String const& idData) const
{
  typedef typename hidden::DataHandlerTraits<DataHandler, Type>::Data DataType;
  typedef MixtureData<DataType> MixtureDataType;
  return static_cast<MixtureDataType const*>(getMixtureData(idData))->dataij_;
}
} // namespace STK


#endif /* STK_IMIXTUREMANAGER_H */
