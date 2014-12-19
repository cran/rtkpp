/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff

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

/** @file STK_IMixtureManager.cpp
 *  @brief In this file we implement the Interface IMixtureManager class.
 **/

#include "DManager/include/STK_IDataHandler.h"

#include "../include/STK_IMixtureData.h"
#include "../include/STK_IMixtureManager.h"

namespace STK
{

/* Default constructor, need an instance of a DataHandler.  */
IMixtureManager::IMixtureManager(IDataHandler const* const p_handler) : p_handler_(p_handler)
{}
/* destructor */
IMixtureManager::~IMixtureManager()
{
  for (DataIterator it = v_data_.begin() ; it != v_data_.end(); ++it)
  { delete (*it);}
}
/* Utility function allowing to find the idModel name from the idData
 *  @param idData the id name of the data we want the idModel
 *  @return the idModel name
 **/
String IMixtureManager::getIdModelName( String const& idData) const
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
IMixture* IMixtureManager::createMixture(String const& idData, int nbCluster)
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
void IMixtureManager::registerMixtureData(IMixtureData* p_data)
{ v_data_.push_back(p_data);}
/* release a data set from v_data_.
 *  @param idData name of the data set to release
 **/
void IMixtureManager::releaseMixtureData(String const& idData)
{
  for (DataIterator it = v_data_.begin(); it != v_data_.end(); ++it)
  { if ((*it)->idData() == idData) {delete (*it); v_data_.erase(it); break;}}
}

/* Utility lookup function allowing to find a MixtureData from its idData
 *  @param idData the id name of the mixture we want to get
 *  @return a pointer on the MixtureData
 **/
IMixtureData* IMixtureManager::getMixtureData( String const& idData) const
{
  for (ConstDataIterator it = v_data_.begin(); it != v_data_.end(); ++it)
  {  if ((*it)->idData() == idData) return (*it);}
  return 0;
}

} // namespace STK

