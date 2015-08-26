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

/** @file STK_GammaMixtureManager.h
 *  @brief In this file we define the GammaMixtureManager class.
 **/


#ifndef STK_GAMMAMIXTUREMANAGER_H
#define STK_GAMMAMIXTUREMANAGER_H

#include "../STK_IMixtureManager.h"
#include "STK_GammaBridge.h"

#define STK_CREATE_MIXTURE(Data, Bridge) \
          Data* p_data = new Data(idData); \
          registerMixtureData(p_data); \
          p_handler()->getData(idData, p_data->dataij_, p_data->nbVariable_ ); \
          p_data->initialize(); \
          Bridge* p_bridge = new Bridge( p_data, idData, nbCluster);  \
          return p_bridge;

namespace STK
{
/** @ingroup Clustering
 *  @brief A mixture manager is a factory class for injection dependency in the
 *  STK++ derived class of the MixtureComposer class.
 *
 *  It allows to handle all the creation and initialization stuff needed by the
 *  (bridged) mixture models of the stkpp library.
 *
 *  @tparam DataHandler is any concrete class from the interface DataHandlerBase
 */
template<class DataHandler>
class GammaMixtureManager : public IMixtureManager<DataHandler>
{
  public:
    typedef IMixtureManager<DataHandler> Base;
    using Base::registerMixtureData;
    using Base::getMixtureData;
    using Base::getIdModel;
    using Base::p_handler;

    // All data handlers will store and return a specific container for
    // the data they handle. The DataHandlerTraits class allow us to know the
    // type of these containers when data is Real and Integer.
    typedef typename hidden::DataHandlerTraits<DataHandler, Real>::Data DataReal;
    // Classes wrapping the Real and Integer containers
    typedef MixtureData<DataReal> MixtureDataReal;

    // All Gamma bridges
    typedef GammaBridge<Clust::Gamma_ajk_bjk_, DataReal> MixtureBridge_ajk_bjk;
    typedef GammaBridge<Clust::Gamma_ajk_bk_,  DataReal> MixtureBridge_ajk_bk;
    typedef GammaBridge<Clust::Gamma_ajk_bj_,  DataReal> MixtureBridge_ajk_bj;
    typedef GammaBridge<Clust::Gamma_ajk_b_,   DataReal> MixtureBridge_ajk_b;
    typedef GammaBridge<Clust::Gamma_ak_bjk_,  DataReal> MixtureBridge_ak_bjk;
    typedef GammaBridge<Clust::Gamma_ak_bk_,   DataReal> MixtureBridge_ak_bk;
    typedef GammaBridge<Clust::Gamma_ak_bj_,   DataReal> MixtureBridge_ak_bj;
    typedef GammaBridge<Clust::Gamma_ak_b_,    DataReal> MixtureBridge_ak_b;
    typedef GammaBridge<Clust::Gamma_aj_bjk_,  DataReal> MixtureBridge_aj_bjk;
    typedef GammaBridge<Clust::Gamma_aj_bk_,   DataReal> MixtureBridge_aj_bk;
    typedef GammaBridge<Clust::Gamma_a_bjk_,   DataReal> MixtureBridge_a_bjk;
    typedef GammaBridge<Clust::Gamma_a_bk_,    DataReal> MixtureBridge_a_bk;

    /** Default constructor, need an instance of a DataHandler.  */
    GammaMixtureManager(DataHandler const& handler): Base(&handler) {}
    /** destructor */
    virtual ~GammaMixtureManager() {}
    /** get the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param param the array to return with the parameters
     **/
    void getParameters(IMixture* p_mixture, ArrayXX& param) const
    {
      Clust::Mixture idModel = getIdModel(p_mixture->idName());
      if (idModel == Clust::unknown_mixture_) return;
      // up-cast... (Yes it's bad....;)...)
      switch (idModel)
      {
        // gamma models
        case Clust::Gamma_ajk_bjk_:
        { static_cast<MixtureBridge_ajk_bjk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ajk_bk_:
        { static_cast<MixtureBridge_ajk_bk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ajk_bj_:
        { static_cast<MixtureBridge_ajk_bj*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ajk_b_:
        { static_cast<MixtureBridge_ajk_b*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ak_bjk_:
        { static_cast<MixtureBridge_ak_bjk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ak_bk_:
        { static_cast<MixtureBridge_ak_bk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ak_bj_:
        { static_cast<MixtureBridge_ak_bj*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ak_b_:
        { static_cast<MixtureBridge_ak_b*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_aj_bjk_:
        { static_cast<MixtureBridge_aj_bjk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_aj_bk_:
        { static_cast<MixtureBridge_aj_bk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_a_bjk_:
        { static_cast<MixtureBridge_a_bjk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_a_bk_:
        { static_cast<MixtureBridge_a_bk*>(p_mixture)->getParameters(param);}
        break;
        default: // idModel is not implemented
        break;
      }
    }
    /** set the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param param the array with the parameters to set
     **/
    virtual void setParameters(IMixture* p_mixture, ArrayXX const& param) const
    {
      Clust::Mixture idModel = getIdModel(p_mixture->idName());
      if (idModel == Clust::unknown_mixture_) return;
      // up-cast... (Yes it's bad....;)...)
      switch (idModel)
      {
        // gamma models
        case Clust::Gamma_ajk_bjk_:
        { static_cast<MixtureBridge_ajk_bjk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ajk_bk_:
        { static_cast<MixtureBridge_ajk_bk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ajk_bj_:
        { static_cast<MixtureBridge_ajk_bj*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ajk_b_:
        { static_cast<MixtureBridge_ajk_b*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ak_bjk_:
        { static_cast<MixtureBridge_ak_bjk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ak_bk_:
        { static_cast<MixtureBridge_ak_bk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ak_bj_:
        { static_cast<MixtureBridge_ak_bj*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ak_b_:
        { static_cast<MixtureBridge_ak_b*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_aj_bjk_:
        { static_cast<MixtureBridge_aj_bjk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_aj_bk_:
        { static_cast<MixtureBridge_aj_bk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_a_bjk_:
        { static_cast<MixtureBridge_a_bjk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_a_bk_:
        { static_cast<MixtureBridge_a_bk*>(p_mixture)->setParameters(param);}
        break;
        default: // idModel is not implemented
        break;
      }
    }

  protected:
    /** create a concrete mixture and initialize it.
     *  @param idModelName, idData Id names of the model and of the data
     *  @param nbCluster number of cluster of the model
     **/
    virtual IMixture* createMixtureImpl(String const&  idModelName, String const& idData, int nbCluster)
    {
      Clust::Mixture idModel = Clust::stringToMixture(idModelName);
      return createMixtureImpl(idModel, idData, nbCluster);
    }
  private:
    /** create a concrete mixture and initialize it.
     *  @param idModel Id name of the model
     *  @param idData Id name of the data
     *  @param nbCluster number of cluster of the model
     **/
    IMixture* createMixtureImpl(Clust::Mixture idModel, String const& idData, int nbCluster)
    {
      switch (idModel)
      {
        // gamma models
        case Clust::Gamma_ajk_bjk_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_ajk_bjk)}
        break;
        case Clust::Gamma_ajk_bk_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_ajk_bk)}
        break;
        case Clust::Gamma_ajk_bj_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_ajk_bj)}
        break;
        case Clust::Gamma_ajk_b_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_ajk_b)}
        break;
        case Clust::Gamma_ak_bjk_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_ak_bjk)}
        break;
        case Clust::Gamma_ak_bk_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_ak_bk)}
        break;
        case Clust::Gamma_ak_bj_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_ak_bj)}
        break;
        case Clust::Gamma_ak_b_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_ak_b)}
        break;
        case Clust::Gamma_aj_bjk_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_aj_bjk)}
        break;
        case Clust::Gamma_aj_bk_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_aj_bk)}
        case Clust::Gamma_a_bjk_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_a_bjk)}
        break;
        case Clust::Gamma_a_bk_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_a_bk)}
        default:
          return 0; // 0 if idModel is not implemented
          break;
      }
      return 0; // 0 if idModel is not a STK++ model
    }
};

} // namespace STK

#undef STK_CREATE_MIXTURE

#endif /* STK_GAMMAMIXTUREMANAGER_H */
