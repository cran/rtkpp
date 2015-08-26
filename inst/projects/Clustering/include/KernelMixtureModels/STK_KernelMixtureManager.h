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

/** @file STK_KernelMixtureManager.h
 *  @brief In this file we define the KernelMixtureManager class.
 **/


#ifndef STK_KERNELMIXTUREMANAGER_H
#define STK_KERNELMIXTUREMANAGER_H

#include "../STK_IMixtureManager.h"
#include "STK_KernelGaussianBridge.h"


#define STK_CREATE_MIXTURE(Data, Bridge) \
          Data* p_data = new Data(idData); \
          registerMixtureData(p_data); \
          p_handler()->getData(idData, p_data->dataij_, p_data->nbVariable_ ); \
          p_data->initialize(); \
          return new Bridge( p_data, idData, nbCluster);

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
class KernelMixtureManager : public IMixtureManager<DataHandler>
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

    // All Kernel bridges
    typedef KernelGaussianBridge<Clust::KernelGaussian_sk_,  DataReal> MixtureKernelGaussianBridge_sk;
    typedef KernelGaussianBridge<Clust::KernelGaussian_s_,  DataReal> MixtureKernelGaussianBridge_s;

    /** Default constructor, need an instance of a DataHandler.  */
    inline KernelMixtureManager(DataHandler const& handler): Base(&handler) {}
    /** destructor */
    inline virtual ~KernelMixtureManager() {}

    /** set the dimension of the kernel mixture model */
    void setDim(IMixture* p_mixture, Real const& dim)
    {
      if (!p_mixture) return;
      Clust::Mixture idModel = getIdModel(p_mixture->idName());
      // up-cast... (Yes it's bad....;)...)
      switch (idModel)
      {
        // Kernel models
        case Clust::KernelGaussian_sk_:
        { static_cast<MixtureKernelGaussianBridge_sk*>(p_mixture)->setDim(dim);}
        break;
        case Clust::KernelGaussian_s_:
        { static_cast<MixtureKernelGaussianBridge_s*>(p_mixture)->setDim(dim);}
        break;
        default: // idModel is not implemented
        break;
      }
    }
    /** set the dimension of the kernel mixture model */
    void setDim(IMixture* p_mixture, PointX const& dim)
    {
      if (!p_mixture) return;
      Clust::Mixture idModel = getIdModel(p_mixture->idName());
      // up-cast... (Yes it's bad....;)...)
      switch (idModel)
      {
        // Kernel models
        case Clust::KernelGaussian_sk_:
        { static_cast<MixtureKernelGaussianBridge_sk*>(p_mixture)->setDim(dim);}
        break;
        case Clust::KernelGaussian_s_:
        { static_cast<MixtureKernelGaussianBridge_s*>(p_mixture)->setDim(dim);}
        break;
        default: // idModel is not implemented
        break;
      }
    }
    /** get the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param param the array to return with the parameters
     **/
    void getParameters(IMixture* p_mixture, ArrayXX& param) const
    {
      if (!p_mixture) return;
      Clust::Mixture idModel = getIdModel(p_mixture->idName());
      // up-cast... (Yes it's bad....;)...)
      switch (idModel)
      {
        // Kernel models
        case Clust::KernelGaussian_sk_:
        { static_cast<MixtureKernelGaussianBridge_sk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::KernelGaussian_s_:
        { static_cast<MixtureKernelGaussianBridge_s*>(p_mixture)->getParameters(param);}
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
      if (!p_mixture) return;
      Clust::Mixture idModel = getIdModel(p_mixture->idName());
      // up-cast... (Yes it's bad....;)...)
      switch (idModel)
      {
        // Kernel models
        case Clust::KernelGaussian_sk_:
        { static_cast<MixtureKernelGaussianBridge_sk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::KernelGaussian_s_:
        { static_cast<MixtureKernelGaussianBridge_s*>(p_mixture)->setParameters(param);}
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
#ifdef STK_MIXTURE_VERY_VERBOSE
      stk_cout << _T("KernelMixtureManager::Entering createMixtureImpl(") <<  idModelName
               << _T(",") << idData << _T(",") << nbCluster << _T(")\n");
#endif
      Clust::Mixture idModel = Clust::stringToMixture(idModelName);
      return createMixtureImpl(idModel, idData, nbCluster);
    }
  private:
    /** create a concrete mixture and initialize it.
     *  @param idModel Id of the model
     *  @param idData Id of the data
     *  @param nbCluster number of cluster of the model
     **/
    IMixture* createMixtureImpl(Clust::Mixture idModel, String const& idData, int nbCluster)
    {
      switch (idModel)
      {
        case Clust::KernelGaussian_sk_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureKernelGaussianBridge_sk)}
        break;
        case Clust::KernelGaussian_s_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureKernelGaussianBridge_s)}
        break;
        default:
          return 0; // 0 if idModel is not implemented
          break;
      }
      return 0; // 0 if idModel is not a STK++ model
    }
};

} // namespace STK

#undef STK_CREATE_MIXTURE

#endif /* STK_KERNELMIXTUREMANAGER_H */
