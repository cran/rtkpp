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

/** @file STK_MixtureManager.h
 *  @brief In this file we define the MixtureManager class.
 **/


#ifndef STK_MIXTUREMANAGER_H
#define STK_MIXTUREMANAGER_H

#include <DManager/include/STK_DataHandlerBase.h>

#include "CategoricalMixtureModels/STK_CategoricalMixtureManager.h"
#include "GammaMixtureModels/STK_GammaMixtureManager.h"
#include "DiagGaussianMixtureModels/STK_DiagGaussianMixtureManager.h"
#include "PoissonMixtureModels/STK_PoissonMixtureManager.h"
#include "KernelMixtureModels/STK_KernelMixtureManager.h"

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
 *  @class MixtureManager
 *  @brief [DEPRECATED] A mixture manager is a factory class for injection dependency in the
 *  STK++ derived class of the MixtureComposer class.
 *
 *  It allows to handle all the creation and initialization stuff needed by the
 *  (bridged) mixture models of the stkpp library.
 *
 *  @tparam DataHandler is any concrete class from the interface DataHandlerBase
 */
template<class DataHandler>
class MixtureManager : public IMixtureManager<DataHandler>
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
    typedef typename hidden::DataHandlerTraits<DataHandler, Integer>::Data DataInt;
    // Classes wrapping the Real and Integer containers
    typedef MixtureData<DataReal> MixtureDataReal;
    typedef MixtureData<DataInt>  MixtureDataInt;

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
    // All Gaussian bridges
    typedef DiagGaussianBridge<Clust::Gaussian_sjk_, DataReal> MixtureBridge_sjk;
    typedef DiagGaussianBridge<Clust::Gaussian_sk_,  DataReal> MixtureBridge_sk;
    typedef DiagGaussianBridge<Clust::Gaussian_sj_,  DataReal> MixtureBridge_sj;
    typedef DiagGaussianBridge<Clust::Gaussian_s_,   DataReal> MixtureBridge_s;
    // All Categorical bridges
    typedef CategoricalBridge<Clust::Categorical_pjk_, DataInt> MixtureBridge_pjk;
    typedef CategoricalBridge<Clust::Categorical_pk_,  DataInt> MixtureBridge_pk;
    // All Poisson bridges
    typedef PoissonBridge<Clust::Poisson_ljk_,  DataInt> MixtureBridge_ljk;
    typedef PoissonBridge<Clust::Poisson_lk_,  DataInt> MixtureBridge_lk;
    typedef PoissonBridge<Clust::Poisson_ljlk_, DataInt> MixtureBridge_ljlk;
    // All Kernel bridges
    typedef KernelGaussianBridge<Clust::KernelGaussian_sk_,  DataReal> MixtureKernelGaussianBridge;

    /** Default constructor, need an instance of a DataHandler.  */
    MixtureManager(DataHandler const& handler) : Base(&handler) {}
    /** destructor */
    virtual ~MixtureManager() {}
    /** Utility function allowing to create and register all the STK++ mixtures
     *  defined in the handler. [DEPRECATED]
     *  @param composer the composer claiming the mixtures
     **/
    void createMixtures(MixtureComposer& composer);
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
        // Gaussian models
        case Clust::Gaussian_sjk_:
        { static_cast<MixtureBridge_sjk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gaussian_sk_:
        { static_cast<MixtureBridge_sk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gaussian_sj_:
        { static_cast<MixtureBridge_sj*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gaussian_s_:
        { static_cast<MixtureBridge_s*>(p_mixture)->getParameters(param);}
        break;
        // Categorical models
        case Clust::Categorical_pjk_:
        { static_cast<MixtureBridge_pjk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Categorical_pk_:
        { static_cast<MixtureBridge_pk*>(p_mixture)->getParameters(param);}
        break;
        // Poisson models
        case Clust::Poisson_ljk_:
        { static_cast<MixtureBridge_ljk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Poisson_lk_:
        { static_cast<MixtureBridge_lk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Poisson_ljlk_:
        { static_cast<MixtureBridge_ljlk*>(p_mixture)->getParameters(param);}
        break;
        // Kernel models
        case Clust::KernelGaussian_sk_:
        { static_cast<MixtureKernelGaussianBridge*>(p_mixture)->getParameters(param);}
        break;
        default: // idModel is not implemented
        break;
      }
    }
    /** set the parameters to an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param param the array to return with the parameters
     **/
    void setParameters(IMixture* p_mixture, ArrayXX const& param) const
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
        // Gaussian models
        case Clust::Gaussian_sjk_:
        { static_cast<MixtureBridge_sjk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gaussian_sk_:
        { static_cast<MixtureBridge_sk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gaussian_sj_:
        { static_cast<MixtureBridge_sj*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gaussian_s_:
        { static_cast<MixtureBridge_s*>(p_mixture)->setParameters(param);}
        break;
        // Categorical models
        case Clust::Categorical_pjk_:
        { static_cast<MixtureBridge_pjk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Categorical_pk_:
        { static_cast<MixtureBridge_pk*>(p_mixture)->setParameters(param);}
        break;
        // Poisson models
        case Clust::Poisson_ljk_:
        { static_cast<MixtureBridge_ljk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Poisson_lk_:
        { static_cast<MixtureBridge_lk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Poisson_ljlk_:
        { static_cast<MixtureBridge_ljlk*>(p_mixture)->setParameters(param);}
        break;
        // Kernel models
        case Clust::KernelGaussian_sk_:
        { static_cast<MixtureKernelGaussianBridge*>(p_mixture)->setParameters(param);}
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
        // Gaussian models
        case Clust::Gaussian_sjk_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_sjk)}
        break;
        case Clust::Gaussian_sk_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_sk)}
        break;
        case Clust::Gaussian_sj_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_sj)}
        break;
        case Clust::Gaussian_s_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureBridge_s)}
        break;
        // Categorical models
        case Clust::Categorical_pjk_:
        { STK_CREATE_MIXTURE(MixtureDataInt, MixtureBridge_pjk)}
        break;
        case Clust::Categorical_pk_:
        { STK_CREATE_MIXTURE(MixtureDataInt, MixtureBridge_pk)}
        break;
        // Poisson models
        case Clust::Poisson_ljk_:
        { STK_CREATE_MIXTURE(MixtureDataInt, MixtureBridge_ljk)}
        break;
        case Clust::Poisson_lk_:
        { STK_CREATE_MIXTURE(MixtureDataInt, MixtureBridge_lk)}
        break;
        case Clust::Poisson_ljlk_:
        { STK_CREATE_MIXTURE(MixtureDataInt, MixtureBridge_ljlk)}
        break;
        // Kernel models
        case Clust::KernelGaussian_sk_:
        { STK_CREATE_MIXTURE(MixtureDataReal, MixtureKernelGaussianBridge)}
        break;
        default:
          return 0; // 0 if idModel is not implemented
          break;
      }
      return 0; // 0 if idModel is not a STK++ model
    }
};

/* Utility function allowing to create and register all the STK++ mixtures
 *  defined in the handler.
 *  @param composer the composer claiming the mixtures
 *  @param nbCluster the number of clusters
 **/
template<class DataHandler>
void MixtureManager<DataHandler>::createMixtures(MixtureComposer& composer)
{
  int nbCluster = composer.nbCluster();
  typedef typename DataHandlerBase<DataHandler>::InfoMap InfoMap;
  for ( typename InfoMap::const_iterator it=p_handler()->info().begin(); it!=p_handler()->info().end(); ++it)
  {
    IMixture* p_mixture = createMixtureImpl(it->second, it->first, nbCluster);
    if (p_mixture) composer.registerMixture(p_mixture);
  }
}

} // namespace STK

#undef STK_CREATE_MIXTURE

#endif /* STK_MIXTUREMANAGER_H */
