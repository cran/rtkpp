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
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_MixtureBridge.h
 *  @brief In this file we define the bridge class between the mixtures and the composer.
 **/

#ifndef STK_MIXTUREBRIDGE_H
#define STK_MIXTUREBRIDGE_H

#include "../STK_IMixture.h"

#include "../GammaMixtureModels/STK_Gamma_ajk_bjk.h"
#include "../GammaMixtureModels/STK_Gamma_ajk_bk.h"
#include "../GammaMixtureModels/STK_Gamma_ajk_bj.h"
#include "../GammaMixtureModels/STK_Gamma_ajk_b.h"
#include "../GammaMixtureModels/STK_Gamma_ak_bjk.h"
#include "../GammaMixtureModels/STK_Gamma_ak_bk.h"
#include "../GammaMixtureModels/STK_Gamma_ak_bj.h"
#include "../GammaMixtureModels/STK_Gamma_ak_b.h"
#include "../GammaMixtureModels/STK_Gamma_aj_bjk.h"
#include "../GammaMixtureModels/STK_Gamma_aj_bk.h"
#include "../GammaMixtureModels/STK_Gamma_a_bjk.h"
#include "../GammaMixtureModels/STK_Gamma_a_bk.h"
#include "../GaussianMixtureModels/STK_Gaussian_sjk.h"
#include "../GaussianMixtureModels/STK_Gaussian_sk.h"
#include "../GaussianMixtureModels/STK_Gaussian_sj.h"
#include "../GaussianMixtureModels/STK_Gaussian_s.h"
#include "../CategoricalMixtureModels/STK_Categorical_pjk.h"
#include "../CategoricalMixtureModels/STK_Categorical_pk.h"

#include "../STK_DataManager.h"

namespace STK
{

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Categorical_pjk model
 **/
template<int Id, class Data> struct MixtureBridgeTraits;
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Categorical_pjk model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Categorical_pjk_, Data>
{
  /** Type of the Mixture model */
  typedef Categorical_pjk<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Categorical_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Categorical_pk model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Categorical_pk_, Data>
{
  /** Type of the mixture model */
  typedef Categorical_pk<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Categorical_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ajk_bjk_ model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gamma_ajk_bjk_, Data>
{
  /** Type of the mixture model */
  typedef Gamma_ajk_bjk<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gamma_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ajk_bk_ model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gamma_ajk_bk_, Data>
{
  /** Type of the mixture model */
  typedef Gamma_ajk_bk<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gamma_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ajk_bj_ model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gamma_ajk_bj_, Data>
{
  /** Type of the mixture model */
  typedef Gamma_ajk_bj<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gamma_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ajk_b_ model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gamma_ajk_b_, Data>
{
  /** Type of the mixture model */
  typedef Gamma_ajk_b<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gamma_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ak_bjk_ model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gamma_ak_bjk_, Data>
{
  /** Type of the mixture model */
  typedef Gamma_ak_bjk<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gamma_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ak_bk_ model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gamma_ak_bk_, Data>
{
  /** Type of the mixture model */
  typedef Gamma_ak_bk<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gamma_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ak_bj_ model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gamma_ak_bj_, Data>
{
  /** Type of the mixture model */
  typedef Gamma_ak_bj<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gamma_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ak_b_ model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gamma_ak_b_, Data>
{
  /** Type of the mixture model */
  typedef Gamma_ak_b<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gamma_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_aj_bjk_ model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gamma_aj_bjk_, Data>
{
  /** Type of the mixture model */
  typedef Gamma_aj_bjk<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gamma_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_aj_bk_ model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gamma_aj_bk_, Data>
{
  /** Type of the mixture model */
  typedef Gamma_aj_bk<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gamma_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_a_bjk_ model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gamma_a_bjk_, Data>
{
  /** Type of the mixture model */
  typedef Gamma_a_bjk<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gamma_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_a_bk_ model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gamma_a_bk_, Data>
{
  /** Type of the mixture model */
  typedef Gamma_a_bk<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gamma_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Gaussian_sjk model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gaussian_sjk_, Data>
{
  /** Type of the mixture model */
  typedef Gaussian_sjk<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gaussian_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Gaussian_sk model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gaussian_sk_, Data>
{
  /** Type of the mixture model */
  typedef Gaussian_sk<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gaussian_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Gaussian_sj model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gaussian_sj_, Data>
{
  /** Type of the mixture model */
  typedef Gaussian_sj<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gaussian_, Data> DataBridge;
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Gaussian_s model
 **/
template<class Data>
struct MixtureBridgeTraits<Clust::Gaussian_s_, Data>
{
  /** Type of the mixture model */
  typedef Gaussian_s<Data> Mixture;
  /** Type of the DataManager */
  typedef DataManager<Clust::Gaussian_, Data> DataBridge;
};

/** @ingroup hidden
 *  Initialize mixture, default implementation */
template<int Id, class Data>
struct InitializeMixtureImpl
{
  typedef typename MixtureBridgeTraits<Id, Data>::Mixture Mixture;
  typedef typename MixtureBridgeTraits<Id, Data>::DataBridge DataBridge;

  // call to IMixtureModel::setData will trigger a call to IMixtureModel::initializeModel
  // and call the specialization initializeModelImpl.
  // Thus it seems that at this level there is no more needs for
  // specialization of this template.
  static void run( Mixture& mixture, DataBridge* p_data)
  { mixture.setData(p_data->m_dataij());}
};

} // namespace hidden

/** @ingroup Clustering
 *  @brief Templated implementation of the IMixture interface allowing
 *  to bridge a STK++ mixture with the composer.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureModel class.
 */
template<int Id, class Data>
class MixtureBridge: public IMixture
{
  public:
    // type of Mixture
    typedef typename hidden::MixtureBridgeTraits<Id, Data>::Mixture Mixture;
    // type of Mixture
    typedef typename hidden::MixtureBridgeTraits<Id, Data>::DataBridge DataBridge;
    // parameters type to get
    typedef typename Clust::MixtureTraits<Mixture>::Param Param;

    /** default constructor.
     *  @param p_data pointer on the DataManager that will be used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    MixtureBridge( DataBridge* p_data, std::string const& idData, int nbCluster)
                 : IMixture( idData, nbCluster)
                 , mixture_( nbCluster)
                 , p_data_(p_data)
    { initializeMixture();}
    /** copy constructor */
    MixtureBridge( MixtureBridge const& mixture)
                 : IMixture(mixture)
                 , mixture_(mixture.mixture_)
                 , p_data_(mixture.p_data_)
    {  mixture_.setData(p_data_->m_dataij());}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual MixtureBridge* clone() const { return new MixtureBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     * @return New instance of class as that of calling object.
     */
    virtual MixtureBridge* create() const
    {
      MixtureBridge* p_bridge = new MixtureBridge( mixture_, idName(), nbCluster());
      p_bridge->p_data_ = p_data_;
      // Bug Fix: set the correct data set
      p_bridge->mixture_.setData(p_bridge->p_data_->m_dataij());
      return p_bridge;
    }
    /** @brief Initialize the mixture model before its use by the composer.
     *  The parameters values are set to their default values if the mixture_
     *  is newly created. if MixtureBridge::initializeStep is used during a
     *  cloning, mixture class have to take care of the existing values of the
     *  parameters.
     **/
    virtual void initializeStep()
    {
      if (!p_composer())
        STKRUNTIME_ERROR_NO_ARG(MixtureBridge::initializeStep,composer is not set);
      mixture_.setMixtureParameters( p_pk(), p_tik(), p_zi());
      if (!mixture_.initializeStep()) throw Clust::initializeStepFail_;
    }
     /** This function must be defined to return the component probability (PDF)
     *  for corresponding sample i and cluster k.
     * @param i,k Sample and Cluster numbers
     * @return the log-component probability
     */
    virtual double lnComponentProbability(int i, int k)
    { return mixture_.lnComponentProbability(i, k);}
    /** This function is equivalent to Mstep and must be defined to update
     * parameters.
     */
    virtual void paramUpdateStep()
    { if (!mixture_.mStep()) throw Clust::mStepFail_;}
    /** @brief This function should be used in order to initialize randomly the
     *  parameters of the mixture.
     */
    virtual void randomInit() { mixture_.randomInit();};
    /** This function should be used for imputation of data.
     *  The default implementation (in the base class) is to do nothing.
     */
    virtual void imputationStep()
    {
      typedef std::vector<std::pair<int,int> >::const_iterator ConstIterator;
      for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
      { p_data_->m_dataij_(it->first, it->second) = mixture_.impute(it->first, it->second);}
    }
    /** This function must be defined for simulation of all the latent variables
     * and/or missing data excluding class labels. The class labels will be
     * simulated by the framework itself because to do so we have to take into
     * account all the mixture laws. do nothing by default.
     */
    virtual void samplingStep()
    {
      typedef std::vector<std::pair<int,int> >::const_iterator ConstIterator;
      for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
      { p_data_->m_dataij_(it->first, it->second) = mixture_.sample(it->first, it->second);}
    }
    /** This function must return the number of free parameters.
     *  @return Number of free parameters
     */
    virtual int nbFreeParameter() const { return mixture_.computeNbFreeParameters();}
    /** This function must return the number of free parameters.
     *  @return Number of free parameters
     */
    virtual int nbVariable() const { return mixture_.nbVariable();}
    /** This function can be used to write summary of parameters to the output stream.
     * @param out Stream where you want to write the summary of parameters.
     */
    virtual void writeParameters(std::ostream& out) const
    { mixture_.writeParameters(out);}
    /** This function can be used in order to the values of the parameters
     *  in an Array2D.
     *  @param param the array with the parameters of the mixture.
     */
    void getParameters(Param& param) const { mixture_.getParameters(param);}

  private:
    /** This function will be used in order to initialize the mixture model
     *  using informations stored by the p_data_ container.
     **/
     void initializeMixture()
     { hidden::InitializeMixtureImpl<Id, Data>::run(mixture_, p_data_);}
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    MixtureBridge( Mixture const& mixture, std::string const& idData, int nbCluster)
                 : IMixture( idData, nbCluster)
                 , mixture_(mixture)
                 , p_data_(0)
    {}
    /** The Mixture to bridge with the composer */
    Mixture mixture_;
    /** Bridge for the data */
    DataBridge* p_data_;
};

} // namespace STK

#endif /* MIXTUREBRIDGE_H */
