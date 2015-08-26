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

/** @file STK_IMixtureBridge.h
 *  @brief In this file we define the interface for all the bridge classes
 *  between the mixtures and the composer.
 **/

#ifndef STK_IMIXTUREBRIDGE_H
#define STK_IMIXTUREBRIDGE_H

#include "STK_IMixture.h"

namespace STK
{

namespace hidden
{
/** @ingroup hidden
 *  MixtureBridgeTraits struct for bridged mixtures
 **/
template<class Derived> struct MixtureBridgeTraits;

}

/** @ingroup Clustering
 *  @brief Interface base class for the bridges of the STK++ mixture.
 *
 *  In this interface we give a default implementation for all
 *  the virtual functions required by the IMixture interface by delegating
 *  to the bridged mixture the computations.
 *
 *  The virtual methods to implement in derived class are
 *  @code
 *  virtual Derived* create() const;
 *  virtual Derived* clone() const;
 *  @endcode
 */
template<class Derived>
class IMixtureBridge: public IMixture, public IRecursiveTemplate<Derived>
{
  public:
    typedef std::vector<std::pair<int,int> >::const_iterator ConstIterator;
    typedef typename hidden::MixtureBridgeTraits<Derived>::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits<Derived>::Data Data;
    typedef typename hidden::MixtureBridgeTraits<Derived>::Parameters Parameters;
    typedef typename hidden::MixtureBridgeTraits<Derived>::ParamHandler ParamHandler;
    enum
    {
      // class of mixture
      idMixtureClass_ = hidden::MixtureBridgeTraits<Derived>::idMixtureClass_
    };

  protected:
    /** default constructor. Remove the missing values from the data set and
     *  initialize the mixture by setting the data set.
     *  @param p_data pointer on the MixtureData that will be used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    IMixtureBridge( MixtureData<Data>* p_data, std::string const& idData, int nbCluster)
                  : IMixture( idData, nbCluster)
                  , mixture_( nbCluster)
                  , p_data_(p_data)
    {}
    /** copy constructor */
    IMixtureBridge( IMixtureBridge const& bridge)
                  : IMixture(bridge)
                  , mixture_(bridge.mixture_)
                  , p_data_(bridge.p_data_)
    {}
    virtual ~IMixtureBridge() {}

  public:
    /** @brief Initialize the mixture model before its use by the composer.
     *  The parameters values are set to their default values if the mixture_
     *  is newly created. if IMixtureBridge::initializeStep is used during a
     *  cloning, mixture class have to take care of the existing values of the
     *  parameters.
     **/
    virtual void initializeStep()
    {
      if (!p_composer())
        STKRUNTIME_ERROR_NO_ARG(IMixtureBridge::initializeStep,composer is not set);
      mixture_.setMixtureParameters( p_pk(), p_nk(), p_tik(), p_zi());
      if (!mixture_.initializeStep()) throw Clust::initializeStepFail_;
    }
     /** This function must be defined to return the component probability (PDF)
     *  for corresponding sample i and cluster k.
     * @param i,k Sample and Cluster numbers
     * @return the log-component probability
     */
    virtual Real lnComponentProbability(int i, int k)
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
    /** This function must return the number of free parameters.
     *  @return Number of free parameters
     */
    virtual int nbFreeParameter() const { return mixture_.computeNbFreeParameters();}
    /** This function must return the number of variables.
     *  @return Number of variables
     */
    virtual int nbVariable() const { return mixture_.nbVariable();}
    /** @brief This function should be used to store any intermediate results
     * during various iterations after the burn-in period.
     * @param iteration Provides the iteration number beginning after the burn-in
     * period.
     */
    virtual void storeIntermediateResults(int iteration)
    { mixture_.storeIntermediateResults(iteration);}
    /**@brief This step can be used to signal to the mixtures that they must
     * release the stored results. This is usually called if the estimation
     * process failed.
     **/
    virtual void releaseIntermediateResults()
    { mixture_.releaseIntermediateResults();}
    /** @brief set the parameters of the model.
     *  This function should be used to set the parameters computed using the
     *  intermediate results.
     **/
    virtual void setParameters() { mixture_.setParameters();}
    /** @brief This step can be used by developer to finalize any thing. It will
     *  be called only once after we finish running the estimation algorithm.
     */
    virtual void finalizeStep() { mixture_.finalizeStep();}
    /** @brief This function should be used for Imputation of data.
     *  The default implementation (in the base class) is to do nothing.
     */
    virtual void imputationStep();
    /** @brief This function must be defined for simulation of all the latent
     * variables and/or missing data excluding class labels. The class labels
     * will be simulated by the framework itself because to do so we have to
     * take into account all the mixture laws.
     */
    virtual void samplingStep();

    /** @return the structure storing the current values of the parameters. */
    ParamHandler const& paramHandler() const { return mixture_.paramHandler();}
    /** set the parameter handler (current values) of the model */
    void setParamHandler(ParamHandler const& param) { mixture_.setParamHandler(param);}
    /** This function is used in order to set the current values of the
     *  parameters to the paramHandler of the mixture.
     *  @param param the array/expression with the parameters of the mixture to
     *  store in the ParamHandler.
     */
    template<class Array>
    void setParameters( Array const& param)
    {
//      ParamHandler handler(this->nbCluster(), param);
      mixture_.setParamHandler(param);
    }

  protected:
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    IMixtureBridge( Mixture const& mixture, std::string const& idData, int nbCluster)
                  : IMixture( idData, nbCluster)
                  , mixture_(mixture)
                  , p_data_(0)
    {}
    /** The Mixture to bridge with the composer */
    Mixture mixture_;
    /** pointer on the data manager */
    MixtureData<Data>* p_data_;
};

// implementation
template< class Derived>
void IMixtureBridge<Derived>::imputationStep()
{
  for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
  { p_data_->dataij_(it->first, it->second) = mixture_.impute(it->first, it->second);}
}
// implementation
template< class Derived>
void IMixtureBridge<Derived>::samplingStep()
{
  for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
  { p_data_->dataij_(it->first, it->second) = mixture_.sample(it->first, it->second);}
}

} // namespace STK

#endif /* IMIXTUREBRIDGE_H */
