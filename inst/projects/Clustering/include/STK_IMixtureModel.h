/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2012  Serge Iovleff

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
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_IMixtureModel.h
 *  @brief In this file we define the main interface class for mixture models.
 **/


#ifndef STK_IMIXTUREMODEL_H
#define STK_IMIXTUREMODEL_H

#include "STK_IMixtureModelBase.h"
#include "Arrays/include/STK_Array1D.h"
#include "Arrays/include/STK_Array2D.h"
#ifdef STK_MIXTURE_DEBUG
#include "Arrays/include/STK_Display.h"
#endif

namespace STK
{
namespace Clust
{
/** Main class for the mixtures traits policy.
 *  The traits struct MixtureTraits must be specialized for any
 *  Mixture deriving from the Interface IMixtureModel.
 **/
template <class Mixture> struct MixtureTraits;

} // namespace Clust


/**@ingroup Clustering
 * @brief Main interface class for mixture models.
 * At this level we create the array of Components.
 *
 * The pseudo virtual methods to implement in derived class are
 * @code
 * Param getParametersImpl() const;
 * // default implementation (do nothing) provided to all these methods
 * void initializeModelImpl();
 * bool initializeStepImpl(); // return true by default
 * void finalizeModelImpl();
 * void setParametersImpl();
 * void storeIntermediateResultsImpl(int iter);
 * void releaseIntermediateResultsImpl();
 * @endcode
 *
 * @sa IMixtureModelBase, IRecursiveTemplate
 **/
template<class Derived>
class IMixtureModel : public IRecursiveTemplate<Derived>, public IMixtureModelBase
{
  public:
    typedef typename Clust::MixtureTraits<Derived>::Array Array;
    typedef typename Clust::MixtureTraits<Derived>::Parameters Parameters;
    typedef typename Clust::MixtureTraits<Derived>::Param Param;

  protected:
    /** Default constructor.
     *  Set the number of cluster and create the parameters struct.
     **/
    IMixtureModel( int nbCluster)
                 : IMixtureModelBase(nbCluster)
                 , p_dataij_(0)
                 , components_(nbCluster, 0)
    {
      for (int k= components_.begin(); k < components_.end(); ++k)
      { components_[k] = new Parameters();}
    }
    /** copy constructor.
     *  Call the clone method of the Parameters class.
     *  The pointer on the data set is copied as-is. Just check if you should not
     *  change it on the copied object.
     *  @param model the model to copy
     **/
    IMixtureModel( IMixtureModel const& model)
                 : IMixtureModelBase(model)
                 , p_dataij_(model.p_dataij_)
                 , components_(model.components_)
    {
      for (int k= components_.begin(); k < components_.end(); ++k)
      { components_[k] = model.components_[k]->clone(); }
    }

  public:
    /** destructor */
    ~IMixtureModel()
    {
      for (int k= components_.begin(); k < components_.end(); ++k)
      { delete components_[k];}
    }
    /** create pattern.  */
    inline IMixtureModel* create() const { return new Derived(this->nbCluster());}
    /** @return a pointer on the current data set */
    inline Array const* p_data() const { return p_dataij_;}
    /** @return a constant reference on the k-th parameter */
    inline Parameters const* p_param(int k) const { return components_[k];}
    /** Set the data set and initialize the model.
     *  @param data the data set to set*/
    void setData(Array const& data)
    {
      p_dataij_ = &data;
      initializeModel();
    }
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k)
    { return components_[k]->computeLnLikelihood(p_dataij_->row(i));}
    /** @brief This function will be called once the model is created and data is set.
     *  @note a stk++ mixture create and initialize all the containers when the data
     *  is set. Thus the default behavior is @c return true.
     */
    bool initializeStep() { return this->asDerived().initializeStepImpl();}
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    {
      for (int k= components_.begin(); k < components_.end(); ++k)
      { components_[k]->storeIntermediateResults(iteration); }
      this->asDerived().storeIntermediateResultsImpl(iteration);
    }
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    {
      for (int k= components_.begin(); k < components_.end(); ++k)
      { components_[k]->releaseIntermediateResults(); }
      // release shared parameters staticitcs if any
      this->asDerived().releaseIntermediateResultsImpl();
    }
    /** call specific model finalization stuff */
    void finalizeStep() { this->asDerived().finalizeModelImpl();}
    /** set the parameters stored in the Stat* structures. */
    void setParameters()
    {
      for (int k= components_.begin(); k < components_.end(); ++k)
      { components_[k]->setParameters(); }
      // set shared parameters if any
      this->asDerived().setParametersImpl();
    }
    /**  @return the parameters of the mixture in the structure Param */
    Param getParameters() const
    { return this->asDerived().getParametersImpl();}

    /** default implementation of initializeModelImpl (do nothing) */
    void initializeModelImpl() {}
    /** default implementation of initializeStepImpl (return true) */
    bool initializeStepImpl() { return true;}
    /** default implementation of finalizeModelImpl (do nothing) */
    void finalizeModelImpl() {}
    /** default implementation of storeIntermediateResultsImpl (do nothing) */
    void storeIntermediateResultsImpl(int iteration) {}
    /** default implementation of setParametersImpl (do nothing) */
    void setParametersImpl() {}
    /** default implementation of releaseIntermediateResultsImpl (do nothing) */
    void releaseIntermediateResultsImpl() {}

  protected:
    /** @brief Initialize the model before its first use.
     * This function is used when the data is set.
     * In this interface, the @c initializeModel() method
     *  - set the number of samples and variables of the mixture model
     *  - resize the parameters of each component with the range of the variables
     *  - call the derived class implemented method
     * @code
     *   initializeModelImpl()
     * @endcode
     * for initialization of the specific model parameters if needed.
     **/
    void initializeModel()
    {
      // set dimensions
      this->setNbSample(p_dataij_->sizeRows());
      this->setNbVariable(p_dataij_->sizeCols());
      // initialize the parameters
      for (int k= components_.begin(); k < components_.end(); ++k)
      { components_[k]->resize(p_dataij_->cols());}
      // call specific model initialization stuff
      this->asDerived().initializeModelImpl();
    }
    /** @return the array with the components */
    inline Array1D<Parameters*> const& components() const { return components_;}
    /** @return the array with the components */
    inline Array1D<Parameters*>& components() { return components_;}
    /** @return a pointer on the k-th parameter */
    inline Parameters* p_param(int k) { return components_[k];}

  private:
    /** pointer on the data set */
    Array const* p_dataij_;
    /** Array of the components of the mixture model */
    Array1D< Parameters* > components_;
};

} // namespace STK

#endif /* STK_IMIXTUREMODEL_H */
