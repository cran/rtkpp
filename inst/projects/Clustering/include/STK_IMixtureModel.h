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
#include <Arrays/include/STK_Array1D.h>
#include <Arrays/include/STK_Array2D.h>
#include <STatistiK/include/STK_Law_Categorical.h>

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

/** Parameters Handler class. All mixture models implemented has an Id defined
 *  in STK::Clust::Mixture enumeration.
 **/
template <int Id> struct ParametersHandler;

/**@ingroup Clustering
 * @brief Main interface class for mixture models.
 * At this level we create the array of Parameters.
 *
 * The pseudo virtual methods to implement in derived class are
 * @code
 * // default implementation (do nothing) provided to all these methods
 * void initializeModelImpl();
 * bool initializeStepImpl(); // return true by default
 * void finalizeStepImpl();
 * void setParametersImpl();
 * void storeIntermediateResultsImpl(int iter);
 * void releaseIntermediateResultsImpl();
 * @endcode
 *
 * @sa IMixtureModelBase, IRecursiveTemplate
 **/
template<class Derived>
class IMixtureModel: public IRecursiveTemplate<Derived>, public IMixtureModelBase
{
  public:
    typedef typename Clust::MixtureTraits<Derived>::Array Array;
    typedef typename Clust::MixtureTraits<Derived>::ParamHandler ParamHandler;

  protected:
    /** Default constructor.
     *  @param nbCluster the number of cluster
     **/
    inline IMixtureModel( int nbCluster)
                        : IMixtureModelBase(nbCluster)
                        , param_(nbCluster)
                        , p_dataij_(0)
    {}
    /** copy constructor.
     *  - The Parameter class is copied using the copy constructor.
     *  - The pointer on the data set is copied as-is. Check if you should not
     *  change it on the copied object.
     *  @param model the model to copy
     **/
    IMixtureModel( IMixtureModel const& model)
                 : IMixtureModelBase(model)
                 , param_(model.param_)
                 , p_dataij_(model.p_dataij_)
    {}

  public:
    /** destructor */
    inline ~IMixtureModel() {}
    /** create pattern.  */
    inline IMixtureModel* create() const { return new Derived(this->nbCluster());}
    /** @return a pointer on the current data set */
    inline Array const* p_data() const { return p_dataij_;}
    /** @return the parameter handler of the model */
    inline ParamHandler const& paramHandler() const { return param_;}

    /** set the parameter handler of the model */
    inline void setParamHandler(ParamHandler const& param) { param_ = param;}
    /** set the parmater handler using an array/expression storing the values */
    template<class Array>
    inline void setParamHandler(ExprBase<Array> const& param) { param_ = param;}
    /** @brief Set the data set.
     *  Setting a (new) data set will trigger the initialization process of the model.
     *  @param data the data set to set
     **/
    inline void setData(Array const& data)
    {
      p_dataij_ = &data;
      initializeModel();
    }
    /** @brief This function will be called once the model is created and data is set.
     *  @note a stk++ mixture create and initialize all the containers when the data
     *  is set. Thus the default behavior is @c return true.
     */
    inline bool initializeStep() { return this->asDerived().initializeStepImpl();}
    /** Store the intermediate results of the Mixture.
      *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    inline void storeIntermediateResults(int iteration)
    {
      param_.storeIntermediateResults(iteration);
      this->asDerived().storeIntermediateResultsImpl(iteration);
    }
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    inline void releaseIntermediateResults()
    {
      param_.releaseIntermediateResults();
      this->asDerived().releaseIntermediateResultsImpl();
    }
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    inline void setParameters()
    {
      param_.setParameters();
      this->asDerived().setParametersImpl();
    }
    /** @brief This function will be called once the model is estimated.
     *  perform specific model finalization stuff */
    inline void finalizeStep() { this->asDerived().finalizeStepImpl();}

    // default implementation of the pseudo-virtual methods
    /** default implementation of initializeModelImpl (do nothing) */
    inline void initializeModelImpl() {}
    /** default implementation of initializeStepImpl (return true) */
    inline bool initializeStepImpl() { return true;}
    /** default implementation of finalizeStepImpl (do nothing) */
    inline void finalizeStepImpl() {}
    /** default implementation of storeIntermediateResultsImpl (do nothing) */
    inline void storeIntermediateResultsImpl(int iteration) {}
    /** default implementation of setParametersImpl (do nothing) */
    inline void setParametersImpl() {}
    /** default implementation of releaseIntermediateResultsImpl (do nothing) */
    inline void releaseIntermediateResultsImpl() {}

    /** @return a simulated value for the jth variable of the ith sample
     *  @param i,j indexes of the data to simulate
     **/
    inline Real sample(int i, int j) const
    { return this->asDerived().rand(i, j, Law::Categorical::rand(p_tik()->row(i)));}

  protected:
    /** @return the parameter handler of the model */
    inline ParamHandler& paramHandler() { return param_;}

    /** @brief Initialize the model before its first use.
     * This function is triggered when data set is set.
     * In this interface, the @c initializeModel() method
     *  - set the number of samples and variables of the mixture model
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
      // call specific model initialization stuff
      this->asDerived().initializeModelImpl();
    }
    /** parameter handler associated with the derived mixture model */
    ParamHandler param_;

  private:
    /** pointer on the data set */
    Array const* p_dataij_;
};

} // namespace STK

#endif /* STK_IMIXTUREMODEL_H */
