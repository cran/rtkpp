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
 *  @brief In this file we define the main class for mixture models.
 **/


#ifndef STK_IMIXTUREMODEL_H
#define STK_IMIXTUREMODEL_H

#include "STK_IMixtureModelBase.h"
#include "STK_MixtureComponent.h"
#include "Arrays/include/STK_Array1D.h"

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
 * At this level we create the array of Components
 *
 * @sa MixtureComponent, IMixtureModelBase, IMultiParameters
 **/
template<class Derived>
class IMixtureModel : public IRecursiveTemplate<Derived>, public IMixtureModelBase
{
  public:
    typedef typename Clust::MixtureTraits<Derived>::Array Array;
    typedef typename Clust::MixtureTraits<Derived>::Component Component;
    typedef typename Clust::MixtureTraits<Derived>::Parameters Parameters;

  protected:
    /** Default constructor.
     * Set the number of cluster and create components with zero pointer on data.
     **/
    IMixtureModel( int nbCluster)
                 : IMixtureModelBase(nbCluster), p_dataij_(0), components_(nbCluster, 0)
    {
      for (int k= components_.begin(); k < components_.end(); ++k)
      { components_[k] = new Component();}
    }
    /** copy constructor.
     *  Call the clone method of the Components class.
     *  The pointer on the data set is copied as-is. Just check if you should no
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
    /** @return the array with the components */
    inline Array1D< Component* > const& components() const { return components_;}
    /** @return a constant reference on the k-th component */
    inline Component const* components(int k) const { return components_[k];}
    /** @return a constant reference on the k-th parameter */
    inline Parameters const* p_param(int k) const { return components_[k]->p_param();}
    /** Set a new data set and initialize the model.
     *  @param data the data set to set*/
    void setData(Array const& data)
    {
      p_dataij_ = &data;
      initializeModel();
    }
    /** default implementation of initializeModelImpl */
    void initializeModelImpl() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k)
    { return components_[k]->computeLnLikelihood(p_dataij_->row(i));}

  protected:
    /** @brief Initialize the model before its first use.
     * This function will be called when the data is set.
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
      if (!p_dataij_)
        STKRUNTIME_ERROR_NO_ARG(IMixtureModel::initializeModel,p_dataij_ is not set);
      // set dimensions
      this->setNbSample(p_dataij_->sizeRows());
      this->setNbVariable(p_dataij_->sizeCols());
      // initialize the parameters
      for (int k= components_.begin(); k < components_.end(); ++k)
      { components_[k]->p_param()->resize(p_dataij_->cols());}
      // call specific model initialization stuff
      this->asDerived().initializeModelImpl();
    }
    /** @return the array with the components */
    inline Array1D<Component*>& components() { return components_;}
    /** @return a pointer on the k-th parameter */
    inline Parameters* p_param(int k) { return components_[k]->p_param();}
    /** pointer on the data set */
    Array const* p_dataij_;
    /** Array of the components of the mixture model */
    Array1D< Component* > components_;
};

} // namespace STK

#endif /* STK_IMIXTUREMODEL_H */
