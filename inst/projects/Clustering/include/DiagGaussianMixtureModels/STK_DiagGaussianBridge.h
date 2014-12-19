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

/** @file STK_DiagGaussianBridge.h
 *  @brief In this file we define the bridge classes between the diagonal
 *  Gaussian mixtures and the composer.
 **/

#ifndef STK_DIAGGAUSSIANBRIDGE_H
#define STK_DIAGGAUSSIANBRIDGE_H

#include "STK_Gaussian_s.h"
#include "STK_Gaussian_sj.h"
#include "STK_Gaussian_sjk.h"
#include "STK_Gaussian_sk.h"

#include "../STK_MixtureData.h"
#include "../STK_IMixtureBridge.h"

namespace STK
{

// forward declaration
template<int Id, class Data> class DiagGaussianBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Gaussian_sjk model
 **/
template<class Data>
struct MixtureBridgeTraits< DiagGaussianBridge< Clust::Gaussian_sjk_, Data> >
{
  /** Type of the mixture model */
  typedef Gaussian_sjk<Data> Mixture;
  enum
  {
    idMixtureClass_ = Clust::Gaussian_
  };
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Gaussian_sk model
 **/
template<class Data>
struct MixtureBridgeTraits< DiagGaussianBridge< Clust::Gaussian_sk_, Data> >
{
  /** Type of the mixture model */
  typedef Gaussian_sk<Data> Mixture;
  enum
  {
    idMixtureClass_ = Clust::Gaussian_
  };
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Gaussian_sj model
 **/
template<class Data>
struct MixtureBridgeTraits< DiagGaussianBridge< Clust::Gaussian_sj_, Data> >
{
  /** Type of the mixture model */
  typedef Gaussian_sj<Data> Mixture;
  enum
  {
    idMixtureClass_ = Clust::Gaussian_
  };
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Gaussian_s model
 **/
template<class Data>
struct MixtureBridgeTraits< DiagGaussianBridge< Clust::Gaussian_s_, Data> >
{
  /** Type of the mixture model */
  typedef Gaussian_s<Data> Mixture;
  enum
  {
    idMixtureClass_ = Clust::Gaussian_
  };
};

} // namespace hidden

} // namespace STK

namespace STK
{
/** @ingroup Clustering
 *  @brief Templated implementation of the IMixture interface allowing
 *  to bridge a STK++ mixture with the composer.
 *
 *  This class inherit from the interface IMixture and delegate almost
 *  all the treatments to the wrapped class.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureModel class.
 */
template<int Id, class Data>
class DiagGaussianBridge: public IMixtureBridge< DiagGaussianBridge<Id,Data> >
{
  public:
    // Base class
    typedef IMixtureBridge< DiagGaussianBridge<Id,Data> > Base;
    // type of Mixture
    typedef typename hidden::MixtureBridgeTraits< DiagGaussianBridge<Id,Data> >::Mixture Mixture;
    // type of data
    typedef typename Data::Type Type;
    // class of mixture
    enum
    {
      idMixtureClass_ = Clust::Gaussian_
    };
    // parameters type to get
    typedef typename Clust::MixtureTraits<Mixture>::Param Param;

    typedef std::vector<std::pair<int,int> >::const_iterator ConstIterator;
    using Base::mixture_;

    /** default constructor. Remove the missing values from the data set and
     *  initialize the mixture by setting the data set.
     *  @param p_data pointer on the MixtureData that will be used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    DiagGaussianBridge( MixtureData<Data>* p_data, std::string const& idData, int nbCluster)
                  : Base( idData, nbCluster)
                  , p_data_(p_data)
    { removeMissing(); initializeMixture();}
    /** copy constructor */
    DiagGaussianBridge( DiagGaussianBridge const& bridge)
                  : Base(bridge)
                  , p_data_(bridge.p_data_)
    { initializeMixture();}
    /** destructor */
    virtual ~DiagGaussianBridge() {}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual DiagGaussianBridge* clone() const { return new DiagGaussianBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual DiagGaussianBridge* create() const
    {
      DiagGaussianBridge* p_bridge = new DiagGaussianBridge( mixture_, this->idName(), this->nbCluster());
      p_bridge->p_data_ = p_data_;
      // Bug Fix: set the correct data set
      p_bridge->mixture_.setData(p_bridge->p_data_->dataij());
      return p_bridge;
    }
    /** This function is used in order to get the current values of the
     *  parameters.
     *  @param params the array with the parameters of the mixture.
     */
    virtual void getParameters(ArrayXX& params) const;

  private:
    /** This function will be used for the imputation of the missing data
     *  at the initialization.
     **/
    void removeMissing();
    /** This function will be used in order to initialize the mixture model
     *  using informations stored by the MixtureData. For example the missing
     *  values in the case of a MixtureData instance.
     **/
    void initializeMixture() { mixture_.setData(p_data_->dataij());}
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    DiagGaussianBridge( Mixture const& mixture, std::string const& idData, int nbCluster)
                 : Base(mixture, idData, nbCluster)
                 , p_data_(0)
    {}
    /** pointer on the data manager */
    MixtureData<Data>* p_data_;
};

// implementation
template<int Id, class Data>
void DiagGaussianBridge<Id, Data>::removeMissing()
{
  Type value = Type();
  int j, old_j = Arithmetic<int>::NA();
  for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
  {
    j = it->second; // get column
    if (j != old_j)
    {
      old_j = j;
      value =  p_data_->dataij_.col(j).safe().mean();
    }
    p_data_->dataij_(it->first, j) = value;
  }
}

template<int Id, class Data>
void DiagGaussianBridge<Id, Data>::getParameters(ArrayXX& params) const
{
  int nbClust = this->nbCluster();
  params.resize(2*nbClust, mixture_.p_data()->cols());
  for (int k= 0; k < nbClust; ++k)
  {
    for (int j= mixture_.p_data()->beginCols();  j< mixture_.p_data()->endCols(); ++j)
    {
      params(baseIdx+2*k  , j) = mixture_.param(baseIdx+k).mean(j);
      params(baseIdx+2*k+1, j) = mixture_.param(baseIdx+k).sigma(j);
    }
  }
}

} // namespace STK

#endif /* STK_DIAGGAUSSIANBRIDGE_H */
