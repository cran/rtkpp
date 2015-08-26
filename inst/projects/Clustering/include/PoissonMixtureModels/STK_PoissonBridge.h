/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015 Serge Iovleff

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

/** @file STK_PoissonBridge.h
 *  @brief In this file we define the bridge class between the Poisson mixture
 *  models and the composer.
 **/

#ifndef STK_POISSONBRIDGE_H
#define STK_POISSONBRIDGE_H

#include "STK_Poisson_ljk.h"
#include "STK_Poisson_ljlk.h"
#include "STK_Poisson_lk.h"

#include "../STK_MixtureData.h"
#include "../STK_IMixtureBridge.h"

namespace STK
{

// forward declaration
template<int Id, class Data> class PoissonBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Poisson_ljk model
 **/
template<class Data_>
struct MixtureBridgeTraits< PoissonBridge<Clust::Poisson_ljk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the Mixture model */
  typedef Poisson_ljk<Data> Mixture;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Poisson_ljk_> ParamHandler;
  /** Structure storing Parameters */
  typedef ArrayXX Parameters;
  enum
  {
    idMixtureClass_ = Clust::Poisson_
  };
};

/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Poisson_lk model
 **/
template<class Data_>
struct MixtureBridgeTraits< PoissonBridge<Clust::Poisson_lk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the Mixture model */
  typedef Poisson_lk<Data> Mixture;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Poisson_lk_> ParamHandler;
  /** Structure storing Parameters */
  typedef ArrayXX Parameters;
  enum
  {
    idMixtureClass_ = Clust::Poisson_
  };
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Poisson_ljlk model
 **/
template<class Data_>
struct MixtureBridgeTraits< PoissonBridge< Clust::Poisson_ljlk_, Data_> >
{
  typedef Data_ Data;
  /** Type of the mixture model */
  typedef Poisson_ljlk<Data> Mixture;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Poisson_ljlk_> ParamHandler;
  /** Structure storing Parameters */
  typedef ArrayXX Parameters;
  enum
  {
    idMixtureClass_ = Clust::Poisson_
  };
};

} // namespace hidden

} // namespace STK

namespace STK
{
/** @ingroup Clustering
 *  @brief Templated implementation of the IMixture interface allowing
 *  to bridge a STK++ Poisson mixture with the composer.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureModel class.
 */
template<int Id, class Data>
class PoissonBridge: public IMixtureBridge< PoissonBridge<Id,Data> >
{
  public:
    // Base class
    typedef IMixtureBridge< PoissonBridge<Id,Data> > Base;
    typedef typename hidden::MixtureBridgeTraits< PoissonBridge<Id,Data> >::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits< PoissonBridge<Id,Data> >::ParamHandler ParamHandler;
    typedef typename hidden::MixtureBridgeTraits< PoissonBridge<Id,Data> >::Parameters Parameters;
    typedef typename Data::Type Type;
    // class of mixture
    enum
    {
      idMixtureClass_ = Clust::Poisson_
    };
    typedef std::vector<std::pair<int,int> >::const_iterator ConstIterator;
    using Base::mixture_;
    using Base::p_data_;

    /** default constructor. Remove the missing values from the data set and
     *  initialize the mixture by setting the data set.
     *  @param p_data pointer on the MixtureData that will be used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    PoissonBridge( MixtureData<Data>* p_data, std::string const& idData, int nbCluster)
                 : Base(p_data, idData, nbCluster)
    { removeMissing(); initializeMixture();}
    /** copy constructor */
    PoissonBridge( PoissonBridge const& bridge): Base(bridge)
    { initializeMixture();}
    /** destructor */
    virtual ~PoissonBridge() {}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual PoissonBridge* clone() const { return new PoissonBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual PoissonBridge* create() const
    {
      PoissonBridge* p_bridge = new PoissonBridge( mixture_, this->idName(), this->nbCluster());
      p_bridge->p_data_ = p_data_;
      // Bug Fix: set the correct data set
      p_bridge->mixture_.setData(p_bridge->p_data_->dataij());
      return p_bridge;
    }
    /** This function is used in order to get the current values of the
     *  parameters.
     *  @param params the array with the parameters of the mixture.
     */
    void getParameters(Parameters& params) const;
    /** This function can be used to write summary of parameters to the output stream.
     *  @param out Stream where you want to write the summary of parameters.
     */
    virtual void writeParameters(std::ostream& out) const;

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
    PoissonBridge( Mixture const& mixture, std::string const& idData, int nbCluster)
                 : Base(mixture, idData, nbCluster)
    {}
    /** @return a safe value for the jth variable
     *  @param  dataij the matrix of the data set
     *  @param j index of the column with the safe value needed */
    static Type safeValue(Data const& dataij, int j)
    {

      int lmin = dataij.col(j).safe().minElt(), lmax = dataij.col(j).safe().maxElt();
      if (lmax -lmin > 10)
      { return Real(dataij.col(j).safe().sum())/dataij.sizeRows();}
      Array2DVector<int> count(Range(lmin, lmax, 0), 0);
      for (int i= dataij.beginRows(); i < dataij.endRows(); ++i)
      {
        if (!Arithmetic<int>::isNA(dataij(i,j)))
          count[dataij(i,j)]++;
      }
      int l; count.maxElt(l);
      return l;
    }

};

// implementation
template<int Id, class Data>
void PoissonBridge<Id, Data>::removeMissing()
{
  Type value = Type();
  int j, old_j = Arithmetic<int>::NA();
  for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
  {
    j = it->second; // get column
    if (j != old_j)
    {
      old_j = j;
      value =  safeValue(p_data_->dataij_, j);
    }
    p_data_->dataij_(it->first, j) = value;
  }
}

template<int Id, class Data>
void PoissonBridge<Id, Data>::getParameters(Parameters& params) const
{
  params.resize(this->nbCluster(), mixture_.p_data()->cols());
  for (int k= params.beginRows(); k < params.endRows(); ++k)
  {
    for (int j= mixture_.p_data()->beginCols();  j < mixture_.p_data()->endCols(); ++j)
    { params(k, j) = mixture_.lambda(k,j);}
  }
}

template<int Id, class Data>
void PoissonBridge<Id, Data>::writeParameters(ostream& os) const
{
  ArrayXX params;
  getParameters(params);
  for (int k= params.beginRows(); k < params.endRows(); ++k)
  {
    os << _T("---> Component ") << k << _T("\n");
    os << _T("lambda = ") << params.row(k);
  }
}

} // namespace STK

#endif /* STK_POISSONBRIDGE_H */
