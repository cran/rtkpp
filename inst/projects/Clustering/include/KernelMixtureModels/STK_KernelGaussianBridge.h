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
 * created on: 23 oct. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_KernelGaussianBridge.h
 *  @brief In this file we define the bridge class between the mixtures and the composer.
 **/

#ifndef STK_KERNELMIXTUREBRIDGE_H
#define STK_KERNELMIXTUREBRIDGE_H

#include "../STK_MixtureData.h"
#include "../STK_IMixtureBridge.h"
#include "STK_KernelGaussian.h"


namespace STK
{
// forward declaration
template<int Id, class Data> class KernelGaussianBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the KernelGaussian_sk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< KernelGaussianBridge< Clust::KernelGaussian_sk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef KernelGaussian_sk Mixture;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::KernelGaussian_sk_> ParamHandler;
  /** Structure storing Parameters */
  typedef ArrayXX Parameters;
  // class of mixture
  enum
  {
    idMixtureClass_ = Clust::Kernel_
  };
};

/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the KernelGaussian_s_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< KernelGaussianBridge< Clust::KernelGaussian_s_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef KernelGaussian_s Mixture;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::KernelGaussian_s_> ParamHandler;
  /** Structure storing Parameters */
  typedef ArrayXX Parameters;
  // class of mixture
  enum
  {
    idMixtureClass_ = Clust::Kernel_
  };
};

} // namespace hidden

} // namespace STK

namespace STK
{
/** @ingroup Clustering
 *  @brief Templated implementation of the IMixture interface allowing
 *  to bridge a STK++ kernel mixture with the composer.
 *
 *  This class inherit from the interface IMixtureBridge.
 *
 *  @note The p_data_ DataManager is wrapping the Gram matrix.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureModel class.
 * @tparam Data is any container storing the data.
 */
template<int Id, class Data>
class KernelGaussianBridge: public IMixtureBridge< KernelGaussianBridge<Id,Data> >
{
  public:
    typedef IMixtureBridge< KernelGaussianBridge<Id,Data> > Base;
    typedef typename hidden::MixtureBridgeTraits< KernelGaussianBridge<Id,Data> >::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits< KernelGaussianBridge<Id,Data> >::ParamHandler ParamHandler;
    typedef typename hidden::MixtureBridgeTraits< KernelGaussianBridge<Id,Data> >::Parameters Parameters;
    typedef typename Data::Type Type;

    // parameters type to get
    using Base::mixture_;
    using Base::p_data_;
    using Base::p_tik;
    using Base::tik;
    using Base::nk;

    /** default constructor. Remove the missing values from the data set and
     *  initialize the mixture by setting the data set.
     *  @param p_data pointer on the MixtureData that will be used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    KernelGaussianBridge( MixtureData<Data>* p_data, std::string const& idData, int nbCluster)
                       : Base( p_data, idData, nbCluster)
                       , p_kii_(&(p_data->dataij_))
                       , dik_(p_kii_->rows(), nbCluster)
    { initializeMixture();}
    /** copy constructor */
    KernelGaussianBridge( KernelGaussianBridge const& bridge)
                       : Base(bridge)
                       , p_kii_(&(p_data_->dataij_))
                       , dik_(bridge.dik_)
    { initializeMixture();}
    /** destructor */
    virtual ~KernelGaussianBridge() {}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual KernelGaussianBridge* clone() const { return new KernelGaussianBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual KernelGaussianBridge* create() const
    {
      KernelGaussianBridge* p_bridge = new KernelGaussianBridge( mixture_, this->idName(), this->nbCluster());
      p_bridge->p_data_ = p_data_;
      p_bridge->p_kii_  = &(p_bridge->p_data_->dataij_);
      p_bridge->initializeMixture();
      return p_bridge;
    }
    /** set the dimension of the kernel mixture model */
    inline void setDim(Real const& dim) { mixture_.setDim(dim);}
    /** This function is equivalent to MStep and must be defined to update
     *  parameters. In a Kernel mixture model, the MStep is defined by
     *  - an update of the distance from the center of the class
     *  - an usual update of the parameters
     */
    virtual void paramUpdateStep()
    {
      compute_dik(); // this update the values needed by the mixture_
      if (!mixture_.mStep()) throw Clust::mStepFail_;
    }
    /** This function can be used in order to the values of the parameters
     *  in an Array2D.
     *  @param param the array with the parameters of the mixture.
     */
    void getParameters(Parameters& param) const;
    /** This function can be used to write summary of parameters to the output stream.
     *  @param out Stream where you want to write the summary of parameters.
     */
    virtual void writeParameters(std::ostream& out) const;

  private:
    /** Compute the intermediate results dik */
    void compute_dik()
    {
      // matrix of size (n,K) with values \sum_{j=1}^n k(x_i,x_j) t_{jk}/t_{.k}
      CArrayXX wik = (*p_kii_ * tik()) / (Const::Vector<Real>(p_kii_->rows()) * nk());
      for (int k= dik_.beginCols(); k<dik_.endCols(); ++k)
      {
        Real scal =tik().col(k).dot(wik.col(k))/nk()[k];
        for (int i= dik_.beginRows(); i<dik_.endRows(); ++i)
        { dik_(i,k) = p_kii_->elt(i,i) - 2. * wik(i,k) + scal  ;}
      }
    }
    /** This function will be used in order to initialize the mixture model
     *  using informations stored by the MixtureData. For example the missing
     *  values in the case of a MixtureData instance.
     *
     *  In the kernel bridge the mixture use the intermediary array @c dik_ computed
     *  by the bridge. The initialization step consist in resizing the array
     *  and to set the pointer to @c mixture_.
     **/
    void initializeMixture()
    {
      dik_.resize(p_kii_->rows(), this->nbCluster());
      mixture_.setData(dik_);
    }
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    KernelGaussianBridge( Mixture const& mixture, std::string const& idData, int nbCluster)
                        : Base( mixture, idData, nbCluster)
                        , p_kii_(0)
                        , dik_()
    {}
    /** constant reference on the gram matrix*/
    Data const* p_kii_;
    /** Array of the intermediate results dik */
    ArrayXX dik_;
};
template<int Id, class Data>
void KernelGaussianBridge<Id, Data>::getParameters(Parameters& param) const
{
  param.resize(this->nbCluster(), 2);
  for (int k= param.beginRows(); k < param.endRows(); ++k)
  {
    param(k, baseIdx  ) = std::sqrt(mixture_.paramHandler().sigma2(k));
    param(k, baseIdx+1) = mixture_.paramHandler().dim(k);
  }
}

/* Write the parameters on the output stream os */
template<int Id, class Data>
void KernelGaussianBridge<Id, Data>::writeParameters(ostream& os) const
{
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    os << _T("---> Component ") << k << _T("\n");
    os << _T("sigma = ") << std::sqrt(mixture_.paramHandler().sigma2(k)) << _T("\n");
    os << _T("dim = ")   << mixture_.paramHandler().dim(k) << _T("\n");
  }
}

} // namespace STK

#endif /* STK_KERNELMIXTUREBRIDGE_H */
