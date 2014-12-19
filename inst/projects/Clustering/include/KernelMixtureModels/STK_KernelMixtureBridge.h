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

/** @file STK_KernelMixtureBridge.h
 *  @brief In this file we define the bridge class between the mixtures and the composer.
 **/

#ifndef STK_KERNELMIXTUREBRIDGE_H
#define STK_KERNELMIXTUREBRIDGE_H

namespace STK
{
// forward declaration
template<int Id, class Data> class KernelMixtureBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the KernelGaussian model
 **/
template<class Data>
struct MixtureBridgeTraits< KernelMixtureBridge< 1, Data> >
{
  /** Type of the mixture model */
  typedef KernelGaussian<Data> Mixture;
};

} // namespace hidden

} // namespace STK

#include <Clustering/include/STK_IMixtureBridge.h>

namespace STK
{
/** @ingroup Clustering
 *  @brief Templated implementation of the IMixture interface allowing
 *  to bridge a STK++ kernel mixture with the composer.
 *
 *  This class inherit from the interface IMixtureBridge.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureModel class.
 */
template<int Id, class Data>
class KernelMixtureBridge: public IMixtureBridge< KernelMixtureBridge<Id,Data> >
{
  public:
    typedef IMixtureBridge< KernelMixtureBridge<Id,Data> > Base;
    typedef typename hidden::MixtureBridgeTraits< KernelMixtureBridge<Id,Data> >::Mixture Mixture;
    typedef typename Clust::MixtureTraits<Mixture>::Param Param;
    typedef typename Data::Type Type;

    // class of mixture
    enum
    {
      idMixtureClass_ = 1
    };
    // parameters type to get
    using Base::mixture_;

    /** default constructor. Remove the missing values from the data set and
     *  initialize the mixture by setting the data set.
     *  @param p_data pointer on the MixtureData that will be used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    KernelMixtureBridge( Data* p_gram, std::string const& idData, int nbCluster)
                       : Base( idData, nbCluster)
                       , p_kii_(p_data)
    { initializeMixture();}
    /** copy constructor */
    KernelMixtureBridge( KernelMixtureBridge const& bridge)
                       : Base(bridge)
                       , p_kii_(bridge.p_kii_)
    {  mixture_.setData(p_data_);}
    /** destructor */
    virtual ~KernelMixtureBridge() {}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual KernelMixtureBridge* clone() const { return new KernelMixtureBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual KernelMixtureBridge* create() const
    {
      KernelMixtureBridge* p_bridge = new KernelMixtureBridge( mixture_, idName(), nbCluster());
      p_bridge->p_kii_ = p_kii_;
      p_bridge->mixture_.setData(p_bridge->p_kii_);
      return p_bridge;
    }
    /** This function is equivalent to MStep and must be defined to update
     *  parameters. In a Kernel mixture model, the MStep is defined by
     *  - an update of the distance from the center of the class
     *  - an usual update of the parameters
     */
    virtual void paramUpdateStep()
    {
      if (!mixture_.mStep()) throw Clust::mStepFail_;
    }
    /** This function can be used in order to the values of the parameters
     *  in an Array2D.
     *  @param param the array with the parameters of the mixture.
     */
    void getParameters(Param& param) const { mixture_.getParameters(param);}

  private:
   /** Compute the */
   void compute_yik()
   {
     // matrix of size (n,K) with values \sum_{j=1}^n k(x_i,x_j) t_{jk}/t_{.k}
     CVectorX wik =  (*p_kii_) * (tik())
                    / Const::Vector<Real>(this->nbSample()) * this->nk();
     // compute \sum_{i=1}^n \sum_{j=1}^n t_{ik} k(x_i,x_j) t_{jk}/(t_{.k})^2
     CSquareX akk = (tik()/nk()).transpose() * (*p_kii_) * (tik()/nk());

     for (int k= yik_.beginCols(); yik_.endCols(); ++k)
     { yik_.col(k) = p_kii_->diagonal() - 2 * wik.col(k) + wik.col(k).dot(p_tik()->col(k));}
   }
    /** This function will be used in order to initialize the mixture model
     *  using informations stored by the MixtureData. For example the missing
     *  values in the case of a MixtureData instance.
     **/
    void initializeMixture(){ mixture_.setData(p_kii_);}
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    KernelMixtureBridge( Mixture const& mixture, std::string const& idData, int nbCluster)
                       : Base( idData, nbCluster)
                       , p_kii_(0)
    {}
    /** reference on the gram matrix*/
    Data const* p_kii_;
    CArrayXX yik_;
};

} // namespace STK

#endif /* STK_KERNELMIXTUREBRIDGE_H */
