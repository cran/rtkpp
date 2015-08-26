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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project: stkpp::Clustering
 * created on: 5 sept. 2013
 * Author:  iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Gamma_ak_bk.h
 *  @brief In this file we define the Gamma_pk_ak_bk and Gamma_p_ak_bk models.
 **/

#ifndef STK_GAMMA_AK_BK_H
#define STK_GAMMA_AK_BK_H

#include "STK_GammaBase.h"
#include <STatistiK/include/STK_Law_Exponential.h>


namespace STK
{
template<class Array>class Gamma_ak_bk;

namespace Clust
{
/** @ingroup Clustering
 * Traits class for the Gamma_ak_bk traits policy
 **/
template<class Array_>
struct MixtureTraits< Gamma_ak_bk<Array_> >
{
  typedef Array_ Array;
  typedef ParametersHandler<Clust::Gamma_ak_bk_> ParamHandler;
};

} // namespace Clust

/** Specialization of the ParametersHandler struct for Gamma_ak_bk model */
template <>
struct ParametersHandler<Clust::Gamma_ak_bk_>: public ParametersHandlerGammaBase<  ParametersHandler<Clust::Gamma_ak_bk_> >
{
  typedef ParametersHandlerGammaBase Base;
  /** shape parameters and statistics */
  MixtureParametersSet<Real> shape_;
  /** scale parameters and statistics */
  MixtureParametersSet<Real> scale_;
  /** @return the shape of the kth cluster and jth variable */
  inline Real const& shapeImpl(int k, int j) const { return shape_[k];}
  /** @return the scale of the kth cluster and jth variable */
  inline Real const& scaleImpl(int k, int j) const { return scale_[k];}
  /** copy operator */
  inline ParametersHandler& operator=( ParametersHandler const& other)
  { Base::operator =(other);
    shape_ = other.shape_; scale_ = other.scale_;
    return *this;
  }
  /** copy operator using an array/expression storing the values */
  template<class Array>
  inline ParametersHandler& operator=( ExprBase<Array> const& param)
  {
    for (int k2= param.beginRows(), k= param.beginRows(); k2 < param.endRows(); k2+=2, k++)
    {
      shape_[k] = param.row(k2).mean();
      scale_[k] = param.row(k2+1).mean();
    }
    return *this;
  }

  /** default constructor */
  ParametersHandler( int nbCluster)
                   : Base(nbCluster), shape_(nbCluster), scale_(nbCluster) {}
  /** copy constructor */
  ParametersHandler( ParametersHandler const& model)
                   : Base(model), shape_(model.shape_), scale_(model.scale_) {}
  /** Initialize the parameters with an array/expression of value */
  template<class Array>
  inline ParametersHandler( int nbCluster, ExprBase<Array> const& param)
                          : Base(nbCluster), shape_(nbCluster), scale_(nbCluster)
  {
    Base::resize(param.cols());
    for (int k2= param.beginRows(), k= param.beginRows(); k2 < param.endRows(); k2+=2, k++)
    {
      shape_[k] = param.row(k2).mean();
      scale_[k] = param.row(k2+1).mean();
    }
  }
  /** destructor */
  inline ~ParametersHandler() {}
  /** Initialize the parameters of the model.
   *  This function initialize the parameters and the statistics.
   **/
  void resize(Range const& range)
  {
    Base::resize(range);
    shape_.initialize(1.);
    scale_.initialize(1.);
  }
  /** Store the intermediate results of the Mixture.
   *  @param iteration Provides the iteration number beginning after the burn-in period.
   **/
  inline void storeIntermediateResults(int iteration)
  { shape_.storeIntermediateResults(iteration); scale_.storeIntermediateResults(iteration);}
  /** Release the stored results. This is usually used if the estimation
   *  process failed.
   **/
  inline void releaseIntermediateResults()
  { shape_.releaseIntermediateResults(); scale_.releaseIntermediateResults();}
  /** set the parameters stored in stat_proba_ and release stat_proba_. */
  inline void setParameters() { shape_.setParameters(); scale_.setParameters();}
};

/** @ingroup Clustering
 *  Gamma_ak_bk is a mixture model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) = \sum_{k=1}^K p_k
 *     \prod_{j=1}^p \left(\frac{x_i^j}{b_{k}}\right)^{a_{k}-1}
 *                   \frac{e^{-x_i^j/b_{k}}}{b_{k} \, \Gamma(a_{k})},
 *      \quad x_i^j>0, \quad i=1,\ldots,n.
 * \f]
 **/
template<class Array>
class Gamma_ak_bk : public GammaBase< Gamma_ak_bk<Array> >
{
  public:
    typedef GammaBase< Gamma_ak_bk<Array> > Base;
    using Base::p_tik; using Base::param_;
    using Base::p_data;
    using Base::meanjk;
    using Base::variancejk;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline Gamma_ak_bk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline Gamma_ak_bk( Gamma_ak_bk const& model): Base(model) {}
    /** destructor */
    inline ~Gamma_ak_bk() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      { sum += Law::Gamma::lpdf(p_data()->elt(i,j), param_.shape_[k], param_.scale_[k]);}
      return sum;
    }
    /** Initialize randomly the parameters of the Gamma mixture. The shape
     *  will be selected randomly using an exponential of parameter mean^2/variance
     *  and the scale will be selected randomly using an exponential of parameter
     *  variance/mean.
     */
    void randomInit();
    /** Compute the mStep. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const { return 2*this->nbCluster();}
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void Gamma_ak_bk<Array>::randomInit()
{
    // compute moments
    this->moments();
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    Real mean = this->meank(k), variance = this->variancek(k);
    param_.shape_[k]= Law::Exponential::rand((mean*mean/variance));
    param_.scale_[k] = Law::Exponential::rand((variance/mean));
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Gamma_ak_bk<Array>::randomInit done\n");
#endif
}

/* Compute the weighted mean and the common variance. */
template<class Array>
bool Gamma_ak_bk<Array>::mStep()
{
  if (!this->moments()) { return false;}
  // estimate a and b
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    // moment estimate and oldest value
    Real x0 = this->meank(k)*this->meank(k)/this->variancek(k);
    Real x1 =  param_.shape_[k];
    if ((x0 <=0.) || (isNA(x0))) return false;

    // get shape
    hidden::invPsiMLog f( (param_.meanLog_[k]-std::log( this->meank(k))).mean() );
    Real a = Algo::findZero(f, x0, x1, 1e-08);
    if (!Arithmetic<Real>::isFinite(a))
    {
#ifdef STK_MIXTURE_DEBUG
        stk_cout << "ML estimation failed in Gamma_ak_bjk::mStep()\n";
        stk_cout << "x0 =" << x0 << _T("\n";);
        stk_cout << "f(x0) =" << f(x0) << _T("\n";);
        stk_cout << "x1 =" << x1 << _T("\n";);
        stk_cout << "f(x1) =" << f(x1) << _T("\n";);
#endif
      a = x0; // use moment estimate
    }
    // set values
    param_.shape_[k]= a;
    param_.scale_[k] = this->meank(k)/a;
  }
  return true;
}


}  // namespace STK

#endif /* STK_GAMMA_AK_BK_H */
