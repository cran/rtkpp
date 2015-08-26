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

/** @file STK_Gamma_aj_bjk.h
 *  @brief In this file we define the Gamma_pk_aj_bjk and Gamma_p_aj_bjk models.
 **/

#ifndef STK_GAMMA_AJ_BJK_H
#define STK_GAMMA_AJ_BJK_H

#include "STK_GammaBase.h"
#include <STatistiK/include/STK_Law_Exponential.h>

namespace STK
{
template<class Array>class Gamma_aj_bjk;

namespace Clust
{
/** @ingroup Clustering
 * Traits class for the Gamma_aj_bjk traits policy
 **/
template<class Array_>
struct MixtureTraits< Gamma_aj_bjk<Array_> >
{
  typedef Array_ Array;
  typedef ParametersHandler<Clust::Gamma_aj_bjk_> ParamHandler;
};

} // namespace Clust

/** Specialization of the ParametersHandler struct for Gamma_aj_bjk model */
template <>
struct ParametersHandler<Clust::Gamma_aj_bjk_>: public ParametersHandlerGammaBase<  ParametersHandler<Clust::Gamma_aj_bjk_> >
{
  typedef ParametersHandlerGammaBase Base;
  /** shape parameters and statistics */
  MixtureParameters<PointX> shape_;
  /** scale parameters and statistics */
  MixtureParametersSet<PointX> scale_;
  /** @return the shape of the kth cluster and jth variable */
  inline Real const& shapeImpl(int k, int j) const { return shape_()[j];}
  /** @return the scale of the kth cluster and jth variable */
  inline Real const& scaleImpl(int k, int j) const { return scale_[k][j];}
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
    int nbCluster = mean_().size();
    for (int j= param.beginCols();  j< param.endCols(); ++j)
    {
      shape_()[j] = 0.;
      for (int k2= param.beginRows(), k= param.beginRows(); k2 < param.endRows(); k2+=2, k++)
      {
        shape_()[j] += param(k2, j);
        scale_[k][j] = param(k2+1, j);
      }
      shape_()[j] /= nbCluster;
    }
    return *this;
  }

  /** default constructor */
  ParametersHandler( int nbCluster)
                   : Base(nbCluster), shape_(), scale_(nbCluster) {}
  /** copy constructor */
  ParametersHandler( ParametersHandler const& model)
                   : Base(model), shape_(model.shape_), scale_(model.scale_) {}
  /** Initialize the parameters with an array/expression of value */
  template<class Array>
  inline ParametersHandler( int nbCluster, ExprBase<Array> const& param)
                          : Base(nbCluster), shape_(), scale_(nbCluster)
  {
    Base::resize(param.cols());
    shape_.resize(param.cols());
    scale_.resize(param.cols());
    for (int j= param.beginCols();  j< param.endCols(); ++j)
    {
      shape_()[j] = 0.;
      for (int k2= param.beginRows(), k= param.beginRows(); k2 < param.endRows(); k2+=2, k++)
      {
        shape_()[j] += param(k2, j);
        scale_[k][j] = param(k2+1, j);
      }
      shape_()[j] /= nbCluster;
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
    shape_.resize(range);
    shape_.initialize(1.);
    scale_.resize(range);
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
 *  Gamma_aj_bjk is a mixture model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) = \sum_{k=1}^K p_k
 *     \prod_{j=1}^p\left(\frac{x_i^j}{b_{jk}}\right)^{a_{j}-1}
 *                   \frac{e^{-x_i^j/b_{jk}}}{b_{jk} \, \Gamma(a_{j})},
 *      \quad x_i^j>0, \quad i=1,\ldots,n.
 * \f]
 **/
template<class Array>
class Gamma_aj_bjk : public GammaBase< Gamma_aj_bjk<Array> >
{
  public:
    typedef GammaBase< Gamma_aj_bjk<Array> > Base;
    using Base::p_tik; using Base::param_;
    using Base::p_nk;
    using Base::p_data;
    using Base::meanjk;
    using Base::variancejk;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline Gamma_aj_bjk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline Gamma_aj_bjk( Gamma_aj_bjk const& model): Base(model) {}
    /** destructor */
    inline ~Gamma_aj_bjk() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      { sum += Law::Gamma::lpdf(p_data()->elt(i,j), param_.shape_()[j], param_.scale_[k][j]);}
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
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*this->nbVariable()+this->nbVariable();}
};

/* Initialize randomly the parameters of the gamma mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void Gamma_aj_bjk<Array>::randomInit()
{
  // compute moments
  this->moments();
  for (int j=p_data()->beginCols(); j < p_data()->endCols(); ++j)
  {
    // random scale for each cluster
    Real value = 0.0;
    for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
    {
      Real mean = meanjk(j,k), variance = variancejk(j,k);
      param_.scale_[k][j] = Law::Exponential::rand((variance/mean));
      value += p_nk()->elt(k) * (mean*mean/variance);
    }
    param_.shape_()[j] = STK::Law::Exponential::rand(value/(this->nbSample()));
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Gamma_aj_bjk<Array>::randomInit done\n");
#endif
}

/* Compute the weighted mean and the common variance. */
template<class Array>
bool Gamma_aj_bjk<Array>::mStep()
{
  if (!this->moments()) { return false;}
  // estimate a and b
  for (int j=p_data()->beginCols(); j < p_data()->endCols(); ++j)
  {
    Real y =0.0, x0 = 0.0;
    for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
    {
      Real mean = meanjk(j,k);
      y  += p_nk()->elt(k) * (param_.meanLog_[k][j]-std::log(mean));
      x0 += p_nk()->elt(k) * (mean*mean/variancejk(j,k));
    }
    // constant, moment estimate and oldest value
    y  /= this->nbSample();
    x0 /= this->nbSample();
    Real x1 = param_.shape_()[j];
    if ((x0 <=0.) || (isNA(x0))) return false;

    // get shape
    hidden::invPsiMLog f(y);
    Real a = Algo::findZero(f, x0, x1, 1e-08);
    if (!Arithmetic<Real>::isFinite(a))
    {
#ifdef STK_MIXTURE_DEBUG
      stk_cout << "ML estimation failed in Gamma_ajk_bjk::mStep()\n";
      stk_cout << "x0 =" << x0 << _T("\n";);
      stk_cout << "f(x0) =" << f(x0) << _T("\n";);
      stk_cout << "x1 =" << x1 << _T("\n";);
      stk_cout << "f(x1) =" << f(x1) << _T("\n";);
#endif
      a = x0; // use moment estimate
    }
    // set values
    param_.shape_()[j] = a;
    // update bjk
    for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
    { param_.scale_[k][j] = meanjk(j, k)/a;}
  }
  return true;
}


}  // namespace STK

#endif /* STK_GAMMA_AJ_BJK_H */
