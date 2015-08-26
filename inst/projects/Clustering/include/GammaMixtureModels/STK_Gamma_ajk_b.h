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
 * created on: 29 août 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Gamma_ajk_b.h
 *  @brief In this file we define the Gamma_pk_ajk_b and Gamma_p_ajk_b mixture models.
 **/


#ifndef STK_GAMMA_AJK_B_H
#define STK_GAMMA_AJK_B_H

#include "STK_GammaBase.h"
#include <STatistiK/include/STK_Law_Exponential.h>

#define MAXITER 400
#define TOL 1e-8


namespace STK
{
template<class Array>class Gamma_ajk_b;

namespace Clust
{
/** @ingroup Clustering
 *  Traits class for the Gamma_ajk_b traits policy. */
template<class Array_>
struct MixtureTraits< Gamma_ajk_b<Array_> >
{
  typedef Array_ Array;
  typedef ParametersHandler<Clust::Gamma_ajk_b_> ParamHandler;
};

} // namespace Clust

/** Specialization of the ParametersHandler struct for Gamma_ajk_b model */
template <>
struct ParametersHandler<Clust::Gamma_ajk_b_>: public ParametersHandlerGammaBase<  ParametersHandler<Clust::Gamma_ajk_b_> >
{
  typedef ParametersHandlerGammaBase Base;
  /** shape parameters and statistics */
  MixtureParametersSet<PointX> shape_;
  /** scale parameters and statistics */
  MixtureParameters<Real> scale_;
  /** @return the shape of the kth cluster and jth variable */
  inline Real const& shapeImpl(int k, int j) const { return shape_[k][j];}
  /** @return the scale of the kth cluster and jth variable */
  inline Real const& scaleImpl(int k, int j) const { return scale_();}
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
    scale_() = 0.;
    for (int k2= param.beginRows(), k= param.beginRows(); k2 < param.endRows(); k2+=2, k++)
    {
      for (int j= param.beginCols();  j< param.endCols(); ++j)
      {
        shape_[k][j] = param(k2, j);
        scale_()    += param(k2+1, j);
      }
    }
    scale_() /= (nbCluster*param.sizeCols());
    return *this;
  }

  /** default constructor */
  ParametersHandler( int nbCluster)
                   : Base(nbCluster), shape_(nbCluster), scale_() {}
  /** copy constructor */
  ParametersHandler( ParametersHandler const& model)
                   : Base(model), shape_(model.shape_), scale_(model.scale_) {}
  /** Initialize the parameters with an array/expression of value */
  template<class Array>
  inline ParametersHandler( int nbCluster, ExprBase<Array> const& param)
                          : Base(nbCluster), shape_(nbCluster), scale_()
  {
    Base::resize(param.cols());
    shape_.resize(param.cols());
    scale_() = 0.;
    for (int k2= param.beginRows(), k= param.beginRows(); k2 < param.endRows(); k2+=2, k++)
    {
      for (int j= param.beginCols();  j< param.endCols(); ++j)
      {
        shape_[k][j] = param(k2, j);
        scale_()    += param(k2+1, j);
      }
    }
    scale_() /= (nbCluster*param.sizeCols());
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
 *  Gamma_ajk_b is a mixture model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) = \sum_{k=1}^K p_k
 *     \prod_{j=1}^p \left(\frac{x_i^j}{b}\right)^{a_{jk}-1}
 *                   \frac{e^{-x_i^j/b}}{b \, \Gamma(a_{jk})},
 *      \quad x_i^j>0, \quad i=1,\ldots,n.
 * \f]
 **/
template<class Array>
class Gamma_ajk_b : public GammaBase<Gamma_ajk_b<Array> >
{
  public:
    typedef GammaBase<Gamma_ajk_b<Array> > Base;
    using Base::p_tik; using Base::param_;
    using Base::p_nk;
    using Base::p_data;
    using Base::meanjk;
    using Base::variancejk;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline Gamma_ajk_b( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline Gamma_ajk_b( Gamma_ajk_b const& model): Base(model) {}
    /** destructor */
    inline ~Gamma_ajk_b() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      { sum += Law::Gamma::lpdf(p_data()->elt(i,j), param_.shape_[k][j], param_.scale_());}
      return sum;
    }
    /** Initialize randomly the parameters of the Gamma mixture. */
    void randomInit();
    /** Compute the weighted mean and the common variance. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*this->nbVariable() + 1;}
};

/* Initialize randomly the parameters of the gamma mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void Gamma_ajk_b<Array>::randomInit()
{
    // compute moments
    this->moments();
  Real value =0.;
  for (int j=p_data()->beginCols(); j < p_data()->endCols(); ++j)
  {
    // random scale for each cluster
    for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
    {
      Real mean = meanjk(j,k), variance = variancejk(j,k);
      param_.shape_[k][j] = Law::Exponential::rand((mean*mean/variance));
      value += p_nk()->elt(k) * variance/mean;
    }
  }
  param_.scale_() = Law::Exponential::rand(value/(this->nbVariable()*this->nbSample()));
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Gamma_ajk_b<Array>::randomInit done\n");
#endif
}

/* Compute the weighted mean and the common variance. */
template<class Array>
bool Gamma_ajk_b<Array>::mStep()
{
  if (!this->moments()) { return false;}
  // start estimations of the ajk and bj
  Real qvalue = this->qValue();
  int iter;
  for(iter=0; iter<MAXITER; ++iter)
  {
    for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
    {
      // compute ajk
      for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
      {
        // moment estimate and oldest value
        Real x0 = meanjk(j,k)*meanjk(j,k)/variancejk(j,k);
        Real x1 = param_.shape_[k][j];
        if ((x0 <=0.) || !Arithmetic<Real>::isFinite(x0)) return false;
        // compute shape
        hidden::invPsi f(param_.meanLog_[k][j] - std::log(param_.scale_()));
        Real a =  Algo::findZero(f, x0, x1, TOL);

        if (!Arithmetic<Real>::isFinite(a))
        {
          param_.shape_[k][j] = x0; // use moment estimate
#ifdef STK_MIXTURE_DEBUG
          stk_cout << _T("ML estimation failed in Gamma_ajk_bj::mStep()\n");
          stk_cout << "x0 =" << x0 << _T("\n";);
          stk_cout << "f(x0) =" << f(x0) << _T("\n";);
          stk_cout << "x1 =" << x1 << _T("\n";);
          stk_cout << "f(x1) =" << f(x1) << _T("\n";);
#endif
        }
        else { param_.shape_[k][j] = a;}
      }
    }
    Real num=0., den = 0.;
    for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
    {
      num += param_.mean_[k].sum()  * p_nk()->elt(k);
      den +=  param_.shape_[k].sum() * p_nk()->elt(k);
    }
    // compute b
    Real b = num/den;
    // divergence
    if (!Arithmetic<Real>::isFinite(b)) { return false;}
    param_.scale_() = b;
    // check convergence
    Real value = this->qValue();
#ifdef STK_MIXTURE_DEBUG
    if (value < qvalue)
    {
      stk_cout << _T("In Gamma_ajk_b::mStep(): mStep diverge\n");
      stk_cout << _T("New value =") << value << _T(", qvalue =") << qvalue << _T("\n");
    }
#endif
    if ((value - qvalue) < TOL) break;
    qvalue = value;
  }
#ifdef STK_MIXTURE_DEBUG
  if (iter == MAXITER)
  {
    stk_cout << _T("In Gamma_ajk_b::mStep(): mStep did not converge\n");
    stk_cout << _T("qvalue =") << qvalue << _T("\n");
  }
#endif
  return true;
}

} // namespace STK

#undef MAXITER
#undef TOL

#endif /* STK_GAMMA_AJK_B_H */
