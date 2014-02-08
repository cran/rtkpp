/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

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

/** @file STK_Gamma_ajk_bjk.h
 *  @brief In this file we define the Gamma_pk_ajk_bjk and Gamma_p_ajk_bjk models.
 **/

#ifndef STK_GAMMA_AJK_BJK_H
#define STK_GAMMA_AJK_BJK_H

#include "STK_GammaBase.h"

#include "../../../STatistiK/include/STK_Law_Exponential.h"

namespace STK
{
template<class Array>class Gamma_ajk_bjk;

namespace Clust
{
/** @ingroup Clustering
 * Traits class for the Gamma_ajk_bjk traits policy
 **/
template<class _Array>
struct MixtureTraits< Gamma_ajk_bjk<_Array> >
{
  typedef _Array Array;
  typedef typename Array::Type Type;
  typedef MixtureComponent<_Array, Gamma_ajk_bjk_Parameters> Component;
  typedef Gamma_ajk_bjk_Parameters        Parameters;
  typedef Array2D<Real>        Param;
};

} // namespace Clust

/** @ingroup Clustering
 *  Gamma_ajk_bjk is a mixture model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) = \sum_{k=1}^K p_k
 *     \prod_{j=1}^p\left(\frac{x_i^j}{b_{jk}}\right)^{a_{jk}-1}
 *                   \frac{e^{-x_i^j/b_{jk}}}{b_{jk} \, \Gamma(a_{jk})},
 *      \quad x_i^j>0, \quad i=1,\ldots,n.
 * \f]
 **/
template<class Array>
class Gamma_ajk_bjk : public GammaBase< Gamma_ajk_bjk<Array> >
{
  public:
    typedef typename Clust::MixtureTraits< Gamma_ajk_bjk<Array> >::Component Component;
    typedef typename Clust::MixtureTraits< Gamma_ajk_bjk<Array> >::Parameters Parameters;
    typedef GammaBase< Gamma_ajk_bjk<Array> > Base;

    using Base::p_tik;
    using Base::p_data;
    using Base::p_param;
    using Base::components;
    using Base::meanjk;
    using Base::variancejk;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline Gamma_ajk_bjk( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline Gamma_ajk_bjk( Gamma_ajk_bjk const& model) : Base(model) {}
    /** destructor */
    inline ~Gamma_ajk_bjk() {}
    /** initialize shape and scale parameters using weighted moment estimates.*/
    inline bool initializeStep() { return mStep();}
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
    { return 2*this->nbCluster()*this->nbVariable();}
};

/* Initialize randomly the parameters of the gamma mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void Gamma_ajk_bjk<Array>::randomInit()
{
    // compute moments
    this->moments();
  for (int j=p_data()->beginCols(); j < p_data()->endCols(); ++j)
  {
    for (int k= baseIdx; k < components().end(); ++k)
    {
      Real mean = meanjk(j,k), variance = variancejk(j,k);
      p_param(k)->shape_[j] = Law::Exponential::rand(mean*mean/variance);
      p_param(k)->scale_[j] = Law::Exponential::rand(variance/mean);
    }
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Gamma_ajk_bjk<Array>::randomInit done\n");
  this->writeParameters(stk_cout);
#endif
}

/* Compute the weighted mean and the common variance. */
template<class Array>
bool Gamma_ajk_bjk<Array>::mStep()
{
  if (!this->moments()) { return false;}
  // estimate a and b
  for (int k= baseIdx; k < p_tik()->endCols(); ++k)
  {
    for (int j=p_data()->beginCols(); j < p_data()->endCols(); ++j)
    {
      // moment estimate and oldest value
      Real x0 = meanjk(j,k)*meanjk(j,k)/variancejk(j,k);
      Real x1 = p_param(k)->shape_[j];
      if ((x0 <=0.) || (Arithmetic<Real>::isNA(x0))) return false;

      // get shape
      hidden::invPsiMLog f(p_param(k)->meanLog_[j]-std::log(p_param(k)->mean_[j]));
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
      p_param(k)->shape_[j] = a;
      p_param(k)->scale_[j] = p_param(k)->mean_[j]/a;
    }
  }
  return true;
}


}  // namespace STK

#endif /* STK_GAMMA_AJK_BJK_H */
