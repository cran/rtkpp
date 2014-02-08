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

/** @file STK_Exponential_bk.h
 *  @brief In this file we define the Exponential_pk_ajk_bk and Exponential_p_ajk_bk models.
 **/

#ifndef STK_EXPONENTIAL_BK_H
#define STK_EXPONENTIAL_BK_H


#include "STK_ExponentialBase.h"

#include "../../../STatistiK/include/STK_Law_Exponential.h"

#define MAXITER 400
#define TOL 1e-8

namespace STK
{
template<class Array>class Exponential_bk;

namespace Clust
{
/** @ingroup Clustering
 * Traits class for the Exponential_bk traits policy
 **/
template<class _Array>
struct MixtureTraits< Exponential_bk<_Array> >
{
  typedef _Array Array;
  typedef typename Array::Type Type;
  typedef Exponential_bk_Parameters  Parameters;
  typedef MixtureComponent<_Array, Parameters> Component;
  typedef Array2D<Real>            Param;
};

} // namespace Clust

/** @ingroup Clustering
 *  Exponential_bk is a mixture model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) = \sum_{k=1}^K p_k
 *     \prod_{j=1}^p\left(\frac{x_i^j}{b_{k}}\right)^{a_{jk}-1}
 *                   \frac{e^{-x_i^j/b_{k}}}{b_{k} \, \Exponential(a_{jk})},
 *      \quad x_i^j>0, \quad i=1,\ldots,n.
 * \f]
 **/
template<class Array>
class Exponential_bk : public ExponentialBase< Exponential_bk<Array> >
{
  public:
    typedef typename Clust::MixtureTraits< Exponential_bk<Array> >::Component Component;
    typedef typename Clust::MixtureTraits< Exponential_bk<Array> >::Parameters Parameters;
    typedef ExponentialBase< Exponential_bk<Array> > Base;

    using Base::p_tik;
    using Base::p_data;
    using Base::p_param;
    using Base::components;
    using Base::meanjk;
    using Base::variancejk;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline Exponential_bk( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline Exponential_bk( Exponential_bk const& model) : Base(model) {}
    /** destructor */
    inline ~Exponential_bk() {}
    /** initialize shape and scale parameters using weighted moment estimates.*/
    inline bool initializeStep() { return mStep();}
    /** Initialize randomly the parameters of the Exponential mixture. The shape
     *  will be selected randomly using an exponential of parameter mean^2/variance
     *  and the scale will be selected randomly using an exponential of parameter
     *  variance/mean.
     */
    void randomInit();
    /** Compute the mStep. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*this->nbVariable()+this->nbCluster();}
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void Exponential_bk<Array>::randomInit()
{
    // compute moments
    this->moments();
  for (int k= baseIdx; k < components().end(); ++k)
  {
    Real value =0;
    for (int j=p_data()->beginCols(); j < p_data()->endCols(); ++j)
    {
      Real mean = meanjk(j,k), variance = variancejk(j,k);
      p_param(k)->shape_[j] = Law::Exponential::rand(mean*mean/variance);
      value += variance/mean;
    }
    p_param(k)->scale_ = Law::Exponential::rand(value/this->nbVariable());
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Exponential_bk<Array>::randomInit done\n");
  this->writeParameters(stk_cout);
#endif
}

/* Compute the weighted mean and the common variance. */
template<class Array>
bool Exponential_bk<Array>::mStep()
{
  if (!this->moments()) { return false;}
  // start estimations of the ajk and bj
  Real qvalue = this->qValue();
  // enter iterative algorithm
  int iter;
  for(iter = 0; iter<MAXITER; ++iter)
  {
    // compute ajk
    for (int k= baseIdx; k < p_tik()->endCols(); ++k)
    {
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      {
        // moment estimate and oldest value
        Real x0 = meanjk(j,k)*meanjk(j,k)/variancejk(j,k);
        Real x1 = p_param(k)->shape_[j];
        if ((x0 <=0.) || !Arithmetic<Real>::isFinite(x0)) return false;
        // compute shape
        hidden::invPsi f(p_param(k)->meanLog_[j] - std::log(p_param(k)->scale_));
        Real a =  Algo::findZero(f, x0, x1, TOL);

        if (!Arithmetic<Real>::isFinite(a))
        {
           p_param(k)->shape_[j] = x0; // use moment estimate
#ifdef STK_MIXTURE_DEBUG
          stk_cout << _T("ML estimation failed in Exponential_bj::mStep()\n");
          stk_cout << "x0 =" << x0 << _T("\n";);
          stk_cout << "f(x0) =" << f(x0) << _T("\n";);
          stk_cout << "x1 =" << x1 << _T("\n";);
          stk_cout << "f(x1) =" << f(x1) << _T("\n";);
#endif
        }
        else { p_param(k)->shape_[j] = a;}
      }
      // compute bk
      p_param(k)->scale_ = p_param(k)->mean_.sum()/p_param(k)->shape_.sum();
    } // end ajk
    // check convergence
    Real value = this->qValue();
#ifdef STK_MIXTURE_VERBOSE
  if (value < qvalue)
  {
    stk_cout << _T("In Exponential_bk::mStep(): mStep diverge\n");
    stk_cout << _T("New value =") << value << _T(", qvalue =") << qvalue << _T("\n");
  }
#endif
    if ((value - qvalue) < TOL) break;
    qvalue = value;
  }
#ifdef STK_MIXTURE_VERBOSE
  if (iter == MAXITER)
  {
    stk_cout << _T("In Exponential_bk::mStep(): mStep did not converge\n");
    stk_cout << _T("qvalue =") << qvalue << _T("\n");
  }
#endif
  return true;
}

}  // namespace STK

#undef MAXITER
#undef TOL

#endif /* STK_EXPONENTIAL_BK_H */
