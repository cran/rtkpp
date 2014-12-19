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
 * created on: Oct 24, 2013
 * Author:   Serge Iovleff
 **/

/** @file STK_Gaussian_s.h
 *  @brief In this file we implement the Gaussian_s class
 **/

#ifndef STK_GAUSSIAN_S_H
#define STK_GAUSSIAN_S_H

#include <Clustering/include/DiagGaussianMixtureModels/STK_DiagGaussianBase.h>

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class Gaussian_s;

namespace Clust
{
/** @ingroup Clust
 *  Traits class for the Gaussian_s traits policy. */
template<class _Array>
struct MixtureTraits< Gaussian_s<_Array> >
{
  typedef _Array Array;
  typedef typename Array::Type Type;
  typedef Gaussian_s_Parameters Parameters;
  typedef Array2D<Real>        Param;
};

} // namespace Clust

/** @ingroup Clustering
 *  The diagonal Gaussian_s mixture model have a density function of the form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d
 *    \frac{1}{\sqrt{2\pi}\sigma} \exp\left\{-\frac{(x^j-\mu^j_k)^2}{2\sigma^2}\right\}.
 * \f]
 **/
template<class Array>
class Gaussian_s : public DiagGaussianBase<Gaussian_s<Array> >
{
  public:
    typedef DiagGaussianBase<Gaussian_s<Array> > Base;
    typedef typename Clust::MixtureTraits< Gaussian_s<Array> >::Parameters Parameters;

    using Base::p_tik;
    using Base::components;
    using Base::p_data;
    using Base::param;


    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline Gaussian_s( int nbCluster) : Base(nbCluster), sigma_(1) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline Gaussian_s( Gaussian_s const& model)
                     : Base(model), sigma_(model.sigma_)
                     , stat_sigma_(model.stat_sigma_)
    {}
    /** destructor */
    inline ~Gaussian_s() {}
    /** Initialize the component of the model.
     *  In this interface, the @c initializeModelImpl()
     *  initialize the shared parameter @c sigma_ for all the components .
     **/
    void initializeModelImpl()
    {
      sigma_ = 1.0;
      for (int k= baseIdx; k < components().end(); ++k)
      { param(k).p_sigma_ = &sigma_;}
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResultsImpl(int iteration)
    { stat_sigma_.update(sigma_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResultsImpl()
    { stat_sigma_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParametersImpl()
    {
      sigma_ = stat_sigma_.param_;
      stat_sigma_.release();
    }
    /** Initialize randomly the parameters of the Gaussian mixture. The centers
     *  will be selected randomly among the data set and the standard-deviation
     *  will be set to 1.
     */
    void randomInit();
    /** Compute the weighted mean and the common standard deviation. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*this->nbVariable()+1;}

  protected:
    /** Common standard deviation */
    Real sigma_;
    /** statistics */
    MixtureStatReal stat_sigma_;
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void Gaussian_s<Array>::randomInit()
{
  this->randomMean();
  // compute the standard deviation
  Real variance = 0.0;
  for (int k= baseIdx; k < components().end(); ++k)
  {
    variance += ( p_tik()->col(k).transpose()
                 * (*p_data() - (Const::Vector<Real>(p_data()->rows()) * param(k).mean_)
                   ).square()
                ).sum();
  }
  sigma_ = ((variance<=0) || !Arithmetic<Real>::isFinite(variance))
           ? 1.
           : std::sqrt(variance/(this->nbSample()*this->nbVariable()));
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Gaussian_s<Array>::randomInit() done\n");
#endif
}

/* Compute the weighted mean and the common standard deviation. */
template<class Array>
bool Gaussian_s<Array>::mStep()
{
  // compute the means
  if (!this->updateMean()) return false;
  // compute the standard deviation
  Real variance = 0.0;
  for (int k= baseIdx; k < components().end(); ++k)
  {
    variance += ( p_tik()->col(k).transpose()
                 * (*p_data() - (Const::Vector<Real>(p_data()->rows()) * param(k).mean_)
                   ).square()
                ).sum();
  }
  if ((variance<=0) || !Arithmetic<Real>::isFinite(variance)) return false;
  sigma_ = std::sqrt(variance/(this->nbSample()*this->nbVariable()));
  return true;
}

} // namespace STK

#endif /* STK_GAUSSIAN_SJK_H */
