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
 * created on: Oct 24, 2013
 * Author:   Serge Iovleff
 **/

/** @file STK_Gaussian_sk.h
 *  @brief In this file we define the Gaussian_sk model
 **/

#ifndef STK_GAUSSIAN_SK_H
#define STK_GAUSSIAN_SK_H

#include "STK_DiagGaussianBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class Gaussian_sk;

namespace Clust
{
/** @ingroup Clustering
 *  Traits class for the Gaussian_sk traits policy. */
template<class _Array>
struct MixtureTraits< Gaussian_sk<_Array> >
{
  typedef _Array Array;

  typedef ParametersHandler<Clust::Gaussian_sk_> ParamHandler;
};

} // namespace hidden

/** Specialization of the ParametersHandler struct for Gaussian_sk model */
template <>
struct ParametersHandler<Clust::Gaussian_sk_>: public DiagGaussianHandlerBase<  ParametersHandler<Clust::Gaussian_sk_> >
{
  /** RowVector and statistics of the means */
  MixtureParametersSet<PointX> mean_;
  /** standard deviation and statistics */
  MixtureParametersSet<Real> sigma_;
  /** @return the mean of the kth cluster and jth variable */
  inline Real const& meanImpl(int k, int j) const { return mean_[k][j];}
  /** @return the standard deviation of the kth cluster and jth variable */
  inline Real const& sigmaImpl(int k, int j) const { return sigma_[k];}
  /** copy operator */
  inline ParametersHandler& operator=( ParametersHandler const& other)
  { mean_ = other.mean_; sigma_ = other.sigma_; return *this;}
  /** copy operator using an array/expression storing the values */
  template<class Array>
  inline ParametersHandler& operator=( ExprBase<Array> const& param)
  {
    for (int k2= param.beginRows(), k= param.beginRows(); k2 < param.endRows(); k2+=2, k++)
    {
      for (int j= param.beginCols();  j< param.endCols(); ++j)
      { mean_[k][j]  = param(k2, j);}
      sigma_[k] = param.row(k2+1).mean();
    }
    return *this;
  }

  /** default constructor */
  ParametersHandler(int nbCluster): mean_(nbCluster), sigma_(nbCluster) {}
  /** copy constructor */
  ParametersHandler(ParametersHandler const& model): mean_(model.mean_), sigma_(model.sigma_) {}
  /** Initialize the parameters with an array/expression of value */
  template<class Array>
  inline ParametersHandler( int nbCluster, ExprBase<Array> const& param)
                          : mean_(nbCluster), sigma_(nbCluster)
  {
    mean_.resize(param.cols());
    for (int k2= param.beginRows(), k= param.beginRows(); k2 < param.endRows(); k2+=2, k++)
    {
      for (int j= param.beginCols();  j< param.endCols(); ++j)
      { mean_[k][j]  = param(k2, j);}
      sigma_[k] = param.row(k2+1).mean();
    }
  }

  /** Initialize the parameters of the model.
   *  This function initialize the parameters and the statistics.
   **/
  void resize(Range const& range)
  {
    mean_.resize(range);
    mean_.initialize(0.);
    sigma_.initialize(1.);
  }
  /** Store the intermediate results of the Mixture.
   *  @param iteration Provides the iteration number beginning after the burn-in period.
   **/
  inline void storeIntermediateResults(int iteration)
  { mean_.storeIntermediateResults(iteration); sigma_.storeIntermediateResults(iteration);}
  /** Release the stored results. This is usually used if the estimation
   *  process failed.
   **/
  inline void releaseIntermediateResults()
  { mean_.releaseIntermediateResults(); sigma_.releaseIntermediateResults();}
  /** set the parameters stored in stat_proba_ and release stat_proba_. */
  inline void setParameters() { mean_.setParameters(); sigma_.setParameters();}
};

/** @ingroup Clustering
 *  The diagonal Gaussian mixture model @c Gaussian_sk assumes an equal standard
 *  deviation in each cluster and has a density function of the form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d
 *    \frac{1}{\sqrt{2\pi}\sigma_{k}} \exp\left\{-\frac{(x^j-\mu^j_k)^2}{2(\sigma_{k})^2}\right\}.
 * \f]
 **/
template<class Array>
class Gaussian_sk : public DiagGaussianBase<Gaussian_sk<Array> >
{
  public:
    typedef DiagGaussianBase<Gaussian_sk<Array> > Base;
    using Base::p_tik;
    using Base::param_;
    using Base::p_nk;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline Gaussian_sk( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline Gaussian_sk( Gaussian_sk const& model) : Base(model) {}
    /** destructor */
    inline ~Gaussian_sk() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      {
        Real mean  = param_.mean_[k][j];
        Real sigma = param_.sigma_[k];
        sum += Law::Normal::lpdf(p_data()->elt(i,j), mean, sigma);
      }
      return sum;
    }
    /** Initialize randomly the parameters of the Gaussian mixture. The centers
     *  will be selected randomly among the data set and the standard-deviations
     *  will be set to 1.
     */
    void randomInit();
    /** Compute the weighted mean and the common standard deviation. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*this->nbVariable() + this->nbCluster();}
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void Gaussian_sk<Array>::randomInit()
{
  this->randomMean();
  // compute the standard deviation
  Real variance;
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    variance = sqrt( ( p_tik()->col(k).transpose()
                     *(*p_data() - (Const::Vector<Real>(p_data()->rows()) * param_.mean_[k])
                      ).square()
                     ).sum() / (p_data()->sizeCols()*p_nk()->elt(k))
                   );
    param_.sigma_[k] = ((variance<=0) || !Arithmetic<Real>::isFinite(variance))
                       ? 1.
                       : std::sqrt(variance/(this->nbSample()*this->nbVariable()));
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Gaussian_sk<Array>::randomInit() done\n");
#endif
}

/* Compute the weighted mean and the common standard deviation. */
template<class Array>
bool Gaussian_sk<Array>::mStep()
{
  // compute the means
  if (!this->updateMean()) return false;
  // compute the standard deviation
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    param_.sigma_[k]
    = sqrt( ( p_tik()->col(k).transpose()
             *(*p_data() - (Const::Vector<Real>(p_data()->rows()) * param_.mean_[k])
              ).square()
            ).sum()
           /(p_data()->sizeCols()*p_nk()->elt(k))
          );
//    if (param(k).sigma_ <= 0.) return false;
  }
  return true;
}

} // namespace STK

#endif /* STK_GAUSSIAN_SK_H */
