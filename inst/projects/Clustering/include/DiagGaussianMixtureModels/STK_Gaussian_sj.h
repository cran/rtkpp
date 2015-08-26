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

/** @file STK_Gaussian_sj.h
 *  @brief In this file we define and implement the Gaussian_sj class
 **/

#ifndef STK_GAUSSIAN_SJ_H
#define STK_GAUSSIAN_SJ_H

#include "STK_DiagGaussianBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class Gaussian_sj;

namespace Clust
{
/** @ingroup hidden
 *  Traits class for the Gaussian_s traits policy. */
template<class _Array>
struct MixtureTraits< Gaussian_sj<_Array> >
{
  typedef _Array Array;

  typedef ParametersHandler<Clust::Gaussian_sj_> ParamHandler;
};

} // namespace hidden

/** Specialization of the ParametersHandler struct for Gaussian_sj model */
template <>
struct ParametersHandler<Clust::Gaussian_sj_>: public DiagGaussianHandlerBase<  ParametersHandler<Clust::Gaussian_sj_> >
{
    /** RowVector and statistics of the means */
    MixtureParametersSet<PointX> mean_;
    /** standard deviation and statistics */
    MixtureParameters<PointX> sigma_;
    /** @return the mean of the kth cluster and jth variable */
    inline Real const& meanImpl(int k, int j) const { return mean_[k][j];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& sigmaImpl(int k, int j) const { return sigma_()[j];}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    { mean_ = other.mean_; sigma_ = other.sigma_; return *this;}
    /** copy operator using an array/expression storing the values */
    template<class Array>
    inline ParametersHandler& operator=( ExprBase<Array> const& param)
    {
      int nbCluster = mean_().size();
      sigma_.initialize(0.);
      for (int j= param.beginCols();  j< param.endCols(); ++j)
      {
        for (int k1= param.beginRows(), k2= param.beginRows(); k2 < param.endRows(); k2+=2, k1++)
        {
          mean_[k1][j]  = param(k2  , j);
          sigma_()[j] += param(k2+1, j);
        }
        sigma_()[j] /= nbCluster;
      }
      return *this;
    }

    /** default constructor */
    ParametersHandler(int nbCluster): mean_(nbCluster), sigma_() {}
    /** copy constructor */
    ParametersHandler(ParametersHandler const& model): mean_(model.mean_), sigma_(model.sigma_) {}
    /** Initialize the parameters with an array/expression of value */
    template<class Array>
    inline ParametersHandler( int nbCluster, ExprBase<Array> const& param)
                            : mean_(nbCluster), sigma_()
    {
      mean_.resize(param.cols());
      sigma_.resize(param.cols());
      sigma_.initialize(0.);
      for (int j= param.beginCols();  j< param.endCols(); ++j)
      {
        for (int k2= param.beginRows(), k= param.beginRows(); k2 < param.endRows(); k2+=2, k++)
        {
          mean_[k][j]  = param(k2  , j);
          sigma_()[j] += param(k2+1, j);
        }
        sigma_()[j] /= nbCluster;
      }
    }

   /** Initialize the parameters of the model.
     *  This function initialize the parameters and the statistics.
     **/
    void resize(Range const& range)
    {
      mean_.resize(range);
      mean_.initialize(0.);
      sigma_.resize(range);
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
 *  The diagonal Gaussian mixture model Gaussian_sj have a density function of the
 *  form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d
 *    \frac{1}{\sqrt{2\pi}\sigma_j} \exp\left\{-\frac{(x^j-\mu^j_k)^2}{2\sigma_j^2}\right\}.
 * \f]
 **/
template<class Array>
class Gaussian_sj : public DiagGaussianBase<Gaussian_sj<Array> >
{
  public:
    typedef DiagGaussianBase<Gaussian_sj<Array> > Base;
    using Base::p_tik;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline Gaussian_sj( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline Gaussian_sj( Gaussian_sj const& model)
                      : Base(model)
    {}
    /** destructor */
    inline ~Gaussian_sj() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      {
        Real mean  = param_.mean_[k][j];
        Real sigma = param_.sigma_()[j];
        sum += Law::Normal::lpdf(p_data()->elt(i,j), mean, sigma);
      }
      return sum;
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
    { return this->nbCluster()*this->nbVariable()+this->nbVariable();}
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void Gaussian_sj<Array>::randomInit()
{
  // compute the initial mean
  this->randomMean();
  // compute the standard deviation
  Array2DPoint<Real> variance(p_data()->cols(), 0.);
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    variance += p_tik()->col(k).transpose()
               *(*p_data() - (Const::Vector<Real>(p_data()->rows()) * param_.mean_[k])
                ).square()
                ;
  }
  // store the standard deviation
  param_.sigma_() = (variance /= this->nbSample()).sqrt();
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Gaussian_sj<Array>::randomInit() done\n");
#endif
}

/* Compute the weighted mean and the common standard deviation. */
template<class Array>
bool Gaussian_sj<Array>::mStep()
{
  // compute the means
  if (!this->updateMean()) return false;
  // compute the standard deviation
  Array2DPoint<Real> variance(p_data()->cols(), 0.);
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    variance += p_tik()->col(k).transpose()
               *(*p_data() - (Const::Vector<Real>(p_data()->rows()) * param_.mean_[k])
                ).square()
                ;
  }
//  if (variance.nbAvailableValues() != this->nbVariable()) return false;
//  if ((variance > 0.).template cast<int>().sum() != this->nbVariable()) return false;
  // compute the standard deviation
  param_.sigma_() = (variance /= this->nbSample()).sqrt();
  return true;
}

} // namespace STK

#endif /* STK_GAUSSIAN_SJ_H */
