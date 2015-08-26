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
 * created on: Oct 24, 2014
 * Author:   Serge Iovleff
 **/

/** @file STK_KernelGaussian.h
 *  @brief In this file we define the KernelGaussian_sk class
 **/

#ifndef STK_KERNELGAUSSIAN_H
#define STK_KERNELGAUSSIAN_H

#include <Arrays/include/STK_Array2DPoint.h>
#include <STatistiK/include/STK_Stat_Online.h>
#include "../STK_IMixtureModel.h"

namespace STK
{

//forward declaration, to allow for recursive template
class KernelGaussian_sk;
class KernelGaussian_s;

namespace Clust
{
/** @ingroup Clustering
 *  Traits class for the KernelGaussian_sk traits policy. */
template<>
struct MixtureTraits< KernelGaussian_sk >
{
  typedef ArrayXX Array;
  typedef typename Array::Type       Type;
  typedef ParametersHandler<KernelGaussian_sk_>  ParamHandler;
};

/** @ingroup Clustering
 *  Traits class for the KernelGaussian_sk traits policy. */
template<>
struct MixtureTraits< KernelGaussian_s >
{
  typedef ArrayXX Array;
  typedef typename Array::Type       Type;
  typedef ParametersHandler<KernelGaussian_s_>  ParamHandler;
};

} // namespace Clust
/** @ingroup Clustering
 *  Base class for the Kernel models
 **/
template<class Derived>
struct KernelHandlerBase: public IRecursiveTemplate<Derived>
{
  /** @return the value of lambda of the kth cluster and jth variable */
  inline Real const& sigma2(int k) const { return this->asDerived().sigma2Impl(k);}
  /** @return the value of lambda of the kth cluster and jth variable */
  inline Real const& dim(int k) const { return this->asDerived().dimImpl(k);}
};


/** Specialization of the ParametersHandler struct for KernelGaussian_sk models*/
template <>
struct ParametersHandler<Clust::KernelGaussian_sk_>: public KernelHandlerBase<  ParametersHandler<Clust::KernelGaussian_sk_> >
{
  /** standard deviation and statistics */
  MixtureParametersSet<Real> sigma2_;
    /** vector of the standard deviations */
  /** standard deviation and statistics */
  PointX dim_;
  /** @return the value of lambda of the kth cluster and jth variable */
  inline Real const& sigma2Impl(int k) const { return sigma2_[k];}
  /** @return the value of lambda of the kth cluster and jth variable */
  inline Real const& dimImpl(int k) const { return dim_[k];}
  /** copy operator */
  inline ParametersHandler& operator=( ParametersHandler const& other)
  {
    sigma2_ = other.sigma2_;
    dim_    = other.dim_;
    return *this;
  }
  /** copy operator using an array/expression storing the values */
  template<class Array>
  inline ParametersHandler& operator=( ExprBase<Array> const& param)
  {
    for (int k= param.beginRows(); k < param.endRows(); k++)
    { sigma2_[k] = param(k, baseIdx);
      dim_[k]    = param(k, baseIdx+1);
    }
    return *this;
  }

  /** default constructor */
  inline ParametersHandler(int nbCluster): sigma2_(nbCluster), dim_(nbCluster)
  {}
  /** copy constructor.
   * @param param the parameters to copy.
   **/
  inline ParametersHandler( ParametersHandler const& param)
                          : sigma2_(param.sigma2_)
                          , dim_(param.dim_)
  {}
  /** Initialize the parameters with an array/expression of value */
  template<class Array>
  inline ParametersHandler( int nbCluster, ExprBase<Array> const& param)
                          : sigma2_(nbCluster)
                          , dim_(nbCluster)
  {
    for (int k= param.beginRows(); k < param.endRows(); k++)
    { sigma2_[k] = param(k, baseIdx);
      dim_[k]    = param(k, baseIdx+1);
    }
  }
  /** destructor */
  inline ~ParametersHandler() {}
  /** Store the intermediate results of the Mixture.
   *  @param iteration Provides the iteration number beginning after the burn-in period.
   **/
  inline void storeIntermediateResults(int iteration)
  { sigma2_.storeIntermediateResults(iteration); }
  /** Release the stored results. This is usually used if the estimation
   *  process failed.
   **/
  inline void releaseIntermediateResults()
  { sigma2_.releaseIntermediateResults(); }
  /** set the parameters stored in stat_proba_ and release stat_proba_. */
  inline void setParameters()
  { sigma2_.setParameters();}
};

/** Specialization of the ParametersHandler struct for KernelGaussian_s models*/
template <>
struct ParametersHandler<Clust::KernelGaussian_s_>: public KernelHandlerBase<  ParametersHandler<Clust::KernelGaussian_s_> >
{
  /** value of the standard deviation */
  Real sigma2_;
  /** Array of the statistics */
  Stat::Online<Real, Real> stat_sigma2_;
  /** vector of the dimensions */
  PointX dim_;
  /** copy operator */
  inline ParametersHandler& operator=( ParametersHandler const& other)
  { sigma2_ = other.sigma2_; dim_ = other.dim_;
    stat_sigma2_ = other.stat_sigma2_;
    return *this;
  }
  /** copy operator using an array/expression storing the values */
  template<class Array>
  inline ParametersHandler& operator=( ExprBase<Array> const& param)
  {
    for (int k= param.beginRows(); k < param.endRows(); k++)
    { sigma2_ = param.col(baseIdx).mean();
      dim_[k] = param(k, baseIdx+1);
    }
    return *this;
  }

  /** @return the value of lambda of the kth cluster and jth variable */
  inline Real const& sigma2Impl(int k) const { return sigma2_;}
  /** @return the value of lambda of the kth cluster and jth variable */
  inline Real const& dimImpl(int k) const { return dim_[k];}

  /** default constructor */
  inline ParametersHandler(int nbCluster): sigma2_(1.), stat_sigma2_(), dim_(nbCluster, 1.)
  {}
  /** copy constructor.
   * @param param the parameters to copy.
   **/
  inline ParametersHandler( ParametersHandler const& param)
                          : sigma2_(param.sigma2_)
                          , stat_sigma2_(param.stat_sigma2_)
                          , dim_(param.dim_)
  {}
  /** Initialize the parameters with an array/expression of value */
  template<class Array>
  inline ParametersHandler( int nbCluster, ExprBase<Array> const& param)
                          : sigma2_(), stat_sigma2_(), dim_(nbCluster)
  {
    for (int k= param.beginRows(); k < param.endRows(); k++)
    { sigma2_ = param.col(baseIdx).mean();
      dim_[k] = param(k, baseIdx+1);
    }
  }
  /** destructor */
  inline ~ParametersHandler() {}
  /** Store the intermediate results of the Mixture.
   *  @param iteration Provides the iteration number beginning after the burn-in period.
   **/
  inline void storeIntermediateResults(int iteration)
  { stat_sigma2_.update(sigma2_); }
  /** Release the stored results. This is usually used if the estimation
   *  process failed.
   **/
  inline void releaseIntermediateResults()
  { stat_sigma2_.release(); }
  /** set the parameters stored in stat_proba_ and release stat_proba_. */
  inline void setParameters()
  {
    sigma2_ = stat_sigma2_.mean_;
    stat_sigma2_.release();
  }
};

/** @ingroup Clustering
 *  The Gaussian mixture model @c KernelGaussian_sk is an isotrope Gaussian
 *  mixture model on a kernel space. It has a density function of the
 *  form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k
 *    \sum_{k=1}^K p_k \left(\frac{1}{\sqrt{2\pi}\sigma_k}\right)^{d_k}
 *    \exp\left\{ -\frac{\|\phi(x)-m_k\|^2}{2\sigma_k^2} \right\}
 * \f]
 * where \f$ \phi \f$ denote a feature mapping from the original space to an RKHS.
 *
 * In a KernelGaussian_sk model, the data set refer to the Gram's matrix.
 **/
class KernelGaussian_sk : public IMixtureModel<KernelGaussian_sk >
{
  public:
    typedef IMixtureModel<KernelGaussian_sk > Base;
    typedef ParametersHandler<Clust::KernelGaussian_sk_> ParamHandler;

    using Base::p_tik;
    using Base::param_;
    using Base::p_nk;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline KernelGaussian_sk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline KernelGaussian_sk( KernelGaussian_sk const& model): Base(model) {}
    /** destructor */
    inline ~KernelGaussian_sk() {}
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const { return this->nbCluster();}
    /** set the dimensions of the kernel mixture model using an unique value */
    inline void setDim(Real const& dim)  { param_.dim_ = dim;}
    /** set the dimension of the kernel mixture model */
    inline void setDim(PointX const& dim)  { param_.dim_ = dim;}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      return(- p_data()->elt(i,k)/(2.*param_.sigma2_[k])
             - (std::log(param_.sigma2_[k])+2.*Const::_LNSQRT2PI_)*param_.dim_[k]/2.);
    }
    /** @return an imputation value for the jth variable of the ith sample */
    inline Real impute(int i, int j) const { return 0.;}
    /** @return a simulated value for the jth variable of the ith sample */
    inline Real rand(int i, int j, int k) const { return 0.;}
    /** Initialize randomly the variances of the Gaussian kernel mixture. */
    void randomInit();
    /** update the variances. */
    bool mStep();
};

/** @ingroup Clustering
 *  The Gaussian mixture model @c KernelGaussian_sk is an isotrope Gaussian
 *  mixture model on a kernel space. It has a density function of the
 *  form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k
 *    \sum_{k=1}^K p_k \left(\frac{1}{\sqrt{2\pi}\sigma_k}\right)^{d_k}
 *    \exp\left\{ -\frac{\|\phi(x)-m_k\|^2}{2\sigma_k^2}  \right\}
 * \f]
 * where \f$ \phi \f$ denote a feature mapping from the original space to an RKHS.
 *
 * In a KernelGaussian_sk model, the data set refer to the Gram's matrix.
 **/
class KernelGaussian_s : public IMixtureModel<KernelGaussian_s >
{
  public:
    typedef IMixtureModel<KernelGaussian_s > Base;
    typedef ParametersHandler<Clust::KernelGaussian_s_> ParamHandler;

    using Base::p_tik;
    using Base::param_;
    using Base::p_nk;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline KernelGaussian_s( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline KernelGaussian_s( KernelGaussian_s const& model): Base(model) {}
    /** destructor */
    inline ~KernelGaussian_s() {}
    /** @return a constant reference on the paremeter handler structure.*/
    inline ParamHandler const& getParameters() const { return param_;}
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const { return 1;}
    /** set the dimensions of the kernel mixture model using an unique value */
    inline void setDim(Real const& dim)  { param_.dim_ = dim;}
    /** set the dimensions of the kernel mixture model */
    inline void setDim(PointX const& dim)  { param_.dim_ = dim;}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      return(- p_data()->elt(i,k)/(2.*param_.sigma2_)
             - (std::log(param_.sigma2_)+2.*Const::_LNSQRT2PI_)*param_.dim_[k]/2.);
    }
    /** @return an imputation value for the jth variable of the ith sample */
    inline Real impute(int i, int j) const { return 0.;}
    /** @return a simulated value for the jth variable of the ith sample */
    inline Real rand(int i, int j, int k) const { return 0.;}
    /** Initialize randomly the variances of the Gaussian kernel mixture. */
    void randomInit();
    /** update the variances. */
    bool mStep();
};

/* Initialize randomly the parameters of the Gaussian mixture. */
inline void KernelGaussian_sk::randomInit()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering KernelGaussian_sk::randomInit()\n");
#endif
  param_.sigma2_() = sum( p_data()->prod(*p_tik()) )/ (*p_nk() * param_.dim_)
                 + PointX(p_tik()->cols()).rand(Law::Normal(0, 0.05)).abs();
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("KernelGaussian_sk::randomInit() done\n");
#ifdef STK_MIXTURE_DEBUG
  stk_cout << param_.sigma2_() << "\n";
#endif
#endif
}

/* Compute the weighted means and the weighted standard deviations. */
inline bool KernelGaussian_sk::mStep()
{
#ifdef STK_MIXTURE_DEBUG
  stk_cout << _T("Entering KernelGaussian_sk::mStep()\n");
#endif
  param_.sigma2_() =  sum( p_data()->prod(*p_tik()) )/ (*p_nk() * param_.dim_);
  //if ((param_.sigma2_() <= 0.).any()) return false; // not work with Array1D
#ifdef STK_MIXTURE_DEBUG
  stk_cout << param_.sigma2_() << "\n";
#endif
  return true;
}

/* Initialize randomly the parameters of the Gaussian mixture. */
inline void KernelGaussian_s::randomInit()
{
  // compute the standard deviation
  param_.sigma2_ = p_data()->prod(*p_tik()).sum()/(this->nbSample() * param_.dim_.sum())
                 + std::abs(Law::generator.randGauss(0, 0.05));
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("KernelGaussian_s::randomInit() done\n");
#endif
}

/* Compute the weighted means and the weighted standard deviations. */
inline bool KernelGaussian_s::mStep()
{
  param_.sigma2_ =  ( p_data()->prod( *p_tik() ) ).sum()/(this->nbSample() * param_.dim_.sum());
  if (param_.sigma2_ <= 0.)  return false;
  return true;
}

} // namespace STK

#endif /* STK_KERNELGAUSSIAN_H */
