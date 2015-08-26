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
 * Project:  stkpp::AAModels
 * Purpose:  Interface base class for AA models.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_GaussianAAModel.h
 *  @brief In this file we declare the class GaussianAAModel for the
 *  auto-associative models.
 **/

#ifndef STK_GAUSSIANAAMODEL_H
#define STK_GAUSSIANAAMODEL_H

#include "STK_IAAModel.h"
#include "Arrays/include/STK_Array2DSquare.h"

#include <StatModels/include/STK_IStatModel.h>
#include <STatistiK/include/STK_MultiLaw_Normal.h>
#include <STatistiK/include/STK_Stat_MultivariateReal.h>

namespace STK
{
/** @ingroup AAModels
 *  @brief Gaussian Auto-Associative models.
 *  A Gaussian Auto-Associative model is a p-dimensional vector \f$\mathbf{y}\f$
 *  with projection function \f$\mathbf{P}\f$ and regression function
 *  \f$\mathbf{R}\f$, if it can be written
 *  \f[
 *    \mathbf{y} = \mathbf{Q}
 *    \left(
 *    \begin{pmatrix}
 *      x_1 \\ \vdots \\
 *      x_d \\
 *      \tilde{r}_{d+1}(\mathbf{x}) \\ \vdots \\
 *      \tilde{r}_p(\mathbf{x})
 *    \end{pmatrix}
 *     + \tilde{\varepsilon}
 *    \right)
 *    = R(\mathbf{x})+ \varepsilon,
 *  \f]
 *  where the \f$\tilde{r}_j(\mathbf{x})\f$, \f$ d+1 \leq j \leq p\f$, are
 *  arbitrary real functions from \f$\mathbb{R}^d\f$ to \f$\mathbb{R}\f$.
 *
 *  The vector \f$\mathbf{x}\f$ is a \f$d\f$-dimensional Gaussian random vector:
 *  \f[
 *    \mathbf{x} \sim \mathcal{N}(\mu_x, \Sigma_x)
 *  \f]
 *  with covariance matrix \f$\Sigma_x = \mathrm{Diag}(\sigma_1^2, \ldots, \sigma_d^2)\f$.
 *
 *  The Gaussian noise \f$\tilde{\varepsilon}\f$ is centered with the following
 *  covariance matrix \f$\Sigma_\varepsilon =
 *  \mathrm{diag}(0,\ldots,0,\sigma^2,\ldots,\sigma^2)\f$.
 *
 *  The GaussianModel class is a factory class which compute the covariance
 *  matrix of x, the residual covariance and the number of free parameters
 *  of the model. It can be sub-classed or used by any class.
 **/
template<class Array>
class GaussianAAModel : public IAAModel<Array>
                      , public IStatModel<Array>
{
  public:
    using IAAModel<Array>::p_regressor_;
    using IAAModel<Array>::p_reducer_;
    using IAAModel<Array>::p_reduced_;
    using IAAModel<Array>::p_predicted_;
    using IAAModel<Array>::p_residuals_;
    using IAAModel<Array>::dim;
    /** Constructor.
     *  @param p_workData a pointer on the data set to process
     **/
    GaussianAAModel( Array* p_workData);
    /** Constructor.
     *  @param workData a reference on the data set to process
     **/
    GaussianAAModel( Array& workData);
    /** virtual destuctor. */
    inline virtual ~GaussianAAModel() {}
    /** @return the ln-likelihood of the projected data set */
    inline Real const& projectedLnLikelihood() const { return projectedLnLikelihood_;}
    /** @return the ln-likelihood of the residuals */
    inline Real const& residualLnLikelihood() const { return residualLnLikelihood_;}
    /** @return the covariance matrix of the projected the data set */
    inline ArraySquareX const& projectedCovariance() const { return projectedCovariance_;}
    /** @return the covairance matrix of the residuals */
    inline ArraySquareX const& residualCovariance() const { return residualCovariance_;}
    /** @return the variance of the residuals */
    inline Real const& residualVariance() const { return residualVariance_;}
    /** Set a new working data set.
     *  @param workData the working data set to use
     **/
    virtual void setWorkData(Array& workData);
    /** compute the covariance matrix of the projected data set.
     *  This method is set public as the projected covariance can be computed
     *  only the first time the data set is projected.
     **/
    void computeProjectedCovariance();
    /** compute the ln-likelihood of the model */
    void computeModelParameters();

  protected:
    /** @brief compute the number of free parameter of the model.
     *  It is given by the number of parameter of the regression function,
     *  the number of variance and covariance of the projected data set
     *  (d * (d+1))/2 and the variance of the residuals.
     **/
    void computeNbFreeParameters();
    /** compute the covariance matrix of the residuals. */
    void computeResidualCovariance();
    /** @brief compute the ln-likelihood of the projected data set
     *  The projected data set is assumed Gaussian with an arbitrary
     *  covariance Array.
     **/
    void computeProjectedLnLikelihood();
    /** @brief compute the ln-likelihood of the projected data set.
     * The residuals are assumed orthogonal to the the projected
     * data set with a single residual variance.
     **/
    void computeResidualLnLikelihood();

  private:
    /** The covariance matrix of the projected data set */
    ArraySquareX projectedCovariance_;
    /** The covariance matrix of the residuals */
    ArraySquareX residualCovariance_;
    /** The total variance of the residuals */
    Real residualVariance_;
    /** likelihood of the projected data set. */
    Real projectedLnLikelihood_;
    /** likelihood of the residuals. */
    Real residualLnLikelihood_;
};

/* Constructor.
 *  @param p_workData a pointer on the data set to process
 **/
template<class Array>
GaussianAAModel<Array>::GaussianAAModel( Array* p_workData)
                                : IAAModel<Array>(p_workData)
                                , IStatModel<Array>(p_workData)
                                , projectedCovariance_()
                                , residualCovariance_()
                                , residualVariance_(0.)
                                , projectedLnLikelihood_(0.)
                                , residualLnLikelihood_(0.)
{ }

// constructor
template<class Array>
GaussianAAModel<Array>::GaussianAAModel( Array& workData)
                                : IAAModel<Array>(workData)
                                , IStatModel<Array>(workData)
                                , projectedCovariance_()
                                , residualCovariance_()
                                , residualVariance_(0.)
                                , projectedLnLikelihood_(0.)
                                , residualLnLikelihood_(0.)
{ }

/* update the container when the data set is modified. **/
template<class Array>
void GaussianAAModel<Array>::setWorkData(Array& workData)
{
  // update data set and flags for the IAAModel part
  IAAModel<Array>::setWorkData(workData);
  // set dimensions to new size for the IStatModel part
  IStatModel<Array>::setData(workData);
}

/* compute the ln-likelihood of the model */
template<class Array>
void GaussianAAModel<Array>::computeModelParameters()
{
#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("GaussianAAModel::computeModelParameters().\n");
#endif
  // compute the number of free parameters
  computeNbFreeParameters();
  // compute the covariance matrix of the residual and the residual variance
  computeResidualCovariance();
  // compute projected lnLikelihood
  computeProjectedLnLikelihood();
  // compute projected lnLikelihood
  computeResidualLnLikelihood();
  // compute complete nLikelihood
  this->setLnLikelihood(projectedLnLikelihood_ + residualLnLikelihood_);

#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("GaussianAAModel::computeModelParameters() done.\n");
#endif
}


/* @brief compute the number of free parameter of the model. **/
template<class Array>
void GaussianAAModel<Array>::computeNbFreeParameters()
{
  this->setNbFreeParameter(p_regressor_->nbFreeParameter() + dim() * (dim()+1)/2 + 1);
}

/* compute the ln-likelihood of the model */
template<class Array>
void GaussianAAModel<Array>::computeProjectedLnLikelihood()
{
#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("GaussianAAModel::computeProjectedLnLikelihood().\n");
#endif
  // range of the column to use
  Range cols = Range(p_reduced_->beginCols(), dim());
  // create a reference with the first columns of the reduced data
  Array reducedData(*p_reduced_, p_reduced_->rows(), cols);
  // create a reference with the first columns of the reduced data
  ArraySquareX reducedCovariance(projectedCovariance(), cols);
  // compute first part of the ln-likelihood
  Point mean(cols, 0.);
  MultiLaw::Normal<Point> normalP(mean, reducedCovariance);
  projectedLnLikelihood_ = normalP.lnLikelihood(reducedData);

#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("GaussianAAModel::computeProjectedLnLikelihood() done.\n");
  stk_cout << _T("projectedLnLikelihood_ = ") << projectedLnLikelihood_ << _T("\n");
#endif
}

/* compute the ln-likelihood of the model */
template<class Array>
void GaussianAAModel<Array>::computeResidualLnLikelihood()
{
#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("GaussianAAModel::computeResidualsLnLikelihood().\n");
#endif
  // compute constant part and determinant part of the log-likelihood
  residualLnLikelihood_ = ( Const::_LNSQRT2PI_ + 0.5*std::log(residualVariance_ ))
                          * (dim() - this->nbVariable()) * this->nbSample();
  // compute second part of the log-likelihood
  for (int i=p_residuals_->beginRows(); i<p_residuals_->endRows(); i++)
  { residualLnLikelihood_ -= p_residuals_->row(i).norm2()/(2.*residualVariance_);}

#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("GaussianAAModel::computeResidualsLnLikelihood() done.\n");
  stk_cout << _T("residualLnLikelihood_ = ") << residualLnLikelihood_ << _T("\n");
#endif
}

/* compute the covariance matrix of the projected data set. */
template<class Array>
void GaussianAAModel<Array>::computeProjectedCovariance()
{
#ifdef STK_AAMODELS_DEBUG
  if (!p_reduced_)
  { STKRUNTIME_ERROR_NO_ARG(GaussianAAModel::computeProjectedCovariance,p_reduced_ is not set);}
#endif
#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("GaussianAAModel::computeProjectedCovariance().\n");
#endif
  Stat::covariance(*p_reduced_, projectedCovariance_);
#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("GaussianAAModel::computeProjectedCovariance() done.\n");
#endif
}
/* compute the covariance matrix of the residuals. */
template<class Array>
void GaussianAAModel<Array>::computeResidualCovariance()
{
#ifdef STK_AAMODELS_DEBUG
  if (!p_residuals_)
  { STKRUNTIME_ERROR_NO_ARG(GaussianAAModel::computeResidualCovariance,p_residuals_ is not set);}
#endif
#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("in GaussianAAModel::computeResidualCovariance().\n");
#endif
  Stat::covariance(*p_residuals_, residualCovariance_);
  residualVariance_ = (residualCovariance_.trace())/Real(this->nbVariable()-dim());
#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("GaussianAAModel::computeResidualCovariance() done.\n");
  stk_cout << _T("residualVariance_ = ") << residualVariance_ << _T("\n");
#endif
}


} // namespace STK

#endif //STK_GAUSSIANAAMODEL_H
