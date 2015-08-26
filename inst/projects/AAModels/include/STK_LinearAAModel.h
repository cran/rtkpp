/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff

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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)                                   */

/*
 * Project:  stkpp::AAModels
 * Purpose:  class for AA Linear models.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_LinearAAModel.h
 *  @brief In this file we define the LinearAAModel class.
 **/

#ifndef STK_LINEARAAMODEL_H
#define STK_LINEARAAMODEL_H

#include "STK_GaussianAAModel.h"
#include <Sdk/include/STK_IRunner.h>
#include <Reduct/include/STK_ILinearReduct.h>
#include <Regress/include/STK_MultidimRegression.h>
#include <Arrays/include/STK_Array2D.h>
#include <Arrays/include/STK_Array2DSquare.h>
#include <Arrays/include/STK_Array2D_Functors.h>
#include <Algebra/include/STK_GramSchmidt.h>
#include <STatistiK/include/STK_Law_Normal.h>

#ifdef STK_AAMODELS_VERBOSE
#include <Arrays/include/STK_Display.h>
#endif


namespace STK
{

namespace Law
{
template<class TYPE>
class IUnivLaw;
}

/** @ingroup AAModels
 *  @brief A Linear AutoAssociative Model (LAAM) is a generalization of
 *  the PCA model for projecting variables on an affine plane of lower
 *  dimension.
 * 
 * A SLAAM is a (probabilistic) model of the form
 *  \f[
 *       Y = \boldsymbol{\mu}
 *                      + \sum_{j=1}^{q} \mathrm{Reg}^j(C^j) + \epsilon
 * \f]
 * in \f$ \mathbf{R}^d \f$.
 *
 *  \f$ Y \f$ is an input d dimensional random variable  and \f$ \epsilon \f$
 *  is the residual random variable.
 *  The (generalized) principal components \f$ C^{j} \f$ are real random
 *  variables estimated by a projection of the centered data on a subspace
 *  \f[
 *      C^j = <a^j,Y>.
 * \f]
 * where the \f$ a^j\f$ form an orthonormal set of vector.
 *
 *  The regression functions
 * \f$ \mathrm{Reg}^j : \mathbb{R} \rightarrow  \mathbb{R}^d \f$
 *  are specifics to the model.
 *
 *  A LAAM is an auto-associative model of the form
 *  \f{eqnarray*} Y_i & = & \mu
 *                       + \sum_{j=1}^{q} <a^j, Y_i - \mu> V^{j} + \epsilon_i
 *  \f}
 *
 *  In this implementation, the axis \f$ a^j \f$ are computed using any
 *  arbitrary derived class of Index.
 *
 *  The principal variables \f$ V^j \f$ are computed using least square
 *  regression. As the regression model is linear there is no necessity
 *  to use the @c setRegresor method. it will be created by this object.
 **/
template<class Array>
class LinearAAModel : public IRunnerUnsupervised<Array, Vector>
                    , public GaussianAAModel<Array>
{
  public :
    typedef IRunnerUnsupervised<Array, Vector> Runner;
    typedef GaussianAAModel<Array> Base;
    using Base::p_workData_;
    using Base::p_regressor_;
    using Base::p_reducer_;
    using Base::p_reduced_;
    using Base::p_predicted_;
    using Base::p_residuals_;
    using Base::dim;
    using Base::regressionStep;
    using Base::reductionStep;
    /** constructor. compute the Linear AA models of the matrix data
     *  using the local variance as criteria.
     *  @param data the data set to modelize
     **/
    LinearAAModel( Array const& data);
    /** copy constructor.
     *  @param model the model to copy
     **/
    LinearAAModel( LinearAAModel const& model);
    /** Virtual destructor */
    virtual ~LinearAAModel();
    /** clone pattern */
    inline virtual LinearAAModel* clone() const { return new LinearAAModel(*this);}
    /** run the estimation of the AA model.
     * The behavior of the estimation is the following :
     * - compute the projected data set and set the result in p_reduced_
     * - regress the workData_ set using the p_reduced_ set as predictors
     * - compute the restored data set using the AA model and set the result
     * in the p_predicted_ container.
     * - compute the residuals
     * The data have to be centered before running the computations.
     * @param dim the dimension of the AA Model
     **/
    bool run( int const& dim);
    /** run the estimation of the weighted AA model.
     * The behavior of the estimation is the following :
     * - compute the projected data set and set the result in p_reduced_
     * - regress the data set using the reduced data set as predictors
     * - compute the restored data set using the AA model and set the result
     * in the p_predicted_ container
     * - compute the residuals and set them in p_residual_
     *
     * The data have to be centered before running the computations.
     *
     * @param dim the dimension of the AA Model
     * @param weights the container of the weights
     **/
    bool run( Vector const& weights, int const& dim);
    /** @brief run the estimation of the AA model.
     * The behavior of the estimation is the following :
     * - compute the projected data set and set the result in p_reduced_
     * - regress the workData_ set using the p_reduced_ set as predictors
     * - compute the restored data set using the AA model and set the result
     * in the p_predicted_ container.
     * - compute the residuals
     * - compute the distribution of the projected data set and the residuals
     *
     * The data have to be centered before running the computations.
     **/
    virtual bool run();
    /** run the estimation of the weighted AA model.
     * The behavior of the estimation is the following :
     * - compute the projected data set and set the result in p_reduced_
     * - regress the data set using the reduced data set as predictors
     * - compute the restored data set using the AA model and set the result
     * in the p_predicted_ container
     * - compute the residuals and set them in p_residual_
     * - compute the distribution of the projected data set and the residuals
     *
     * The data have to be centered before running the computations.
     *
     * @param weights the container of the weights
     **/
    virtual bool run(Vector const& weights);
    /** Simulate a centered auto-associative linear model in
     * \f$ \mathbb{R}^p \f$ of the form
     * \f[
     *    X = X.P.P' + \epsilon
     * \f]
     * with \f$ P'P = I_d \f$ and d < p.
     * @param law the law to use in order to simulate the data.
     * @param mu the position parameter of the AA model
     * @param std the standard deviation of the gaussian noise
     * @param proj the simulated projection matrix. The dimension of the
     * container give the dimension of the AA model.
     * @param data the data to simulate. The dimension of the container
     * give the number of the samples and variables.
     **/
    static void simul( const Law::IUnivLaw<Real>& law
                     , Vector const& mu
                     , Real const& std
                     , Array& proj
                     , Array& data
                     );
    protected:
      /** working data set;*/
    Array workData_;
};

/* Constructors.                                                             */
/* Default constructor : compute the Linear AA models of the matrix X
 *  using the local variance as criteria.
**/
template<class Array>
LinearAAModel<Array>::LinearAAModel( Array const& data)
                                   : Runner(), Base(0), workData_(data)
{
  Runner::setData(workData_);
  Base::setWorkData(workData_);
  p_regressor_ = new MultidimRegression<Array, Vector>();
}

/** copy constructor.
 *  @param model the model to copy
 **/
template<class Array>
LinearAAModel<Array>::LinearAAModel( LinearAAModel const& model)
                                   : Runner(), Base(0), workData_(model.workData_)
{
  Runner::setData(workData_);
  Base::setWorkData(workData_);
  p_regressor_ = new MultidimRegression<Array, Vector>();
}



/* Virtual destructor */
template<class Array>
LinearAAModel<Array>::~LinearAAModel()
{ if (p_regressor_) delete p_regressor_;}

/* run the probabilistic AAModel */
template<class Array>
bool LinearAAModel<Array>::run( int const& dim)
{
  this->setDimension(dim);
  return run();
}


/* run the probabilistic AAModel */
template<class Array>
bool LinearAAModel<Array>::run()
{
  // compute AAM
  try
  {
    if (!p_reducer_)
    {
      this->msg_error_ = STKERROR_NO_ARG(LinearAAModel::run(),reducer is not set);
      return false;
    }
    if (!p_regressor_)
    {
      this->msg_error_ = STKERROR_NO_ARG(LinearAAModel::run(),regressor is not be set);
      return false;
    }
    // set p_workData to the reducer
    p_reducer_->setData(*p_workData_);
    // compute the projected data set
    reductionStep();
    // compute the projected covariance
    this->computeProjectedCovariance();
#ifdef STK_AAMODELS_VERBOSE
    stk_cout << _T("In LinearAAModel::run(), reduction done.\n");
#endif
    // set data
    p_regressor_->setY(p_workData_);
    p_regressor_->setX(p_reduced_);
    // compute the regression function
    regressionStep();
#ifdef STK_AAMODELS_VERBOSE
    stk_cout << _T("In LinearAAModel::run(), regression done.\n");
#endif
    this->computeModelParameters();
    // check if data have been standardized or centered
    if (this->isStandardized()) { this->destandardizeResults();}
    else if (this->isCentered()){ this->decenterResults();}
  }
  catch (Exception const& error)
  {
    this->msg_error_ = _T("Error in LinearAAModel::run():\nWhat: ");
    this->msg_error_ += error.error();
    return false;
  }
  return true;
}

/* run the probabilistic AAModel */
template<class Array>
bool LinearAAModel<Array>::run( Vector const& weights, int const& dim)
{
  this->setDimension(dim);
  return run(weights);
}
template<class Array>
bool LinearAAModel<Array>::run( Vector const& weights)
{
  try
  {
    if (!p_reducer_) STKRUNTIME_ERROR_NO_ARG(LinearAAModel::run,reducer is not set);
    if (!p_regressor_) STKRUNTIME_ERROR_NO_ARG(LinearAAModel::run,regressor is not set);
    // set p_workData to the reducer
    p_reducer_->setData(*p_workData_);
    // compute the weighted reduced data set
    reductionStep(weights);
    // compute the projected covariance
    this->computeProjectedCovariance();
#ifdef STK_AAMODELS_VERBOSE
    stk_cout << _T("In LinearAAModel::run(weights), reduction done.\n");
#endif
    // set data
    p_regressor_->setY(p_workData_);
    p_regressor_->setX(p_reduced_);
    // compute the weighted regression vectors
    regressionStep(weights);
#ifdef STK_AAMODELS_VERBOSE
    stk_cout << _T("In LinearAAModel::run(weights), regression done.\n");
#endif
    this->computeModelParameters();
    // check if data have been standardized or centered
    if (this->isStandardized()) { this->destandardizeResults();}
    else if (this->isCentered()){ this->decenterResults();}
  }
  catch (Exception const& error)
  {
    this->msg_error_ = _T("Error in LinearAAModel::run(weights): ");
    this->msg_error_ += error.error();
    return false;
  }
  return true;
}

/* Simulate a centered auto-associative linear model of the form
 * \f[
 *    X = X.P.P' + \epsilon
 * \f]
 * with \f$ P'P = I_d \f$ and d < p.
 * @param law, the law to use in order to simulate the data.
 * @param std the standard deviation of the gaussian noise
 * @param data the data to simulate. The dimension of the container
 * give the number of the samples and variables.
 * @param proj the simulated projection matrix. The dimension of the
 * container give the dimension of the AA model.
 **/
template<class Array>
void LinearAAModel<Array>::simul( const Law::IUnivLaw<Real>& law
                         , Vector const& mu, Real const& std
                         , Array& proj, Array& data
                         )
{
  if (data.cols() != mu.range())
    STKRUNTIME_ERROR_NO_ARG(LinearAAModel::simul,data.cols()!=mu.range());
  // simul AA model
  data.rand(law);
#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("data simulated\n");
#endif
  // choose arbitrary coefficients for proj
  proj.randGauss();
#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("proj simulated\n");
#endif
  // orthonormalize proj
  gramSchmidt(proj);
#ifdef STK_AAMODELS_VERBOSE
  stk_cout << _T("proj orthonormalized\n");
#endif

#ifdef STK_AAMODELS_DEBUG
  stk_cout << _T("proj =") << proj;
#endif
  // compute data
  data = data * (proj * proj.transpose());

  // add noise to the model
  for (int j= data.beginCols(); j< data.endCols(); j++)
  {
    Law::Normal noise(mu[j], std);
    for (int i= data.beginRows(); i< data.endRows(); i++)
      data(i, j) += noise.rand();
  }
}

} // namespace STK

#endif // STK_LINEARAAMODEL_H
