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
 * Project: stkpp::AAModels
 * Purpose: Interface bas class for all AA Models.
 * Author:  iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_IAAModel.h
 *  @brief  In this file we define the Interface Base class for all AutoAssociative Models
 **/

#ifndef STK_IAAMODEL_H
#define STK_IAAMODEL_H

#include <Arrays/include/STK_Array2DPoint.h>
#include <Arrays/include/STK_Array2DVector.h>

#include <STatistiK/include/STK_Stat_MultivariateReal.h>
#include <STatistiK/include/STK_Stat_Transform.h>

#include "Reduct/include/STK_IReducer.h"
#include "Regress/include/STK_IRegression.h"

#ifdef STK_AAMODELS_VERBOSE
#include "Arrays/include/STK_Display.h"
#endif


namespace STK
{
/** @ingroup @ingroup AAModels
 * @brief Interface base class for all Auto-Associative models.
 *
 * A function \f$ g \f$ is an auto-associative function
 * of dimension d if it is a map from \f$ \mathbb{R}^p \f$ to
 * \f$ \mathbb{R}^p \f$ that can be written \f$ g=R\circ P \f$ where
 * P (the ``Projecttion'' or the ``Reduction'') is a map from \f$\mathbb{R}^p\f$
 * to \f$ \mathbb{R}^d \f$ and R (the ``Regression'') is a map from
 * \f$\mathbb{R}^d\f$ to \f$\mathbb{R}^p\f$ .
 *
 * The IAAModel class is an abstract base class for all these AAM and a
 * factory method that propose some tools that can be used in derived classes.
 * End user can have to set a reducer and a regressor in order to
 * implement this interface.
 */
template<class Array>
class IAAModel
{
  private:
    /** regression type */
    typedef IRegression<Array, Array, Vector > Regressor;
    /** reducer type */
    typedef IReducer<Array, Vector > Reducer;

  protected:
    /** Constructor.
     *  @param p_workData A pointer on the the working data set
     **/
    IAAModel( Array* p_workData);
    /** Constructor.
     *  @param workData the working data set
     **/
    IAAModel( Array& workData);
  public:
    /** destructor. */
    ~IAAModel();
    /** @return the working data set */
    inline Array const& workData() const { return *p_workData_;}
    /**  @return a reference on the Index object */
    inline Reducer* const& p_reducer() const { return p_reducer_;}
    /** @return a reference on the regresor object */
    inline Regressor* const& p_regressor() const { return p_regressor_;}
    /** @return the reduced data set */
    inline Array* const& p_reduced() const { return p_reduced_;}
    /** @return a pointer on the predicted values */
    inline Array* const& p_predicted() const { return p_predicted_;}
    /** @return A ptr on the residuals */
    inline Array* const& p_residuals() const { return p_residuals_;}
    /** @return the dimension of the model */
    inline int dim() const { return dim_;}
    /** @return @c true if the data set is centered, @c false otherwise */
    inline bool isCentered() const { return isCentered_;}
    /** @return @c true if the data set is standardized, @c false otherwise */
    inline bool isStandardized() const { return isStandardized_;}
    /** @return the mean of the data set */
    inline PointX const& mean() const { return mean_;}
    /** @return the standard deviation of the data set */
    inline PointX const& std() const { return std_;}
    /** @param dim the dimension of the model to set */
    void setDimension( int const& dim);
    /** @param workData the working data set to treat */
    void setWorkData( Array& workData);
    /** @param p_reducer a pointer on the reduction dimension method to use */
    void setReducer( Reducer* p_reducer);
    /** @param p_regressor a pointer on the regression method to use */
    void setRegressor( Regressor* p_regressor);
    /** delete the reducer set to this model by the method @c setReducer. */
    void freeReducer();
    /** delete the regressor set to this model by the method @c setRegressor.*/
    void freeRegressor();
    /** center the data set workData_. */
    void center();
    /** weighted centering of the data set.
     *  @param weights the weights of the samples
     **/
    void center( Vector const& weights);
    /** standardize the data set. */
    void standardize();
    /** weighted standardization the data set.
     *  @param weights the weights of the samples
     **/
    void standardize( Vector const& weights);
    /** compute the reduction of the data set and store the result in
     *  the @c p_reduced_ container. The reducer p_reducer have to be set.
     **/
    void reductionStep();
    /** compute the weighted dimension reduction of the data set and store the
     *  result in the @c p_reduced_ container.
     **/
    void reductionStep( Vector const& weights);
    /** compute the regression of the original data set and set the results in
     *  @c p_predicted and @c p_residuals. **/
    void regressionStep();
    /** compute the weighted regression  of the original data set using the
     *  reduced data set as predictor and set the results in @c p_predicted and
     *  @c p_residuals.
     **/
    void regressionStep( Vector const& weights);

    /** decenter the predicted data set. This will invalidate the results
     *  of p_regressor_.  */
    void decenterResults();
    /** destandardize the predicted data set and the residuals.  This will
     *  invalidate the results of p_regressor_. */
    void destandardizeResults();

  protected:
    /** pointer on the regression method. The regression method will
     *  be chosen in the derived class.
     **/
    Regressor* p_regressor_;
    /** pointer on the reeducer. The dimension reduction method will be chosen
     *  in the derived class.
     **/
    Reducer* p_reducer_;
    /** Array of the local data set. */
    Array* p_workData_;
    /** Array of the reduced data set : the data set is shared with
     * @c p_reducer and set when the @c regression method is call.
     **/
    Array* p_reduced_;
    /** Array of the predicted data set: the data set is shared with
     * @c p_regressor and set when the @c regression method is call. **/
    Array* p_predicted_;
    /** Array of the residuals: the data set is shared with @c p_regressor
     *  and set when the @c regression method is call.
     **/
    Array* p_residuals_;

  private:
    /** The dimension of the AA Model. */
    int dim_;
    /** vector of the means of the input data set. */
    PointX mean_;
    /** vector of the standard deviation of the input data set. */
    PointX std_;
    /** a boolean @c true if the working data set is centered, @c false
     *  otherwise. */
    bool isCentered_;
    /** a boolean @c true if the working data set is standardized, @c false
     *  otherwise */
    bool isStandardized_;
};

/* Constructor.
 *  @param p_workData A pointer on the the working data set
 **/
template<class Array>
IAAModel<Array>::IAAModel( Array* p_workData)
                        : p_regressor_(0)
                        , p_reducer_(0)
                        , p_workData_(p_workData)
                        , p_reduced_(0)
                        , p_predicted_(0)
                        , p_residuals_(0)
                        , dim_(0)
                        , mean_()
                        , std_()
                        , isCentered_(false)
                        , isStandardized_(false)
{}

template<class Array>
IAAModel<Array>::IAAModel( Array& workData)
                        : p_regressor_(0)
                        , p_reducer_(0)
                        , p_workData_(&workData)
                        , p_reduced_(0)
                        , p_predicted_(0)
                        , p_residuals_(0)
                        , dim_(0)
                        , mean_()
                        , std_()
                        , isCentered_(false)
                        , isStandardized_(false)
{ }

/* destructor */
template<class Array>
IAAModel<Array>::~IAAModel()
{}

/* set the dimension of the model */
template<class Array>
void IAAModel<Array>::setDimension( const int& dim)
{ dim_ = dim;}

/* set the working set with the Data to treat.
 * @param p_reducer a pointer on the reduction dimension method to use
 */
template<class Array>
void IAAModel<Array>::setWorkData( Array& workData)
{
  p_workData_ = &workData;
  isCentered_     = false;
  isStandardized_ = false;
}
/* set the reduction dimension method.
 * @param p_reducer a pointer on the reduction dimension method to use
 */
template<class Array>
void IAAModel<Array>::setReducer( Reducer* p_reducer)
{ p_reducer_ = p_reducer;}
/* set the regression method.
 * @param p_regressor a pointer on the regresssion method to use
 */
template<class Array>
void IAAModel<Array>::setRegressor( Regressor* p_regressor)
{ p_regressor_ = p_regressor;}

/** delete the reducer allocated set to this model. */
template<class Array>
void IAAModel<Array>::freeReducer()
{
  if (p_reducer_) delete p_reducer_;
  p_reducer_ = 0;
}
/** delete the regressor allocated to this model.*/
template<class Array>
void IAAModel<Array>::freeRegressor()
{
  if (p_regressor_) delete p_regressor_;
  p_regressor_ = 0;
}

/* standardize the local data set */
template<class Array>
void IAAModel<Array>::center()
{
#ifdef STK_AAMODELS_DEBUG
  if (!p_workData_)
    STKRUNTIME_ERROR_NO_ARG(IAAModel<Array>::center,workData is not initialized);
#endif
  // we have to decenter in case they have been centered with weights
  if (isCentered_)
  {
    Stat::decenter(*p_workData_, mean_);
    isCentered_ = false;
  }
  Stat::center(*p_workData_, mean_);
  isCentered_ = true;
}

/* center the local data set */
template<class Array>
void IAAModel<Array>::center(Vector const& weights)
{
#ifdef STK_AAMODELS_DEBUG
  if (!p_workData_)
    STKRUNTIME_ERROR_NO_ARG(IAAModel<Array>::center,workData is not initialized);
#endif
  // we have to decenter in case of this is not the same weights
  if (isCentered_)
  {
    Stat::decenter(*p_workData_, mean_);
    isCentered_ = false;
  }
  // center
  Stat::center(*p_workData_, weights, mean_);
  isCentered_ = true;
}

/* standardize the local data set */
template<class Array>
void IAAModel<Array>::standardize()
{
#ifdef STK_AAMODELS_DEBUG
  if (!p_workData_)
    STKRUNTIME_ERROR_NO_ARG(IAAModel<Array>::standardize,workData is not initialized);
#endif
  // we have to destandardize in case they have been standardized with weights
  if (isStandardized_)
  {
    Stat::destandardize(*p_workData_, mean_, std_);
    isStandardized_ = false;
  }
  Stat::standardize(*p_workData_, mean_, std_);
  isStandardized_ = true;
}

/* standardize the local data set */
template<class Array>
void IAAModel<Array>::standardize(Vector const& weights)
{
#ifdef STK_AAMODELS_DEBUG
  if (!p_workData_)
    STKRUNTIME_ERROR_NO_ARG(IAAModel<Array>::standardize,workData is not initialized);
#endif
  // we have to destandardize in case of this is not the same weights
  if (isStandardized_)
  {
    Stat::destandardize(*p_workData_, mean_, std_);
    isStandardized_ = false;
  }
  // standardize
  Stat::standardize(*p_workData_, weights, mean_, std_);
  isStandardized_ = true;
}


/* compute the dimension reduction **/
template<class Array>
void IAAModel<Array>::reductionStep()
{
#ifdef STK_AAMODELS_DEBUG
  if (!p_reducer_)
    STKRUNTIME_ERROR_NO_ARG(IAAModel<Array>::reduction,reducer is not set);
#endif
  // compute axis
  p_reducer_->setDimension(dim_);
  p_reducer_->run();
  // compute matrix multiplication
  p_reduced_= p_reducer_->p_reduced();
}

/* compute the weighted dimension reduction **/
template<class Array>
void IAAModel<Array>::reductionStep( Vector const& weights)
{
#ifdef STK_AAMODELS_DEBUG
  if (!p_reducer_)
    STKRUNTIME_ERROR_NO_ARG(IAAModel<Array>::reduction,reducer is not set);
#endif
  // compute axis
  p_reducer_->setDimension(dim_);
  p_reducer_->run(weights);
  // get the reduced data set
  p_reduced_= p_reducer_->p_reduced();
}

/* compute the regression **/
template<class Array>
void IAAModel<Array>::regressionStep()
{
#ifdef STK_AAMODELS_DEBUG
  if (!p_regressor_)
    STKRUNTIME_ERROR_NO_ARG(IAAModel<Array>::regressionStep,regressor is not set);
#endif
  // compute regression
  p_regressor_->run();
  // get results
  p_predicted_ = p_regressor_->p_predicted();
  p_residuals_ = p_regressor_->p_residuals();
}
/* compute the weighted regression **/
template<class Array>
void IAAModel<Array>::regressionStep( Vector const& weights)
{
#ifdef STK_AAMODELS_DEBUG
  if (!p_regressor_)
    STKRUNTIME_ERROR_NO_ARG(IAAModel<Array>::regressionStep,regressor is not set);
#endif
  p_regressor_->run(weights);
  // get results
  p_predicted_ = p_regressor_->p_predicted();
  p_residuals_ = p_regressor_->p_residuals();
}

/* destandardize the predicted result and residuals */
template<class Array>
void IAAModel<Array>::decenterResults()
{
#ifdef STK_AAMODELS_DEBUG
  if (!p_predicted_)
    STKRUNTIME_ERROR_NO_ARG(IAAModel<Array>::decenterResults,predictions are not computed);
#endif
  Stat::decenter(*p_predicted_, mean_);
}
/* destandardize the predicted result and residuals */
template<class Array>
void IAAModel<Array>::destandardizeResults()
{
#ifdef STK_AAMODELS_DEBUG
  if (!p_predicted_)
    STKRUNTIME_ERROR_NO_ARG(IAAModel<Array>::destandardizeResults,predictions are not computed);
  if (!p_residuals_)
    STKRUNTIME_ERROR_NO_ARG(IAAModel<Array>::destandardizeResults,residuals are not computed);
#endif
  Stat::destandardize(*p_predicted_, mean_, std_);
  Stat::destandardize(*p_residuals_, std_);
}


} // namespace STK

#endif /* STK_IAAMODEL_H */
