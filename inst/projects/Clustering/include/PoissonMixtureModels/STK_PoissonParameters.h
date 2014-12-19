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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 22 juil. 2013
 * Purpose: define the gamma parameters.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_PoissonParameters.h
 *  @brief In this file we define the parameters class for the Poisson mixture models
 **/

#ifndef STK_POISSONPARAMETERS_H
#define STK_POISSONPARAMETERS_H

#include <cmath>

#include "Arrays/include/STK_Const_Arrays.h"
#include "Arrays/include/STK_Array2DPoint.h"
#include "Arrays/include/STK_Display.h"

#include "STatistiK/include/STK_Law_Poisson.h"

namespace STK
{
/** @ingroup Clustering
 *  @brief Interface base class for the parameters of a poisson model.
  */
template<class Parameters>
class PoissonParametersBase : public IRecursiveTemplate<Parameters>
{
  protected:
    typedef IRecursiveTemplate<Parameters> Base;
    /** default constructor.*/
    inline PoissonParametersBase(): Base(), tk_(0) {}
    /** copy constructor.*/
    inline PoissonParametersBase( PoissonParametersBase const& param)
                                : tk_(param.tk_)
    {}
    /** Destructor */
    inline ~PoissonParametersBase() {}
  public:
    /** resize the parameters (default implementation).
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range) {}
    /** @return the j-th scale value */
    inline Real lambda(int j) const {return this->asDerived().lambdaImpl(j);}
    /** Value of the sum of the k-th column of the tik */
    Real tk_;
    /** compute the log Likelihood of an observation.
     *  @param rowData the observation
     **/
    template<class RowVector>
    Real computeLnLikelihood( RowVector const& rowData) const
    {
      Real sum =0.;
      for (Integer j= rowData.begin(); j < rowData.end(); ++j)
      { sum += Law::Poisson::lpdf(rowData[j], lambda(j));}
      return sum;
    }
};


/** @ingroup Clustering
 *  Structure encapsulating the parameters of the poisson_ljk model.
 **/
class Poisson_ljk_Parameters: public PoissonParametersBase<Poisson_ljk_Parameters>
{
  public:
    typedef PoissonParametersBase<Poisson_ljk_Parameters> Base;
    /** default constructor */
    inline Poisson_ljk_Parameters() : Base(), lambda_() {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Poisson_ljk_Parameters( Poisson_ljk_Parameters const& param)
                                 : Base(param), lambda_(param.lambda_)
                                 , stat_lambda_(param.stat_lambda_)
    {}
    /** destructor */
    inline ~Poisson_ljk_Parameters() {}
    /** @return the j-th scale value */
    inline Real lambdaImpl(int j) const {return lambda_[j];}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      lambda_.resize(range); lambda_ = 1.;
      stat_lambda_.initialize(range);
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_lambda_.update(lambda_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_lambda_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      lambda_ = stat_lambda_.param_;
      stat_lambda_.release();
    }
    /** vector of the lambda */
    Array2DPoint<Real> lambda_;
    /** Array of the statistics */
    MixtureStatVector stat_lambda_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the poisson_ljk model.
 **/
class Poisson_lk_Parameters: public PoissonParametersBase<Poisson_lk_Parameters>
{
  public:
    typedef PoissonParametersBase<Poisson_lk_Parameters> Base;
    /** default constructor */
    inline Poisson_lk_Parameters() : Base(), lambda_() {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Poisson_lk_Parameters( Poisson_lk_Parameters const& param)
                                : Base(param), lambda_(param.lambda_)
                                , stat_lambda_(param.stat_lambda_)
    {}
    /** destructor */
    inline ~Poisson_lk_Parameters() {}
    /** @return the j-th scale value */
    inline Real lambdaImpl(int j) const {return lambda_;}
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration) { stat_lambda_.update(lambda_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults() { stat_lambda_.release();}
    /** set the parameters stored in stat_lambda_ and release stat_lambda_. */
    void setParameters()
    {
      lambda_ = stat_lambda_.param_;
      stat_lambda_.release();
    }
    /** lambda value */
    Real lambda_;
    /** Array of the statistics */
    MixtureStatReal stat_lambda_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the poisson_ljlk model.
 **/
class Poisson_ljlk_Parameters: public PoissonParametersBase<Poisson_ljlk_Parameters>
{
  public:
    typedef PoissonParametersBase<Poisson_ljlk_Parameters> Base;
    /** default constructor */
    inline Poisson_ljlk_Parameters() : Base(), lambdak_(1.), p_lambdaj_(0) {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Poisson_ljlk_Parameters( Poisson_ljlk_Parameters const& param)
                                  : Base(param)
                                  , lambdak_(param.lambdak_)
                                  , p_lambdaj_(param.p_lambdaj_)
                                  , stat_lambdak_(param.stat_lambdak_)

    {}
    /** destructor */
    inline ~Poisson_ljlk_Parameters() {}
    /** @return the j-th lambda value */
    inline Real lambdaImpl(int j) const {return lambdak_ * p_lambdaj_->elt(j);}
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_lambdak_.update(lambdak_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_lambdak_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      lambdak_ = stat_lambdak_.param_;
      stat_lambdak_.release();
    }
    /** value of lambdak */
    Real lambdak_;
    /** vector of the lambdaj */
    PointX* p_lambdaj_;
    /** statistics */
    MixtureStatReal stat_lambdak_;
};


} // namespace STK

#endif /* STK_POISSONPARAMETERS_H */
