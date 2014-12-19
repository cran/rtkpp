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
 * Purpose: define the categorical parameters.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_CategoricalParameters.h
 *  @brief In this file we define the parameters class for the Categorical
 *  mixture models.
 **/

#ifndef STK_CATEGORICALPARAMETERS_H
#define STK_CATEGORICALPARAMETERS_H

#include <cmath>

#include "Arrays/include/STK_Const_Arrays.h"
#include "Arrays/include/STK_Display.h"

#include "../STK_MixtureParamStat.h"

namespace STK
{

/** @ingroup Clustering
 *  @brief Interface base class for the parameters of a diagonal categorical
 *  multivariate model.
  */
template<class Parameters>
class CategoricalParametersBase : public IRecursiveTemplate<Parameters>
{
  protected:
    typedef IRecursiveTemplate<Parameters> Base;
    /** default constructor.*/
    inline CategoricalParametersBase() {}
    /** constructor with specified range
     *  @param range the range of the variables
     **/
    inline CategoricalParametersBase( Range const& range) {}
    /** copy constructor.*/
    inline CategoricalParametersBase( CategoricalParametersBase const& param): Base(param) {}
    /** destructor */
    inline ~CategoricalParametersBase() {}

  public:
    /** @return the j-th probability value of the l-th modality */
    inline Real proba(int j, int l) const {return this->asDerived().probaImpl(j, l);}
    /** @return the j-th probability distribution */
    inline Array2DVector<Real> const& proba(int j) const {return this->asDerived().probaImpl(j);}

    /** compute the log Likelihood of an observation.
     *  @param rowData the observation
     **/
    template<class RowVector>
    Real computeLnLikelihood( RowVector const& rowData) const
    {
      Real sum =0.;
      for (Integer j= rowData.begin(); j < rowData.end(); ++j)
      {
        Real prob = proba(j, rowData[j]);
        if (prob <= 0.) return -Arithmetic<Real>::infinity();
        sum += std::log(prob);
      }
      return sum;
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the categorical_pjk model.
 **/
class Categorical_pjkParameters: public CategoricalParametersBase<Categorical_pjkParameters>
{
  public:
    typedef CategoricalParametersBase<Categorical_pjkParameters> Base;
    /** default constructor */
    inline Categorical_pjkParameters() : Base(), proba_() {}
    /** constructor with specified range
     *  @param range the range of the variables
     **/
    inline Categorical_pjkParameters( Range const& range): Base(range), proba_(range)
    { stat_proba_.resize(range);}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Categorical_pjkParameters( Categorical_pjkParameters const& param)
                                    : Base(param), proba_(param.proba_)
    {}
    /** destructor */
    inline ~Categorical_pjkParameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Categorical_pjkParameters& operator=( Categorical_pjkParameters const& other)
    {
      proba_ = other.proba_;
      stat_proba_ = other.stat_proba_;
      return *this;
    }
    /** @return the j-th probability value of the l-th modality */
    inline Real probaImpl(int j, int l) const {return proba_[j][l];}
    /** @return the j-th probability distribution */
    inline Array2DVector<Real> const& probaImpl(int j) const { return proba_[j];}
    /** resize the set of parameter
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      proba_.resize(range);
      stat_proba_.resize(range);
    }
    /** utility function allowing to resize the probability vector with a
     *  given Range for the modalities.
     *  @param rangeMod the range of the modalities of the categorical distribution
     **/
    inline void initializeParameters(Range const& rangeMod)
    {
      for(int j=proba_.begin(); j< proba_.end(); j++)
      {
        proba_[j].resize(rangeMod);
        proba_[j] = 1./rangeMod.size();
        stat_proba_[j].initialize(rangeMod);
      }
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    {
      for(int j=proba_.begin(); j< proba_.end(); j++)
      { stat_proba_[j].update(proba_[j]);}
    }
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    {
      for(int j=proba_.begin(); j< proba_.end(); j++)
      { stat_proba_[j].release();}
    }
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      for(int j=proba_.begin(); j< proba_.end(); j++)
      { proba_[j] = stat_proba_[j].param_;
        stat_proba_[j].release();
      }
    }
    /** Array of the probabilities */
    Array1D< VectorX > proba_;
    /** Array of the statistics */
    Array1D< MixtureStatVector > stat_proba_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the categorical_pk model.
 **/
class Categorical_pkParameters: public CategoricalParametersBase<Categorical_pkParameters>
{
  public:
    typedef CategoricalParametersBase<Categorical_pkParameters> Base;
    /** default constructor */
    inline Categorical_pkParameters() : Base(), proba_(),stat_proba_() {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Categorical_pkParameters( Categorical_pkParameters const& param)
                                   : Base(param), proba_(param.proba_), stat_proba_(param.stat_proba_)
    {}
    /** destructor */
    inline ~Categorical_pkParameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Categorical_pkParameters& operator=( Categorical_pkParameters const& other)
    {
      proba_ = other.proba_;
      stat_proba_ = other.stat_proba_;
      return *this;
    }
    /** @return the j-th probability value of the l-th modality (does not depend of j) */
    inline Real probaImpl(int j, int l) const { return proba_[l];}
    /** @return the j-th probability distribution (does not depend of j) */
    inline Array2DVector<Real> const& probaImpl(int j) const {return proba_;}
    /** resize the set of parameter
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range) {}
    /** utility function allowing to resize the probability vector with a
     *  given Range for the modalities.
     *  @param rangeMod the range of modalities of the categorical distribution
     **/
    inline void initializeParameters(Range const& rangeMod)
    {
      proba_.resize(rangeMod);
      proba_ = 1./rangeMod.size();
      stat_proba_.initialize(rangeMod);
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_proba_.update(proba_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_proba_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      proba_ = stat_proba_.param_;
      stat_proba_.release();
    }
    /** probabilities of each modalities */
    Array2DVector<Real> proba_;
    /** Array of the statistics */
    MixtureStatVector stat_proba_;
};

} // namespace STK

#endif /* STK_CATEGORICALPARAMETERS_H */
