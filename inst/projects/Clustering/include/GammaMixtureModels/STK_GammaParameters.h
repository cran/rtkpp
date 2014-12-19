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
 * Purpose: define the gamma parameters structures.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_GammaParameters.h
 *  @brief In this file we define the parameters class for the Gamma mixture models
 **/

#ifndef STK_GAMMAPARAMETERS_H
#define STK_GAMMAPARAMETERS_H

#include <cmath>

#include "../STK_MixtureParamStat.h"
#include "Arrays/include/STK_Const_Arrays.h"
#include "Arrays/include/STK_Array2DPoint.h"
#include "Arrays/include/STK_Display.h"

#include "STatistiK/include/STK_Law_Gamma.h"


namespace STK
{

/** @ingroup Clustering
 *  @brief Interface base class for the parameters of a multivariate model.
  */
template<class Parameters>
class GammaParametersBase : public IRecursiveTemplate<Parameters>
{
  protected:
    typedef IRecursiveTemplate<Parameters> Base;
    /** default constructor.*/
    inline GammaParametersBase(): tk_(0) {}
    /** copy constructor.*/
    inline GammaParametersBase( GammaParametersBase const& param)
                              : tk_(param.tk_), mean_(param.mean_), meanLog_(param.meanLog_), variance_(param.variance_)
    {}
    /** Destructor */
    inline ~GammaParametersBase() {}
  public:
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline GammaParametersBase& operator=( GammaParametersBase const& other)
    {
      tk_ = other.tk_;
      mean_ = other.mean_;
      meanLog_ = other.meanLog_;
      variance_ = other.variance_;
      return *this;
    }
    /** @return the j-th shape value */
    inline Real shape(int j) const {return this->asDerived().shapeImpl(j);}
    /** @return the j-th scale value */
    inline Real scale(int j) const {return this->asDerived().scaleImpl(j);}
    /** Value of the sum of the k-th column of the tik */
    Real tk_;
    /** vector of the mean of the observations */
    Array2DPoint<Real> mean_;
    /** vector of the mean log of the observations */
    Array2DPoint<Real> meanLog_;
    /** vector of the variance of the observations */
    Array2DPoint<Real> variance_;
    /** compute the log Likelihood of an observation.
     *  @param rowData the observation
     **/
    template<class RowVector>
    Real computeLnLikelihood( RowVector const& rowData) const
    {
      Real sum =0.;
      for (Integer j= rowData.begin(); j < rowData.end(); ++j)
      { sum += Law::Gamma::lpdf(rowData[j], shape(j), scale(j));}
      return sum;
    }
};


/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gamma_ajk_bj model.
 **/
class Gamma_ajk_bjk_Parameters: public GammaParametersBase<Gamma_ajk_bjk_Parameters>
{
  public:
    typedef GammaParametersBase<Gamma_ajk_bjk_Parameters> Base;
    /** default constructor */
    inline Gamma_ajk_bjk_Parameters()
                      : Base(), shape_(), scale_(),stat_shape_(), stat_scale_() {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Gamma_ajk_bjk_Parameters( Gamma_ajk_bjk_Parameters const& param)
                                   : Base(param)
                                   , stat_shape_(param.stat_shape_)
                                   , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~Gamma_ajk_bjk_Parameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Gamma_ajk_bjk_Parameters& operator=( Gamma_ajk_bjk_Parameters const& other)
    {
      Base::operator =(other);
      shape_ = other.shape_;
      scale_ = other.scale_;
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** @return the j-th shape value */
    inline Real shapeImpl(int j) const {return shape_[j];}
    /** @return the j-th scale value */
    inline Real scaleImpl(int j) const {return scale_[j];}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      shape_.resize(range); shape_ = 1.;
      scale_.resize(range); scale_ = 1.;
      mean_.resize(range); mean_ = 1.;
      meanLog_.resize(range); meanLog_ = 0.;
      variance_.resize(range); variance_ = 1.;
      stat_shape_.initialize(range);
      stat_scale_.initialize(range);
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_shape_.update(shape_); stat_scale_.update(scale_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_shape_.release(); stat_scale_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      shape_ = stat_shape_.param_;
      stat_shape_.release();
      scale_ = stat_scale_.param_;
      stat_scale_.release();
    }
    /** vector of the shape */
    Array2DPoint<Real> shape_;
    /** vector of the scale */
    Array2DPoint<Real> scale_;
    /** Array of the statistics */
    MixtureStatVector stat_shape_;
    /** Array of the statistics */
    MixtureStatVector stat_scale_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gamma_ajk_bk model.
 **/
class Gamma_ajk_bk_Parameters: public GammaParametersBase<Gamma_ajk_bk_Parameters>
{
  public:
    typedef GammaParametersBase<Gamma_ajk_bk_Parameters> Base;
    /** default constructor */
    inline Gamma_ajk_bk_Parameters()
                       : Base(), shape_(), scale_(1.), stat_shape_(), stat_scale_() {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Gamma_ajk_bk_Parameters( Gamma_ajk_bk_Parameters const& param)
                                  : Base(param)
                                  , shape_(param.shape_)
                                  , scale_(param.scale_)
                                  , stat_shape_(param.stat_shape_)
                                  , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~Gamma_ajk_bk_Parameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Gamma_ajk_bk_Parameters& operator=( Gamma_ajk_bk_Parameters const& other)
    {
      Base::operator =(other);
      shape_ = other.shape_;
      scale_ = other.scale_;
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** @return the j-th shape value */
    inline Real shapeImpl(int j) const {return shape_[j];}
    /** @return the j-th scale value */
    inline Real scaleImpl(int j) const {return scale_;}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      shape_.resize(range);    shape_ = 1.;
      mean_.resize(range);     mean_  = 1.;
      meanLog_.resize(range);  meanLog_ = 0.;
      variance_.resize(range); variance_ = 1.;
      stat_shape_.initialize(range);
      stat_scale_.initialize();
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_shape_.update(shape_); stat_scale_.update(scale_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_shape_.release(); stat_scale_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      shape_ = stat_shape_.param_;
      stat_shape_.release();
      scale_ = stat_scale_.param_;
      stat_scale_.release();
    }
    /** vector of the shape */
    Array2DPoint<Real> shape_;
    /** vector of the scale */
    Real scale_;
    /** Array of the statistics */
    MixtureStatVector stat_shape_;
    /** Array of the statistics */
    MixtureStatReal stat_scale_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gamma_ajk_bj model.
 **/
class Gamma_ajk_bj_Parameters: public GammaParametersBase<Gamma_ajk_bj_Parameters>
{
  public:
    typedef GammaParametersBase<Gamma_ajk_bj_Parameters> Base;
    /** default constructor */
    inline Gamma_ajk_bj_Parameters() : Base(), shape_(), p_scale_(0) {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Gamma_ajk_bj_Parameters( Gamma_ajk_bj_Parameters const& param)
                                  : Base(param)
                                  , shape_(param.shape_)
                                  , p_scale_(param.p_scale_)
                                  , stat_shape_(param.stat_shape_)

    {}
    /** destructor */
    inline ~Gamma_ajk_bj_Parameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Gamma_ajk_bj_Parameters& operator=( Gamma_ajk_bj_Parameters const& other)
    {
      Base::operator =(other);
      shape_ = other.shape_;
      p_scale_ = other.p_scale_;
      stat_shape_ = other.stat_shape_;
      return *this;
    }
    /** @return the j-th shape value */
    inline Real shapeImpl(int j) const {return shape_[j];}
    /** @return the j-th scale value */
    inline Real scaleImpl(int j) const {return p_scale_->elt(j);}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      shape_.resize(range); shape_ = 1.;
      mean_.resize(range); mean_ = 1.;
      meanLog_.resize(range); meanLog_ = 0.;
      variance_.resize(range); variance_ = 1.;
      stat_shape_.initialize(range);
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_shape_.update(shape_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_shape_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      shape_ = stat_shape_.param_;
      stat_shape_.release();
    }
    /** vector of the shape */
    Array2DPoint<Real> shape_;
    /** vector of the scale */
    Array2DPoint<Real>* p_scale_;
    /** Array of the statistics */
    MixtureStatVector stat_shape_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gamma_ajk_b model.
 **/
class Gamma_ajk_b_Parameters: public GammaParametersBase<Gamma_ajk_b_Parameters>
{
  public:
    typedef GammaParametersBase<Gamma_ajk_b_Parameters> Base;
    /** default constructor */
    inline Gamma_ajk_b_Parameters() : Base(), shape_(), p_scale_(0) {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Gamma_ajk_b_Parameters( Gamma_ajk_b_Parameters const& param)
                                 : Base(param)
                                 , shape_(param.shape_)
                                 , p_scale_(param.p_scale_)
                                 , stat_shape_(param.stat_shape_)
    {}
    /** destructor */
    inline~Gamma_ajk_b_Parameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Gamma_ajk_b_Parameters& operator=( Gamma_ajk_b_Parameters const& other)
    {
      Base::operator =(other);
      shape_ = other.shape_;
      p_scale_ = other.p_scale_;
      stat_shape_ = other.stat_shape_;
      return *this;
    }
    /** @return the j-th shape value */
    inline Real shapeImpl(int j) const {return shape_[j];}
    /** @return the j-th scale value */
    inline Real scaleImpl(int j) const {return *p_scale_;}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      shape_.resize(range); shape_ = 1.;
      mean_.resize(range); mean_ = 1.;
      meanLog_.resize(range); meanLog_ = 0.;
      variance_.resize(range); variance_ = 1.;
      stat_shape_.initialize(range);
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_shape_.update(shape_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_shape_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      shape_ = stat_shape_.param_;
      stat_shape_.release();
    }
    /** vector of the shape */
    Array2DPoint<Real> shape_;
    /** pointer on the scale */
    Real* p_scale_;
    /** Array of the statistics */
    MixtureStatVector stat_shape_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gamma_ak_bjk model.
 **/
class Gamma_ak_bjk_Parameters: public GammaParametersBase<Gamma_ak_bjk_Parameters>
{
  public:
    typedef GammaParametersBase<Gamma_ak_bjk_Parameters> Base;
    /** default constructor */
    inline Gamma_ak_bjk_Parameters() : Base(), shape_(1), scale_() {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Gamma_ak_bjk_Parameters( Gamma_ak_bjk_Parameters const& param)
                                  : Base(param)
                                  , shape_(param.shape_)
                                  , scale_(param.scale_)
                                  , stat_shape_(param.stat_shape_)
                                  , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~Gamma_ak_bjk_Parameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Gamma_ak_bjk_Parameters& operator=( Gamma_ak_bjk_Parameters const& other)
    {
      Base::operator =(other);
      shape_ = other.shape_;
      scale_ = other.scale_;
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** @return the j-th shape value */
    inline Real shapeImpl(int j) const {return shape_;}
    /** @return the j-th scale value */
    inline Real scaleImpl(int j) const {return scale_[j];}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      scale_.resize(range);   scale_ = 1.;
      mean_.resize(range);    mean_ = 1.;
      meanLog_.resize(range); meanLog_ = 0.;
      variance_.resize(range);variance_ = 1.;
      stat_scale_.initialize(range);
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_shape_.update(shape_); stat_scale_.update(scale_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_shape_.release(); stat_scale_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      shape_ = stat_shape_.param_;
      stat_shape_.release();
      scale_ = stat_scale_.param_;
      stat_scale_.release();
    }
    /** vector of the shape */
    Real shape_;
    /** vector of the scale */
    Array2DPoint<Real> scale_;
    /** Array of the statistics */
    MixtureStatReal stat_shape_;
    /** Array of the statistics */
    MixtureStatVector stat_scale_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gamma_ak_bk model.
 **/
class Gamma_ak_bk_Parameters: public GammaParametersBase<Gamma_ak_bk_Parameters>
{
  public:
    typedef GammaParametersBase<Gamma_ak_bk_Parameters> Base;
    /** default constructor */
    inline Gamma_ak_bk_Parameters() : Base(), shape_(1), scale_(1) {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Gamma_ak_bk_Parameters( Gamma_ak_bk_Parameters const& param)
                                  : Base(param)
                                  , shape_(param.shape_)
                                  , scale_(param.scale_)
                                  , stat_shape_(param.stat_shape_)
                                  , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~Gamma_ak_bk_Parameters() {}
    /** @return the j-th shape value */
    inline Real shapeImpl(int j) const {return shape_;}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Gamma_ak_bk_Parameters& operator=( Gamma_ak_bk_Parameters const& other)
    {
      Base::operator =(other);
      shape_ = other.shape_;
      scale_ = other.scale_;
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** @return the j-th scale value */
    inline Real scaleImpl(int j) const {return scale_;}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      mean_.resize(range);    mean_ = 1.;
      meanLog_.resize(range); meanLog_ = 0.;
      variance_.resize(range);variance_ = 1.;
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_shape_.update(shape_); stat_scale_.update(scale_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_shape_.release(); stat_scale_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      shape_ = stat_shape_.param_;
      stat_shape_.release();
      scale_ = stat_scale_.param_;
      stat_scale_.release();
    }
   /** vector of the shape */
    Real shape_;
    /** vector of the scale */
    Real scale_;
    /** Array of the statistics */
    MixtureStatReal stat_shape_;
    /** Array of the statistics */
    MixtureStatReal stat_scale_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gamma_ak_bj model.
 **/
class Gamma_ak_bj_Parameters: public GammaParametersBase<Gamma_ak_bj_Parameters>
{
  public:
    typedef GammaParametersBase<Gamma_ak_bj_Parameters> Base;
    /** default constructor */
    inline Gamma_ak_bj_Parameters() : Base(), shape_(1), p_scale_(0) {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Gamma_ak_bj_Parameters( Gamma_ak_bj_Parameters const& param)
                                  : Base(param)
                                  , shape_(param.shape_)
                                  , p_scale_(param.p_scale_)
                                  , stat_shape_(param.stat_shape_)
    {}
    /** destructor */
    inline ~Gamma_ak_bj_Parameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Gamma_ak_bj_Parameters& operator=( Gamma_ak_bj_Parameters const& other)
    {
      Base::operator =(other);
      shape_ = other.shape_;
      p_scale_ = other.p_scale_;
      stat_shape_ = other.stat_shape_;
      return *this;
    }
    /** @return the j-th shape value */
    inline Real shapeImpl(int j) const {return shape_;}
    /** @return the j-th scale value */
    inline Real scaleImpl(int j) const {return p_scale_->elt(j);}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      mean_.resize(range);    mean_ = 1.;
      meanLog_.resize(range); meanLog_ = 0.;
      variance_.resize(range);variance_ = 1.;
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_shape_.update(shape_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_shape_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      shape_ = stat_shape_.param_;
      stat_shape_.release();
    }
    /** vector of the shape */
    Real shape_;
    /** vector of the scale */
    Array2DPoint<Real>* p_scale_;
    /** Array of the statistics */
    MixtureStatReal stat_shape_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gamma_ak_b model.
 **/
class Gamma_ak_b_Parameters: public GammaParametersBase<Gamma_ak_b_Parameters>
{
  public:
    typedef GammaParametersBase<Gamma_ak_b_Parameters> Base;
    /** default constructor */
    inline Gamma_ak_b_Parameters() : Base(), shape_(1), p_scale_(0) {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Gamma_ak_b_Parameters( Gamma_ak_b_Parameters const& param)
                                : Base(param)
                                , shape_(param.shape_)
                                , p_scale_(param.p_scale_)
                                , stat_shape_(param.stat_shape_)
    {}
    /** destructor */
    inline ~Gamma_ak_b_Parameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Gamma_ak_b_Parameters& operator=( Gamma_ak_b_Parameters const& other)
    {
      Base::operator =(other);
      shape_ = other.shape_;
      p_scale_ = other.p_scale_;
      stat_shape_ = other.stat_shape_;
      return *this;
    }
    /** @return the j-th shape value */
    inline Real shapeImpl(int j) const {return shape_;}
    /** @return the j-th scale value */
    inline Real scaleImpl(int j) const {return *p_scale_;}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      mean_.resize(range);    mean_ = 1.;
      meanLog_.resize(range); meanLog_ = 0.;
      variance_.resize(range);variance_ = 1.;
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_shape_.update(shape_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_shape_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      shape_ = stat_shape_.param_;
      stat_shape_.release();
    }
    /** vector of the shape */
    Real shape_;
    /** vector of the scale */
    Real* p_scale_;
    /** Array of the statistics */
    MixtureStatReal stat_shape_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gamma_aj_bjk model.
 **/
class Gamma_aj_bjk_Parameters: public GammaParametersBase<Gamma_aj_bjk_Parameters>
{
  public:
    typedef GammaParametersBase<Gamma_aj_bjk_Parameters> Base;
    /** default constructor */
    inline Gamma_aj_bjk_Parameters() : Base(), p_shape_(0), scale_() {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Gamma_aj_bjk_Parameters( Gamma_aj_bjk_Parameters const& param)
                                  : Base(param)
                                  , p_shape_(param.p_shape_)
                                  , scale_(param.scale_)
                                  , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline~Gamma_aj_bjk_Parameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Gamma_aj_bjk_Parameters& operator=( Gamma_aj_bjk_Parameters const& other)
    {
      Base::operator =(other);
      p_shape_ = other.p_shape_;
      scale_ = other.scale_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** @return the j-th shape value */
    inline Real shapeImpl(int j) const {return p_shape_->elt(j);}
    /** @return the j-th scale value */
    inline Real scaleImpl(int j) const {return scale_[j];}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      scale_.resize(range);    scale_ = 1.;
      mean_.resize(range);     mean_ = 1.;
      meanLog_.resize(range);  meanLog_ = 0.;
      variance_.resize(range); variance_ = 1.;
      stat_scale_.initialize(range);
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_scale_.update(scale_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_scale_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      scale_ = stat_scale_.param_;
      stat_scale_.release();
    }
    /** vector of the shape */
    Array2DPoint<Real>* p_shape_;
    /** vector of the scale */
    Array2DPoint<Real> scale_;
    /** Array of the statistics */
    MixtureStatVector stat_scale_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gamma_aj_bk model.
 **/
class Gamma_aj_bk_Parameters: public GammaParametersBase<Gamma_aj_bk_Parameters>
{
  public:
    typedef GammaParametersBase<Gamma_aj_bk_Parameters> Base;
    /** default constructor */
    inline Gamma_aj_bk_Parameters() : Base(), p_shape_(0), scale_(1.) {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Gamma_aj_bk_Parameters( Gamma_aj_bk_Parameters const& param)
                                 : Base(param)
                                 , p_shape_(param.p_shape_)
                                 , scale_(param.scale_)
                                , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~Gamma_aj_bk_Parameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Gamma_aj_bk_Parameters& operator=( Gamma_aj_bk_Parameters const& other)
    {
      Base::operator =(other);
      p_shape_ = other.p_shape_;
      scale_ = other.scale_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** @return the j-th shape value */
    inline Real shapeImpl(int j) const {return p_shape_->elt(j);}
    /** @return the j-th scale value */
    inline Real scaleImpl(int j) const {return scale_;}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      mean_.resize(range); mean_ = 1.;
      meanLog_.resize(range); meanLog_ = 0.;
      variance_.resize(range); variance_ = 1.;
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_scale_.update(scale_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_scale_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      scale_ = stat_scale_.param_;
      stat_scale_.release();
    }
    /** vector of the shape */
    Array2DPoint<Real>* p_shape_;
    /** vector of the scale */
    Real scale_;
    /** Array of the statistics */
    MixtureStatReal stat_scale_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gamma_a_bjk model.
 **/
class Gamma_a_bjk_Parameters: public GammaParametersBase<Gamma_a_bjk_Parameters>
{
  public:
    typedef GammaParametersBase<Gamma_a_bjk_Parameters> Base;
    /** default constructor */
    inline Gamma_a_bjk_Parameters() : Base(), p_shape_(0), scale_() {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Gamma_a_bjk_Parameters( Gamma_a_bjk_Parameters const& param)
                                  : Base(param)
                                  , p_shape_(param.p_shape_)
                                  , scale_(param.scale_)
                                  , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline~Gamma_a_bjk_Parameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Gamma_a_bjk_Parameters& operator=( Gamma_a_bjk_Parameters const& other)
    {
      Base::operator =(other);
      p_shape_ = other.p_shape_;
      scale_ = other.scale_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** @return the j-th shape value */
    inline Real shapeImpl(int j) const {return *p_shape_;}
    /** @return the j-th scale value */
    inline Real scaleImpl(int j) const {return scale_[j];}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      scale_.resize(range);    scale_ = 1.;
      mean_.resize(range);     mean_  = 1.;
      meanLog_.resize(range);  meanLog_ = 0.;
      variance_.resize(range); variance_ = 1.;
      stat_scale_.initialize(range);
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_scale_.update(scale_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_scale_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      scale_ = stat_scale_.param_;
      stat_scale_.release();
    }
    /** vector of the shape */
    Real* p_shape_;
    /** vector of the scale */
    Array2DPoint<Real> scale_;
    /** Array of the statistics */
    MixtureStatVector stat_scale_;
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gamma_a_bk model.
 **/
class Gamma_a_bk_Parameters: public GammaParametersBase<Gamma_a_bk_Parameters>
{
  public:
    typedef GammaParametersBase<Gamma_a_bk_Parameters> Base;
    /** default constructor */
    inline Gamma_a_bk_Parameters() : Base(), p_shape_(0), scale_(1.) {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline Gamma_a_bk_Parameters( Gamma_a_bk_Parameters const& param)
                                 : Base(param)
                                 , p_shape_(param.p_shape_)
                                 , scale_(param.scale_)
                                 , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~Gamma_a_bk_Parameters() {}
    /** overwrite the parameters with other.
     *  @param other the parameters to copy
     **/
    inline Gamma_a_bk_Parameters& operator=( Gamma_a_bk_Parameters const& other)
    {
      Base::operator =(other);
      p_shape_ = other.p_shape_;
      scale_ = other.scale_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** @return the j-th shape value */
    inline Real shapeImpl(int j) const {return *p_shape_;}
    /** @return the j-th scale value */
    inline Real scaleImpl(int j) const {return scale_;}
    /** resize the parameters.
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      mean_.resize(range);    mean_ = 1.;
      meanLog_.resize(range); meanLog_ = 0.;
      variance_.resize(range); variance_ = 1.;
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_scale_.update(scale_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_scale_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      scale_ = stat_scale_.param_;
      stat_scale_.release();
    }
    /** vector of the shape */
    Real* p_shape_;
    /** vector of the scale */
    Real scale_;
    /** Array of the statistics */
    MixtureStatReal stat_scale_;
};

} // namespace STK

#endif /* STK_GAMMAPARAMETERS_H */
