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

/** @file STK_DiagGaussianParameters.h
 *  @brief In this file we define the parameters class for the Gaussian mixture models
 **/

#ifndef STK_DIAGGAUSSIANPARAMETERS_H
#define STK_DIAGGAUSSIANPARAMETERS_H

#include <cmath>

#include "Arrays/include/STK_Const_Arrays.h"
#include "Arrays/include/STK_Array2DPoint.h"
#include "Arrays/include/STK_Display.h"

#include "STatistiK/include/STK_Law_Normal.h"

#include "../STK_MixtureParamStat.h"

namespace STK
{

/** @ingroup Clustering
 *  Structure encapsulating the parameters of the gaussian_sjk model.
 **/
class KernelGaussian_Parameters
{
  public:
    /** default constructor */
    inline KernelGaussian_Parameters(): sigma_() {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline KernelGaussian_Parameters( KernelGaussian_Parameters const& param)
                                    : sigma_(param.sigma_)
                                    , stat_sigma_(param.stat_sigma_)
    {}
    /** destructor */
    inline ~KernelGaussian_Parameters() {}
    /** @return the j-th sigma value */
    inline Real sigmaImpl(int j) const {return sigma_[j];}
    /** resize the set of parameter
     *  @param range range of the parameters
     **/
    inline void resize(Range const& range)
    {
      sigma_ = 1.;
      stat_sigma_.initialize();
    }
    /** Store the intermediate results of the Mixture.
     *  @param iteration Provides the iteration number beginning after the burn-in period.
     **/
    void storeIntermediateResults(int iteration)
    { stat_sigma_.update(sigma_);}
    /** Release the stored results. This is usually used if the estimation
     *  process failed.
     **/
    void releaseIntermediateResults()
    { stat_sigma_.release();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    void setParameters()
    {
      sigma_ = stat_sigma_.param_;
      stat_sigma_.release();
    }
    /** vector of the standard deviation */
    Real sigma_;
    /** Array of the statistics */
    MixtureStatReal stat_sigma_;
};

} // namespace STK

#endif /* STK_DIAGGAUSSIANPARAMETERS_H */
