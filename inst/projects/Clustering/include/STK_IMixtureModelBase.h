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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_IMixtureModelBase.h
 *  @brief In this file we define the interface base class for mixture models.
 **/

#ifndef STK_IMIXTUREMODELBASE_H
#define STK_IMIXTUREMODELBASE_H

#include "STKernel/include/STK_Range.h"
#include "STK_Clust_Util.h"
#include "Arrays/include/STK_CArrayPoint.h"
#include "Arrays/include/STK_CArrayVector.h"
#include "Arrays/include/STK_CArray.h"

namespace STK
{
/** @ingroup Clustering
 *  @brief Base class for Mixture model.
 *
 * All the mixture parameters: pk, tik, nk and zi are accessed by pointer and
 * set to this class using the method
 * @code
 *   void setMixtureParameters(CPointX const* p_pk, CPointX const* p_nk, ArrayXX const* p_tik, CVectorXi const* p_zi);
 * @endcode
 * so that they can be used in a composed model.
 * @sa MixtureComposer
 */
class IMixtureModelBase
{
  protected:
    /** default constructor
     * @param nbCluster number of cluster */
    IMixtureModelBase( int nbCluster);
    /** copy constructor.
     *  @note the pointer are initialized to 0.
     *  @param model the model to clone
     **/
    IMixtureModelBase( IMixtureModelBase const& model);

  public:
    /** destructor */
    ~IMixtureModelBase();
    /** @return the number of cluster */
    inline int nbCluster() const { return nbCluster_;}
    /** @return the total available observations */
    inline int nbSample() const { return nbSample_;}
    /** @return the Log of the total available observations */
    inline Real lnNbSample() const
    { return (nbSample_ <= 0) ? -Arithmetic<Real>::infinity() : std::log((Real)nbSample_);}
    /** @return the total available variables */
    inline int nbVariable() const { return nbVariable_;}

    /** @return the proportions of each mixtures */
    inline CPointX const* p_pk() const { return p_pk_;};
    /** @return the estimated number of individuals of each mixtures */
    inline CPointX const* p_nk() const { return p_nk_;};
    /** @return the tik probabilities */
    inline CArrayXX const* p_tik() const { return p_tik_;};
    /** @return  the zi class label */
    inline CVectorXi const* p_zi() const { return p_zi_;};
    /** Set the parameters of the  mixture model using external pk, tik and zi.
     *  @param p_pk pointer on the proportions of the mixture model
     *  @param p_nk pointer on the estimated number of individuals of each mixtures
     *  @param p_tik pointer on the posterior probabilities
     *  @param p_zi pointer on the class labels
     * */
    void setMixtureParameters( CPointX const* p_pk
                             , CPointX const* p_nk
                             , CArrayXX const* p_tik
                             , CVectorXi const* p_zi
                             );

  protected:
    /** The proportions of each mixtures */
    CPointX const* p_pk_;
    /** The estimated number of individuals of each mixtures */
    CPointX const* p_nk_;
    /** The tik probabilities */
    CArrayXX const* p_tik_;
    /** The zik class label */
    CVectorXi const* p_zi_;
    /** Set the number of sample of the model
     *  @param nbSample number of sample of the model
     * */
    inline void setNbSample( int nbSample) { nbSample_ = nbSample;}
    /** Set the number of variables of the model
     *  @param nbVariable number of variables of the model
     * */
    inline void setNbVariable( int nbVariable) { nbVariable_ = nbVariable;}

  private:
    /** number of cluster. */
    int nbCluster_;
    /** total available samples */
    int nbSample_;
    /** total available variables */
    int nbVariable_;
};

} // namespace STK

#endif /* STK_IMIXTUREMODELBASE_H */
