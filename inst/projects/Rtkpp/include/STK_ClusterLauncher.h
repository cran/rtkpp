/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  rtkpp
 * created on: 4 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_ClusterLauncher.h
 *  @brief In this file we define the ClusterLauncher which
 *  construct properly a mixture model.
 **/


#ifndef STK_CLUSTERLAUNCHER_H
#define STK_CLUSTERLAUNCHER_H


namespace STK
{

/**The ClusterLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
class ClusterLauncher : public IRunnerBase
{
  public:
    /** constructor.
     * @param model a reference on the current model
     * @param strategy the strategy defined in R
     **/
    ClusterLauncher( SEXP model, SEXP nbCluster, SEXP modelNames, SEXP strategy, SEXP critName );
    /** constructor with a list of component.
     * @param model a reference on the current model
     * @param strategy the strategy defined in R
     **/
    ClusterLauncher( SEXP model, SEXP nbCluster, SEXP strategy, SEXP critName );
    /** destructor. */
    virtual ~ClusterLauncher();
    /** run the estimation */
    bool run();
    /** @return the model */
    inline Rcpp::S4 const& s4_model() const { return s4_model_;}

  protected:
    /** strategy from the R side */
    Rcpp::S4              s4_model_;
    /** strategy from the R side */
    Rcpp::S4              s4_strategy_;
    /** vector with the number of cluster to try */
    Rcpp::IntegerVector   v_nbCluster_;
    /** vector with the model names to try */
    Rcpp::CharacterVector v_modelNames_;
    /** character string with the model selection criterion name */
    std::string           critName_;

  private:
    /** get the parameters */
    void getParameters(Rcpp::S4& s4_component, std::string const& idData);
    /** get the diagonal Gaussian parameters */
    void getDiagGaussianParameters(Rcpp::S4& s4_component, std::string const& idData);
    /** get the gamma parameters */
    void getGammaParameters(Rcpp::S4& s4_component, std::string const& idData);
    /** get the gamma parameters */
    void getCategoricalParameters(Rcpp::S4& s4_component, std::string const& idData);

    /** data handler */
    RDataHandler handler_;
    /** manager */
    MixtureManager<RDataHandler> manager_;
    /** pointer on the main composer */
    IMixtureComposer* p_composer_;
    /** Is the model with heterogeneous data ? */
    bool isHeterogeneous_;
    /** Select the best model among the modelNames and nbCluster given.
     *  @return the value of the best criteria.
     **/
    Real selectSingleBestModel();
    /** Select the best model among the modelNames and nbCluster given.
     *  @return the value of the best criteria.
     **/
    Real selectHeteroBestModel();
};


}  // namespace STK

#endif /* STK_CLUSTERLAUNCHER_H */
