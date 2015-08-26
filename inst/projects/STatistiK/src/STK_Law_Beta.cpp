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
 * Project:  stkpp::STatistiK::Law
 * created on: 23 janv. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Beta.cpp
 *  @brief In this file we implement the Beta distibution law.
 **/


#ifndef IS_RTKPP_LIB
#include "../include/STK_Law_Beta.h"
#include "../include/STK_Law_Gamma.h"
#include <Analysis/include/STK_Funct_betaRatio.h>
#endif

//
namespace STK
{
namespace Law
{

#ifndef IS_RTKPP_LIB

/*  Generate a pseudo Beta random variate. */
Real Beta::rand() const
{
  Real g1 = Law::Gamma::rand(alpha_, 1.);
  return g1/(g1+Law::Gamma::rand(beta_, 1.));
}

/*
 *  Give the value of the pdf at x.
 */
Real Beta::pdf( Real const& x) const
{
  // trivial cases
  if (!Arithmetic<Real>::isFinite(x)||(x<0.)||(x>1)) return 0.0;
  // compute result
  return 0.;
}

/*
 * Give the value of the log-pdf at x.
 */
Real Beta::lpdf( Real const& x) const
{
  // check NA value
  if (isNA(x)) return Arithmetic<Real>::NA();
  // check parameter
  if (Arithmetic<Real>::isInfinite(x))
    return -Arithmetic<Real>::infinity();
  // compute result
  return 0.;
}

/*
 * The cumulative distribution function at t.
 */
Real Beta::cdf( Real const& t) const
{
  // check NA value
  if (isNA(t)) return Arithmetic<Real>::NA();
  // compute result
  return (Arithmetic<Real>::isInfinite(t)) ? (t < 0.) ? 0.0 : 1.0
                                           :  Funct::betaRatio(alpha_,beta_,t, false);
}
    
/*
 * The inverse cumulative distribution function at p.
 */
Real Beta::icdf( Real const& p) const
{
  // check parameter
  if ((p > 1.) || (p < 0.))
    STKDOMAIN_ERROR_1ARG(Beta::icdf,p,argument outside [0;1]);
  // trivial cases
  if (p == 0.) return 0.;
  if (p == 1.) return 1.;
  // result 
  return 0.;
}

/*  Generate a pseudo Beta random variate with the specified parameters.
 *  (static)
 */
Real Beta::rand( Real const& alpha, Real const& beta)
{
  Real g1 = Law::Gamma::rand(alpha, 1.);
  return g1/(g1+Law::Gamma::rand(beta, 1.));
}

Real Beta::pdf(const Real& x, const Real& alpha, const Real& beta)
{
  // trivial cases
  if (!Arithmetic<Real>::isFinite(x)||(x<0.)||(x>1)) return 0.0;
  // compute result
  return 0.;
}

Real Beta::lpdf(const Real& x, const Real& alpha, const Real& beta)
{
  // check NA value
  if (isNA(x)) return Arithmetic<Real>::NA();
  // check parameter
  if (Arithmetic<Real>::isInfinite(x))
    return -Arithmetic<Real>::infinity();
  // compute result
  return 0.;
}

Real Beta::cdf(const Real& t, const Real& alpha, const Real& beta)
{
  // check NA value
  if (isNA(t)) return Arithmetic<Real>::NA();
  // compute result
  return (Arithmetic<Real>::isInfinite(t)) ? (t < 0.) ? 0.0 : 1.0
                                           :  Funct::betaRatio(alpha,beta,t, false);
}

Real Beta::icdf(const Real& p, const Real& alpha, const Real& beta)
{
  // check parameter
  if ((p > 1.) || (p < 0.))
    STKDOMAIN_ERROR_1ARG(Beta::icdf,p,argument outside [0;1]);
  // trivial cases
  if (p == 0.) return 0.;
  if (p == 1.) return 1.;
  // result
  return 0.;
}

#endif /* !IS_RTKPP_LIB */

} // namespace Law

} // namespace STK

