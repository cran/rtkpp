/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015 Serge Iovleff

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
 * Project: stkpp::STatistiK::Law
 * Purpose: Logistic probability distribution.
 * Author:  Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Logistic.h
 *  @brief In this file we define the Logistic probability law class.
 **/

#ifndef IS_RTKPP_LIB
#include "../include/STK_Law_Logistic.h"
#endif

namespace STK
{
namespace Law
{

#ifndef IS_RTKPP_LIB

/* @brief Generate a pseudo logisticized Logistic random variate.
 *
 *  Generate a pseudo logisticized Logistic random variate
 *  with location parameter @c mu_ and scale @c scale_.
 *  @return a pseudo logistic random variate
 **/
Real Logistic::rand() const
{
  return 0;
}
/* @param x a real value
 *  @return the value of the logistic pdf at @c x
 **/
Real Logistic::pdf( Real const& x) const
{
  return 0;
}
/* @return Give the value of the log-pdf at x.
 *  @param x a real value
 **/
Real Logistic::lpdf( Real const& x) const
{
  return 0;
}
/* @brief Compute the cumulative distribution function at t of
 *  the standard logistic distribution.
 *
 *  The cumulative distribution function of the logistic distribution is
 *  also a scaled version of the Hyperbolic function.
 *  \f[
 *   F(t; \mu, s) = \frac{1}{1+e^{-\frac{t-\mu}{s}}}
 *   = \frac{1}{2} + \frac{1}{2} \;\operatorname{tanh}\!\left(\frac{t-\mu}{2s}\right).
 *   \f]
 * 
 *  @param t a real value
 *  @return the cumulative distribution function value at t
 **/
Real Logistic::cdf( Real const& t) const
{
  return 0;
}

/* @brief Compute the inverse cumulative distribution function at p
 *  of the standard logistic distribution.
 *
 *  The inverse cumulative distribution function (quantile function) of the
 *  logistic distribution is a generalization of the logit function.
 *  It is defined as follows:
 *  \f[
 *      Q(p;\mu,s) = \mu + s\,\ln\left(\frac{p}{1-p}\right).
 *  \f]
 *  @param p a probability number.
 *  @return the inverse cumulative distribution function value at p.
 **/
Real Logistic::icdf( Real const& p) const
{
  return 0;
}

/* @brief Generate a pseudo Logistic random variate.
 *
 *  Generate a pseudo Logistic random variate with location @c mu and
 *  scale @c scale parameters.
 *  @param mu mean of the Logistic distribution
 *  @param scale scale of the Logistic distribution
 *  @return a pseudo logistic random variate, centered in @c mu and with
 *  scale @c scale
 **/
Real Logistic::rand( Real const& mu, Real const& scale)
{
  return 0;
}
/* @param x a real value
 *  @param mu mean of the logistic law
 *  @param scale scale of the logistic law
 *  @return the value of the logistic pdf at @c x
 **/
Real Logistic::pdf( Real const& x, Real const& mu, Real const& scale)
{
  return 0;
}
/* @return Give the value of the log-pdf at x.
 *  @param x a real value
 *  @param mu mean of the logistic law
 *  @param scale scale of the logistic law
 **/
Real Logistic::lpdf( Real const& x, Real const& mu, Real const& scale)
{
  return 0;
}
/* @brief Compute the cumulative distribution function at t of
 *  the standard logistic distribution.
 *
 *  @param t a real value
 *  @return the cumulative distribution function value at t
 **/
Real Logistic::cdf( Real const& t, Real const& mu, Real const& scale)
{
  return 0;
}

/* @brief Compute the inverse cumulative distribution function at p
 *  of the standard logistic distribution.
 *
 *  @param p a probability number.
 *  @return the inverse cumulative distribution function value at p.
 **/
Real Logistic::icdf( Real const& p, Real const& mu, Real const& scale)
{
  return 0;
}

#endif

} // namespace Law

} // namespace STK

