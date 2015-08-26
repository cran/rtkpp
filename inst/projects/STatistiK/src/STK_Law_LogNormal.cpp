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
 * Purpose: LogNormal probability distribution.
 * Author:  Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_LogNormal.cpp
 *  @brief In this file we implement the LogNormal probability law class.
 **/


#ifndef IS_RTKPP_LIB
#include "../include/STK_Law_LogNormal.h"
#endif

namespace STK
{

namespace Law
{

#ifndef IS_RTKPP_LIB

/* @brief Generate a pseudo log-normalized LogNormal random variate.
 *
 *  Generate a pseudo log-normalized LogNormal random variate
 *  with location parameter @c mu_ and sigma @c sigma_.
 *  @return a pseudo log-normal random variate
 **/
Real LogNormal::rand() const
{
  return 0;
}
/* @param x a real value
 *  @return the value of the log-normal pdf at @c x
 **/
Real LogNormal::pdf( Real const& x) const
{
  return 0;
}
/* @return Give the value of the log-pdf at x.
 *  @param x a real value
 **/
Real LogNormal::lpdf( Real const& x) const
{
  return 0;
}
/* @brief Compute the cumulative distribution function at t of
 *  the standard log-normal distribution.
 *  @param t a real value
 *  @return the cumulative distribution function value at t
 **/
Real LogNormal::cdf( Real const& t) const
{
  return 0;
}
/* @brief Compute the inverse cumulative distribution function at p
 *  of the standard log-normal distribution.
 *
 *  @param p a probability number.
 *  @return the inverse cumulative distribution function value at p.
 **/
Real LogNormal::icdf( Real const& p) const
{
  return 0;
}

/* @brief Generate a pseudo LogNormal random variate.
 *
 *  Generate a pseudo LogNormal random variate with location @c mu and
 *  sigma @c sigma parameters.
 *  @param mu location of the LogNormal distribution
 *  @param sigma sigma of the LogNormal distribution
 *  @return a pseudo log-normal random variate, centered in @c mu and with
 *  sigma @c sigma
 **/
Real LogNormal::rand( Real const& mu, Real const& sigma)
{
  return 0;
}
/* @param x a real value
 *  @param mu location of the log-normal law
 *  @param sigma sigma of the log-normal law
 *  @return the value of the log-normal pdf at @c x
 **/
Real LogNormal::pdf( Real const& x, Real const& mu, Real const& sigma)
{
  return 0;
}
/* @return Give the value of the log-pdf at x.
 *  @param x a real value
 *  @param mu location of the log-normal law
 *  @param sigma sigma of the log-normal law
 **/
Real LogNormal::lpdf( Real const& x, Real const& mu, Real const& sigma)
{
  return 0;
}
/* @brief Compute the cumulative distribution function at t of
 *  the standard log-normal distribution.
 *  @param t a real value
 *  @param mu, sigma location and sigma of the log-normal law
 *  @return the cumulative distribution function value at t
 **/
Real LogNormal::cdf( Real const& t, Real const& mu, Real const& sigma)
{
  return 0;
}
/* @brief Compute the inverse cumulative distribution function at p
 *  of the standard log-normal distribution.
 *
 *  @param p a probability number.
 *  @param mu, sigma location and sigma of the log-normal law
 *  @return the inverse cumulative distribution function value at p.
 **/
Real LogNormal::icdf( Real const& p, Real const& mu, Real const& sigma)
{
  return 0;
}

#endif


} // namespace Law

} // namespace STK
