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

/** @file STK_Law_HyperGeometric.cpp
 *  @brief In this file we implement the HyperGeometric distribution.
 **/


#ifndef IS_RTKPP_LIB
#include "../include/STK_Law_HyperGeometric.h"
#endif

namespace STK
{

namespace Law
{

#ifndef IS_RTKPP_LIB

/* @return a random hypergeometric variate. */
Integer HyperGeometric::rand() const
{
  return 0;
}
/* @brief compute the probability distribution function (density)
 *  Give the value of the pdf at the point x.
 *  @param x an Integer value
 *  @return the value of the pdf
 **/
Real HyperGeometric::pdf(Integer const& x) const
{
  return 0;
}
/* @brief compute the log probability distribution function
 *  Give the value of the log-pdf at the point x.
 *  @param x an Integer value
 *  @return the value of the log-pdf
 **/
Real HyperGeometric::lpdf(Integer const& x) const
{
  return 0;
}
/* @brief compute the cumulative distribution function
 *  Give the probability that a HyperGeometric random variate is less or equal
 *  to t.
 *  @param t a real value
 *  @return the value of the cdf
 **/
Real HyperGeometric::cdf(Real const& t) const
{
  return 0;
}
/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param p a probability number
 **/
Integer HyperGeometric::icdf(Real const& p) const
{
  return 0;
}

/* @brief random hypergeometric variate generation.
 *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
 *  @return a Integer random variate.
 **/
Integer HyperGeometric::rand( int nbSuccesses, int nbFailures, int nbDraws)
{
  return 0;
}
/* @brief compute the probability distribution function.
 *  Give the value of the pdf at the point x.
 *  @param x an Integer value
 *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
 *  @return the value of the pdf
 **/
Real HyperGeometric::pdf(Integer x, int nbSuccesses, int nbFailures, int nbDraws)
{
  return 0;
}
/* @brief compute the log probability distribution function.
 *  @param x an Integer value
 *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
 *  @return the value of the log-pdf
 **/
Real HyperGeometric::lpdf(Integer x, int nbSuccesses, int nbFailures, int nbDraws)
{
  return 0;
}
/* @brief compute the cumulative distribution function
 *  Give the probability that a HyperGeometric random variate is less or equal
 *  to t.
 *  @param t a real value
 *  @return the value of the cdf
 **/
Real HyperGeometric::cdf(Real const& t, int nbSuccesses, int nbFailures, int nbDraws)
{
  return 0;
}
/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param p a probability number
 **/
Integer HyperGeometric::icdf(Real const& p, int nbSuccesses, int nbFailures, int nbDraws)
{
  return 0;
}

#endif

} // namespace Law

} // namespace STK


