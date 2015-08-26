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
 * Project: STatistiK
 * Purpose: Implementation of the Cauchy Distribution
 * Author:  Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Cauchy.cpp
 *  @brief In this file we implement the Cauchy distribution.
 **/

#ifndef IS_RTKPP_LIB
#include "../include/STK_Law_Cauchy.h"
#include "../include/STK_Law_Util.h"
#endif


namespace STK
{

namespace Law
{

#ifndef IS_RTKPP_LIB

/*  Generate a pseudo Cauchy random variate. */
Real Cauchy::rand() const
{ return mu_ + scale_ * Real(std::tan( double(Const::_PI_ * generator.randUnif())));}

/*
 *  Generate a pseudo Cauchy random variate with the specified parameters.
 *  (static)
 */
Real Cauchy::rand( Real const& mu, Real const& scale)
{
#ifdef STK_DEBUG
  // check parameters
  if (  !Arithmetic<Real>::isFinite(mu) || !Arithmetic<Real>::isFinite(scale) || scale <= 0)
    STKDOMAIN_ERROR_2ARG(Cauchy::Cauchy,mu, scale,argument error);
#endif
  return mu + scale * Real(std::tan( double(Const::_PI_ * generator.randUnif()) ));
}

/*  Give the value of the pdf at x.*/
Real Cauchy::pdf( Real const& x) const
{
  // check NA value
  if (isNA(x)) return Arithmetic<Real>::NA();
  // trivial case
  if (Arithmetic<Real>::isInfinite(x)) return 0.0;

  // general case
  Real y = (x - mu_) / scale_;
  return 1. / (Const::_PI_ * scale_ * (1. + y * y));
}

/*  Give the value of the pdf at x.*/
Real Cauchy::pdf( Real const& x, Real const& mu, Real const& scale)
{
#ifdef STK_DEBUG
  // check parameters
  if (  !Arithmetic<Real>::isFinite(mu) || !Arithmetic<Real>::isFinite(scale) || scale <= 0)
    STKDOMAIN_ERROR_2ARG(Cauchy::pdf,mu, scale,argument error);
#endif
  // check NA value
  if (isNA(x)) return Arithmetic<Real>::NA();
  // trivial case
  if (Arithmetic<Real>::isInfinite(x)) return 0.0;

  // general case
  Real y = (x - mu) / scale;
  return 1. / (Const::_PI_ * scale * (1. + y * y));
}

/* Give the value of the log-pdf at x. */
Real Cauchy::lpdf( Real const& x) const
{
  // check NA value
  if (isNA(x)) return Arithmetic<Real>::NA();
  // trivial case
  if (Arithmetic<Real>::isInfinite(x)) return -Arithmetic<Real>::infinity();

  // general case
  Real y = (x - mu_) / scale_;
  return -Real(std::log( double(Const::_PI_ * scale_ * (1. + y * y)) ));
}

/* Give the value of the log-pdf at x. */
Real Cauchy::lpdf( Real const& x, Real const& mu, Real const& scale)
{
#ifdef STK_DEBUG
  // check parameters
  if (  !Arithmetic<Real>::isFinite(mu) || !Arithmetic<Real>::isFinite(scale) || scale <= 0)
    STKDOMAIN_ERROR_2ARG(Cauchy::lpdf,mu, scale,argument error);
#endif
  // check NA value
  if (isNA(x)) return Arithmetic<Real>::NA();
  // trivial case
  if (Arithmetic<Real>::isInfinite(x))
    return -Arithmetic<Real>::infinity();

  // general case
  Real y = (x - mu) / scale;
  return - Real(std::log( double(Const::_PI_ * scale * (1. + y * y)) ));
}

/* The cumulative distribution function at t.
 */
Real Cauchy::cdf( Real const& t) const
{
  // check NA value
  if (isNA(t)) return Arithmetic<Real>::NA();
  // check parameter
  if (Arithmetic<Real>::isInfinite(t))
   return (t < 0.) ? 0.0 : 1.0;

  /* http://www.faqs.org/faqs/fr/maths/maths-faq-3/
   * arctan on [0, 1[:
   * if x<0 atan(x)= -atan(-x)
   * elseif x>1 atan(x)= Pi/2-atan(1/x).
   */
  Real td = (t - mu_)/scale_;
  if (std::abs(td) > 1)
  {
    Real y = Real(std::atan( 1./double(td))) / Const::_PI_;
    return (td > 0) ? (1. - y) : (-y);
  }
  return 0.5 + Real(std::atan( double(td))) / Const::_PI_;
}

/*
 * The inverse cumulative distribution function at p.
 */
Real Cauchy::icdf( Real const& p) const
{
  // check NA value
  if (isNA(p)) return Arithmetic<Real>::NA();
  // check parameter
  if ((p > 1.) || (p < 0.))
    STKDOMAIN_ERROR_1ARG(Cauchy::icdf,p,p must be a probability);
  // trivial cases
  if (p == 0.)  return -Arithmetic<Real>::infinity();
  if (p == 1.)  return  Arithmetic<Real>::infinity();

  // general case
  // tan(pi * (p - 1/2)) = -cot(pi * p) = -1/tan(pi * p)
  return mu_ - scale_ / Real(std::tan( double(Const::_PI_ * p) ));
}

/* The cumulative distribution function at t.
 */
Real Cauchy::cdf( Real const& t, Real const& mu, Real const& scale)
{
  // check NA value
  if (isNA(t)) return Arithmetic<Real>::NA();
  // check parameter
  if (Arithmetic<Real>::isInfinite(t))
   return (t < 0.) ? 0.0 : 1.0;

  /* http://www.faqs.org/faqs/fr/maths/maths-faq-3/
   * arctan on [0, 1[:
   * if x<0 atan(x)= -atan(-x)
   * elseif x>1 atan(x)= Pi/2-atan(1/x).
   */
  Real td = (t - mu)/scale;
  if (std::abs(td) > 1)
  {
    Real y = Real(std::atan( 1./double(td))) / Const::_PI_;
    return (td > 0) ? (1. - y) : (-y);
  }
  return 0.5 + Real(std::atan( double(td))) / Const::_PI_;
}
    
/*
 * The inverse cumulative distribution function at p.
 */
Real Cauchy::icdf( Real const& p, Real const& mu, Real const& scale)
{
  // check NA value
  if (isNA(p)) return Arithmetic<Real>::NA();
  // check parameter
  if ((p > 1.) || (p < 0.))
    STKDOMAIN_ERROR_1ARG(Cauchy::icdf,p,p must be a probability);
 // trivial cases
 if (p == 0.)  return -Arithmetic<Real>::infinity();
 if (p == 1.)  return  Arithmetic<Real>::infinity();

  // general case
  // tan(pi * (p - 1/2)) = -cot(pi * p) = -1/tan(pi * p) 
  return mu - scale / Real(std::tan( double(Const::_PI_ * p) ));
}

#endif

} // namespace Law

} // namespace STK

