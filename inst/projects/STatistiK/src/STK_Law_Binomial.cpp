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

/** @file STK_Law_Binomial.cpp
 *  @brief In this file we implement the Binomial distribution.
 **/

#ifndef IS_RTKPP_LIB
#include "../include/STK_Law_Binomial.h"
#include "Analysis/include/STK_Funct_raw.h"
#endif


namespace STK
{

namespace Law
{

#ifndef IS_RTKPP_LIB

int Binomial::rand() const
{
  return 0;
}

Real Binomial::pdf(int const& x) const
{
  // trivial cases
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();
  // compute result
  return Funct::binomial_pdf_raw(x, n_, prob_);
}

Real Binomial::lpdf(int const& x) const
{
  // trivial cases
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();
  // compute result
  return Funct::binomial_lpdf_raw(x, n_, prob_);
}
Real Binomial::cdf(Real const& t) const
{
  // trivial cases
  if (Arithmetic<Real>::isNA(t)) return Arithmetic<Real>::NA();
  return 0.;
}

int Binomial::icdf(Real const& p) const
{
  // trivial cases
  if (Arithmetic<Real>::isNA(p)) return Arithmetic<Real>::NA();
  return 0.;
}

int Binomial::rand(int n, Real const& prob)
{
  return 0;
}

Real Binomial::pdf(int x, int n, Real const& prob)
{
  // trivial cases
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();
  // compute result
  return Funct::binomial_pdf_raw(x, n, prob);
}

Real Binomial::lpdf(int x, int n, Real const& prob)
{
  // trivial cases
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();
  // compute result
  return Funct::binomial_lpdf_raw(x, n, prob);
}
Real Binomial::cdf(Real const& t, int n, Real const& prob)
{
  return 0.;
}

int Binomial::icdf(Real const& p, int n, Real const& prob)
{
  return 0.;
}

#endif

} // namespace Law

} // namespace STK

