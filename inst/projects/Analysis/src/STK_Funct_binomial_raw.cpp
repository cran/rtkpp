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
 * Project:  Analysis
 * Purpose:  implementation of the binomial function
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Funct_binomial_raw.cpp
 *  @brief In this file we implement the raw functions around the binomial function.
 *  density.
 **/

#include <cmath>

#include "../include/STK_Const_Math.h"
#include "../include/STK_Funct_gamma.h"
#include "../include/STK_Funct_Util.h"
#include "../include/STK_Funct_raw.h"

namespace STK
{

namespace Funct
{
/* @ingroup Analysis
 *  @brief Compute the binomial probability mass function.
 *
 *  Compute the function
 *  \f[ B(n,p,x) = \frac{ p^{x} (1-p)^{n-x}}{\binom{n}{x}} \f].
 *  The function assume that there is no NA or infinite values and
 *  that <em>0<=p<=1, n>=0, 0<=x<=n</em>
 *  @param x,n,p parameters of the binomial density
 */
Real binomial_pdf_raw(Real const& x, Real const& n, Real const& p)
{
  // trivial and limit cases
  if (p == 0) return((x == 0) ? 1 : 0);
  if (p == 1) return((x == n) ? 1 : 0);
  if (x == 0) return((n == 0) ? 1 : std::exp(n*log1p(-p)) );
  if (x == n) return std::pow(p,n);
  // other cases
  Real y = n-x;
  return ( std::exp(- Const::_LNSQRT2PI_ - 0.5 *( std::log(x) + log1p(-x/n))
                    + gammaLnStirlingError(n)-gammaLnStirlingError(x)-gammaLnStirlingError(y)
                    - dev0(x, n*p) - dev0(y, n- n*p)
                   )
        );
}

/* @ingroup Analysis
 *  @brief Compute the binomial probability mass function.
 *
 *  Compute the function
 *  \f[ B(n,p,x) = p^{x} (1-p)^{n-x} \binom{n}{x} \f].
 *  The function assume that there is no NA or infinite values and
 *  that <em>0<=p<=1, n>=0, 0<=x<=n</em>
 *  @param x,n,p parameters of the binomial density
 */
Real binomial_pdf_raw(int x, int n, Real const& p)
{
  // trivial and limit cases
  if (p == 0) return((x == 0) ? 1 : 0);
  if (p == 1) return((x == n) ? 1 : 0);
  if (x == 0) return((n == 0) ? 1 : std::exp(n*log1p(-p)) );
  if (x == n) return std::pow(p,n);
  // other cases
  int y = n-x;
  return ( std::exp(- Const::_LNSQRT2PI_ - 0.5 *( std::log(x) + log1p(-x/(Real)n))
                    + gammaLnStirlingError(n)-gammaLnStirlingError(x)-gammaLnStirlingError(y)
                    - dev0(x, n*p) - dev0(y, n - n*p)
                   )
        );
}
/* @ingroup Analysis
 *  @brief Compute the binomial probability mass function.
 *
 *  Compute the function
 *  \f[ B(n,p,x) = \frac{ p^{x} (1-p)^{n-x}}{\binom{n}{x}} \f].
 *  The function assume that there is no NA or infinite values and
 *  that <em>0<=p<=1, n>=0, 0<=x<=n</em>
 *  @param x,n,p parameters of the binomial density
 */
Real binomial_lpdf_raw(Real const& x, Real const& n, Real const& p)
{
  // trivial and limit cases
  if (p == 0) return((x == 0) ? 0 : -Arithmetic<Real>::infinity());
  if (p == 1) return((x == n) ? 0 : -Arithmetic<Real>::infinity());
  if (x == 0) return((n == 0) ? 0 : n*log1p(-p));
  if (x == n) return(n*std::log(p));
  // other cases
  Real y = n-x;
  return ( - Const::_LNSQRT2PI_ - 0.5 *( std::log(x) + log1p(-x/n))
           + gammaLnStirlingError(n)-gammaLnStirlingError(x)-gammaLnStirlingError(y)
           - dev0(x, n*p) - dev0(y, n- n*p)
        );
}

/* @ingroup Analysis
 *  @brief Compute the binomial probability mass function.
 *
 *  Compute the function
 *  \f[ B(n,p,x) = p^{x} (1-p)^{n-x} \binom{n}{x} \f].
 *  The function assume that there is no NA or infinite values and
 *  that <em>0<=p<=1, n>=0, 0<=x<=n</em>
 *  @param x,n,p parameters of the binomial density
 */
Real binomial_lpdf_raw(int x, int n, Real const& p)
{
  // trivial and limit cases
  if (p == 0) return((x == 0) ? 0 : -Arithmetic<Real>::infinity());
  if (p == 1) return((x == n) ? 0 : -Arithmetic<Real>::infinity());
  if (x == 0) return((n == 0) ? 0 : n*log1p(-p));
  if (x == n) return(n*std::log(p));
  // other cases
  int y = n-x;
  return( - Const::_LNSQRT2PI_ - 0.5 *( std::log(x) + log1p(-x/(Real)n))
          + gammaLnStirlingError(n)-gammaLnStirlingError(x)-gammaLnStirlingError(y)
          - dev0(x, n*p) - dev0(y, n - n*p)
        );
}

} // namespace Funct
} // namespace STK
