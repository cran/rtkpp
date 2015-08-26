/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2007  Serge Iovleff

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
 * Purpose:  implementation of the beta function
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Funct_beta_raw.cpp
 *  @brief In this file we implement the beta_raw function.
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
 *  @brief Compute the beta density function.
 *
 *  Compute the function (beta pdf)
 *  \f[ B(a,b,x) = \frac{ x^{a-1} (1-x)^{b-1}}{B(a,b)} \f]
 *  using  the partial deviance \f$ (a+b) * (p*log(x/p)+q*log((1-x)/q)) \f$.
 *  @param a,b,x parameters of the beta density
 */
Real beta_pdf_raw(Real const& x, Real const& a, Real const& b)
{
  // trivial cases
  if (x < 0 || x > 1) return(0);
  if (x == 0)
  {
    if(a > 1) return(0);
    if(a < 1) return(Arithmetic<Real>::infinity());
    return(b); // a == 1
  }
  if (x == 1)
  {
    if(b > 1) return(0);
    if(b < 1) return(Arithmetic<Real>::infinity());
    return(a); // b == 1
  }
  // limit cases with x \in (0,1)
  if (a == 0 || b == 0) { return(0);}
  // general case
  Real s = a+b, sx = s*x, sy = s*(1.-x);
  return ( std::exp(- Const::_LNSQRT2PI_
                    + 0.5 *( a<b ? std::log(a) + log1p(-a/s) : std::log(b) + log1p(-b/s))
                    + gammaLnStirlingError(s)-gammaLnStirlingError(a)-gammaLnStirlingError(b)
                    - dev0(a, sx) - dev0(b, sy)
                    - std::log(x) - log1p(-x)
                   )
         );
}


} // namespace Funct

} // namespace STK
