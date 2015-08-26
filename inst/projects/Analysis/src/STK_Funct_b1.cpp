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
 * Project:  stkpp::Analysis::Funct
 * Purpose:  implementation of the b1 function
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Funct_b1.cpp
 *  @brief In this file we implement the b1 function.
 **/

#include <cmath>

#include "../include/STK_Funct_Util.h"
#include "../include/STK_Const_Math.h"
#include "../include/STK_Funct_gamma.h"

namespace STK
{
  
namespace Funct
{
/* @ingroup Analysis
 *  Compute the function
 *  \f[ B_1(a,b,x) = \frac{ x^{a} (1-x)^{b}}{B(a,b)} \f]
 *  using  the partial deviance \f$ (a+b) * (p*log(x/p)+q*log((1-x)/q)) \f$.
 *  @param a,b,x parameters of the beta density
 *  @param xm1 true if @e x is to be taken as @e 1-x
 **/
Real b1(Real const& a, Real const& b, Real const& x, bool xm1)
{
  if (x == 0) return 0;
  if (x == 1) return 0;
  Real s = a+b, sx = s*x, sy = s*(1.-x);
  return ( std::exp(- Const::_LNSQRT2PI_
                    + 0.5 * ( a<b ? std::log(a) + log1p(-a/s) : std::log(b) + log1p(-b/s))
                    + (gammaLnStirlingError(s)-gammaLnStirlingError(a)-gammaLnStirlingError(b))
                    - (xm1 ? dev0(a, sy)+dev0(b, sx) : dev0(a, sx)+dev0(b, sy))
                   )
        );
}

} // namespace Funct

} // namespace STK
