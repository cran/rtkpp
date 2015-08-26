/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2015  Serge Iovleff
    
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
 * Purpose:  raw mathematical functions
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Funct_raw.h
 *  @brief In this file we declare raw the functions.
 * 
 *  raw functions are generic functions that can be used in various part
 *  of the STatistiK project. No test is done about the arguments.
 **/

#ifndef STK_FUNCT_RAW_H
#define STK_FUNCT_RAW_H

#include "STKernel/include/STK_Integer.h"
#include "STKernel/include/STK_Real.h"

namespace STK
{

namespace Funct
{

/** @ingroup Analysis
 *  @brief Compute the beta density function.
 *
 *  Compute the function (beta pdf)
 *  \f[ B(a,b,x) = \frac{ x^{a-1} (1-x)^{b-1}}{B(a,b)}. \f]
 *  @param a,b,x parameters of the beta density
 */
Real beta_pdf_raw(Real const& x, Real const& a, Real const& b);

/** @ingroup Analysis
 *  @brief Compute the generalized binomial probability mass function.
 *  Compute the function
 *  \f[ B(n,p,x) = p^{x} (1-p)^{n-x} \frac{\Gamma(n+1)}{\Gamma(x+1)\Gamma(n-x+1)} \f].
 *  The function assume that there is no NA or infinite values and
 *  that <em>0<=p<=1, n>=0, 0<=x<=n</em>
 *  @param x,n,p parameters of the binomial density
 */
Real binomial_pdf_raw(Real const& x, Real const& n, Real const& p);

/** @ingroup Analysis
 *  @brief Compute the binomial probability mass function.
 *
 *  Compute the function
 *  \f[ B(n,p,x) = p^{x} (1-p)^{n-x} \binom{n}{x} \f].
 *  The function assume that there is no NA or infinite values and
 *  that <em>0<=p<=1, n>=0, 0<=x<=n</em>
 *  @param x,n,p parameters of the binomial density
 */
Real binomial_pdf_raw(int x, int n, Real const& p);

/** @ingroup Analysis
 *  @brief Compute the generalized binomial log-probability mass function.
 *  Compute the function
 *  \f[ B(n,p,x) = x\log(p) (n-x)\log(1-p)
 *      \log\left(\frac{\Gamma(n+1)}{\Gamma(x+1)\Gamma(n-x+1)}\right)
 *  \f].
 *  The function assume that there is no NA or infinite values and
 *  that <em>0<=p<=1, n>=0, 0<=x<=n</em>
 *  @param x,n,p parameters of the binomial density
 */
Real binomial_lpdf_raw(Real const& x, Real const& n, Real const& p);

/** @ingroup Analysis
 *  @brief Compute the binomial log-probability mass function.
 *
 *  Compute the function
 *  \f[ B(n,p,x) = x\log(p) (n-x)\log(1-p)
 *      \log\left(\frac{\Gamma(n+1)}{\Gamma(x+1)\Gamma(n-x+1)}\right)
 *  \f].
 *  The function assume that there is no NA or infinite values and
 *  that <em>0<=p<=1, n>=0, 0<=x<=n</em>
 *  @param x,n,p parameters of the binomial density
 */
Real binomial_lpdf_raw(int x, int n, Real const& p);

/** @ingroup Analysis
 *  @brief Compute the Poisson density.
 *  Compute the function:
 *  \f[
 *    p(x, \lambda) = e^{-\lambda} \frac{\lambda^x}{\Gamma(x+1)}
 *  \f]
 *  with good accuracy using the partial deviance.
 *  This is the version for x Real.
 *
 *  @see http://www.herine.net/stat/software/dbinom.html
 *
 *  @param x value to evaluate the function
 *  @param lambda value of the parameter
 */
Real poisson_pdf_raw(Real const& x, Real const& lambda);

/** @ingroup Analysis
 *  @brief Compute the poisson density with integer value.
 *  Compute the function:
 *  \f[
 *    p(x, \lambda) = e^{-\lambda} \frac{\lambda^x}{x!}
 *  \f]
 *  with good accuracy using the partial deviance.
 *  This is the version for x integer.
 *
 *  @see http://www.herine.net/stat/software/dbinom.html
 *
 *  @param x value to evaluate the function
 *  @param lambda value of the parameter
 */
Real poisson_pdf_raw(int const& x, Real const& lambda);

/** @ingroup Analysis
 *  @brief Compute the log-poisson density.
 *  Compute the function:
 *  \f[
 *    p(x, \lambda) = -\lambda + x \log(\lambda) -\log(\Gamma(x+1))
 *  \f]
 *  with good accuracy using the partial deviance.
 *  This is the version for x Real.
 *
 *  @see http://www.herine.net/stat/software/dbinom.html
 *
 *  @param x value to evaluate the function
 *  @param lambda value of the parameter
 */
Real poisson_lpdf_raw(Real const& x, Real const& lambda);

/** @ingroup Analysis
 *  @brief Compute the log-poisson density with integer value.
 *  Compute the function:
 *  \f[
 *    p(x, \lambda) =  -\lambda + x \log(\lambda) -\log(x!)
 *  \f]
 *  with good accuracy using the partial deviance.
 *  This is the version for x integer.
 *
 *  @see http://www.herine.net/stat/software/dbinom.html
 *
 *  @param x value to evaluate the function
 *  @param lambda value of the parameter
 */
Real poisson_lpdf_raw(int const& x, Real const& lambda);

/** @ingroup Analysis
 *  @brief Compute the psi function psi(x).
 */
Real psi_raw(Real x);

/** @ingroup Analysis
 *  @brief Compute the error function erf(a)
 */
Real erf_raw(Real const& a);

/** @ingroup Analysis
 *  @brief Compute the complementary error function erfc(a)
 */
Real erfc_raw(Real const& a);

/** @ingroup Analysis
 *  @brief Compute the cumulative distribution function of
 *  the normal density
 */
Real normal_cdf_raw(Real const& x);

/** @ingroup Analysis
 *  @brief compute the probability distribution function of
 *  the normal density
 */
Real normal_pdf_raw(Real const& x);

} // namespace Funct

} // namespace STK

#endif // STK_FUNCT_RAW_H
