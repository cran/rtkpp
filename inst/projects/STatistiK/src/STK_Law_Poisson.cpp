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
 * created on: 8 dec. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/* @file STK_Law_Poisson.cpp
 *  @brief In this file we implement the Poisson distribution.
 **/



#ifndef IS_RTKPP_LIB
#include "../include/STK_Law_Poisson.h"
#include "../include/STK_Law_Util.h"
#include <Analysis/include/STK_Funct_gammaRatio.h>
#include <Analysis/include/STK_Funct_gamma.h>
#include <Analysis/include/STK_Funct_raw.h>
#include <Analysis/include/STK_Funct_Util.h>
#include <Analysis/include/STK_Const_Math.h>
#endif


namespace STK
{

namespace Law
{

#ifndef IS_RTKPP_LIB

/* The inverse cumulative distribution function at p.*/
static Real gauss_icdf_fast(Real const& p)
{
 // trivial cases
 if (p == 0.5) return  0.;

 const Real a[6] =
 {
   -3.969683028665376e+01
 ,  2.209460984245205e+02
 , -2.759285104469687e+02
 ,  1.383577518672690e+02
 , -3.066479806614716e+01
 ,  2.506628277459239e+00
 };
 const Real b[5] =
 {
   -5.447609879822406e+01
 ,  1.615858368580409e+02
 , -1.556989798598866e+02
 ,  6.680131188771972e+01
 , -1.328068155288572e+01
 };
 const Real c[6] =
 {
   -7.784894002430293e-03
 , -3.223964580411365e-01
 , -2.400758277161838e+00
 , -2.549732539343734e+00
 ,  4.374664141464968e+00
 ,  2.938163982698783e+00
 };
 const Real d[4] =
 {
   7.784695709041462e-03
 , 3.224671290700398e-01
 , 2.445134137142996e+00
 , 3.754408661907416e+00
 };

 Real t, u, q = std::min(p, 1-p);
 if (q > 0.02425)
 {
  /* Rational approximation for central region. */
  u = q-0.5;
  t = u*u;
  u = u*(((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4])*t+a[5])
       /(((((b[0]*t+b[1])*t+b[2])*t+b[3])*t+b[4])*t+1.);
 }
 else
 {
  /* Rational approximation for tail region. */
  t = sqrt(-2.*log(q));
  u = (((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5])
     /((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1);
 }
 return (p > 0.5 ? - u : u);
}

/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param p a probability number
 **/
static int icdf_poisson_raw(Real const& p, Real const& lambda)
{
  // check values 0 and 1
  Real q=1, el = std::exp(-lambda), s = p/el;
  if (s<q) return 0;
  q += lambda;
  if (s<q) return 1;
  q += lambda*lambda/2;
  if (s<q) return 2;
  // otherwise search
  q = gauss_icdf_fast(p);
  Real sqlambda = std::sqrt(lambda);
  int k = round(lambda + q * sqlambda + (q-1.)*(q+1.)*(1.-q/(12.*sqlambda))/6);
  if (Funct::gammaRatioQ(k, lambda) >= p)
  { /* decreasing search */
    while(1)
    {
      if ( Funct::gammaRatioQ(k-1, lambda) < p) return k;
      k--;
    }
  }
  else
  { /* increasing search */
    while(1)
    {
      k++;
      if( Funct::gammaRatioQ(k, lambda) >= p) return k;
    }
  }
  return k; // avoid warning at compilation
}


/* @return a Poisson random variate . */
int Poisson::rand() const
{
 return icdf_poisson_raw(Law::generator.randUnif(), lambda_);
}
/* @brief compute the probability distribution function.
 *  Give the value of the pdf at the point x.
 *  @param x a binary value
 *  @return the value of the pdf
 **/
Real Poisson::pdf(int const& x) const
{
  // check trivial values
  if (x<0) return( 0. );
  // if lambda is 0, we have P(X=0) = 1
  if (lambda_==0.) return( (x==0) ? 1. : 0. );
  // special value
  if (x==0) return( Real(std::exp(-lambda_)) );
  // stirling approximation and deviance
  return( std::exp(-Funct::gammaLnStirlingError((Real)x)-Funct::dev0((Real)x, lambda_))
          /(Const::_SQRT2PI_*std::sqrt(x))
        );
}

/* @brief compute the log probability distribution function.
 *  Give the value of the log-pdf at the point x.
 *  @param x an integer value
 *  @return the value of the log-pdf
 **/
Real Poisson::lpdf(int const& x) const
{
  // check trivial values
  if (x<0) return( -Arithmetic<Real>::infinity() );
  // if lambda is 0, we have P(X=0) = 1
  if (lambda_==0.) return( (x==0) ? 0 : -Arithmetic<Real>::infinity() );
  // special value
  if (x==0) return( -lambda_ );
  // stirling approximation and deviance
  return( -Funct::gammaLnStirlingError(x)-Funct::dev0(x, lambda_)
          -Const::_LNSQRT2PI_-std::log((Real)x)/2.
        );
}
/* @brief compute the cumulative distribution function
 *  Give the probability that a Poisson random variate is less or equal
 *  to t.
 *  @param t a real value
 *  @return the value of the cdf
 **/
Real Poisson::cdf(Real const& t) const
{
  // check trivial values
  if (t < 0)    return 0;
  if (lambda_ == 0.) return( 1. );
  // use gamma ratio function
  return Funct::gammaRatioQ(std::floor(t) + 1, lambda_);
}
/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param p a probability number
 **/
int Poisson::icdf(Real const& p) const
{
  // check trivial values
  if (p == 0.) return 0;
  if (p == 1.) return Arithmetic<int>::max(); // infty does not exists for int
  // check values 0 and 1
  Real el = std::exp(-lambda_);
  if (p<el) return 0;
  if (p<(1+lambda_)*el) return 1;
  // othewise search
  Real q = gauss_icdf_fast(p), sqlambda = std::sqrt(lambda_);
  int k = std::max(0., round(lambda_ + q * sqlambda + (q-1.)*(q+1.)*(1.-q/(12.*sqlambda))/6));
  if (cdf(k) >= p)
  { /* decreasing search */
    while(1)
    {
      if ( cdf(k - 1) < p) return k;
      k--;
    }
  }
  else
  { /* increasing search */
    while(1)
    {
      k++;
      if( cdf(k) >= p) return k;
    }
  }
  return k; // avoid warning at compilation
}

/* @param lambda the mean
 *  @return a int random variate.
 **/
int Poisson::rand(Real const& lambda)
{ return icdf_poisson_raw(Law::generator.randUnif(), lambda);}
/* @brief compute the probability distribution function
 *  Give the value of the pdf at the point x.
 *  @param x a binary value
 *  @param lambda the mean
 *  @return the value of the pdf
 **/
Real Poisson::pdf(int const& x, Real const& lambda)
{
  // check trivial values
  if (x<0) return( 0. );
  // if lambda is 0, we have P(X=0) = 1
  if (lambda==0.) return( (x==0) ? 1. : 0. );
  // special value
  if (x==0) return( Real(std::exp(-lambda)) );
  // stirling approximation and deviance
  return( std::exp(-Funct::gammaLnStirlingError((Real)x)-Funct::dev0((Real)x, lambda))
          /(Const::_SQRT2PI_*std::sqrt(x))
        );
}

/* @brief compute the log probability distribution function
 *  Give the value of the log-pdf at the point x.
 *  @param x a binary value
 *  @param lambda the mean
 *  @return the value of the log-pdf
 **/
Real Poisson::lpdf(int const& x, Real const& lambda)
{
  // check trivial values
  if (x<0) return( -Arithmetic<Real>::infinity() );
  // if lambda is 0, we have P(X=0) = 1
  if (lambda==0.) return( (x==0) ? 0 : -Arithmetic<Real>::infinity() );
  // special value
  if (x==0) return( -lambda );
  // stirling approximation and deviance
  return( -Funct::gammaLnStirlingError((Real)x)-Funct::dev0((Real)x, lambda)
          -Const::_LNSQRT2PI_-std::log((Real)x)/2.
        );
}

/* @brief compute the cumulative distribution function
 *  Give the probability that a Poisson random variate is less or equal
 *  to t.
 *  @param t a real value
 *  @return the value of the cdf
 **/
Real Poisson::cdf(Real const& t, Real const& lambda)
{
  // check trivial values
  if (t < 0) return 0;
  if (lambda == 0.) return( 1. );
  // use gamma ratio function
  return Funct::gammaRatioQ(std::floor(t) + 1, lambda);
}
/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param p a probability number
 **/
int Poisson::icdf(Real const& p, Real const& lambda)
{
  // check trivial values
  if (p == 0.) return 0;
  if (p == 1.) return Arithmetic<int>::max();
  // check values 0 and 1
  Real el = std::exp(-lambda);
  if (p<el) return 0;
  if (p< (1+lambda)*el) return 1;
  // othewise search
  Real q = gauss_icdf_fast(p), sqlambda = std::sqrt(lambda);
  int k = round(lambda + q * sqlambda + (q-1.)*(q+1.)*(1.-q/(12.*sqlambda))/6);
  if (cdf(k, lambda) >= p)
  { /* decreasing search */
    while(1)
    {
      if ( cdf(k - 1,lambda) < p) return k;
      k--;
    }
  }
  else
  { /* increasing search */
    while(1)
    {
      k++;
      if( cdf(k, lambda) >= p) return k;
    }
  }
  return k; // avoid warning at compilation
}

#endif

} // namespace Law

} // namespace STK


