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

/** @file STK_Law_Uniform.h
 *  @brief In this file we implement the (continuous) uniform distribution law.
 **/

#ifndef STK_LAW_UNIFORM_H
#define STK_LAW_UNIFORM_H

#include "STK_Law_IUnivLaw.h"
#include "../include/STK_Law_Util.h"
#include <STKernel/include/STK_Real.h>
#include <Sdk/include/STK_Macros.h>

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief class for the Uniform law distribution.
 * 
 *  In probability theory and statistics, the <em>continuous uniform distribution</em>
 *  or rectangular distribution is a family of symmetric probability distributions
 *  such that for each member of the family, all intervals of the same length on
 *  the distribution's support are equally probable. The support is defined by
 *  the two parameters, @e a and @e b, which are its minimum and maximum values.
 *
 *  The probability density function of the continuous uniform distribution is:
 *  \f[
 *    f(x; a, b) = \frac{1}{b-a} 1_{ a \leq x \leq b}.
 *  \f]
**/
class Uniform : public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** constructor.
     *  @param a,b the lower and upper bounds
     **/
    inline Uniform( Real const& a =0., Real const& b =1.)
                  : Base(_T("Uniform")), a_(a), b_(b), range_(b_ - a_)
    {
      if (range_ <= 0.)
        STKINVALIDARGUMENT_ERROR_2ARG(Uniform::Uniform, a_, b_,invalid parameters);
    }
    /** copy constructor.
     *  @param law the law to copy
     **/
    inline Uniform( Uniform const& law): Base(law), a_(law.a_), b_(law.b_), range_(law.range_)
    {};
    /** destructor. */
	  inline virtual ~Uniform() {}
    /** @return the lower bound */
    inline Real const& a() const { return a_;}
    /** @return the upper bound */
    inline Real const& b() const { return b_;}
    /** @return the value b-a */
    inline Real const& range() const { return range_;}
    /** @param a set the lower bound */
    inline void setA(Real const& a) { a_ =a;}
    /** @param b set the upper bound */
    inline void setB(Real const& b){ b_ =b;}

    /** Generate a pseudo Uniform random variate. */
    virtual Real rand() const;
    /** Give the value of the pdf at x.
     *  @param x a real value
     **/
    virtual Real pdf( Real const& x) const;
    /** Give the value of the log-pdf at x.
     *  @param x a real value
     **/
    virtual Real lpdf( Real const& x) const;
    /** The cumulative distribution function is
     * \f[
     *  F(t; a,b)= \frac{t - a}{b-a}
     * \f]
     *  @param t a real value
     **/
    virtual Real cdf( Real const& t) const;
    /** The inverse cumulative distribution function is
     * \f[
     * F^{-1}(p; \lambda) = p (b-a) + a.
     * \f]
     *  @param p a probability
     **/
    virtual Real icdf( Real const& p) const;

    /** Generate a pseudo Uniform random variate.
     *  @param a,b the lower and upper bounds
     **/
    static Real rand( Real const& a, Real const& b);
    /** Give the value of the pdf at x.
     *  @param x a real value
     *  @param a,b the lower and upper bounds
     **/
    static Real pdf( Real const& x, Real const& a, Real const& b);
    /** Give the value of the log-pdf at x.
     *  @param p a probablility
     *  @param a,b the lower and upper bounds
     **/
    static Real lpdf( Real const& p, Real const& a, Real const& b);
    /** Give the value of the cdf at t.
     *  @param t a real value
     *  @param a,b the lower and upper bounds
     **/
    static Real cdf( Real const& t, Real const& a, Real const& b);
    /** Give the value of the quantile at @e p.
     *  @param p a probability
     *  @param a,b the lower and upper bounds
     **/
    static Real icdf( Real const& p, Real const& a, Real const& b);

  protected:
    /** The lower bound. */
    Real a_;
    /** The upper bound. */
    Real b_;

  private:
    Real range_;
};

/* Generate a pseudo Uniform random variate. */
inline Real Uniform::rand() const
{
  return ((range_ <= 1.) || (a_ <=1.)) ? a_ + range_ * generator.randUnif()
                                       : a_ * (1. + generator.randUnif()*range_/a_) ;
}
/* Give the value of the pdf at x.
 *  @param x a real value
 **/
inline Real Uniform::pdf( Real const& x) const
{
  if (!Arithmetic<Real>::isFinite(x) ) return x;
  if ((x < a_)||(x > b_)) return 0.;
  return 1./range_;
}
/* Give the value of the log-pdf at x.
 *  @param x a real value
 **/
inline Real Uniform::lpdf( Real const& x) const
{
  if (!Arithmetic<Real>::isFinite(x) ) return x;
  if ((x < a_)||(x > b_)) return -Arithmetic<Real>::infinity();
  return -std::log(range_);
}
/* The cumulative distribution function is
 * \f[
 *  F(t; a,b)= \frac{t - a}{b-a}
 * \f]
 *  @param t a real value
 **/
inline Real Uniform::cdf( Real const& t) const
{
  if (!Arithmetic<Real>::isFinite(t) ) return t;
  if (t <= a_) return 0.;
  if (t >= b_) return 1.;
  return (b_ - t)/range_;
}

/* The inverse cumulative distribution function is
 * \f[
 * F^{-1}(p; \lambda) = p (b-a) + a.
 * \f]
 *  @param p a probability
 **/
inline Real Uniform::icdf( Real const& p) const
{
  // check parameter
  if ((p > 1.) || (p < 0.))
    STKDOMAIN_ERROR_1ARG(Exponential::icdf,p,invalid argument);

  if (!Arithmetic<Real>::isFinite(p) ) return p;
  if (p == 1.) return b_;
  if (p == 0.) return a_;
  return a_ + p * range_;
}

/* Generate a pseudo Uniform random variate with the specified
 *  parameter.
 *  @param scale the scale of the distribution
 **/
inline Real Uniform::rand( Real const& a, Real const& b)
{
  return( (b-a <= 1.) ? a + (b-a) * generator.randUnif()
                      : a + generator.randDblExc(b-a));
}
/* Give the value of the pdf at x.
 *  @param x a real value
 *  @param scale the scale of the distribution
 **/
inline Real Uniform::pdf( Real const& x, Real const& a, Real const& b)
{
  if (!Arithmetic<Real>::isFinite(x) ) return x;
  if ((x < a)||(x > b)) return 0.;
  return 1./(b-a);
}
/* Give the value of the log-pdf at x.
 *  @param x a real value
 *  @param scale the scale of the distribution
 **/
inline Real Uniform::lpdf( Real const& x, Real const& a, Real const& b)
{
  if (!Arithmetic<Real>::isFinite(x) ) return x;
  if ((x < a)||(x > b)) return -Arithmetic<Real>::infinity();
  return -std::log(b-a);
}

inline Real Uniform::cdf(const Real& t, const Real& a, const Real& b)
{ return (b - t)/(b-a);}


inline Real Uniform::icdf(const Real& p, const Real& a, const Real& b)
{ return std::max(a,std::min((1.-p) * a + p * b, b));}

} // namespace Law

} // namespace STK

#endif /*STK_LAW_UNIFORM_H*/
