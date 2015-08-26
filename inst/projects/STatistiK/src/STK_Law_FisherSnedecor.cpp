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
 * Purpose:  FisherSnedecor probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_FisherSnedecor.cpp
 *  @brief In this file we implement the FisherSnedecor probability distribution.
 **/


#ifndef IS_RTKPP_LIB
#include "../include/STK_Law_FisherSnedecor.h"
#endif


namespace STK
{

namespace Law
{

#ifndef IS_RTKPP_LIB

/* @return a pseudo FisherSnedecor random variate. */
Real FisherSnedecor::rand() const
{
  return 0;
}
/* @return the value of the pdf
 *  @param x a positive real value
 **/
Real FisherSnedecor::pdf(Real const& x) const
{
  return 0;
}
/* @return the value of the log-pdf
 *  @param x a positive real value
 **/
Real FisherSnedecor::lpdf(Real const& x) const
{
  return 0;
}
/* @return the cumulative distribution function
 *  @param t a positive real value
 **/
Real FisherSnedecor::cdf(Real const& t) const
{
  return 0;
}
/* @return the inverse cumulative distribution function
 *  @param p a probability number
 **/
Real FisherSnedecor::icdf(Real const& p) const
{
  return 0;
}

/* @return a pseudo FisherSnedecor random variate with the specified parameters.
 *  @param df1, df2 degree of freedom parameters
 **/
Real FisherSnedecor::rand( int df1, int df2)
{
  return 0;
}
/* @return the value of the pdf
 *  @param x a positive real value
 *  @param df1, df2 degree of freedom parameters
 **/
Real FisherSnedecor::pdf(Real const& x, int df1, int df2)
{
  return 0;
}
/* @return the value of the log-pdf
 *  @param x a positive real value
 *  @param df1, df2 degree of freedom parameters
 **/
Real FisherSnedecor::lpdf(Real const& x, int df1, int df2)
{
  return 0;
}
/* @return the cumulative distribution function
 *  @param t a positive real value
 *  @param df1, df2 degree of freedom parameters
 **/
Real FisherSnedecor::cdf(Real const& t, int df1, int df2)
{
  return 0;
}

/* @return the inverse cumulative distribution function
 *  @param p a probability number
 *  @param df1, df2 degree of freedom parameters
 **/
Real FisherSnedecor::icdf(Real const& p, int df1, int df2)
{
  return 0;
}

#endif

} // namespace Law

} // namespace STK

