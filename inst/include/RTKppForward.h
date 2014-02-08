/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  rtkpp
 * created on: 22 août 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file RTKppForward.h
 *  @brief In this file we include the header files needed for the integration of stkpp to Rcpp.
 **/


#ifndef RTKPPFORWARD_H
#define RTKPPFORWARD_H

#include <RcppCommon.h>
#include "STKpp.h"

/* forward declarations */
namespace STK
{
template <typename Type> class RcppVector;
template <typename Type> class RcppMatrix;
}


/* Rcpp integration */
namespace Rcpp
{
  /* support for wrap */
  template<typename Derived>
  SEXP wrap(STK::ITContainer<Derived, STK::hidden::Traits<Derived>::structure_> const& obj);

  namespace traits
  {
    /* support for as */
    template<typename Type> class Exporter< STK::RcppVector<Type> >;
    template<typename Type> class Exporter< STK::RcppMatrix<Type> >;
  } // namespace traits

} // namespace Rcpp

#endif /* RTKPPFORWARD_H */
