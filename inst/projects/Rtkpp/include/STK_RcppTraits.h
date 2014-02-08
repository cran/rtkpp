/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

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
 * created on: 28 juil. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file RcppTraits.h
 *  @brief In this file we (re)implement the Traits class of the Rcpp package.
 **/


#ifndef RCPPTRAITS_H
#define RCPPTRAITS_H

#ifndef R_INTERNALS_H_

#define NILSXP       0    /* nil = NULL */
#define SYMSXP       1    /* symbols */
#define LISTSXP      2    /* lists of dotted pairs */
#define CLOSXP       3    /* closures */
#define ENVSXP       4    /* environments */
#define PROMSXP      5    /* promises: [un]evaluated closure arguments */
#define LANGSXP      6    /* language constructs (special lists) */
#define SPECIALSXP   7    /* special forms */
#define BUILTINSXP   8    /* builtin non-special forms */
#define CHARSXP      9    /* "scalar" string type (internal only)*/
#define LGLSXP      10    /* logical vectors */
/* 11 and 12 were factors and ordered factors in the 1990s */
#define INTSXP      13    /* integer vectors */
#define REALSXP     14    /* real variables */
#define CPLXSXP     15    /* complex variables */
#define STRSXP      16    /* string vectors */
#define DOTSXP      17    /* dot-dot-dot object */
#define ANYSXP      18    /* make "any" args work.
           Used in specifying types for symbol
           registration to mean anything is okay  */
#define VECSXP      19    /* generic vectors */
#define EXPRSXP     20    /* expressions vectors */
#define BCODESXP    21    /* byte code */
#define EXTPTRSXP   22    /* external pointer */
#define WEAKREFSXP  23    /* weak reference */
#define RAWSXP      24    /* raw bytes */
#define S4SXP       25    /* S4, non-vector */

/* used for detecting PROTECT issues in memory.c */
#define NEWSXP      30    /* fresh node creaed in new page */
#define FREESXP     31    /* node released by GC */

#define FUNSXP      99    /* Closure or Builtin or Special */

#endif

namespace STK
{

namespace hidden
{
/**
 * template that returns the SEXP type that is appropriate for
 * the type T, this is allways VECSXP (lists) unless it is specialized
 */
template <typename T> struct r_sexptype_traits{ enum{ rtype = VECSXP }; } ;
template<> struct r_sexptype_traits<int>{ enum{ rtype = INTSXP } ; } ;
template<> struct r_sexptype_traits<const int>{ enum{ rtype = INTSXP } ; } ;
template<> struct r_sexptype_traits<double>{ enum{ rtype = REALSXP } ; } ;
template<> struct r_sexptype_traits<const double>{ enum{ rtype = REALSXP } ; } ;
template<> struct r_sexptype_traits<bool>{ enum{ rtype = LGLSXP } ; } ;
template<> struct r_sexptype_traits<std::string>{ enum{ rtype = STRSXP } ; } ;
//template<> struct r_sexptype_traits<Rcomplex>{ enum{ rtype = CPLXSXP } ; } ;
//template<> struct r_sexptype_traits<Rbyte>{ enum{ rtype = RAWSXP } ; } ;


template<> struct r_sexptype_traits<unsigned int>{ enum{ rtype = REALSXP } ; } ;
template<> struct r_sexptype_traits<float>{ enum{ rtype = REALSXP } ; } ;

/* long are represented as numeric vectors which allows more precision
   to be preserved than int */
template<> struct r_sexptype_traits<long>{ enum{ rtype = REALSXP } ; } ;
template<> struct r_sexptype_traits<unsigned long>{ enum{ rtype = REALSXP } ; } ;

/* long double are represented as numeric vectors because we cannot do better in R
   some precision will be lost though
*/
template<> struct r_sexptype_traits<long double>{ enum{ rtype = REALSXP } ; } ;

/* short are represented as integer vector in R
*/
template<> struct r_sexptype_traits<short>{ enum{ rtype = INTSXP } ; } ;
template<> struct r_sexptype_traits<unsigned short>{ enum{ rtype = INTSXP } ; } ;

/* std::complex */
//template<> struct r_sexptype_traits< std::complex<double> >{ enum{ rtype = CPLXSXP } ; } ;
//template<> struct r_sexptype_traits< std::complex<float> >{ enum{ rtype = CPLXSXP } ; } ;

/** @ingroup hidden
 * Indicates the storage type associated with a SEXP type
 * for example: RcppTraits<int>::rtype is a INTSXP
 * The default is VECSXP.
 */
template<typename Type> struct RcppTraits
{
 enum
  { Rtype_ = r_sexptype_traits<Type>::rtype };
};


} // namespace hidden

} // namespace STK

#ifndef R_INTERNALS_H_

#undef NILSXP
#undef SYMSXP
#undef LISTSXP
#undef CLOSXP
#undef ENVSXP
#undef PROMSXP
#undef LANGSXP
#undef SPECIALSXP
#undef BUILTINSXP
#undef CHARSXP
#undef LGLSXP

#undef INTSXP
#undef REALSXP
#undef CPLXSXP
#undef STRSXP
#undef DOTSXP
#undef ANYSXP

#undef VECSXP
#undef EXPRSXP
#undef BCODESXP
#undef EXTPTRSXP
#undef WEAKREFSXP
#undef RAWSXP
#undef S4SXP

#undef NEWSXP
#undef FREESXP

#undef FUNSXP

#endif

#endif /* RCPPTRAITS_H */
