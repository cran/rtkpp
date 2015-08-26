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
 * Project:  stkpp::Reduct
 * created on: 23 juin 2011
 * Purpose:  enum and other utilities for the Reduct project.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Reduct_Util.h
 *  @brief In this file we define utilities enum and functions for the Reduct
 *  project.
 **/

#ifndef STK_REDUCT_UTIL_H_
#define STK_REDUCT_UTIL_H_

#include "STKernel/include/STK_String.h"

namespace STK
{

namespace Reduct
{

/** @ingroup Reduct
 *  dimension reduction we can apply to the original data set.
 **/
enum TypeReduction
{
  /** unknown reduction*/
  unknown_reduction_ =0
  /** local projected variance */
  , localVariance_
  /** total projected variance (pca)*/
  , totalVariance_
  /** multidimensional scaling */
  , mds_
};

/** @ingroup Reduct
 *  convert a String to a TypeReduction.
 *  @param type the String we want to convert
 *  @return the TypeReduction represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_ type is returned.
 **/
TypeReduction stringToTypeReduction( String const& type);

/** @ingroup Reduct
 *  convert a TypeReduction to a String.
 *  @param type the type of reduction we want to convert
 *  @return the string associated to this type.
 **/
String typeReductionToString( TypeReduction const& type);

/** Type of proximity graph to used in order to compute the local variance:
 * - prim_ the minimal spanning tree
 * - distance_ the first neighbors
 * - unknown_ unknown type of graph
 */
enum TypeGraph { unknown_graph_, prim_, distance_ };
/** convert a String to a TypeGraph.
 *  @param type the type of graph in a string
 *  @return the TypeGraph represented by the String @c type. If the string
 *  does not match any known name, the @c unknown_ type is returned.
 **/
TypeGraph stringToTypeGraph( String const& type);
/** convert a TypeGraph to a String.
 *  @param type the type of graph we want to convert to a string
 *  @return the string associated to this type of graph
 **/
String typeGraphToString( TypeGraph const& type);
/** Constructor. the TypeGraph and the number of neighbors are
 *  given by the user and are not modified.
 *  @param p_data the data set to process
 *  @param type type of proximity graph to build
 *  @param nbNeighbor number of neighbors to use in the proximity graph
 */

} // namespace Reduct

} // namespace STK

#endif /* STK_REDUCT_UTIL_H_ */
