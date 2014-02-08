/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

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
 * Project:  stkpp::Algebra
 * created on: 20 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_lapack_SymEigen.cpp
 *  @brief In this file we implement the SymEigen class.
 **/

#include<cmath>

#include "../include/STK_ISymEigen.h"

namespace STK
{

/* @brief Constructor
 *  @param data reference on a symmetric square matrix
 *  @param ref @c true if we overwrite the data set, @c false otherwise
 */
ISymEigen::ISymEigen( CArraySquareXX const& data, bool ref)
         : Base()
         , norm_(0.), rank_(0), det_(0.)
         , eigenVectors_(data, ref)
         , eigenValues_(data.size(), 0.)
         , SupportEigenVectors_(2*data.size(), 0)
{}
/* Copy constructor.
 *  @param S the EigenValue to copy
 **/
ISymEigen::ISymEigen( ISymEigen const& eigen)
                    : Base(eigen)
                    , norm_(eigen.norm_), rank_(eigen.rank_), det_(eigen.det_)
                    , eigenVectors_(eigen.eigenVectors_)
                    , eigenValues_(eigen.eigenValues_)
                    , SupportEigenVectors_(eigen.SupportEigenVectors_)
{}

/* Operator = : overwrite the ISymEigen with S.
 *  @param S ISymEigen to copy
 *  @return a reference on this
 **/
ISymEigen& ISymEigen::operator=( ISymEigen const& eigen)
{
  norm_   = eigen.norm_;     // norm of the matrix
  rank_   = eigen.rank_;     // rank of the matrix
  det_    = eigen.det_;      // determinant of the matrix
  eigenVectors_ = eigen.eigenVectors_;
  eigenValues_  = eigen.eigenValues_;
  SupportEigenVectors_ = eigen.SupportEigenVectors_;
  return *this;
}

/* finalize the computation by computing the rank, the trace norm and the
 * determinant of the matrix.
 **/
void ISymEigen::finalizeStep()
{
  norm_ = 0;
  rank_ = eigenValues_.size();
  det_  = 0;
  // compute trace norm sign of the determinant
  int s = 1;
  for (int i=eigenValues_.begin(); i< eigenValues_.end(); ++i )
  {
    Real value = eigenValues_.elt(i);
    norm_ += value;
    s     *= sign(value);
    if (std::abs(value) < Arithmetic<Real>::epsilon()) { rank_--;}
  }
  // compute the determinant for full rank matrices
  if (rank_ == eigenValues_.size())
  {
    Real sum = 0.;
    for (int i=eigenValues_.begin(); i< eigenValues_.end(); ++i )
    { sum += std::log(std::abs(eigenValues_.elt(i)));}
    det_ = s* std::exp(sum);
  }
}

} // namespace STK


