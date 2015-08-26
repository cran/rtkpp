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
 * Project:  stkpp::Algebra
 * Purpose:  Define The Interface ISvd Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_ISvd.h
 *  @brief In this file we define the interface class ISvd.
 **/
 
#ifndef STK_ISVD_H
#define STK_ISVD_H

#include <Sdk/include/STK_IRunner.h>
#include <Sdk/include/STK_IRecursiveTemplate.h>
#include <STKernel/include/STK_Real.h>

#include "STK_Algebra_Util.h"

namespace STK
{
/** @ingroup Algebra
 *  @brief Compute the Singular Value Decomposition of an array.
 * 
 *  The method take as:
 *  - input: A matrix A(nrow,ncol)
 *  - output:
 *    -# U Array (nrow,ncolU).
 *    -# D Vector (ncol)
 *    -# V Array (ncol,ncol).
 *  and perform the decomposition: 
 *  - A = UDV'
 *  U can have more cols than A,
 *  and it is possible to compute some (all) vectors of Ker(A).
 **/
template<class Derived>
class ISvd  : public IRunnerBase, public IRecursiveTemplate<Derived>
{
  protected:
    typedef typename hidden::AlgebraTraits<Derived>::ArrayU ArrayU;
    typedef typename hidden::AlgebraTraits<Derived>::ArrayD ArrayD;
    typedef typename hidden::AlgebraTraits<Derived>::ArrayV ArrayV;
    typedef TransposeOperator<ArrayV> ArrayVT;
    /** Default constructor
     *  @param A the matrix to decompose.
     *  @param ref if true, U_ is a reference of A.
     *  @param withU if @c true save the left housolder transforms in @c U_.
     *  @param withV if @c true save the right housolder transforms in @c V_.
     **/
    inline ISvd( ArrayU const& A, bool ref, bool withU = true, bool withV = true)
        : U_(A, ref), V_(), VT_(V_), D_()
        , withU_(withU), withV_(withV), norm_(0), rank_(0)
    {}
    /** constructor with other kind of array/expression
     *  @param A the matrix/expression to decompose.
     *  @param withU if @c true save the left housolder transforms in @c U_.
     *  @param withV if @c true save the right housolder transforms in @c V_.
     */
    template<class OtherDerived>
    inline ISvd( ArrayBase<OtherDerived> const& A, bool withU = true, bool withV = true)
               : U_(A), V_(), VT_(V_), D_()
               , withU_(withU), withV_(withV), norm_(0), rank_(0)
    {}
    /** Copy Constructor
     *  @param S the Svd to copy
     **/
    inline ISvd( ISvd const& S)
        : U_(S.U_, S.U_.isRef()), V_(S.V_), VT_(V_), D_(S.D_)
        , withU_(S.withU_), withV_(S.withV_)
        , norm_(S.norm_), rank_(S.rank_)
    {}
    /** destructor. */
    inline virtual ~ISvd() {}
    /** Operator = : overwrite the ISvd with S.
     *  @param S the Svd to copy
     **/
    ISvd& operator=(const ISvd &S)
    {
      U_ =S.U_; V_ = S.V_;  D_ = S.D_;
      withU_ =  S.withU_; withV_ = S.withV_;
      norm_ = S.norm_; rank_ = S.rank_;
      return *this;
    }
    /** Finalize any operations that have to be done after the computation
     *  of the decomposition
     **/
    virtual void finalize()
    {
      // Compute the true max norm
      norm_ = D_.front();
      // Compute the rank
      rank_ = 0;
      for (int i=D_.begin(); i<D_.end(); i++)
        if (norm_+D_[i] != norm_) { rank_++;}
        else break;

    }
  public:
    /// @return the number of rows of U_
    inline int nrowU() const { return U_.sizeRows();}
    /// @return the number of columns of U_
    inline int ncolU() const { return U_.sizeCols();}
    /// @return the number of columns of D_
    inline int nrowD() const { return D_.sizeRows();}
    /// @return the number of columns of D_
    inline int ncolD() const { return D_.sizeCols();}
    /// @return the number of rows of V_
    inline int nrowV() const { return V_.sizeRows();}
    /// @return the number of columns of V_
    inline int ncolV() const { return V_.sizeCols();}
    /// @return the norm of the matrix
    inline Real normSup()  const { return norm_;}
    /// @return the rank of the matrix
    inline int rank()  const { return rank_;}
    /// @return U
    inline ArrayU const& getU() const { return U_;}
    /// @return  V
    inline ArrayV const& getV() const { return V_;}
    /// @return  V
    inline ArrayVT const& getVT() const { return VT_;}
    /// @return D
    inline ArrayD const&  getD() const { return D_;}
    /** implement the run method */
    virtual bool run()
    {
      if (U_.empty()) { return true;}
      // compute Svd decomposition
      if (!this->asDerived().runImpl()) return false;
      finalize();
      return true;
    }
    /** Compute the svd of the Array A and copy the data
     *  see the corresponding constructor Take care that if U_ was previously
     *  a reference, it cannot be modified.
     *  @param A is the matrix to decompose.
     *  @param withU if true, we save the left housolder transforms
     *  in U_.
     *  @param withV if true, we save the right housolder transforms
     *  in V_.
     **/
    template<class OtherArray>
    void setData( OtherArray const& A, bool withU = true, bool withV = true)
    {
      U_ = A;           // Copy A in U_
      withU_ = withU;   // copy withU_ value
      withV_ = withV;   // copy withV_ value
      V_.resize(0,0), D_.resize(0);
    }

  protected:
    /// U_ matrix
    ArrayU U_;
    /// V_ matrix
    ArrayV V_;
    /// transposed V_
    ArrayVT VT_;
    /// Diagonal array of the singular values
    ArrayD D_;
    /// Compute U_ ?
    bool withU_;
    /// Compute V_ ?
    bool withV_;
    /// norm of the matrix (largest singular value)
    Real norm_;
    /// rank of the matrix
    int  rank_;
};

} // namespace STK

#endif /* STK_ISVD_H */
