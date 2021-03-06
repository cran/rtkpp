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
 * Project:  stkpp::Arrays
 * created on: 13 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ExprBase.h
 *  @brief In this file we define the base class for Arrays and Expressions
 **/

#ifndef STK_EXPRBASE_H
#define STK_EXPRBASE_H

#include "STK_ITContainer.h"

/// utility macro allowing to construct binary operators
#define MAKE_BINARY_OPERATOR(OPERATOR, FUNCTOR) \
  template<typename Rhs> \
  inline BinaryOperator< FUNCTOR<Type, typename hidden::Traits<Rhs>::Type>, Derived, Rhs> const \
  (OPERATOR)( ExprBase<Rhs> const& other) const \
  { return BinaryOperator<FUNCTOR<Type, typename hidden::Traits<Rhs>::Type>, Derived, Rhs>(this->asDerived(), other.asDerived()) ;}

/// utility macro allowing to construct unary operators
#define MAKE_UNARY_OPERATOR_NOARG(FUNCTION, FUNCTOR) \
  inline UnaryOperator<FUNCTOR<Type>, Derived> const FUNCTION() const \
  { return UnaryOperator<FUNCTOR<Type>, Derived>(this->asDerived()); }

/// utility macro allowing to construct unary operators
#define MAKE_UNARY_OPERATOR_1ARG(FUNCTION, FUNCTOR) \
  inline UnaryOperator<FUNCTOR<Type>, Derived> const FUNCTION(Type const number) const \
  { return UnaryOperator<FUNCTOR<Type>, Derived>(this->asDerived(), FUNCTOR<Type>(number)); }


// forward declarations
namespace STK
{
template<class Derived> class ExprBase;
template<class Derived> class ArrayBase;
template<class Derived, class Rhs> struct  ProductReturnType;
template<class Derived> class  ArrayInitializer;
} // namespace STK

#include "STKernel/include/STK_Functors.h"

#include "products/STK_ProductOperators.h"
#include "operators/STK_TransposeOperator.h"
#include "operators/STK_DiagonalOperator.h"
#include "operators/STK_UnaryOperators.h"
#include "operators/STK_BinaryOperators.h"
#include "operators/STK_DotOperators.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief base class for templated evaluation expressions and visitors.
 *
 * This class is the base that is inherited by all objects (matrix, vector,
 * point and expression). The common API for these objects is contained in
 * this class.
 *
 *  @tparam Derived is the derived type, e.g., a matrix, vector, type or
 *  an expression.
 **/

template<class Derived>
class ExprBase : public ITContainer<Derived, hidden::Traits<Derived>::structure_>
{
  protected:
    /** Default constructor */
    ExprBase() : Base() {}
    /** destructor */
    ~ExprBase() {}

  public:
    typedef ITContainer<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    //--------------
    // Start Visitors
    /** @brief Visit the container using a constant visitor
     *  @param visitor the visitor to run
     **/
    template<typename Visitor>
    typename Visitor::ReturnType visit(Visitor& visitor) const;
    /** @brief compute the number of non-zero element in an expression.
     *  For example
     *  @code
     *  (a > 0).count();
     *  @endcode
     *  compute the number of positive element of the array @c a.
     *  @return the number of non-zero element in the expression.*/
    int count() const;
    /** @brief check if there is any non-zero element in an expression.
     * For example
     *  @code
     *  (a > 0).any();
     *  @endcode
     *  will return @c true if there exists positive numbers in the expression @c a.
     *  @return @c true if at least one element is not zero in the expression, @c false otherwise.*/
    bool const any() const;
    /** @brief check if all the elements in an expression are not zero.
     *  For example
     *  @code
     *  (a > 0).all();
     *  @endcode
     *  will return @c true if all the elements in the expression @c a are positive.
     *   @return @c true if all the elements are not zero in the expression, @c false otherwise.*/
    bool const all() const;
    /** @return the number of available values in the array (not count NA values).*/
    int nbAvailableValues() const;

    /** @return the minimum of all elements of this using a Visitor
      * and puts in (row, col) its location.
      * @sa maxElt(int,int), visitor(), minElt()
      */
    Type const minElt( int& row, int& col) const;
    /** @return the minimum of all elements of @c *this which are not NA values
     *  using a Visitor and puts in (row, col) its location.
      * @sa maxEltSafe(int, int), visitor(), minElt()
      */
    Type const minEltSafe( int& row, int& col) const;
    /** @return the maximum of all elements of this using a Visitor
      * and puts in (row, col) its location.
      * @sa minElt(), visitor()
      */
    Type const maxElt( int& row, int& col) const;
    /** @return the maximum of all elements of this using a Visitor
      * and puts in (row, col) its location.
      * @sa minElt(), visitor()
      */
    Type const maxEltSafe( int& row, int& col) const;
    /** @return the minimum of all elements of this using a Visitor
      * and puts in  pos its location.
      * @note Have to be used for col-vector, row-vector or diagonal matrix only.
      * @sa maxElt(), visitor()
      */
    Type const minElt( int& pos) const;
    /** @return the minimum of all elements of this using a Visitor
      * and puts in  pos its location.
      * @note Have to be used for col-vector, row-vector or diagonal matrix only.
      * @sa maxElt(), visitor()
      */
    Type const minEltSafe( int& pos) const;
    /** @return the maximum of all elements of this using a Visitor
      * and puts in pos its location.
      * @note Have to be used for col-vector, row-vector or diagonal matrix only.
      * @sa minElt(), visitor()
      */
    Type const maxElt( int& pos) const;
    /** @return the maximum of all elements of this using a Visitor
      * and puts in pos its location.
      * @note Have to be used for col-vector, row-vector or diagonal matrix only.
      * @sa minElt(), visitor()
      */
    Type const maxEltSafe( int& pos) const;
    /** @return the minimum of all elements of this using a Visitor.
      * @sa maxElt(row, col), visitor()
      */
    Type const minElt() const;
    /** @return the minimum of all elements of this using a Visitor
      * @sa maxElt(row, pos), visitor()
      */
    Type const minEltSafe() const;
    /** @return the maximum of all elements of this using a Visitor.
      * @sa minElt(row, col), visitor()
      */
    Type const maxElt() const;
    /** @return the maximum of all elements of this using a Visitor
      * @sa minElt(row, pos), visitor()
      */
    Type const maxEltSafe() const;

    /** @return the sum of all the elements of this using a Visitor*/
    Type const sum() const;
    /** @return safely the sum of all the elements of this */
    Type const sumSafe() const;
    /** @return the norm of this*/
    Type const norm() const;
    /** @return the norm of this*/
    Type const normSafe() const;
    /** @return the square norm of this*/
    Type const norm2() const;
    /** @return the square norm of this*/
    Type const norm2Safe() const;
    /** @return the infinite norm of this*/
    Type const normInf() const;

    /** @return the mean of all the elements of this*/
    Type const mean() const;
    /** @return safely the mean of all the elements of this*/
    Type const meanSafe() const;
    /** @return the variance of all the elements of this*/
    Type const variance() const;
    /** @return the variance of all the elements of this*/
    Type const varianceSafe() const;
    /** @return the variance with given mean of all the elements of this*/
    Type const variance(Type const& mean) const;
    /** @return safely the variance with given mean of all the elements of this*/
    Type const varianceSafe(Type const& mean) const;

    /** @return the weighted sum of all the elements of this using a Visitor
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wsum(ExprBase<Rhs> const& weights) const;
    /** @return the safely weighted sum of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wsumSafe(ExprBase<Rhs> const& weights) const;
    /** @return the weighted norm of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wnorm(ExprBase<Rhs> const& weights) const;
    /** @return the weighted norm of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wnormSafe(ExprBase<Rhs> const& weights) const;
    /** @return the weighted square norm of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wnorm2(ExprBase<Rhs> const& weights) const;
    /** @return the weighted square norm of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wnorm2Safe(ExprBase<Rhs> const& weights) const;
    /** @return the weighted mean of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wmean(ExprBase<Rhs> const& weights) const;
    /** @return the safely weighted mean of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wmeanSafe(ExprBase<Rhs> const& weights) const;
    /** @return the weighted variance of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wvariance(ExprBase<Rhs> const& weights) const;
    /** @return the weighted variance of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wvarianceSafe(ExprBase<Rhs> const& weights) const;
    /** @return the variance with given mean of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wvariance(Type const& mean, ExprBase<Rhs> const& weights) const;
    /** @return safely the variance with given mean of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wvarianceSafe(Type const& mean, ExprBase<Rhs> const& weights) const;
    // Visitors terminated
    //--------------------

    /** @return an expression from the difference of this and  other. */
    MAKE_BINARY_OPERATOR(operator-,DifferenceOp)
    /** @return an expression with the addition of this and other. */
    MAKE_BINARY_OPERATOR(operator+,SumOp)
    /** @return an expression with the quotient of this and other. */
    MAKE_BINARY_OPERATOR(operator/,DivOp)
    /** @return an expression with the product (element by element) of this and other.*/
    MAKE_BINARY_OPERATOR(prod,ProductOp)
    /** @return an expression with the min of this and other.*/
    MAKE_BINARY_OPERATOR(min,MinOp)
    /** @return an expression from the max of this and other.*/
    MAKE_BINARY_OPERATOR(max,MaxOp)
    /** @return an expression from the \< operator of this and other*/
    MAKE_BINARY_OPERATOR(operator<,LessOp)
    /** @return an expression from the \<= operator of this and other.*/
    MAKE_BINARY_OPERATOR(operator<=,LeqOp)
    /** @return an expression from the \> operator of this and other.*/
    MAKE_BINARY_OPERATOR(operator>,GreaterOp)
    /** @return an expression from the \>= operator of this and other.*/
    MAKE_BINARY_OPERATOR(operator>=,GeqOp)
    /** @return an expression from the == operator of this and other.*/
    MAKE_BINARY_OPERATOR(operator==,EqualOp)
    /** @return an expression from the != operator of this and other.*/
    MAKE_BINARY_OPERATOR(operator!=,NotEqualOp)

    /** @return an expression of the opposite of this */
    MAKE_UNARY_OPERATOR_NOARG(operator-, OppositeOp)
    /** @return which values of this is a NA value */
    MAKE_UNARY_OPERATOR_NOARG(isNA, IsNaOp)
    /** @return not of the values */
    MAKE_UNARY_OPERATOR_NOARG(neg, NegOp)
    /** @return which values of this are finite value */
    MAKE_UNARY_OPERATOR_NOARG(isFinite, IsFiniteOp)
    /** @return which values of this are infinite value */
    MAKE_UNARY_OPERATOR_NOARG(isInfinite, IsInfiniteOp)
    /** @return an expression of the absolute value of this */
    MAKE_UNARY_OPERATOR_NOARG(abs, AbsOp)
    /** @return an expression of the exponential of this. */
    MAKE_UNARY_OPERATOR_NOARG(exp, ExpOp)
    /** @return an expression of the logarithm of this. */
    MAKE_UNARY_OPERATOR_NOARG(log, LogOp)
    /** @return an expression of the square root of this. */
    MAKE_UNARY_OPERATOR_NOARG(sqrt, SqrtOp)
    /** @return an expression of the cosine of this. */
    MAKE_UNARY_OPERATOR_NOARG(cos, CosOp)
    /** @return an expression of the sine of this.  */
    MAKE_UNARY_OPERATOR_NOARG(sin, SinOp)
    /** @return an expression of the arc cosine of this. */
    MAKE_UNARY_OPERATOR_NOARG(acos, AcosOp)
    /** @return an expression of the arc sine of this. */
    MAKE_UNARY_OPERATOR_NOARG(asin, AsinOp)
    /** @return an expression of the tan of this. */
    MAKE_UNARY_OPERATOR_NOARG(tan, TanOp)
    /** @return an expression of the inverse of this.  */
    MAKE_UNARY_OPERATOR_NOARG(inverse, InverseOp)
    /** @return an expression of the square of this. */
    MAKE_UNARY_OPERATOR_NOARG(square, SquareOp)
    /** @return an expression of the cube of this. */
    MAKE_UNARY_OPERATOR_NOARG(cube, CubeOp)

    /** @return an expression of the minimum of this and a number. */
    MAKE_UNARY_OPERATOR_1ARG(min, MinimumOp)
    /** @return an expression of the maximum of this and a number. */
    MAKE_UNARY_OPERATOR_1ARG(max, MaximumOp)
    /** @return an expression of this with each elements incremented by
     *  the constant number */
    MAKE_UNARY_OPERATOR_1ARG(operator+, AddOp)
    /** @return an expression of this scaled by the number factor number */
    MAKE_UNARY_OPERATOR_1ARG(operator*, MultipleOp)
    /** @return an expression of the power of this. */
    MAKE_UNARY_OPERATOR_1ARG(pow, PowOp)
    /** @return an expression of this divided by the number value number */
    MAKE_UNARY_OPERATOR_1ARG(safeInverse, SafeInverseOp)
    /** @return an expression of this divided by the number value number */
    MAKE_UNARY_OPERATOR_1ARG(operator/, QuotientOp)
    // boolean operations
    /** @return an expression of *this < number. */
    MAKE_UNARY_OPERATOR_1ARG(operator<, LessThanOp)
    /** @return an expression of *this <= number. */
    MAKE_UNARY_OPERATOR_1ARG(operator<=, LeqThanOp)
    /** @return an expression of *this > number. */
    MAKE_UNARY_OPERATOR_1ARG(operator>, GreaterThanOp)
    /** @return the expression of *this >= number. */
    MAKE_UNARY_OPERATOR_1ARG(operator>=, GeqThanOp)
    /** @return an expression of *this == number. */
    MAKE_UNARY_OPERATOR_1ARG(operator==, EqualThanOp)
    /** @return an expression of *this != number. */
    MAKE_UNARY_OPERATOR_1ARG(operator!=, NotEqualThanOp)

    // handle the case number + expression
    /** @return an expression of number + this */
    inline friend UnaryOperator<AddOp<Type>, Derived> const
    operator+(Type const number, ExprBase<Derived> const& other)
    { return other.asDerived() + number;}
    /** @return an expression of number - this */
    inline UnaryOperator<AddOp<Type>, Derived> const
    operator-(Type const number) const
    { return UnaryOperator<AddOp<Type>, Derived>(this->asDerived(), AddOp<Type>(-number));}
    /** @return a safe value of this */
    inline UnaryOperator<SafeOp<Type>, Derived> const safe(Type const number = Type()) const
    { return UnaryOperator<SafeOp<Type>, Derived>(this->asDerived(), SafeOp<Type>(number)); }
    // friends
    // handle the case number - expression
    inline friend UnaryOperator<AddOppositeOp<Type>, Derived > const
    operator-(Type const number, ExprBase<Derived> const& other)
    { return UnaryOperator<AddOppositeOp<Type>, Derived>(other.asDerived(), AddOppositeOp<Type>(number));}
    // handle the case number * expression
    inline friend UnaryOperator< MultipleOp<Type>, Derived> const
    operator*(Type const number, ExprBase<Derived> const& other)
    { return other.asDerived()*number; }
    // templated
    /** @return an expression of *this with the  Type type casted to  OtherType. */
    template<typename CastedType>
    inline UnaryOperator<CastOp<Type, CastedType>, Derived> const cast() const
    { return UnaryOperator<CastOp<Type, CastedType>, Derived>(this->asDerived());}
    /** @return an expression of funct0 to this. */
    template< template<typename> class OtherOperator>
    inline UnaryOperator<OtherOperator<Type>, Derived> const funct0() const
    { return UnaryOperator<OtherOperator<Type>, Derived>(this->asDerived());}
    /** @return an expression of funct1 to this. */
    template< template<typename> class OtherOperator>
    inline UnaryOperator<OtherOperator<Type>, Derived> const funct1(Type const number) const
    { return UnaryOperator<OtherOperator<Type>, Derived>(this->asDerived(), OtherOperator<Type>(number));}

    /** @return the transposed expression of this. */
    inline TransposeOperator<Derived> transpose() const
    { return TransposeOperator<Derived> (this->asDerived());}

    /** @return the diagonal expression of this (work only for 1D expressions). */
    inline DiagonalizeOperator<Derived> diagonalize() const
    { return DiagonalizeOperator<Derived> (this->asDerived());}
    /** @return the diagonal of this expression (work only for square expressions). */
    inline DiagonalOperator<Derived> diagonal() const
    { return DiagonalOperator<Derived> (this->asDerived());}

    /** @return the j-th column of this. */
    inline ColOperator<Derived> col(int j) const
    { return ColOperator<Derived> (this->asDerived(), j);}
    /** @return the i-th row of this. */
    inline RowOperator<Derived> row(int i) const
    { return RowOperator<Derived> (this->asDerived(), i);}
    /** @return the sub-vector(I) of this. */
    inline SubOperator<Derived> sub(Range I) const
    { return SubOperator<Derived> (this->asDerived(), I);}

    /** @returns the dot product of this with other.
      * @sa squaredNorm(), norm(), DotProduct
      */
    template<class Rhs>
    typename hidden::Promote<Type, typename Rhs::Type>::result_type const
    dot(ExprBase<Rhs> const& other) const;
    /** @returns the safe dot product of this with other.
      * @sa squaredNorm(), norm(), DotProduct
      */
    template<class Rhs>
    typename hidden::Promote<Type, typename Rhs::Type>::result_type const
    dotSafe(ExprBase<Rhs> const& other) const;

    /** @return the matrix multiplication of this with other.*/
    template<typename Rhs>
    typename ProductReturnType<Derived, Rhs>::ReturnType const
    operator*( ExprBase<Rhs> const& other) const;
};

#undef MAKE_BINARY_OPERATOR
#undef MAKE_UNARY_OPERATOR_NOARG
#undef MAKE_UNARY_OPERATOR_1ARG

} // namespace STK

#endif /* STK_EXPRBASE_H_ */
