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
 * Project:  stkpp::AAModels
 * Purpose:  Implementation of the ILinearReduct interface using the local variance.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_LocalVariance.h In this file we define the LocalVariance class.
 *  @brief A LocalVariance is an implementation of the Interface class Index.
 **/

#ifndef STK_LOCALVARIANCE_H
#define STK_LOCALVARIANCE_H

#include "STK_Reduct_Util.h"
#include "STK_ILinearReduct.h"

#include <Arrays/include/STK_Array2D.h>
#include <STatistiK/include/STK_Stat_MultivariateReal.h>
#include <Algebra/include/STK_SymEigen.h>

namespace STK
{
/** @ingroup Reduct
 *  @brief A LocalVariance is an implementation of the abstract
 *  @c ILinearReduct class.
 *
 *  A LocalVariance is an Index which maximize the projected local
 *  variance of the data set.
 *  The class can use the algorithm of Prim or the minimal distance
 *  in order to compute the proximity graph defining the local variance.
 *
 *  This class derive from ILinearReduct which derive itself from IRunnerUnsupervised.
 *  The @c run() and @c run(weights) methods have been implemented in the
 *  ILinearReduct interface using the pure virtual methods @c maximizeStep()
 *  and @c MaximizeCriteria(weights).
 **/
template<class Array>
class LocalVariance : public ILinearReduct<Array, Vector>
{
  public:
    typedef typename Array::Row RowVector;
    typedef ILinearReduct<Array, Vector> Base;
    using Base::p_data_;
    using Base::p_reduced_;
    using Base::axis_;
    using Base::idx_values_;
    /** Constructor. the TypeGraph and the number of neighbors are
     *  given by the user and cannot modified.
     *  @param p_data pointyer on the data set to process
     *  @param type type of proximity graph to build
     *  @param nbNeighbor number of neighbors to use in the proximity graph
     */
    LocalVariance( Array const* p_data
                 , Reduct::TypeGraph type = Reduct::distance_
                 , int nbNeighbor =1);
    /** Constructor. the TypeGraph and the number of neighbors are
     *  given by the user and cannot modified.
     *  @param data data set to process
     *  @param type type of proximity graph to build
     *  @param nbNeighbor number of neighbors to use in the proximity graph
     */
    LocalVariance( Array const& data
                 , Reduct::TypeGraph type = Reduct::distance_
                 , int nbNeighbor =1);
    /** copy Constructor.
     *  @param reducer the reducer to copy
     */
    LocalVariance( LocalVariance const& reducer);

    /** clone pattern */
    inline virtual LocalVariance* clone() const { return new LocalVariance(*this);}
    /** Destructor */
    virtual ~LocalVariance();
    /** @return the number of neighbors used in the local covariance.*/
    inline int nbNeighbor() const { return nbNeighbor_;}
    /**@return an array with the index of the neighbors of an individual */
    inline ArrayXXi const& pred() const { return neighbors_;}
    /** @return the covariance matrix of the data set */
    inline ArraySquareX const& covariance() const { return covariance_;}
    /**@return the local covariance matrix computed using the proximity graph. */
    inline ArraySquareX const& localCovariance() const { return localCovariance_;}

  protected:
    /**
     * Compute the proximity graph if the data set is modified.
     */
    virtual void update();

    /** number of neighbors */
    Reduct::TypeGraph type_;
    /** number of neighbors */
    int nbNeighbor_;
    /** Predecessor array : store the spanning tree or the minimal distance
     * to the neighbors */
    ArrayXXi neighbors_;
    /** distance matrix : store the minimal distance to the neighbors */
    ArrayXX dist_;

    /** the local covariance Array */
    ArraySquareX localCovariance_;
    /** the covariance Array */
    ArraySquareX covariance_;
    /** Compute the axis by maximizing the ratio of the local variance on the
     *  total variance of the data set.
     */
    virtual void maximizeStep();
    /** Compute the axis by maximizing the ratio of the weighted local variance
     *  on the weighted total variance of the data set.
     *  @param weights the weights to use
     */
    virtual void maximizeStep( Vector const& weights);

    /** compute the minimal spanning tree */
    void prim();
    /** compute the minimal distance graph */
    void minimalDistance();

    /** compute the covariances matrices of the data set */
    void computeCovarianceMatrices();
    /** compute the weighted covariances matrices of the data set */
    void computeCovarianceMatrices( Vector const& weights);

  private:
    /** compute the axis using the first eigenvectors of the matrix
     * of covariance and local covariance
     **/
    void computeAxis();
    /** Default constructor. The TypeGraph and the number of neighbors are
     *  given by the user and are not modified.
     *  @param type type of proximity graph to build
     *  @param nbNeighbor number of neighbors to use in the proximity graph
     */
    LocalVariance( Reduct::TypeGraph type = Reduct::distance_, int nbNeighbor =1);
};

/*
 * Constructor.
 * @param data the input data set
 */
template<class Array>
LocalVariance<Array>::LocalVariance( Array const* p_data
                                   , Reduct::TypeGraph type, int nbNeighbor)
                                   : Base(p_data)
                                   , type_(type)
                                   , nbNeighbor_(nbNeighbor)
                                   , neighbors_()
                                   , dist_()

{
  if (!p_data) return;
  neighbors_.resize(p_data->rows(), Range(1,nbNeighbor_));
  dist_.resize(p_data->rows(), Range(1,nbNeighbor_));
  dist_ = Arithmetic<Real>::max();
  // compute minimal proximity graph of the data set
  switch (type_)
  {
    case Reduct::prim_:
      prim();
      break;
    case Reduct::distance_:
      minimalDistance();
      break;
    case Reduct::unknown_graph_:
      STKRUNTIME_ERROR_NO_ARG(LocalVariance::LocalVariance(data, type, nbNeighbor),unknown proximity graph);
      break;
  };
}

/*
 * Constructor.
 * @param data the input data set
 */
template<class Array>
LocalVariance<Array>::LocalVariance( Array const& data
                                   , Reduct::TypeGraph type, int nbNeighbor)
                                  : Base(data)
                                  , type_(type)
                                  , nbNeighbor_(nbNeighbor)
                                  , neighbors_(data.rows(), Range(1,nbNeighbor_))
                                  , dist_(data.rows(), Range(1,nbNeighbor_), Arithmetic<Real>::max())
{
  // compute minimal proximity graph of the data set
  switch (type_)
  {
    case Reduct::prim_:
      prim();
      break;
    case Reduct::distance_:
      minimalDistance();
      break;
    case Reduct::unknown_graph_:
      STKRUNTIME_ERROR_NO_ARG(LocalVariance::LocalVariance(data, type, nbNeighbor),unknown proximity graph);
      break;
  };
}

/** copy Constructor.
 *  @param reducer the reducer to copy
 */
template<class Array>
LocalVariance<Array>::LocalVariance( LocalVariance const& reducer)
                                  : Base(reducer)
                                  , type_(reducer.type_)
                                  , nbNeighbor_(reducer.nbNeighbor_)
                                  , neighbors_(reducer.neighbors_)
                                  , dist_(reducer.dist_)
{}

/*
 * Destructor
 */
template<class Array>
LocalVariance<Array>::~LocalVariance() {}

/*
 * set the data set to use.
 */
template<class Array>
void LocalVariance<Array>::update()
{
#ifdef STK_REDUCT_DEBUG
  if (!p_data_)
  { STKRUNTIME_ERROR_NO_ARG(LocalVariance::update,data is not set);}
#endif
  // update dimensions of the containers for the proximity graph
  neighbors_.resize(p_data_->rows(), Range(1,nbNeighbor_));
  // compute minimal proximity graph of the data set
  switch (type_)
  {
    case Reduct::prim_:
      prim();
      break;
    case Reduct::distance_:
      minimalDistance();
      break;
    case Reduct::unknown_graph_:
      STKRUNTIME_ERROR_NO_ARG(LocalVariance::update,unknown proximity graph);
      break;
  };
}

/*
 * Compute the axis by maximizing the ratio of the local variance on the
 * total variance of the data set.
 */
template<class Array>
void LocalVariance<Array>::maximizeStep()
{
#ifdef STK_REDUCT_DEBUG
  if (!p_data_)
  { STKRUNTIME_ERROR_NO_ARG(LocalVariance::maximizeStep,data is not set);}
#endif
  // compute covariance matrices
  computeCovarianceMatrices();
  // compute the axis
  computeAxis();
}

/*
 * Compute the axis by maximizing the ratio of the local variance on the
 * total variance of the data set.
 */
template<class Array>
void LocalVariance<Array>::maximizeStep(Vector const& weights)
{
#ifdef STK_REDUCT_DEBUG
  if (!p_data_)
  { STKRUNTIME_ERROR_NO_ARG(LocalVariance::maximizeStep,data is not set);}
#endif
  // compute covariance matrices
  computeCovarianceMatrices(weights);
  // compute the axis
  computeAxis();
}

/* compute the covariances matrices of the data set
 *  @param nbNeighbor number of neighbors to look at
 **/
template<class Array>
void LocalVariance<Array>::computeCovarianceMatrices()
{
  Stat::Multivariate<Array, Real> stats_;
  // compute the covariance matrix
  stats_.setData(p_data_);
  stats_.run();
  covariance_.move(stats_.covariance());
  // constants
  const Real pond = 2* nbNeighbor_ * p_data_->sizeRows();

  // compute local covariance matrix
  localCovariance_.resize(p_data_->cols());
  for (int j=p_data_->beginCols(); j<p_data_->endCols(); j++)
  {
    for (int k=p_data_->beginCols(); k<p_data_->endCols(); k++)
    {
      Real sum = 0.0;
      for (int i=p_data_->beginRows(); i<p_data_->endRows(); i++)
      {
        for (int l = 1; l <= nbNeighbor_; ++l)
        {
          sum += ((*p_data_)(i, j) - (*p_data_)(neighbors_(i, l), j))
               * ((*p_data_)(i, k) - (*p_data_)(neighbors_(i, l), k));
        }
      }
      localCovariance_(j, k) = sum/pond;
    }
  }
}

/* compute the weighted covariances matrices of the data set */
template<class Array>
void LocalVariance<Array>::computeCovarianceMatrices( Vector const& weights)
{
  Stat::Multivariate<Array, Real> stats_;
  // compute the covariance matrix
  stats_.setData(p_data_);
  stats_.run(weights);
  covariance_.move(stats_.covariance());

  // get dimensions
  const Real pond = 2* nbNeighbor_ * p_data_->sizeRows() ;
  // compute weighted local covariance matrix
  localCovariance_.resize(p_data_->cols());
  for (int i=p_data_->beginCols(); i<p_data_->endCols(); i++)
  {
    for (int j=p_data_->beginCols(); j<p_data_->endCols(); j++)
    {
      Real sum = 0.0;
      for (int k=p_data_->beginRows(); k<p_data_->endRows(); k++)
      {
        for (int l = 1; l <= nbNeighbor_; ++l)
        {
          sum += (weights[k] * weights[neighbors_(k, l)])
               * ((*p_data_)(k, i) - (*p_data_)(neighbors_(k, l), i))
               * ((*p_data_)(k, j) - (*p_data_)(neighbors_(k, l), j));
        }
      }
      localCovariance_(i, j) = sum / pond;
    }
  }
}

/* compute the axis
 **/
template<class Array>
void LocalVariance<Array>::computeAxis()
{
  // compute the number of axis, cannot be greater than the dimension of the data
  Range range(p_data_->beginCols(), std::min(this->dim_, p_data_->sizeCols()));
  // compute the eigenvalues decomposition of the local covariance
  SymEigen<ArraySquareX>* decomp = new SymEigen<ArraySquareX>(localCovariance_);
  decomp->run();
  // compute the generalized inverse of the local covariance
  ArraySquareX inv_local;
  decomp->ginv(inv_local);
  // compute the product
  ArraySquareX prod = inv_local * covariance_;
  // compute the eigenvalues decomposition of the product
  decomp->setData(prod);
  decomp->run();
  // save axis and index values
  axis_.resize(p_data_->cols(), range);
  idx_values_.resize(range);
  axis_       = decomp->rotation().col(range);
  idx_values_ = decomp->eigenValues()[range];
  // we can safely remove the decomposition
  delete decomp;
}

template<class Array>
void LocalVariance<Array>::prim()
{
  // get dimensions
  const int begin_ind = p_data_->beginRows();
  const int last_ind = p_data_->lastIdxRows();
  /* value vector : store the minimal value. */
  Vector value(p_data_->rows(), Arithmetic<Real>::max());
  /* position of the points */
  Array1D<int> ipos(p_data_->rows());
  // Initialization the position array
  for (int i=begin_ind; i<=last_ind; i++) ipos[i] = i;

  // start Prim algorithm
  //Initialization of the root
  value[begin_ind] = 0.0;               // the root have value 0.0
  neighbors_(begin_ind, 1) = begin_ind;          // and have itself as predecessor
  int imin = begin_ind;             // the index of the current minimal value
  Real    kmin = 0.0;                   // the current minimal value
  // begin iterations
  for (int iter = last_ind; iter>=begin_ind; iter--)
  {
    // put the minimal key at the end of the array key_
    value.swap(imin, iter);  // put the minimal value to the end
    ipos.swap(imin, iter);   // update the position of current minimal point
    // Update the value for the neighbors points and find minimal value
    imin = begin_ind;
    kmin = value[begin_ind];
    // reference on the current point
    int icur = ipos[iter];
    RowVector P(p_data_->row(icur), true);
    // update distance of the neighbors point
    for (int i=begin_ind; i<iter; i++)
    {
      // check if we have a better distance for the neighbors
      Real d=dist(P, p_data_->row(ipos[i]));
      if (d < value[i])
      {
        value[i] = d;
        neighbors_(ipos[i], 1) = icur;
      }
      // minimal key
      if (kmin>value[i]) { imin=i; kmin = value[i];}
    }
  }
}

template<class Array>
void LocalVariance<Array>::minimalDistance()
{
  dist_.resize(p_data_->rows(), Range(1,nbNeighbor_));
  dist_ = Arithmetic<Real>::max();
  // get dimensions
  const int begin_ind = p_data_->beginRows();
  const int last_ind = p_data_->lastIdxRows();
  // start minimal distance algorithm
  for (int j = begin_ind; j<last_ind; j++)
  {
    // reference on the current point
    RowVector curPoint(p_data_->row(j), true);
    // update distance of the neighbors point
    for (int i=j+1; i<=last_ind; i++)
    {
      // compute distance between point i and point j
      Real d=dist(curPoint, p_data_->row(i));
      // check if we get a better distance for the point j
      if (dist_(i, nbNeighbor_) > d )
      {
        // check if we get a better distance for the point i
        int pos = nbNeighbor_;
        while (dist_(i, pos) > d && pos-- > 1) {}
        pos++;
        // shift values
        for (int k= nbNeighbor_ -1; k>= pos; k--)
        {
          dist_(i, k+1) = dist_(i, k);
          neighbors_(i, k+1) = neighbors_(i, k);
        }
        // set minimal distance in place
        dist_(i, pos) = d;
        neighbors_(i, pos) = j;
      }
      // check if we get a better distance for the point j
      if (dist_(j, nbNeighbor_) > d )
      {
        // insertion sorting algorihtm
        int pos = nbNeighbor_;
        while (dist_(j, pos) > d && pos-- > 1) {}
        pos++;
        // shift valuesconst
        for (int k= nbNeighbor_ -1; k>= pos; k--)
        {
          dist_(j, k+1) = dist_(j, k);
          neighbors_(j, k+1) = neighbors_(j, k);
        }
        // set minimal distance in place
        dist_(j, pos) = d;
        neighbors_(j, pos) = i;
      }
    }
  }
}


} // namespace STK
#endif /* STK_LOCALVARIANCE_H */
