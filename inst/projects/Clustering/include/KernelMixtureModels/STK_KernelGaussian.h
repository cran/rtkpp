/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013 Serge Iovleff

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: Oct 24, 2013
 * Author:   Serge Iovleff
 **/

/** @file STK_KernelGaussian.h
 *  @brief In this file we implement the KernelGaussian class
 **/

#ifndef STK_KERNELGAUSSIAN_H
#define STK_KERNELGAUSSIAN_H

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class KernelGaussian;

namespace Clust
{
/** @ingroup Clustering
 *  Traits class for the KernelGaussian traits policy. */
template<class _Array>
struct MixtureTraits< KernelGaussian<_Array> >
{
  typedef _Array Array;
  typedef typename Array::Type       Type;
  typedef KernelGaussian_Parameters  Parameters;
  typedef Array2D<Real>              Param;
};

} // namespace hidden

/** @ingroup Clustering
 *  The diagonal Gaussian mixture model @c KernelGaussian is
 *  an isotrope Gaussian model on a kernel space. It as a density function of the
 *  form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d
 *    \frac{1}{\sqrt{2\pi}\sigma^j_{k}} \exp\left\{-\frac{(\phi(x^j_)-\mu^j_{ik})^2}{2(\sigma^j_{k})^2}\right\}.
 * \f]
 **/
template<class Array>
class KernelGaussian : public IMixtureModel<KernelGaussian<Array> >
{
  public:
    typedef IMixtureModel<KernelGaussian<Array> > Base;

    typedef typename Clust::MixtureTraits< KernelGaussian<Array> >::Parameters Parameters;

    using Base::p_tik;
    using Base::p_nk;
    using Base::components;
    using Base::p_data;
    using Base::param;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    KernelGaussian( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    KernelGaussian( KernelGaussian const& model) : Base(model) {}
    /** destructor */
    ~KernelGaussian() {}
    /** Initialize randomly the parameters of the Gaussian mixture. The
     *  standard-deviation will be set to 1.
     */
    void randomInit();
    /** Compute the weighted mean and the common standard deviation. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster();}
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void KernelGaussian<Array>::randomInit()
{
  // compute the standard deviation
  for (int k= baseIdx; k < components().end(); ++k)
  {
    param(k).sigma_ = p_data()->col(k).dot(p_tik()->col(k))/p_nk()->elt(k);
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("KernelGaussian<Array>::randomInit() done\n");
#endif
}

/* Compute the weighted means and the weighted standard deviations. */
template<class Array>
bool KernelGaussian<Array>::mStep()
{
  // compute the standard deviation
  for (int k= baseIdx; k < components().end(); ++k)
  {
    param(k).sigma_ = p_data()->col(k).dot(p_tik()->col(k))/p_nk()->elt(k);
    if (param(k).sigma_ <= 0.) return false;
  }
  return true;
}

} // namespace STK

#endif /* STK_KERNELGAUSSIAN_H */
