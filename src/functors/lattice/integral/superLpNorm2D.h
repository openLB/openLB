/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef SUPER_LP_NORM_2D_H
#define SUPER_LP_NORM_2D_H

#include <memory>
#include <vector>

#include "utilities/functorPtr.h"

namespace olb {

template<typename T, typename W> class SuperF2D;
template<typename T> class SuperIndicatorF2D;
template<typename T> class SuperGeometry2D;

/// Functor that returns the Lp norm over omega of the the euklid norm of the input functor
/**
 * Maintains block level BlockLpNorm2D functors as required.
 *
 * P == 0: inf norm
 * P >= 1: p norm
 */
template <typename T, typename W, int P>
class SuperLpNorm2D : public SuperF2D<T,W> {
private:
  FunctorPtr<SuperF2D<T,W>>        _f;
  FunctorPtr<SuperIndicatorF2D<T>> _indicatorF;
public:
  /**
   * \param f          (non-)owning pointer or reference to SuperF2D<T,W>
   * \param indicatorF (non-)owning pointer or reference to SuperIndicatorF2D<T>.
   *                   Describes the subset to be integrated.
   **/
  SuperLpNorm2D(FunctorPtr<SuperF2D<T,W>>&&        f,
                FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF);

  /**
   * Legacy constructor accepting super geometry reference.
   *
   * \param f          (non-)owning pointer or reference to SuperF2D<T,W>
   * \param indicatorF (non-)owning pointer or reference to SuperIndicatorF2D<T>.
   *                   Describes the subset to be integrated.
   **/
  SuperLpNorm2D(FunctorPtr<SuperF2D<T,W>>&& f,
                SuperGeometry2D<T>&,
                FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF);

  /**
   * \param f         (non-)owning pointer or reference to SuperF2D<T,W>
   * \param geometry  super geometry required to construct SuperIndicatorMaterial2D using materials
   * \param materials vector of material numbers to be included in the Lp norm
   **/
  SuperLpNorm2D(FunctorPtr<SuperF2D<T,W>>&& f,
                SuperGeometry2D<T>& geometry,
                std::vector<int>    materials);

  /**
   * \param f        (non-)owning pointer or reference to SuperF2D<T,W>
   * \param geometry super geometry required to construct SuperIndicatorMaterial2D using material
   * \param material single material number to be included in the Lp norm
   **/
  SuperLpNorm2D(FunctorPtr<SuperF2D<T,W>>&& f,
                SuperGeometry2D<T>& geometry,
                int                 material);

  bool operator() (W output[], const int input[]) override;
};


/// Functor that returns the L1 norm over omega of the the euklid norm of the input functor
template <typename T, typename W = T>
using SuperL1Norm2D = SuperLpNorm2D<T,W,1>;

/// Functor that returns the L2 norm over omega of the the euklid norm of the input functor
template <typename T, typename W = T>
using SuperL2Norm2D = SuperLpNorm2D<T,W,2>;

/// Functor that returns the Linf norm over omega of the the euklid norm of the input functor
template <typename T, typename W = T>
using SuperLinfNorm2D = SuperLpNorm2D<T,W,0>;

}

#endif
