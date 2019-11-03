/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef ERROR_NORM_3D_HH
#define ERROR_NORM_3D_HH

#include "superErrorNorm3D.h"
#include "utilities/functorDsl3D.h"

namespace olb {


template <typename T, typename W, int P>
template <typename DESCRIPTOR>
SuperRelativeErrorLpNorm3D<T,W,P>::SuperRelativeErrorLpNorm3D(
  SuperLattice3D<T,DESCRIPTOR>&      sLattice,
  FunctorPtr<SuperF3D<T,W>>&&        f,
  FunctorPtr<AnalyticalF3D<T,W>>&&   wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperIdentity3D<T,W>([&]()
{
  using namespace functor_dsl;

  auto wantedLatticeF = restrict(wantedF.toShared(), sLattice);

  return norm<P>(wantedLatticeF - f.toShared(), indicatorF.toShared())
         / norm<P>(wantedLatticeF, indicatorF.toShared());
}())
{
  this->getName() = "relErrorNormL" + std::to_string(P);
}

template <typename T, typename W, int P>
template <typename DESCRIPTOR>
SuperRelativeErrorLpNorm3D<T,W,P>::SuperRelativeErrorLpNorm3D(
  SuperLatticeF3D<T,DESCRIPTOR>&     f,
  FunctorPtr<AnalyticalF3D<T,W>>&&   wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperRelativeErrorLpNorm3D(f.getSuperLattice(),
                               f,
                               std::forward<decltype(wantedF)>(wantedF),
                               std::forward<decltype(indicatorF)>(indicatorF))
{ }


template <typename T, typename W, int P>
template <typename DESCRIPTOR>
SuperAbsoluteErrorLpNorm3D<T,W,P>::SuperAbsoluteErrorLpNorm3D(
  SuperLattice3D<T,DESCRIPTOR>&      sLattice,
  FunctorPtr<SuperF3D<T,W>>&&        f,
  FunctorPtr<AnalyticalF3D<T,W>>&&   wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperIdentity3D<T,W>([&]()
{
  using namespace functor_dsl;

  return norm<P>(restrict(wantedF.toShared(), sLattice) - f.toShared(),
                 indicatorF.toShared());
}())
{
  this->getName() = "absErrorNormL" + std::to_string(P);
}

template <typename T, typename W, int P>
template <typename DESCRIPTOR>
SuperAbsoluteErrorLpNorm3D<T,W,P>::SuperAbsoluteErrorLpNorm3D(
  SuperLatticeF3D<T,DESCRIPTOR>&     f,
  FunctorPtr<AnalyticalF3D<T,W>>&&   wantedF,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperAbsoluteErrorLpNorm3D(f.getSuperLattice(),
                               f,
                               std::forward<decltype(wantedF)>(wantedF),
                               std::forward<decltype(indicatorF)>(indicatorF))
{ }


}

#endif
