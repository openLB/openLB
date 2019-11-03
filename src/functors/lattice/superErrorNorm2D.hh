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

#ifndef ERROR_NORM_2D_HH
#define ERROR_NORM_2D_HH

#include "superErrorNorm2D.h"
#include "utilities/functorDsl2D.h"

namespace olb {


template <typename T, typename W, int P>
template <typename DESCRIPTOR>
SuperRelativeErrorLpNorm2D<T,W,P>::SuperRelativeErrorLpNorm2D(
  SuperLattice2D<T,DESCRIPTOR>&      sLattice,
  FunctorPtr<SuperF2D<T,W>>&&        f,
  FunctorPtr<AnalyticalF2D<T,W>>&&   wantedF,
  FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF)
  : SuperIdentity2D<T,W>([&]()
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
SuperRelativeErrorLpNorm2D<T,W,P>::SuperRelativeErrorLpNorm2D(
  SuperLatticeF2D<T,DESCRIPTOR>&     f,
  FunctorPtr<AnalyticalF2D<T,W>>&&   wantedF,
  FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF)
  : SuperRelativeErrorLpNorm2D(f.getSuperLattice(),
                               f,
                               std::forward<decltype(wantedF)>(wantedF),
                               std::forward<decltype(indicatorF)>(indicatorF))
{ }


template <typename T, typename W, int P>
template <typename DESCRIPTOR>
SuperAbsoluteErrorLpNorm2D<T,W,P>::SuperAbsoluteErrorLpNorm2D(
  SuperLattice2D<T,DESCRIPTOR>&      sLattice,
  FunctorPtr<SuperF2D<T,W>>&&        f,
  FunctorPtr<AnalyticalF2D<T,W>>&&   wantedF,
  FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF)
  : SuperIdentity2D<T,W>([&]()
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
SuperAbsoluteErrorLpNorm2D<T,W,P>::SuperAbsoluteErrorLpNorm2D(
  SuperLatticeF2D<T,DESCRIPTOR>&     f,
  FunctorPtr<AnalyticalF2D<T,W>>&&   wantedF,
  FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF)
  : SuperAbsoluteErrorLpNorm2D(f.getSuperLattice(),
                               f,
                               std::forward<decltype(wantedF)>(wantedF),
                               std::forward<decltype(indicatorF)>(indicatorF))
{ }


}

#endif
