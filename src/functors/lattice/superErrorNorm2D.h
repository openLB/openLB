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

#ifndef ERROR_NORM_2D_H
#define ERROR_NORM_2D_H

#include "superBaseF2D.h"

namespace olb {


/// Relative error norm functor
/**
 * Calculates "LpNorm(wantedF - f) / LpNorm(wantedF)"
 **/
template <typename T, typename W, int P>
class SuperRelativeErrorLpNorm2D : public SuperIdentity2D<T,W> {
public:
  template <typename DESCRIPTOR>
  SuperRelativeErrorLpNorm2D(
    SuperLattice2D<T,DESCRIPTOR>&      sLattice,
    FunctorPtr<SuperF2D<T,W>>&&        f,
    FunctorPtr<AnalyticalF2D<T,W>>&&   wantedF,
    FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF);
  template <typename DESCRIPTOR>
  SuperRelativeErrorLpNorm2D(
    SuperLatticeF2D<T,DESCRIPTOR>&     f,
    FunctorPtr<AnalyticalF2D<T,W>>&&   wantedF,
    FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF);
};

template <typename T, typename W=T>
using SuperRelativeErrorL1Norm2D = SuperRelativeErrorLpNorm2D<T,W,1>;

template <typename T, typename W=T>
using SuperRelativeErrorL2Norm2D = SuperRelativeErrorLpNorm2D<T,W,2>;

template <typename T, typename W=T>
using SuperRelativeErrorLinfNorm2D = SuperRelativeErrorLpNorm2D<T,W,0>;


/// Absolute error norm functor
/**
 * Calculates "LpNorm(wantedF - f)"
 **/
template <typename T, typename W, int P>
class SuperAbsoluteErrorLpNorm2D : public SuperIdentity2D<T,W> {
public:
  template <typename DESCRIPTOR>
  SuperAbsoluteErrorLpNorm2D(
    SuperLattice2D<T,DESCRIPTOR>&      sLattice,
    FunctorPtr<SuperF2D<T,W>>&&        f,
    FunctorPtr<AnalyticalF2D<T,W>>&&   wantedF,
    FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF);
  template <typename DESCRIPTOR>
  SuperAbsoluteErrorLpNorm2D(
    SuperLatticeF2D<T,DESCRIPTOR>&     f,
    FunctorPtr<AnalyticalF2D<T,W>>&&   wantedF,
    FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF);
};

template <typename T, typename W=T>
using SuperAbsoluteErrorL1Norm2D = SuperAbsoluteErrorLpNorm2D<T,W,1>;

template <typename T, typename W=T>
using SuperAbsoluteErrorL2Norm2D = SuperAbsoluteErrorLpNorm2D<T,W,2>;

template <typename T, typename W=T>
using SuperAbsoluteErrorLinfNorm2D = SuperAbsoluteErrorLpNorm2D<T,W,0>;


}

#endif
