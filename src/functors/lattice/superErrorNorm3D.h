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

#ifndef ERROR_NORM_3D_H
#define ERROR_NORM_3D_H

#include "superBaseF3D.h"

namespace olb {


/// Relative error norm functor
/**
 * Calculates "LpNorm(wantedF - f) / LpNorm(wantedF)"
 **/
template <typename T, typename W, int P>
class SuperRelativeErrorLpNorm3D : public SuperIdentity3D<T,W> {
public:
  template <typename DESCRIPTOR>
  SuperRelativeErrorLpNorm3D(
    SuperLattice3D<T,DESCRIPTOR>&      sLattice,
    FunctorPtr<SuperF3D<T,W>>&&        f,
    FunctorPtr<AnalyticalF3D<T,W>>&&   wantedF,
    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  template <typename DESCRIPTOR>
  SuperRelativeErrorLpNorm3D(
    SuperLatticeF3D<T,DESCRIPTOR>&     f,
    FunctorPtr<AnalyticalF3D<T,W>>&&   wantedF,
    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
};

template <typename T, typename W=T>
using SuperRelativeErrorL1Norm3D = SuperRelativeErrorLpNorm3D<T,W,1>;

template <typename T, typename W=T>
using SuperRelativeErrorL2Norm3D = SuperRelativeErrorLpNorm3D<T,W,2>;

template <typename T, typename W=T>
using SuperRelativeErrorLinfNorm3D = SuperRelativeErrorLpNorm3D<T,W,0>;


/// Absolute error norm functor
/**
 * Calculates "LpNorm(wantedF - f)"
 **/
template <typename T, typename W, int P>
class SuperAbsoluteErrorLpNorm3D : public SuperIdentity3D<T,W> {
public:
  template <typename DESCRIPTOR>
  SuperAbsoluteErrorLpNorm3D(
    SuperLattice3D<T,DESCRIPTOR>&      sLattice,
    FunctorPtr<SuperF3D<T,W>>&&        f,
    FunctorPtr<AnalyticalF3D<T,W>>&&   wantedF,
    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  template <typename DESCRIPTOR>
  SuperAbsoluteErrorLpNorm3D(
    SuperLatticeF3D<T,DESCRIPTOR>&     f,
    FunctorPtr<AnalyticalF3D<T,W>>&&   wantedF,
    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
};

template <typename T, typename W=T>
using SuperAbsoluteErrorL1Norm3D = SuperAbsoluteErrorLpNorm3D<T,W,1>;

template <typename T, typename W=T>
using SuperAbsoluteErrorL2Norm3D = SuperAbsoluteErrorLpNorm3D<T,W,2>;

template <typename T, typename W=T>
using SuperAbsoluteErrorLinfNorm3D = SuperAbsoluteErrorLpNorm3D<T,W,0>;


}

#endif
