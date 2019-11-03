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

#ifndef SUPER_PLANE_INTEGRAL_FLUX_MASS_2D_HH
#define SUPER_PLANE_INTEGRAL_FLUX_MASS_2D_HH

#include "superPlaneIntegralFluxMass2D.h"

#include "io/ostreamManager.h"
#include "utilities/vectorHelpers.h"
#include "utilities/functorPtr.hh"
#include "functors/lattice/indicator/indicator2D.hh"

#include "functors/lattice/superLatticeLocalF2D.h"
#include "functors/lattice/superCalcF2D.h"
#include "functors/lattice/superCalcF2D.hh"

namespace olb {


template<typename T>
SuperPlaneIntegralFluxMass2D<T>::SuperPlaneIntegralFluxMass2D(
  FunctorPtr<SuperF2D<T>>&& velocityF,
  FunctorPtr<SuperF2D<T>>&& densityF,
  SuperGeometry2D<T>&       geometry,
  T conversationFactorMass,
  T conversationFactorTime,
  const HyperplaneLattice2D<T>& hyperplaneLattice,
  FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
  FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF2D<T>(
      (*densityF) * (*velocityF),
      geometry,
      hyperplaneLattice,
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      std::forward<decltype(subplaneIndicator)>(subplaneIndicator),
      mode),
    _velocityF(std::move(velocityF)),
    _densityF(std::move(densityF)),
    _conversationFactorMass(conversationFactorMass),
    _conversationFactorTime(conversationFactorTime)
{
  this->getName() = "SuperPlaneIntegralFluxMass2D";
}

template<typename T>
SuperPlaneIntegralFluxMass2D<T>::SuperPlaneIntegralFluxMass2D(
  FunctorPtr<SuperF2D<T>>&& velocityF,
  FunctorPtr<SuperF2D<T>>&& densityF,
  SuperGeometry2D<T>&       geometry,
  T conversationFactorMass,
  T conversationFactorTime,
  const Hyperplane2D<T>& hyperplane,
  FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
  FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF2D<T>(
      (*densityF) * (*velocityF),
      geometry,
      hyperplane,
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      std::forward<decltype(subplaneIndicator)>(subplaneIndicator),
      mode),
    _velocityF(std::move(velocityF)),
    _densityF(std::move(densityF)),
    _conversationFactorMass(conversationFactorMass),
    _conversationFactorTime(conversationFactorTime)
{
  this->getName() = "SuperPlaneIntegralFluxMass2D";
}

template<typename T>
SuperPlaneIntegralFluxMass2D<T>::SuperPlaneIntegralFluxMass2D(
  FunctorPtr<SuperF2D<T>>&& velocityF,
  FunctorPtr<SuperF2D<T>>&& densityF,
  SuperGeometry2D<T>&       geometry,
  T conversationFactorMass,
  T conversationFactorTime,
  const Hyperplane2D<T>& hyperplane,
  FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF2D<T>(
      (*densityF) * (*velocityF),
      geometry,
      hyperplane,
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      mode),
    _velocityF(std::move(velocityF)),
    _densityF(std::move(densityF)),
    _conversationFactorMass(conversationFactorMass),
    _conversationFactorTime(conversationFactorTime)
{
  this->getName() = "SuperPlaneIntegralFluxMass2D";
}

template<typename T>
SuperPlaneIntegralFluxMass2D<T>::SuperPlaneIntegralFluxMass2D(
  FunctorPtr<SuperF2D<T>>&& velocityF,
  FunctorPtr<SuperF2D<T>>&& densityF,
  SuperGeometry2D<T>&       geometry,
  T conversationFactorMass,
  T conversationFactorTime,
  const Vector<T,2>& origin, const Vector<T,2>& u,
  std::vector<int> materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF2D<T>(
      (*densityF) * (*velocityF),
      geometry,
      origin, u,
      std::move(materials),
      mode),
    _velocityF(std::move(velocityF)),
    _densityF(std::move(densityF)),
    _conversationFactorMass(conversationFactorMass),
    _conversationFactorTime(conversationFactorTime)
{
  this->getName() = "SuperPlaneIntegralFluxMass2D";
}

template<typename T>
SuperPlaneIntegralFluxMass2D<T>::SuperPlaneIntegralFluxMass2D(
  FunctorPtr<SuperF2D<T>>&& velocityF,
  FunctorPtr<SuperF2D<T>>&& densityF,
  SuperGeometry2D<T>&       geometry,
  T conversationFactorMass,
  T conversationFactorTime,
  const Vector<T,2>& origin, const Vector<T,2>& u,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF2D<T>(
      (*densityF) * (*velocityF),
      geometry,
      origin, u,
      mode),
    _velocityF(std::move(velocityF)),
    _densityF(std::move(densityF)),
    _conversationFactorMass(conversationFactorMass),
    _conversationFactorTime(conversationFactorTime)
{
  this->getName() = "SuperPlaneIntegralFluxMass2D";
}

template<typename T>
bool SuperPlaneIntegralFluxMass2D<T>::operator()(T output[], const int input[])
{
  const bool result = SuperPlaneIntegralF2D<T>::operator()(output, input);

  output[0] = output[0] * _conversationFactorMass / _conversationFactorTime;

  return result;
}

template<typename T>
void SuperPlaneIntegralFluxMass2D<T>::print(
  std::string regionName, std::string massFluxSiScaleName)
{
  OstreamManager clout("SuperPlaneIntegralFluxMass2D");
  int input[1] = { };
  T output[this->getTargetDim()]; // = 5
  operator()(output, input);
  if (regionName != "") {
    clout << "regionName=" << regionName << "; regionSize[m]=" << output[1]
          << std::flush;
  } else {
    clout << "regionSize[m]=" << output[1] << std::flush;
  }
  if (singleton::mpi().isMainProcessor()) {
    if (massFluxSiScaleName == "mcg/s") { // milli gramm
      std::cout << "; massFlowRate[mcg/s]=" << output[0] * T(1.e6)
                << std::endl;
    } else if (massFluxSiScaleName == "mg/s") { // micro gramm
      std::cout << "; massFlowRate[mg/s]=" << output[0] * T(1.e3)
                << std::endl;
    } else {
      std::cout << "; massFlowRate[kg/s]=" << output[0] << std::endl;
    }
  }
}


}

#endif
