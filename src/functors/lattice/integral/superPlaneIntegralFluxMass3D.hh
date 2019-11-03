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

#ifndef SUPER_PLANE_INTEGRAL_FLUX_MASS_3D_HH
#define SUPER_PLANE_INTEGRAL_FLUX_MASS_3D_HH

#include "superPlaneIntegralFluxMass3D.h"

#include "io/ostreamManager.h"
#include "utilities/vectorHelpers.h"
#include "utilities/functorPtr.hh"
#include "functors/lattice/indicator/indicator2D.hh"

#include "functors/lattice/superLatticeLocalF3D.h"
#include "functors/lattice/superCalcF3D.h"
#include "functors/lattice/superCalcF3D.hh"

namespace olb {


template<typename T>
SuperPlaneIntegralFluxMass3D<T>::SuperPlaneIntegralFluxMass3D(
  FunctorPtr<SuperF3D<T>>&& velocityF,
  FunctorPtr<SuperF3D<T>>&& densityF,
  SuperGeometry3D<T>&       geometry,
  T conversationFactorMass,
  T conversationFactorTime,
  const HyperplaneLattice3D<T>& hyperplaneLattice,
  FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
  FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF3D<T>(
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
  this->getName() = "SuperPlaneIntegralFluxMass3D";
}

template<typename T>
SuperPlaneIntegralFluxMass3D<T>::SuperPlaneIntegralFluxMass3D(
  FunctorPtr<SuperF3D<T>>&& velocityF,
  FunctorPtr<SuperF3D<T>>&& densityF,
  SuperGeometry3D<T>&       geometry,
  T conversationFactorMass,
  T conversationFactorTime,
  const Hyperplane3D<T>& hyperplane,
  FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
  FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF3D<T>(
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
  this->getName() = "SuperPlaneIntegralFluxMass3D";
}

template<typename T>
SuperPlaneIntegralFluxMass3D<T>::SuperPlaneIntegralFluxMass3D(
  FunctorPtr<SuperF3D<T>>&& velocityF,
  FunctorPtr<SuperF3D<T>>&& densityF,
  SuperGeometry3D<T>&       geometry,
  T conversationFactorMass,
  T conversationFactorTime,
  const Hyperplane3D<T>& hyperplane,
  FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF3D<T>(
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
  this->getName() = "SuperPlaneIntegralFluxMass3D";
}

template<typename T>
SuperPlaneIntegralFluxMass3D<T>::SuperPlaneIntegralFluxMass3D(
  FunctorPtr<SuperF3D<T>>&& velocityF,
  FunctorPtr<SuperF3D<T>>&& densityF,
  SuperGeometry3D<T>&       geometry,
  T conversationFactorMass,
  T conversationFactorTime,
  const Vector<T,3>& origin, const Vector<T,3>& u, const Vector<T,3>& v,
  std::vector<int> materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D<T>(
      (*densityF) * (*velocityF),
      geometry,
      origin, u, v,
      std::move(materials),
      mode),
    _velocityF(std::move(velocityF)),
    _densityF(std::move(densityF)),
    _conversationFactorMass(conversationFactorMass),
    _conversationFactorTime(conversationFactorTime)
{
  this->getName() = "SuperPlaneIntegralFluxMass3D";
}

template<typename T>
SuperPlaneIntegralFluxMass3D<T>::SuperPlaneIntegralFluxMass3D(
  FunctorPtr<SuperF3D<T>>&& velocityF,
  FunctorPtr<SuperF3D<T>>&& densityF,
  SuperGeometry3D<T>&       geometry,
  T conversationFactorMass,
  T conversationFactorTime,
  const Vector<T,3>& origin, const Vector<T,3>& u, const Vector<T,3>& v,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D<T>(
      (*densityF) * (*velocityF),
      geometry,
      origin, u, v,
      mode),
    _velocityF(std::move(velocityF)),
    _densityF(std::move(densityF)),
    _conversationFactorMass(conversationFactorMass),
    _conversationFactorTime(conversationFactorTime)
{
  this->getName() = "SuperPlaneIntegralFluxMass3D";
}

template<typename T>
SuperPlaneIntegralFluxMass3D<T>::SuperPlaneIntegralFluxMass3D(
  FunctorPtr<SuperF3D<T>>&& velocityF,
  FunctorPtr<SuperF3D<T>>&& densityF,
  SuperGeometry3D<T>&       geometry,
  T conversationFactorMass,
  T conversationFactorTime,
  IndicatorCircle3D<T>& circle,
  std::vector<int> materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D<T>(
      (*densityF) * (*velocityF),
      geometry,
      circle,
      std::move(materials),
      mode),
    _velocityF(std::move(velocityF)),
    _densityF(std::move(densityF)),
    _conversationFactorMass(conversationFactorMass),
    _conversationFactorTime(conversationFactorTime)
{
  this->getName() = "SuperPlaneIntegralFluxMass3D";
}

template<typename T>
SuperPlaneIntegralFluxMass3D<T>::SuperPlaneIntegralFluxMass3D(
  FunctorPtr<SuperF3D<T>>&& velocityF,
  FunctorPtr<SuperF3D<T>>&& densityF,
  SuperGeometry3D<T>&       geometry,
  T conversationFactorMass,
  T conversationFactorTime,
  IndicatorCircle3D<T>& circle,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D<T>(
      (*densityF) * (*velocityF),
      geometry,
      circle,
      mode),
    _velocityF(std::move(velocityF)),
    _densityF(std::move(densityF)),
    _conversationFactorMass(conversationFactorMass),
    _conversationFactorTime(conversationFactorTime)
{
  this->getName() = "SuperPlaneIntegralFluxMass3D";
}

template<typename T>
bool SuperPlaneIntegralFluxMass3D<T>::operator()(T output[], const int input[])
{
  const bool result = SuperPlaneIntegralF3D<T>::operator()(output, input);

  output[0] = output[0] * _conversationFactorMass / _conversationFactorTime;

  return result;
}

template<typename T>
void SuperPlaneIntegralFluxMass3D<T>::print(
  std::string regionName, std::string massFluxSiScaleName)
{
  OstreamManager clout("SuperPlaneIntegralFluxMass3D");
  int input[1] = { };
  T output[this->getTargetDim()]; // = 5
  operator()(output, input);
  if (regionName != "") {
    clout << "regionName=" << regionName << "; regionSize[m^2]=" << output[1]
          << std::flush;
  } else {
    clout << "regionSize[m^2]=" << output[1] << std::flush;
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
