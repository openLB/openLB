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

#ifndef SUPER_PLANE_INTEGRAL_FLUX_F_2D_HH
#define SUPER_PLANE_INTEGRAL_FLUX_F_2D_HH

#include "superPlaneIntegralFluxF2D.h"

#include "io/ostreamManager.h"
#include "utilities/vectorHelpers.h"
#include "utilities/functorPtr.hh"
#include "functors/lattice/indicator/indicator2D.hh"

namespace olb {


template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF2D<T, F>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry2D<T>&           geometry,
  const HyperplaneLattice2D<T>& hyperplaneLattice,
  FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
  FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF2D<T>(
      std::unique_ptr<SuperF2D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      hyperplaneLattice,
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      std::forward<decltype(subplaneIndicator)>(subplaneIndicator),
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF2D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF2D<T, F>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry2D<T>&    geometry,
  const Hyperplane2D<T>& hyperplane,
  FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
  FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF2D<T>(
      std::unique_ptr<SuperF2D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      hyperplane,
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      std::forward<decltype(subplaneIndicator)>(subplaneIndicator),
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF2D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF2D<T, F>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry2D<T>&    geometry,
  const Hyperplane2D<T>& hyperplane,
  FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF2D<T>(
      std::unique_ptr<SuperF2D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      hyperplane,
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF2D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF2D<T, F>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry2D<T>& geometry,
  const Vector<T,2>& origin, const Vector<T,2>& u,
  std::vector<int> materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF2D<T>(
      std::unique_ptr<SuperF2D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      origin, u,
      std::move(materials),
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF2D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF2D<T, F>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry2D<T>& geometry,
  const Vector<T,2>& origin, const Vector<T,2>& u,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF2D<T>(
      std::unique_ptr<SuperF2D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      origin, u,
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF2D";
}


template<typename T>
void SuperPlaneIntegralFluxPressure2D<T>::print(
  std::string regionName, std::string fluxSiScaleName, std::string meanSiScaleName)
{
  OstreamManager clout("SuperPlaneIntegralFluxPressure2D");
  int input[1] = { };
  T output[this->getTargetDim()];
  this->operator()(output, input);
  if ( regionName != "" ) {
    clout << "regionName=" << regionName << "; regionSize[m]=" << output[1] << std::flush;
  } else {
    clout << "regionSize[m]=" << output[1] << std::flush;
  }
  if ( singleton::mpi().isMainProcessor() ) {
    if ( fluxSiScaleName == "MN" ) {
      std::cout << "; force[MN]=" << output[0]/T(1.e6) << std::flush;
    } else if ( fluxSiScaleName == "kN") {
      std::cout << "; force[kN]=" << output[0]/T(1.e3) << std::flush;
    } else {
      std::cout << "; force[N]=" << output[0] << std::flush;
    }
    if ( meanSiScaleName == "mmHg" ) {
      std::cout << "; meanPressure[mmHg]=" << fabs(output[0])/output[1]/T(133.322) << std::endl;
    } else {
      std::cout << "; meanPressure[Pa]=" << fabs(output[0])/output[1] << std::endl;
    }
  }
}


template<typename T>
void SuperPlaneIntegralFluxVelocity2D<T>::print(
  std::string regionName, std::string fluxSiScaleName, std::string meanSiScaleName)
{
  OstreamManager clout("SuperPlaneIntegralFluxVelocity2D");
  int input[1] = { };
  T output[this->getTargetDim()];
  this->operator()(output, input);
  if ( regionName != "" ) {
    clout << "regionName=" << regionName << "; regionSize[m]=" << output[1] << std::flush;
  } else {
    clout << "regionSize[m]=" << output[1] << std::flush;
  }
  if ( singleton::mpi().isMainProcessor() ) {
    std::cout << "; flowRate[m^2/s]=" << output[0] << std::flush;
    if ( meanSiScaleName == "mm/s" ) {
      std::cout << "; meanVelocity[mm/s]=" << output[0]/output[1]*T(1.e3) << std::endl;
    } else {
      std::cout << "; meanVelocity[m/s]=" << output[0]/output[1] << std::endl;
    }
  }
}


}

#endif
