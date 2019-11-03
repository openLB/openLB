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

#ifndef SUPER_PLANE_INTEGRAL_FLUX_F_3D_HH
#define SUPER_PLANE_INTEGRAL_FLUX_F_3D_HH

#include "superPlaneIntegralFluxF3D.h"

#include "io/ostreamManager.h"
#include "utilities/vectorHelpers.h"
#include "utilities/functorPtr.hh"
#include "functors/lattice/indicator/indicator2D.hh"

namespace olb {


template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF3D<T, F>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry3D<T>&           geometry,
  const HyperplaneLattice3D<T>& hyperplaneLattice,
  FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
  FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF3D<T>(
      std::unique_ptr<SuperF3D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      hyperplaneLattice,
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      std::forward<decltype(subplaneIndicator)>(subplaneIndicator),
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF3D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF3D<T, F>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry3D<T>&    geometry,
  const Hyperplane3D<T>& hyperplane,
  FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
  FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF3D<T>(
      std::unique_ptr<SuperF3D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      hyperplane,
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      std::forward<decltype(subplaneIndicator)>(subplaneIndicator),
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF3D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF3D<T, F>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry3D<T>&    geometry,
  const Hyperplane3D<T>& hyperplane,
  FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF3D<T>(
      std::unique_ptr<SuperF3D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      hyperplane,
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF3D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF3D<T, F>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry3D<T>& geometry,
  const Vector<T,3>& origin, const Vector<T,3>& u, const Vector<T,3>& v,
  std::vector<int> materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D<T>(
      std::unique_ptr<SuperF3D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      origin, u, v,
      std::move(materials),
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF3D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF3D<T, F>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry3D<T>& geometry,
  const Vector<T,3>& origin, const Vector<T,3>& u, const Vector<T,3>& v,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D<T>(
      std::unique_ptr<SuperF3D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      origin, u, v,
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF3D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF3D<T, F>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry3D<T>& geometry,
  const Vector<T,3>& origin, const Vector<T,3>& normal,
  std::vector<int> materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D<T>(
      std::unique_ptr<SuperF3D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      origin, normal,
      std::move(materials),
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF3D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF3D<T, F>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry3D<T>& geometry,
  const Vector<T,3>& origin, const Vector<T,3>& normal,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D<T>(
      std::unique_ptr<SuperF3D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      origin, normal,
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF3D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF3D<T, F>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry3D<T>& geometry,
  const Vector<T,3>& normal,
  std::vector<int> materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D<T>(
      std::unique_ptr<SuperF3D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      normal,
      std::move(materials),
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF3D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF3D<T, F>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry3D<T>& geometry,
  const Vector<T,3>& normal,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D<T>(
      std::unique_ptr<SuperF3D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      normal,
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF3D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF3D<T, F>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry3D<T>&   geometry,
  IndicatorCircle3D<T>& circle,
  std::vector<int> materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D<T>(
      std::unique_ptr<SuperF3D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      circle,
      std::move(materials),
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF3D";
}

template<typename T, template<typename, typename> class F>
template<typename DESCRIPTOR>
SuperPlaneIntegralFluxF3D<T, F>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry3D<T>&   geometry,
  IndicatorCircle3D<T>& circle,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D<T>(
      std::unique_ptr<SuperF3D<T>>(new F<T, DESCRIPTOR>(sLattice, converter)),
      geometry,
      circle,
      mode)
{
  this->getName() = "SuperPlaneIntegralFluxF3D";
}

template<typename T>
void SuperPlaneIntegralFluxPressure3D<T>::print(
  std::string regionName, std::string fluxSiScaleName, std::string meanSiScaleName)
{
  OstreamManager clout("SuperPlaneIntegralFluxPressure3D");
  int input[1] = { };
  T output[this->getTargetDim()];
  this->operator()(output, input);
  if ( regionName != "" ) {
    clout << "regionName=" << regionName << "; regionSize[m^2]=" << output[1]
          << std::flush;
  }
  else {
    clout << "regionSize[m^2]=" << output[1] << std::flush;
  }
  if ( singleton::mpi().isMainProcessor() ) {
    if ( fluxSiScaleName == "MN" ) {
      std::cout << "; force[MN]=" << output[0] / T(1.e6) << std::flush;
    }
    else if ( fluxSiScaleName == "kN" ) {
      std::cout << "; force[kN]=" << output[0] / T(1.e3) << std::flush;
    }
    else {
      std::cout << "; force[N]=" << output[0] << std::flush;
    }
    if ( meanSiScaleName == "mmHg" ) {
      std::cout << "; meanPressure[mmHg]=" << std::abs(output[0])/output[1]/T(133.322) << std::endl;
    }
    else {
      std::cout << "; meanPressure[Pa]=" << std::abs(output[0])/output[1] << std::endl;
    }
  }
}


template<typename T>
void SuperPlaneIntegralFluxVelocity3D<T>::print(
  std::string regionName, std::string fluxSiScaleName, std::string meanSiScaleName)
{
  OstreamManager clout("SuperPlaneIntegralFluxVelocity3D");
  int input[1] = { };
  T output[this->getTargetDim()];
  this->operator()(output, input);
  if ( regionName != "" ) {
    clout << "regionName=" << regionName << "; regionSize[m^2]=" << output[1]
          << std::flush;
  }
  else {
    clout << "regionSize[m^2]=" << output[1] << std::flush;
  }
  if ( singleton::mpi().isMainProcessor() ) {
    if ( fluxSiScaleName == "ml/s" ) {
      std::cout << "; volumetricFlowRate[ml/s]=" << output[0] * T(1.e6)
                << std::flush;
    }
    else if ( fluxSiScaleName == "l/s" ) {
      std::cout << "; volumetricFlowRate[l/s]=" << output[0] * T(1.e3)
                << std::flush;
    }
    else {
      std::cout << "; volumetricFlowRate[m^3/s]=" << output[0] << std::flush;
    }
    if ( meanSiScaleName == "mm/s" ) {
      std::cout << "; meanVelocity[mm/s]=" << output[0] / output[1] * T(1.e3)
                << std::endl;
    }
    else {
      std::cout << "; meanVelocity[m/s]=" << output[0] / output[1] << std::endl;
    }
  }
}


}

#endif
