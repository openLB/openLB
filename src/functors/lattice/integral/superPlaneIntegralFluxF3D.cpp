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

#include "superPlaneIntegralFluxF3D.h"
#include "superPlaneIntegralFluxF3D.hh"
#include "dynamics/latticeDescriptors.h"

namespace olb {

template class SuperPlaneIntegralFluxF3D<double, SuperLatticePhysPressure3D>;
template class SuperPlaneIntegralFluxF3D<double, SuperLatticePhysVelocity3D>;

template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysPressure3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const HyperplaneLattice3D<double>&,
  FunctorPtr<SuperIndicatorF3D<double>>&&,
  FunctorPtr<IndicatorF2D<double>>&&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysPressure3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Hyperplane3D<double>&,
  FunctorPtr<SuperIndicatorF3D<double>>&&,
  FunctorPtr<IndicatorF2D<double>>&&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysPressure3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Hyperplane3D<double>&,
  FunctorPtr<SuperIndicatorF3D<double>>&&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysPressure3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Vector<double,3>&, const Vector<double,3>&, const Vector<double,3>&,
  std::vector<int>,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysPressure3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Vector<double,3>&, const Vector<double,3>&, const Vector<double,3>&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysPressure3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Vector<double,3>&, const Vector<double,3>&,
  std::vector<int>,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysPressure3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Vector<double,3>&, const Vector<double,3>&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysPressure3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Vector<double,3>&,
  std::vector<int>,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysPressure3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Vector<double,3>&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysPressure3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  IndicatorCircle3D<double>&,
  std::vector<int>,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysPressure3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  IndicatorCircle3D<double>&,
  BlockDataReductionMode);

template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysVelocity3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const HyperplaneLattice3D<double>&,
  FunctorPtr<SuperIndicatorF3D<double>>&&,
  FunctorPtr<IndicatorF2D<double>>&&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysVelocity3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Hyperplane3D<double>&,
  FunctorPtr<SuperIndicatorF3D<double>>&&,
  FunctorPtr<IndicatorF2D<double>>&&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysVelocity3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Hyperplane3D<double>&,
  FunctorPtr<SuperIndicatorF3D<double>>&&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysVelocity3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Vector<double,3>&, const Vector<double,3>&, const Vector<double,3>&,
  std::vector<int>,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysVelocity3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Vector<double,3>&, const Vector<double,3>&, const Vector<double,3>&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysVelocity3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Vector<double,3>&, const Vector<double,3>&,
  std::vector<int>,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysVelocity3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Vector<double,3>&, const Vector<double,3>&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysVelocity3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Vector<double,3>&,
  std::vector<int>,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysVelocity3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  const Vector<double,3>&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysVelocity3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  IndicatorCircle3D<double>&,
  std::vector<int>,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF3D<double,SuperLatticePhysVelocity3D>::SuperPlaneIntegralFluxF3D(
  SuperLattice3D<double,descriptors::D3Q19<>>&, const UnitConverter<double,descriptors::D3Q19<>>&, SuperGeometry3D<double>&,
  IndicatorCircle3D<double>&,
  BlockDataReductionMode);

template class SuperPlaneIntegralFluxPressure3D<double>;
template class SuperPlaneIntegralFluxVelocity3D<double>;

}
