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

#include "superPlaneIntegralFluxF2D.h"
#include "superPlaneIntegralFluxF2D.hh"
#include "dynamics/latticeDescriptors.h"

namespace olb {

template class SuperPlaneIntegralFluxF2D<double, SuperLatticePhysPressure2D>;
template class SuperPlaneIntegralFluxF2D<double, SuperLatticePhysVelocity2D>;

template SuperPlaneIntegralFluxF2D<double,SuperLatticePhysPressure2D>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<double,descriptors::D2Q9<>>&, const UnitConverter<double,descriptors::D2Q9<>>&, SuperGeometry2D<double>&,
  const HyperplaneLattice2D<double>&,
  FunctorPtr<SuperIndicatorF2D<double>>&&,
  FunctorPtr<IndicatorF2D<double>>&&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF2D<double,SuperLatticePhysPressure2D>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<double,descriptors::D2Q9<>>&, const UnitConverter<double,descriptors::D2Q9<>>&, SuperGeometry2D<double>&,
  const Hyperplane2D<double>&,
  FunctorPtr<SuperIndicatorF2D<double>>&&,
  FunctorPtr<IndicatorF2D<double>>&&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF2D<double,SuperLatticePhysPressure2D>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<double,descriptors::D2Q9<>>&, const UnitConverter<double,descriptors::D2Q9<>>&, SuperGeometry2D<double>&,
  const Hyperplane2D<double>&,
  FunctorPtr<SuperIndicatorF2D<double>>&&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF2D<double,SuperLatticePhysPressure2D>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<double,descriptors::D2Q9<>>&, const UnitConverter<double,descriptors::D2Q9<>>&, SuperGeometry2D<double>&,
  const Vector<double,2>&, const Vector<double,2>&,
  std::vector<int>,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF2D<double,SuperLatticePhysPressure2D>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<double,descriptors::D2Q9<>>&, const UnitConverter<double,descriptors::D2Q9<>>&, SuperGeometry2D<double>&,
  const Vector<double,2>&, const Vector<double,2>&,
  BlockDataReductionMode);

template SuperPlaneIntegralFluxF2D<double,SuperLatticePhysVelocity2D>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<double,descriptors::D2Q9<>>&, const UnitConverter<double,descriptors::D2Q9<>>&, SuperGeometry2D<double>&,
  const HyperplaneLattice2D<double>&,
  FunctorPtr<SuperIndicatorF2D<double>>&&,
  FunctorPtr<IndicatorF2D<double>>&&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF2D<double,SuperLatticePhysVelocity2D>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<double,descriptors::D2Q9<>>&, const UnitConverter<double,descriptors::D2Q9<>>&, SuperGeometry2D<double>&,
  const Hyperplane2D<double>&,
  FunctorPtr<SuperIndicatorF2D<double>>&&,
  FunctorPtr<IndicatorF2D<double>>&&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF2D<double,SuperLatticePhysVelocity2D>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<double,descriptors::D2Q9<>>&, const UnitConverter<double,descriptors::D2Q9<>>&, SuperGeometry2D<double>&,
  const Hyperplane2D<double>&,
  FunctorPtr<SuperIndicatorF2D<double>>&&,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF2D<double,SuperLatticePhysVelocity2D>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<double,descriptors::D2Q9<>>&, const UnitConverter<double,descriptors::D2Q9<>>&, SuperGeometry2D<double>&,
  const Vector<double,2>&, const Vector<double,2>&,
  std::vector<int>,
  BlockDataReductionMode);
template SuperPlaneIntegralFluxF2D<double,SuperLatticePhysVelocity2D>::SuperPlaneIntegralFluxF2D(
  SuperLattice2D<double,descriptors::D2Q9<>>&, const UnitConverter<double,descriptors::D2Q9<>>&, SuperGeometry2D<double>&,
  const Vector<double,2>&, const Vector<double,2>&,
  BlockDataReductionMode);

template class SuperPlaneIntegralFluxPressure2D<double>;
template class SuperPlaneIntegralFluxVelocity2D<double>;

}
