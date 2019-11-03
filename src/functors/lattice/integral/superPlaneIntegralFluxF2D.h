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

#ifndef SUPER_PLANE_INTEGRAL_FLUX_F_2D_H
#define SUPER_PLANE_INTEGRAL_FLUX_F_2D_H

#include "superPlaneIntegralF2D.h"
#include "functors/lattice/superLatticeLocalF2D.h"

namespace olb {


/// Template class for building flux integrals based on SuperLatticePhysF2D functors
/**
 * e.g. phys pressure flux is derived from SuperPlaneIntegralFluxF2D<T,SuperLatticePhysPressure2D> and only adds a print method.
 * All required constructors are provided by this class and need only be inherited.
 *
 * The constructors correspond to SuperPlaneIntegralF2D's with the difference
 * that they accept a super lattice and a unit converter reference instead of
 * a functor.
 * The appropropriate SuperLatticePhysF2D derived functor is then instantiated
 * internally as defined by the template argument F.
 **/
template<typename T, template<typename, typename> class F>
class SuperPlaneIntegralFluxF2D : public SuperPlaneIntegralF2D<T> {
public:
  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF2D(SuperLattice2D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry2D<T>& geometry,
                            const HyperplaneLattice2D<T>& hyperplaneLattice,
                            FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
                            FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF2D(SuperLattice2D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry2D<T>&    geometry,
                            const Hyperplane2D<T>& hyperplane,
                            FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
                            FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF2D(SuperLattice2D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry2D<T>&     geometry,
                            const Hyperplane2D<T>&  hyperplane,
                            FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF2D(SuperLattice2D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry2D<T>& geometry,
                            const Vector<T,2>& origin,
                            const Vector<T,2>& u,
                            std::vector<int> materials,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF2D(SuperLattice2D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry2D<T>& geometry,
                            const Vector<T,2>& origin,
                            const Vector<T,2>& u,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

};

/// Pressure flux line integral
/**
 * Only adds a print method.
 * Calculation is implemented in SuperPlaneIntegralF2D which is constructed around
 * SuperLatticePhysPressure2D by SuperPlaneIntegralFluxF2D.
 **/
template<typename T>
class SuperPlaneIntegralFluxPressure2D final
    : public SuperPlaneIntegralFluxF2D<T, SuperLatticePhysPressure2D> {
public:
  using SuperPlaneIntegralFluxF2D<T, SuperLatticePhysPressure2D>::SuperPlaneIntegralFluxF2D;

  void print(std::string regionName      = "",
             std::string fluxSiScaleName = "N",
             std::string meanSiScaleName = "Pa");
};

/// Velocity flux line integral
/**
 * Only adds a print method.
 * Calculation is implemented in SuperPlaneIntegralF2D which is constructed around
 * SuperLatticePhysVelocity2D by SuperPlaneIntegralFluxF2D.
 **/
template<typename T>
class SuperPlaneIntegralFluxVelocity2D final
    : public SuperPlaneIntegralFluxF2D<T, SuperLatticePhysVelocity2D> {
public:
  using SuperPlaneIntegralFluxF2D<T, SuperLatticePhysVelocity2D>::SuperPlaneIntegralFluxF2D;

  void print(std::string regionName      = "",
             std::string fluxSiScaleName = "m^2/s",
             std::string meanSiScaleName = "m/s");
};


}

#endif
