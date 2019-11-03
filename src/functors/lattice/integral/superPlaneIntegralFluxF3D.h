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

#ifndef SUPER_PLANE_INTEGRAL_FLUX_F_3D_H
#define SUPER_PLANE_INTEGRAL_FLUX_F_3D_H

#include "superPlaneIntegralF3D.h"
#include "functors/lattice/superLatticeLocalF3D.h"

namespace olb {


/// Template class for building flux integrals based on SuperLatticePhysF3D functors
/**
 * e.g. phys pressure flux is derived from SuperPlaneIntegralFluxF3D<T,SuperLatticePhysPressure3D> and only adds a print method.
 * All required constructors are provided by this class and need only be inherited.
 *
 * The constructors correspond to SuperPlaneIntegralF3D's with the difference
 * that they accept a super lattice and a unit converter reference instead of
 * a functor.
 * The appropropriate SuperLatticePhysF3D derived functor is then instantiated
 * internally as defined by the template argument F.
 *
 * See SuperPlaneIntegralF3D for further documentation.
 **/
template<typename T, template<typename, typename> class F>
class SuperPlaneIntegralFluxF3D : public SuperPlaneIntegralF3D<T> {
public:
  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF3D(SuperLattice3D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry3D<T>& geometry,
                            const HyperplaneLattice3D<T>& hyperplaneLattice,
                            FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
                            FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF3D(SuperLattice3D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry3D<T>&    geometry,
                            const Hyperplane3D<T>& hyperplane,
                            FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
                            FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF3D(SuperLattice3D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry3D<T>&     geometry,
                            const Hyperplane3D<T>&  hyperplane,
                            FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF3D(SuperLattice3D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry3D<T>& geometry,
                            const Vector<T,3>& origin,
                            const Vector<T,3>& u, const Vector<T,3>& v,
                            std::vector<int> materials,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF3D(SuperLattice3D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry3D<T>& geometry,
                            const Vector<T,3>& origin,
                            const Vector<T,3>& u, const Vector<T,3>& v,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF3D(SuperLattice3D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry3D<T>& geometry,
                            const Vector<T,3>& origin,
                            const Vector<T,3>& normal,
                            std::vector<int> materials,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF3D(SuperLattice3D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry3D<T>& geometry,
                            const Vector<T,3>& origin,
                            const Vector<T,3>& normal,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF3D(SuperLattice3D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry3D<T>& geometry,
                            const Vector<T,3>& normal,
                            std::vector<int> materials,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF3D(SuperLattice3D<T, DESCRIPTOR>&     sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry3D<T>& geometry,
                            const Vector<T,3>& normal,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry3D<T>& geometry,
                            IndicatorCircle3D<T>& circle,
                            std::vector<int> materials,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  template<typename DESCRIPTOR>
  SuperPlaneIntegralFluxF3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter,
                            SuperGeometry3D<T>& geometry,
                            IndicatorCircle3D<T>& circle,
                            BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

};

/// Pressure flux plane integral
/**
 * Only adds a print method.
 * Calculation is implemented in SuperPlaneIntegralF3D which is constructed around
 * SuperLatticePhysPressure3D by SuperPlaneIntegralFluxF3D.
 **/
template<typename T>
class SuperPlaneIntegralFluxPressure3D final
    : public SuperPlaneIntegralFluxF3D<T, SuperLatticePhysPressure3D> {
public:
  using SuperPlaneIntegralFluxF3D<T, SuperLatticePhysPressure3D>::SuperPlaneIntegralFluxF3D;

  void print(std::string regionName      = "",
             std::string fluxSiScaleName = "N",
             std::string meanSiScaleName = "Pa");
};

/// Velocity flux plane integral
/**
 * Only adds a print method.
 * Calculation is implemented in SuperPlaneIntegralF3D which is constructed around
 * SuperLatticePhysVelocity3D by SuperPlaneIntegralFluxF3D.
 **/
template<typename T>
class SuperPlaneIntegralFluxVelocity3D final
    : public SuperPlaneIntegralFluxF3D<T, SuperLatticePhysVelocity3D> {
public:
  using SuperPlaneIntegralFluxF3D<T, SuperLatticePhysVelocity3D>::SuperPlaneIntegralFluxF3D;

  void print(std::string regionName      = "",
             std::string fluxSiScaleName = "m^3/s",
             std::string meanSiScaleName = "m/s");
};


}

#endif
