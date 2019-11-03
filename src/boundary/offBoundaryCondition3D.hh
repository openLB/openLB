/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Kratzke, Mathias J. Krause
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

/** \file
 * A helper for initialising 3D boundaries -- generic implementation.
 */
#ifndef OFF_BOUNDARY_CONDITION_3D_HH
#define OFF_BOUNDARY_CONDITION_3D_HH

#include "offBoundaryCondition3D.h"
#include "offBoundaryInstantiator3D.h"
#include "offBoundaryPostProcessors3D.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"

namespace olb {

/**
* Boundary Managers provide specific Boundary Processors by creating them
*/

template<typename T, typename DESCRIPTOR, class MixinDynamics>
class BouzidiBoundaryManager3D {
public:

  static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getOnePointZeroVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist);
  static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getTwoPointZeroVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist);
  static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getOnePointVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist);
  static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getTwoPointVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist);
  static Dynamics<T,DESCRIPTOR>*
  getOffDynamics(T location[DESCRIPTOR::d]);
  static Dynamics<T,DESCRIPTOR>*
  getOffDynamics(T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q]);
};

////////// BouzidiBoundaryManager3D /////////////////////////////////////////

template<typename T, typename DESCRIPTOR, class MixinDynamics>
PostProcessorGenerator3D<T,DESCRIPTOR>*
BouzidiBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getOnePointZeroVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist)
{
  return new ZeroVelocityBounceBackPostProcessorGenerator3D
         <T, DESCRIPTOR>(x, y, z, iPop, dist);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
PostProcessorGenerator3D<T,DESCRIPTOR>*
BouzidiBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getTwoPointZeroVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist)
{
  return new ZeroVelocityBouzidiLinearPostProcessorGenerator3D
         <T, DESCRIPTOR>(x, y, z, iPop, dist);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
PostProcessorGenerator3D<T,DESCRIPTOR>*
BouzidiBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getOnePointVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist)
{
  return new VelocityBounceBackPostProcessorGenerator3D
         <T, DESCRIPTOR>(x, y, z, iPop, dist);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
PostProcessorGenerator3D<T,DESCRIPTOR>*
BouzidiBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getTwoPointVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist)
{
  return new VelocityBouzidiLinearPostProcessorGenerator3D
         <T, DESCRIPTOR>(x, y, z, iPop, dist);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
Dynamics<T,DESCRIPTOR>*
BouzidiBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getOffDynamics(T location[DESCRIPTOR::d])
{
  return new OffDynamics<T, DESCRIPTOR>(location);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
Dynamics<T,DESCRIPTOR>*
BouzidiBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getOffDynamics(T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q])
{
  return new OffDynamics<T, DESCRIPTOR>(location, distances);
}

////////// Convenience wrappers for boundary functions ////////////////////////

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition3D<T, DESCRIPTOR>::addOffDynamics(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addOffDynamics(indicator);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition3D<T, DESCRIPTOR>::addZeroVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  IndicatorF3D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial3D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  BlockIndicatorMaterial3D<T> boundaryIndicator(blockGeometryStructure, material);
  addZeroVelocityBoundary(boundaryIndicator, bulkIndicator, geometryIndicator);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition3D<T, DESCRIPTOR>::addVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  IndicatorF3D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial3D<T> boundaryIndicator(blockGeometryStructure, material);
  BlockIndicatorMaterial3D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  addVelocityBoundary(boundaryIndicator, bulkIndicator, geometryIndicator);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition3D<T, DESCRIPTOR>::defineU(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  AnalyticalF3D<T,T>& u, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  BlockIndicatorMaterial3D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  defineU(indicator, bulkIndicator, u);
}

////////// Factory functions //////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
OffLatticeBoundaryCondition3D<T,DESCRIPTOR>*
createBouzidiBoundaryCondition3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block)
{
  return new OffBoundaryConditionInstantiator3D <
         T, DESCRIPTOR,
         BouzidiBoundaryManager3D<T,DESCRIPTOR, MixinDynamics> > (block);
}

}  // namespace olb

#endif
