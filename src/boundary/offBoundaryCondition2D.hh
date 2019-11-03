/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2016 Jonas Kratzke, Mathias J. Krause
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
 * A helper for initialising 2D boundaries -- generic implementation.
 */
#ifndef OFF_BOUNDARY_CONDITION_2D_HH
#define OFF_BOUNDARY_CONDITION_2D_HH

#include "offBoundaryCondition2D.h"
#include "offBoundaryInstantiator2D.h"
#include "offBoundaryPostProcessors2D.h"

namespace olb {

/**
* Boundary Managers provide specific Boundary Processors by creating them
*/

////////// BouzidiBoundaryManager2D /////////////////////////////////////////

template<typename T, typename DESCRIPTOR, class MixinDynamics>
class BouzidiBoundaryManager2D {
public:

  static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getOnePointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getTwoPointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getOnePointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getTwoPointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static Dynamics<T,DESCRIPTOR>*
  getOffDynamics(T location[DESCRIPTOR::d]);
  static Dynamics<T,DESCRIPTOR>*
  getOffDynamics(T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q]);
};


template<typename T, typename DESCRIPTOR, class MixinDynamics>
PostProcessorGenerator2D<T,DESCRIPTOR>*
BouzidiBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getOnePointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new ZeroVelocityBounceBackPostProcessorGenerator2D
         <T, DESCRIPTOR>(iX, iY, iPop, dist);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
PostProcessorGenerator2D<T,DESCRIPTOR>*
BouzidiBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getTwoPointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new ZeroVelocityBouzidiLinearPostProcessorGenerator2D
         <T, DESCRIPTOR>(iX, iY, iPop, dist);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
PostProcessorGenerator2D<T,DESCRIPTOR>*
BouzidiBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getOnePointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new VelocityBounceBackPostProcessorGenerator2D
         <T, DESCRIPTOR>(iX, iY, iPop, dist);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
PostProcessorGenerator2D<T,DESCRIPTOR>*
BouzidiBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getTwoPointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new VelocityBouzidiLinearPostProcessorGenerator2D
         <T, DESCRIPTOR>(iX, iY, iPop, dist);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
Dynamics<T,DESCRIPTOR>*
BouzidiBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getOffDynamics(T location[DESCRIPTOR::d])
{
  return new OffDynamics<T, DESCRIPTOR>(location);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
Dynamics<T,DESCRIPTOR>*
BouzidiBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getOffDynamics(T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q])
{
  return new OffDynamics<T, DESCRIPTOR>(location, distances);
}

////////// BounceBackBoundaryManager2D /////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
class BounceBackBoundaryManager2D {
public:

  static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getOnePointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getTwoPointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getOnePointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getTwoPointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static Dynamics<T,DESCRIPTOR>*
  getOffDynamics(T location[DESCRIPTOR::d]);
  static Dynamics<T,DESCRIPTOR>*
  getOffDynamics(T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q]);
};

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
BounceBackBoundaryManager2D<T,DESCRIPTOR>::
getOnePointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new ZeroVelocityBounceBackPostProcessorGenerator2D
         <T, DESCRIPTOR>(iX, iY, iPop, dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
BounceBackBoundaryManager2D<T,DESCRIPTOR>::
getTwoPointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new ZeroVelocityBouzidiLinearPostProcessorGenerator2D
         <T, DESCRIPTOR>(iX, iY, iPop, dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
BounceBackBoundaryManager2D<T,DESCRIPTOR>::
getOnePointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new VelocityBounceBackPostProcessorGenerator2D
         <T, DESCRIPTOR>(iX, iY, iPop, dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
BounceBackBoundaryManager2D<T,DESCRIPTOR>::
getTwoPointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new VelocityBouzidiLinearPostProcessorGenerator2D
         <T, DESCRIPTOR>(iX, iY, iPop, dist);
}

template<typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR>*
BounceBackBoundaryManager2D<T,DESCRIPTOR>::
getOffDynamics(T location[DESCRIPTOR::d])
{
  return new OffDynamics<T, DESCRIPTOR>(location);
}

template<typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR>*
BounceBackBoundaryManager2D<T,DESCRIPTOR>::
getOffDynamics(T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q])
{
  return new OffDynamics<T, DESCRIPTOR>(location, distances);
}

////////// Convenience wrappers for boundary functions ////////////////////////

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition2D<T, DESCRIPTOR>::addOffDynamics(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addOffDynamics(indicator);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition2D<T, DESCRIPTOR>::addZeroVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY,
  IndicatorF2D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial2D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  addZeroVelocityBoundary(blockGeometryStructure, iX, iY,
                          bulkIndicator, geometryIndicator);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition2D<T, DESCRIPTOR>::addZeroVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  IndicatorF2D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial2D<T> boundaryIndicator(blockGeometryStructure, material);
  BlockIndicatorMaterial2D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  addZeroVelocityBoundary(boundaryIndicator,
                          bulkIndicator,
                          geometryIndicator);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition2D<T, DESCRIPTOR>::addZeroVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial2D<T> boundaryIndicator(blockGeometryStructure, material);
  BlockIndicatorMaterial2D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  addZeroVelocityBoundary(boundaryIndicator, bulkIndicator);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition2D<T, DESCRIPTOR>::addVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY,
  IndicatorF2D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial2D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  addVelocityBoundary(blockGeometryStructure, iX, iY,
                      bulkIndicator, geometryIndicator);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition2D<T, DESCRIPTOR>::addVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial2D<T> boundaryIndicator(blockGeometryStructure, material);
  BlockIndicatorMaterial2D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  addVelocityBoundary(boundaryIndicator, bulkIndicator, geometryIndicator);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition2D<T, DESCRIPTOR>::addVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial2D<T> boundaryIndicator(blockGeometryStructure, material);
  BlockIndicatorMaterial2D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  addVelocityBoundary(boundaryIndicator, bulkIndicator);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPressureBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  IndicatorF2D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial2D<T> boundaryIndicator(blockGeometryStructure, material);
  BlockIndicatorMaterial2D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  addPressureBoundary(boundaryIndicator, bulkIndicator, geometryIndicator);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPressureBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial2D<T> boundaryIndicator(blockGeometryStructure, material);
  BlockIndicatorMaterial2D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  addPressureBoundary(boundaryIndicator, bulkIndicator);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition2D<T, DESCRIPTOR>::
defineU(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
        AnalyticalF2D<T,T>& u, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  BlockIndicatorMaterial2D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  defineU(indicator, bulkIndicator, u);
}

template<typename T, typename DESCRIPTOR>
void OffLatticeBoundaryCondition2D<T, DESCRIPTOR>::
defineRho(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
          AnalyticalF2D<T,T>& rho, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  BlockIndicatorMaterial2D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  defineRho(indicator, bulkIndicator, rho);
}



////////// Factory functions //////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
OffLatticeBoundaryCondition2D<T,DESCRIPTOR>*
createBouzidiBoundaryCondition2D(BlockLatticeStructure2D<T,DESCRIPTOR>& block)
{
  return new OffBoundaryConditionInstantiator2D <
         T, DESCRIPTOR,
         BouzidiBoundaryManager2D<T,DESCRIPTOR, MixinDynamics> > (block);
}


}  // namespace olb

#endif
