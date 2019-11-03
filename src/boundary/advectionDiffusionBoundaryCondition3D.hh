/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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
#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_HH
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_HH

#include "advectionDiffusionBoundaries.h"
#include "advectionDiffusionBoundaries.hh"
#include "advectionDiffusionBoundaryCondition3D.h"
#include "advectionDiffusionBoundaryInstantiator3D.h"
#include "advectionDiffusionBoundaryInstantiator3D.hh"
#include "momentaOnBoundaries.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "advectionDiffusionMomentaOnBoundaries.h"
#include "advectionDiffusionMomentaOnBoundaries.hh"


namespace olb {

template<typename T, typename DESCRIPTOR, class MixinDynamics>
class AdvectionDiffusionBoundaryManager3D {
public:
  template<int direction, int orientation>
  static Momenta<T,DESCRIPTOR>* getTemperatureBoundaryMomenta();
  template<int direction, int orientation>
  static Dynamics<T,DESCRIPTOR>* getTemperatureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation>
  static PostProcessorGenerator3D<T,DESCRIPTOR>* getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int plane, int normal1, int normal2>
  static Momenta<T,DESCRIPTOR>* getTemperatureBoundaryEdgeMomenta();
  template<int plane, int normal1, int normal2>
  static Dynamics<T,DESCRIPTOR>* getTemperatureBoundaryEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int plane, int normal1, int normal2>
  static PostProcessorGenerator3D<T,DESCRIPTOR>* getTemperatureBoundaryEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int normalX, int normalY, int normalZ>
  static Momenta<T,DESCRIPTOR>* getTemperatureBoundaryCornerMomenta();
  template<int normalX, int normalY, int normalZ>
  static Dynamics<T,DESCRIPTOR>* getTemperatureBoundaryCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int normalX, int normalY, int normalZ>
  static PostProcessorGenerator3D<T,DESCRIPTOR>* getTemperatureBoundaryCornerProcessor(int x, int y, int z);

};

//==================  Flat ================================
template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureBoundaryMomenta()
{
  return new EquilibriumBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new AdvectionDiffusionBoundariesDynamics<T,DESCRIPTOR,MixinDynamics,direction,orientation>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

//==================  Edges ================================
template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureBoundaryEdgeMomenta()
{
  return new EquilibriumBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureBoundaryEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new AdvectionDiffusionEdgesDynamics<T,DESCRIPTOR,MixinDynamics,plane,normal1,normal2>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureBoundaryEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

//==================  Corners ================================
template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureBoundaryCornerMomenta()
{
  return new EquilibriumBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureBoundaryCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new AdvectionDiffusionCornerDynamics3D<T,DESCRIPTOR,MixinDynamics,xNormal,yNormal,zNormal>(omega, momenta);
}



template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureBoundaryCornerProcessor(int x, int y, int z)
{
  return nullptr;
}

//================== Convenience wrappers for boundary functions ======

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>::addTemperatureBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addTemperatureBoundary(indicator, x0, x1, y0, y1, z0, z1, omega);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>::addTemperatureBoundary(
  BlockIndicatorF3D<T>& indicator, T omega, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addTemperatureBoundary(indicator,
                         margin, blockGeometryStructure.getNx()-1 -margin,
                         margin, blockGeometryStructure.getNy()-1 -margin,
                         margin, blockGeometryStructure.getNz()-1 -margin,
                         omega);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>::addTemperatureBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material, T omega, bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addTemperatureBoundary(indicator, omega, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>::addConvectionBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, int z0, int z1)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addConvectionBoundary(indicator, x0, x1, y0, y1, z0, z1);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>::addConvectionBoundary(
  BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addConvectionBoundary(indicator,
                        margin, blockGeometryStructure.getNx()-1 -margin,
                        margin, blockGeometryStructure.getNy()-1 -margin,
                        margin, blockGeometryStructure.getNz()-1 -margin);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>::addConvectionBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material, bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addConvectionBoundary(indicator, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>::addExtFieldBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int offset, int x0, int x1, int y0, int y1, int z0, int z1)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addExtFieldBoundary(indicator, offset, x0, x1, y0, y1, z0, z1);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>::addExtFieldBoundary(
  BlockIndicatorF3D<T>& indicator, int offset, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addExtFieldBoundary(indicator, offset,
                      margin, blockGeometryStructure.getNx()-1 -margin,
                      margin, blockGeometryStructure.getNy()-1 -margin,
                      margin, blockGeometryStructure.getNz()-1 -margin);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>::addExtFieldBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int offset, bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addExtFieldBoundary(indicator, offset, includeOuterCells);
}


template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>::addZeroDistributionBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, int z0, int z1)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addZeroDistributionBoundary(indicator, x0, x1, y0, y1, z0, z1);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>::addZeroDistributionBoundary(
  BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addZeroDistributionBoundary(indicator,
                              margin, blockGeometryStructure.getNx()-1 -margin,
                              margin, blockGeometryStructure.getNy()-1 -margin,
                              margin, blockGeometryStructure.getNz()-1 -margin);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>::addZeroDistributionBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material, bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addZeroDistributionBoundary(indicator, includeOuterCells);
}


//================== Factory functions ================================
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR>*
createAdvectionDiffusionBoundaryCondition3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block)
{
  return new AdvectionDiffusionBoundaryConditionInstantiator3D<T, DESCRIPTOR,
         AdvectionDiffusionBoundaryManager3D<T,DESCRIPTOR, MixinDynamics> > (block);
}

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void createAdvectionDiffusionBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& sBC)
{
  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>* blockBC =
      createAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR,MixinDynamics>(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getADblockBCs().push_back(blockBC);
  }
}


}  // namespace olb

#endif
