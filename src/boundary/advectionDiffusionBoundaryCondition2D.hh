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
 * A helper for initialising 2D boundaries -- generic implementation.
 */
#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_2D_HH
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_2D_HH

#include "advectionDiffusionBoundaries.h"
#include "advectionDiffusionBoundaries.hh"
#include "advectionDiffusionMomentaOnBoundaries.h"
#include "advectionDiffusionMomentaOnBoundaries.hh"
#include "advectionDiffusionBoundaryCondition2D.h"
#include "advectionDiffusionBoundaryInstantiator2D.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"


namespace olb {


////////// AdvectionDiffusionBoundaryManager2D /////////////////////////////////////////

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getTemperatureBoundaryMomenta()
{
  return new EquilibriumBM<T,DESCRIPTOR>;
  // return new RegularizedTemperatureBM<T,DESCRIPTOR,direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new AdvectionDiffusionBoundariesDynamics<T,DESCRIPTOR,MixinDynamics,direction,orientation>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return nullptr;
}

//==================  Corners ================================

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,DESCRIPTOR>*
AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getTemperatureCornerMomenta()
{
  return new EquilibriumBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new AdvectionDiffusionCornerDynamics2D<T,DESCRIPTOR,MixinDynamics,xNormal,yNormal>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,DESCRIPTOR>*
AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getTemperatureCornerProcessor(int x, int y)
{
  return nullptr;
}



template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getRegularizedTemperatureBoundaryMomenta()
{
  return new RegularizedTemperatureBM<T,DESCRIPTOR,direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getRegularizedTemperatureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedAdvectionDiffusionRLBdynamics<T,DESCRIPTOR,MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getRegularizedTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return nullptr;
}

//==================  Corners ================================

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,DESCRIPTOR>*
AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getRegularizedTemperatureCornerMomenta()
{
  return new EquilibriumBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getRegularizedTemperatureCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new AdvectionDiffusionCornerDynamics2D<T,DESCRIPTOR,MixinDynamics,xNormal,yNormal>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,DESCRIPTOR>*
AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getRegularizedTemperatureCornerProcessor(int x, int y)
{
  return nullptr;
}

////////// Convenience wrappers for boundary functions ////////////////////////

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition2D<T, DESCRIPTOR>::addTemperatureBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, T omega)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addTemperatureBoundary(indicator, x0, x1, y0, y1, omega);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition2D<T, DESCRIPTOR>::addTemperatureBoundary(
  BlockIndicatorF2D<T>& indicator, T omega, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addTemperatureBoundary(indicator,
                         margin, blockGeometryStructure.getNx()-1 -margin,
                         margin, blockGeometryStructure.getNy()-1 -margin,
                         omega);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeAdvectionDiffusionBoundaryCondition2D<T, DESCRIPTOR>::addTemperatureBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  T omega, bool includeOuterCells)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addTemperatureBoundary(indicator, omega, includeOuterCells);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getRegularizedHeatFluxBoundaryMomenta(T *heatFlux)
{
  return new RegularizedHeatFluxBM<T,DESCRIPTOR,direction,orientation>(heatFlux);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getRegularizedHeatFluxBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedAdvectionDiffusionRLBdynamics<T,DESCRIPTOR,MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getRegularizedHeatFluxBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return nullptr;
}

//==================  Corners ================================

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,DESCRIPTOR>*
AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getRegularizedHeatFluxCornerMomenta()
{
  return new EquilibriumBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,DESCRIPTOR>* AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getRegularizedHeatFluxCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new AdvectionDiffusionCornerDynamics2D<T,DESCRIPTOR,MixinDynamics,xNormal,yNormal>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,DESCRIPTOR>*
AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getRegularizedHeatFluxCornerProcessor(int x, int y)
{
  return nullptr;
}


////////// Factory functions //////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
OnLatticeAdvectionDiffusionBoundaryCondition2D<T,DESCRIPTOR>*
createAdvectionDiffusionBoundaryCondition2D(BlockLatticeStructure2D<T,DESCRIPTOR>& block)
{
  return new AdvectionDiffusionBoundaryConditionInstantiator2D<T, DESCRIPTOR,
         AdvectionDiffusionBoundaryManager2D<T,DESCRIPTOR, MixinDynamics> > (block);
}

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void createAdvectionDiffusionBoundaryCondition2D(sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>& sBC)
{
  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, DESCRIPTOR>* blockBC =
      createAdvectionDiffusionBoundaryCondition2D<T,DESCRIPTOR,MixinDynamics>(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getADblockBCs().push_back(blockBC);
  }
}


}  // namespace olb

#endif
