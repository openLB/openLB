/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Albert Mink
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

#ifndef RTLBM_BOUNDARY_HH
#define RTLBM_BOUNDARY_HH

#include "rtlbmBoundary.h"
#include "rtlbmBoundaryDynamics.h"
#include "momentaOnBoundaries.h"
#include "advectionDiffusionBoundaryInstantiator3D.h" // detects orientation of boundary
#include "advectionDiffusionBoundaries.h"

namespace olb {


template<typename T, typename DESCRIPTOR>
class RtlbmDiffuseBoundaryManager3D {
public:
  template<int direction, int orientation>
  static Momenta<T,DESCRIPTOR>* getTemperatureBoundaryMomenta()
  {
    return new EquilibriumBM<T,DESCRIPTOR>;
  }
  template<int direction, int orientation>
  static Dynamics<T,DESCRIPTOR>* getTemperatureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
  {
    return new RtlbmDiffuseBoundaryDynamics<T,DESCRIPTOR,direction,orientation>(omega, momenta);
  }
  template<int direction, int orientation>
  static PostProcessorGenerator3D<T,DESCRIPTOR>* getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
  {
    return nullptr;
  }

  template<int plane, int normal1, int normal2>
  static Momenta<T,DESCRIPTOR>* getTemperatureBoundaryEdgeMomenta()
  {
    return new EquilibriumBM<T,DESCRIPTOR>;
  }
  template<int plane, int normal1, int normal2>
  static Dynamics<T,DESCRIPTOR>* getTemperatureBoundaryEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
  {
    return new RtlbmDiffuseEdgeBoundaryDynamics<T,DESCRIPTOR,plane,normal1,normal2>(omega, momenta);
  }
  template<int plane, int normal1, int normal2>
  static PostProcessorGenerator3D<T,DESCRIPTOR>* getTemperatureBoundaryEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
  {
    return nullptr;
  }

  template<int normalX, int normalY, int normalZ>
  static Momenta<T,DESCRIPTOR>* getTemperatureBoundaryCornerMomenta()
  {
    return new EquilibriumBM<T,DESCRIPTOR>;
  }
  template<int normalX, int normalY, int normalZ>
  static Dynamics<T,DESCRIPTOR>* getTemperatureBoundaryCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
  {
    return new RtlbmDiffuseCornerBoundaryDynamics<T,DESCRIPTOR,normalX,normalY,normalZ>(omega, momenta);
  }
  template<int normalX, int normalY, int normalZ>
  static PostProcessorGenerator3D<T,DESCRIPTOR>* getTemperatureBoundaryCornerProcessor(int x, int y, int z)
  {
    return nullptr;
  }
};


//================== creator function ================================
// blockLattice creator
template<typename T, typename DESCRIPTOR>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR>* createRtlbmDiffuseBoundaryCondition3D( BlockLatticeStructure3D<T,DESCRIPTOR>& block )
{
  return new AdvectionDiffusionBoundaryConditionInstantiator3D<T, DESCRIPTOR, RtlbmDiffuseBoundaryManager3D<T,DESCRIPTOR> > (block);
}

// superLattice creator
template<typename T, typename DESCRIPTOR>
void createRtlbmDiffuseBoundaryCondition3D( sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sBC )
{
  sBC.setOverlap(1);
  for (int iC = 0; iC < sBC.getSuperLattice().getLoadBalancer().size(); iC++) {
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR>* blockBC =
      createRtlbmDiffuseBoundaryCondition3D<T,DESCRIPTOR>( sBC.getSuperLattice().getExtendedBlockLattice(iC) );
    sBC.getADblockBCs().push_back(blockBC);
  }
}


template<typename T, typename DESCRIPTOR>
class RtlbmDiffuseConstBoundaryManager3D {
public:
  template<int direction, int orientation>
  static Momenta<T,DESCRIPTOR>* getTemperatureBoundaryMomenta()
  {
    return new EquilibriumBM<T,DESCRIPTOR>;
  }
  template<int direction, int orientation>
  static Dynamics<T,DESCRIPTOR>* getTemperatureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
  {
    return new RtlbmDiffuseConstBoundaryDynamics<T,DESCRIPTOR,direction,orientation>(omega, momenta);
  }
  template<int direction, int orientation>
  static PostProcessorGenerator3D<T,DESCRIPTOR>* getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
  {
    return nullptr;
  }

  template<int plane, int normal1, int normal2>
  static Momenta<T,DESCRIPTOR>* getTemperatureBoundaryEdgeMomenta()
  {
    return new EquilibriumBM<T,DESCRIPTOR>;
  }
  template<int plane, int normal1, int normal2>
  static Dynamics<T,DESCRIPTOR>* getTemperatureBoundaryEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
  {
    return new RtlbmDiffuseConstEdgeBoundaryDynamics<T,DESCRIPTOR,plane,normal1,normal2>(omega, momenta);
  }
  template<int plane, int normal1, int normal2>
  static PostProcessorGenerator3D<T,DESCRIPTOR>* getTemperatureBoundaryEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
  {
    return nullptr;
  }

  template<int normalX, int normalY, int normalZ>
  static Momenta<T,DESCRIPTOR>* getTemperatureBoundaryCornerMomenta()
  {
    return new EquilibriumBM<T,DESCRIPTOR>;
  }
  template<int normalX, int normalY, int normalZ>
  static Dynamics<T,DESCRIPTOR>* getTemperatureBoundaryCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
  {
    return new RtlbmDiffuseConstCornerBoundaryDynamics<T,DESCRIPTOR,normalX,normalY,normalZ>(omega, momenta);
  }
  template<int normalX, int normalY, int normalZ>
  static PostProcessorGenerator3D<T,DESCRIPTOR>* getTemperatureBoundaryCornerProcessor(int x, int y, int z)
  {
    return nullptr;
  }
};

//================== creator function ================================
// blockLattice creator
template<typename T, typename DESCRIPTOR>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR>* createRtlbmDiffuseConstBoundaryCondition3D( BlockLatticeStructure3D<T,DESCRIPTOR>& block )
{
  return new AdvectionDiffusionBoundaryConditionInstantiator3D<T, DESCRIPTOR, RtlbmDiffuseConstBoundaryManager3D<T,DESCRIPTOR> > (block);  // TODO AM mark as placeholder
}

// superLattice creator
template<typename T, typename DESCRIPTOR>
void createRtlbmDiffuseConstBoundaryCondition3D( sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sBC )
{
  sBC.setOverlap(1);
  for (int iC = 0; iC < sBC.getSuperLattice().getLoadBalancer().size(); iC++) {
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR>* blockBC =
      createRtlbmDiffuseConstBoundaryCondition3D<T,DESCRIPTOR>( sBC.getSuperLattice().getExtendedBlockLattice(iC) );
    sBC.getADblockBCs().push_back(blockBC);
  }
}


template<typename T, typename DESCRIPTOR>
class RtlbmDirectedBoundaryManager3D {
public:
  // ====== flats
  template<int direction, int orientation>
  static Momenta<T,DESCRIPTOR>* getTemperatureBoundaryMomenta()
  {
    return new EquilibriumBM<T,DESCRIPTOR>;
  }
  template<int direction, int orientation>
  static Dynamics<T,DESCRIPTOR>* getTemperatureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
  {
    return new RtlbmDirectedBoundaryDynamics<T,DESCRIPTOR,direction,orientation>(omega, momenta);
  }
  template<int direction, int orientation>
  static PostProcessorGenerator3D<T,DESCRIPTOR>* getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
  {
    return nullptr;
  }

  // ====== edged
  template<int plane, int normal1, int normal2>
  static Momenta<T,DESCRIPTOR>* getTemperatureBoundaryEdgeMomenta()
  {
    return new EquilibriumBM<T,DESCRIPTOR>;
  }
  template<int plane, int normal1, int normal2>
  static Dynamics<T,DESCRIPTOR>* getTemperatureBoundaryEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
  {
    return new AdvectionDiffusionEdgesDynamics<T,DESCRIPTOR,BGKdynamics<T,DESCRIPTOR>,plane,normal1,normal2>(omega, momenta); // placeholder
  }
  template<int plane, int normal1, int normal2>
  static PostProcessorGenerator3D<T,DESCRIPTOR>* getTemperatureBoundaryEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
  {
    return nullptr;
  }

  // ====== corners
  template<int normalX, int normalY, int normalZ>
  static Momenta<T,DESCRIPTOR>* getTemperatureBoundaryCornerMomenta()
  {
    return new EquilibriumBM<T,DESCRIPTOR>;
  }
  template<int normalX, int normalY, int normalZ>
  static Dynamics<T,DESCRIPTOR>* getTemperatureBoundaryCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
  {
    return new AdvectionDiffusionCornerDynamics3D<T,DESCRIPTOR,BGKdynamics<T,DESCRIPTOR>,normalX,normalY,normalZ>(omega, momenta); // placeholder
  }
  template<int normalX, int normalY, int normalZ>
  static PostProcessorGenerator3D<T,DESCRIPTOR>* getTemperatureBoundaryCornerProcessor(int x, int y, int z)
  {
    return nullptr;
  }
};

//================== creator function ================================
// blockLattice creator
template<typename T, typename DESCRIPTOR>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR>* createRtlbmDirectedBoundaryCondition3D( BlockLatticeStructure3D<T,DESCRIPTOR>& block )
{
  return new AdvectionDiffusionBoundaryConditionInstantiator3D<T, DESCRIPTOR, RtlbmDirectedBoundaryManager3D<T,DESCRIPTOR> > (block);
}

// superLattice creator
template<typename T, typename DESCRIPTOR>
void createRtlbmDirectedBoundaryCondition3D( sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sBC )
{
  sBC.setOverlap(1);
  for (int iC = 0; iC < sBC.getSuperLattice().getLoadBalancer().size(); iC++) {
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR>* blockBC =
      createRtlbmDirectedBoundaryCondition3D<T,DESCRIPTOR>( sBC.getSuperLattice().getExtendedBlockLattice(iC) );
    sBC.getADblockBCs().push_back(blockBC);
  }
}

}  // namespace olb

#endif
