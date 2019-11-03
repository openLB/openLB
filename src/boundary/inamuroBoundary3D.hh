/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006,2007 Orestis Malaspinas and Jonas Latt
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

#ifndef INAMURO_BOUNDARY_3D_HH
#define INAMURO_BOUNDARY_3D_HH

#include "inamuroBoundary3D.h"
#include "inamuroNewtonRaphsonDynamics.h"
#include "inamuroNewtonRaphsonDynamics.hh"
#include "boundaryInstantiator3D.h"
#include "momentaOnBoundaries3D.h"

namespace olb {

template<typename T, typename DESCRIPTOR, class MixinDynamics>
class InamuroBoundaryManager3D {
public:
  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getVelocityBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getPressureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getPressureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int direction, int orientation> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv = NULL);

  template<int plane, int normal1, int normal2> static Momenta<T,DESCRIPTOR>*
  getExternalVelocityEdgeMomenta();
  template<int plane, int normal1, int normal2> static Dynamics<T,DESCRIPTOR>*
  getExternalVelocityEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int plane, int normal1, int normal2> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getExternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int plane, int normal1, int normal2> static Momenta<T,DESCRIPTOR>*
  getInternalVelocityEdgeMomenta();
  template<int plane, int normal1, int normal2> static Dynamics<T,DESCRIPTOR>*
  getInternalVelocityEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int plane, int normal1, int normal2> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getInternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int xNormal, int yNormal, int zNormal> static Momenta<T,DESCRIPTOR>*
  getExternalVelocityCornerMomenta();
  template<int xNormal, int yNormal, int zNormal> static Dynamics<T,DESCRIPTOR>*
  getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal, int zNormal> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getExternalVelocityCornerProcessor(int x, int y, int z);

  template<int xNormal, int yNormal, int zNormal> static Momenta<T,DESCRIPTOR>*
  getInternalVelocityCornerMomenta();
  template<int xNormal, int yNormal, int zNormal> static Dynamics<T,DESCRIPTOR>*
  getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal, int zNormal> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getInternalVelocityCornerProcessor(int x, int y, int z);
};

////////// InamuroBoundaryManager3D /////////////////////////////////////////

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new BasicDirichletBM<T,DESCRIPTOR, VelocityBM, direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new InamuroNewtonRaphsonDynamics<T,DESCRIPTOR, MixinDynamics, direction, orientation>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new BasicDirichletBM<T,DESCRIPTOR, PressureBM, direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new InamuroNewtonRaphsonDynamics<T,DESCRIPTOR, MixinDynamics, direction, orientation>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityEdgeMomenta()
{
  return new FixedVelocityBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, plane,normal1,normal2>(x0,x1, y0,y1, z0,z1);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityEdgeMomenta()
{
  return new InnerEdgeVelBM3D<T,DESCRIPTOR, plane,normal1,normal2>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerProcessor(int x, int y, int z)
{
  return new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, xNormal,yNormal,zNormal> (x,y,z);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM3D<T,DESCRIPTOR, xNormal,yNormal,zNormal>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InamuroBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityCornerProcessor(int x, int y, int z)
{
  return nullptr;
}


////////// Factory functions //////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* createInamuroBoundaryCondition3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block)
{
  return new BoundaryConditionInstantiator3D <
         T, DESCRIPTOR,
         InamuroBoundaryManager3D<T,DESCRIPTOR, MixinDynamics> > (block);
}

}  // namespace olb

#endif
