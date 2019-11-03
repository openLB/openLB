/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
#ifndef BOUNDARY_CONDITION_3D_HH
#define BOUNDARY_CONDITION_3D_HH

#include "boundaryCondition3D.h"
#include "boundaryInstantiator3D.h"
#include "momentaOnBoundaries3D.h"

namespace olb {

template<typename T, typename DESCRIPTOR, class MixinDynamics>
class RegularizedBoundaryManager3D {
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
  getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv=NULL);

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

template<typename T, typename DESCRIPTOR, class MixinDynamics>
class InterpolationBoundaryManager3D {
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
  getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv=NULL);

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


////////// RegularizedBoundaryManager3D /////////////////////////////////////////

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new RegularizedVelocityBM<T,DESCRIPTOR, direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new RegularizedPressureBM<T,DESCRIPTOR, direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv)
{
  return nullptr;
}


template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityEdgeMomenta()
{
  return new FixedVelocityBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, plane,normal1,normal2>(x0,x1, y0,y1, z0,z1);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityEdgeMomenta()
{
  return new InnerEdgeVelBM3D<T,DESCRIPTOR, plane,normal1,normal2>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerProcessor(int x, int y, int z)
{
  return new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, xNormal,yNormal,zNormal> (x,y,z);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM3D<T,DESCRIPTOR, xNormal,yNormal,zNormal>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,DESCRIPTOR>*
RegularizedBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityCornerProcessor(int x, int y, int z)
{
  return nullptr;
}


////////// InterpolationBoundaryManager3D /////////////////////////////////////////

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new BasicDirichletBM<T,DESCRIPTOR,VelocityBM, direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new PlaneFdBoundaryProcessorGenerator3D
         <T,DESCRIPTOR, direction,orientation>(x0,x1, y0,y1, z0,z1);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new BasicDirichletBM<T,DESCRIPTOR,PressureBM, direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new PlaneFdBoundaryProcessorGenerator3D
         <T,DESCRIPTOR, direction,orientation>(x0, x1, y0, y1, z0, z1);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv)
{
  return new StraightConvectionBoundaryProcessorGenerator3D
         <T,DESCRIPTOR,direction,orientation>(x0, x1, y0, y1, z0, z1, uAv);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityEdgeMomenta()
{
  return new FixedVelocityBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, plane,normal1,normal2>(x0,x1, y0,y1, z0,z1);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityEdgeMomenta()
{
  return new InnerEdgeVelBM3D<T,DESCRIPTOR, plane,normal1,normal2>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerProcessor(int x, int y, int z)
{
  return new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, xNormal,yNormal,zNormal> (x,y,z);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM3D<T,DESCRIPTOR, xNormal,yNormal,zNormal>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,DESCRIPTOR>*
InterpolationBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityCornerProcessor(int x, int y, int z)
{
  return nullptr;
}

////////// Convenience wrappers for boundary functions ////////////////////////

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, int z0, int z1,
  T omega)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addVelocityBoundary(indicator, x0, x1, y0, y1, z0, z1, omega);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addVelocityBoundary(
  BlockIndicatorF3D<T>& indicator, T omega, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addVelocityBoundary(indicator,
                      margin, blockGeometryStructure.getNx()-1 -margin,
                      margin, blockGeometryStructure.getNy()-1 -margin,
                      margin, blockGeometryStructure.getNz()-1 -margin,
                      omega);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  T omega, bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addVelocityBoundary(indicator, omega, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addSlipBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, int z0, int z1)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addSlipBoundary(indicator, x0, x1, y0, y1, z0, z1);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addSlipBoundary(
  BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addSlipBoundary(indicator,
                  margin, blockGeometryStructure.getNx()-1 -margin,
                  margin, blockGeometryStructure.getNy()-1 -margin,
                  margin, blockGeometryStructure.getNz()-1 -margin);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addSlipBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material, bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addSlipBoundary(indicator, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addPartialSlipBoundary(
  T tuner, BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, int z0, int z1)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addPartialSlipBoundary(tuner, indicator, x0, x1, y0, y1, z0, z1);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addPartialSlipBoundary(
  T tuner, BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addPartialSlipBoundary(tuner, indicator,
                         margin, blockGeometryStructure.getNx()-1 -margin,
                         margin, blockGeometryStructure.getNy()-1 -margin,
                         margin, blockGeometryStructure.getNz()-1 -margin);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addPartialSlipBoundary(
  T tuner, BlockGeometryStructure3D<T>& blockGeometryStructure, int material, bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addPartialSlipBoundary(tuner, indicator, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addWallFunctionBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, int z0, int z1,
  UnitConverter<T, DESCRIPTOR> const& converter,
  wallFunctionParam<T> const& wallFunctionParam,
  IndicatorF3D<T>*            geoIndicator)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addWallFunctionBoundary(indicator,
                          x0, y1, y0, y1, z0, z1,
                          converter, wallFunctionParam, geoIndicator);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addWallFunctionBoundary(
  BlockIndicatorF3D<T>& indicator,
  UnitConverter<T, DESCRIPTOR> const& converter,
  wallFunctionParam<T> const& wallFunctionParam,
  IndicatorF3D<T>*            geoIndicator,
  bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addWallFunctionBoundary(indicator,
                          margin, blockGeometryStructure.getNx()-1 -margin,
                          margin, blockGeometryStructure.getNy()-1 -margin,
                          margin, blockGeometryStructure.getNz()-1 -margin,
                          converter, wallFunctionParam, geoIndicator);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addWallFunctionBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  UnitConverter<T, DESCRIPTOR> const& converter,
  wallFunctionParam<T> const& wallFunctionParam,
  IndicatorF3D<T>*            geoIndicator,
  bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addWallFunctionBoundary(indicator, converter, wallFunctionParam, geoIndicator, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addPressureBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, int z0, int z1,
  T omega)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addPressureBoundary(indicator, x0, x1, y0, y1, z0, z1, omega);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addPressureBoundary(
  BlockIndicatorF3D<T>& indicator, T omega, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addPressureBoundary(indicator,
                      margin, blockGeometryStructure.getNx()-1 -margin,
                      margin, blockGeometryStructure.getNy()-1 -margin,
                      margin, blockGeometryStructure.getNz()-1 -margin,
                      omega);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addPressureBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  T omega, bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addPressureBoundary(indicator, omega, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addConvectionBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, int z0, int z1,
  T omega, T* uAv)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addConvectionBoundary(indicator, x0, x1, y0, y1, z0, z1, omega, uAv);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addConvectionBoundary(
  BlockIndicatorF3D<T>& indicator, T omega, T* uAv, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addConvectionBoundary(indicator,
                        margin, blockGeometryStructure.getNx()-1 -margin,
                        margin, blockGeometryStructure.getNy()-1 -margin,
                        margin, blockGeometryStructure.getNz()-1 -margin,
                        omega, uAv);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addConvectionBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material, T omega, T* uAv, bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addConvectionBoundary(indicator, omega, uAv, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, int z0, int z1, T addend, int latticeNumber)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addFreeEnergyWallBoundary(indicator, x0, x1, y0, y1, z0, z1, addend, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  BlockIndicatorF3D<T>& indicator, T addend, int latticeNumber, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addFreeEnergyWallBoundary(indicator,
                      margin, blockGeometryStructure.getNx()-1 -margin,
                      margin, blockGeometryStructure.getNy()-1 -margin,
                      margin, blockGeometryStructure.getNz()-1 -margin,
                      addend, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material, T addend,
  int latticeNumber, bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addFreeEnergyWallBoundary(indicator, addend, latticeNumber, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyInletBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, int z0, int z1, T omega, std::string type, int latticeNumber)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addFreeEnergyInletBoundary(indicator, x0, x1, y0, y1, z0, z1, omega, type, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyInletBoundary(
  BlockIndicatorF3D<T>& indicator, T omega, std::string type,
  int latticeNumber, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addFreeEnergyInletBoundary(indicator,
                      margin, blockGeometryStructure.getNx()-1 -margin,
                      margin, blockGeometryStructure.getNy()-1 -margin,
                      margin, blockGeometryStructure.getNz()-1 -margin,
                      omega, type, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyInletBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material, 
  T omega, std::string type, int latticeNumber, bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addFreeEnergyInletBoundary(indicator, omega, type, latticeNumber, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyOutletBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, int z0, int z1, T omega, std::string type, int latticeNumber)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addFreeEnergyOutletBoundary(indicator, x0, x1, y0, y1, z0, z1, omega, type, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyOutletBoundary(
  BlockIndicatorF3D<T>& indicator, T omega, std::string type,
  int latticeNumber, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addFreeEnergyOutletBoundary(indicator,
                      margin, blockGeometryStructure.getNx()-1 -margin,
                      margin, blockGeometryStructure.getNy()-1 -margin,
                      margin, blockGeometryStructure.getNz()-1 -margin,
                      omega, type, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyOutletBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
  T omega, std::string type, int latticeNumber, bool includeOuterCells)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometryStructure, material);
  addFreeEnergyOutletBoundary(indicator, omega, type, latticeNumber, includeOuterCells);
}

////////// Factory functions //////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* createLocalBoundaryCondition3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block)
{
  return new BoundaryConditionInstantiator3D <
         T, DESCRIPTOR,
         RegularizedBoundaryManager3D<T,DESCRIPTOR, MixinDynamics> > (block);
}

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* createInterpBoundaryCondition3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block)
{
  return new BoundaryConditionInstantiator3D <
         T, DESCRIPTOR,
         InterpolationBoundaryManager3D<T,DESCRIPTOR, MixinDynamics> > (block);
}

}  // namespace olb

#endif
