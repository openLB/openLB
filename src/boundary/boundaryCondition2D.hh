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
 * A helper for initialising 2D boundaries -- generic implementation.
 */
#ifndef BOUNDARY_CONDITION_2D_HH
#define BOUNDARY_CONDITION_2D_HH

#include "boundaryCondition2D.h"
#include "boundaryInstantiator2D.h"


namespace olb {

template<typename T, typename DESCRIPTOR, class MixinDynamics>
class RegularizedBoundaryManager2D {
public:
  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getVelocityBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getPressureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getPressureBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, T* uAv=NULL);

  template<int xNormal, int yNormal> static Momenta<T,DESCRIPTOR>*
  getExternalVelocityCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,DESCRIPTOR>*
  getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getExternalVelocityCornerProcessor(int x, int y);

  template<int xNormal, int yNormal> static Momenta<T,DESCRIPTOR>*
  getInternalVelocityCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,DESCRIPTOR>*
  getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getInternalVelocityCornerProcessor(int x, int y);
};

template<typename T, typename DESCRIPTOR, class MixinDynamics>
class InterpolationBoundaryManager2D {
public:
  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getVelocityBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getPressureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getPressureBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, T* uAv=NULL);

  template<int xNormal, int yNormal> static Momenta<T,DESCRIPTOR>*
  getExternalVelocityCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,DESCRIPTOR>*
  getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getExternalVelocityCornerProcessor(int x, int y);

  template<int xNormal, int yNormal> static Momenta<T,DESCRIPTOR>*
  getInternalVelocityCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,DESCRIPTOR>*
  getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getInternalVelocityCornerProcessor(int x, int y);
};


////////// RegularizedBoundaryManager2D /////////////////////////////////////////

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new RegularizedVelocityBM<T,DESCRIPTOR, direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new RegularizedPressureBM<T,DESCRIPTOR, direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, T* uAv)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,DESCRIPTOR>*
RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,DESCRIPTOR>* RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,DESCRIPTOR>*
RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerProcessor(int x, int y)
{
  return new OuterVelocityCornerProcessorGenerator2D<T,DESCRIPTOR, xNormal,yNormal> (x,y);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,DESCRIPTOR>*
RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM2D<T,DESCRIPTOR, xNormal,yNormal>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,DESCRIPTOR>* RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,DESCRIPTOR>*
RegularizedBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityCornerProcessor
(int x, int y)
{
  return nullptr;
}


////////// InterpolationBoundaryManager2D /////////////////////////////////////////

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new BasicDirichletBM<T,DESCRIPTOR,VelocityBM,direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return new StraightFdBoundaryProcessorGenerator2D
         < T,DESCRIPTOR,direction,orientation  >  (x0, x1, y0, y1);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new BasicDirichletBM<T,DESCRIPTOR,PressureBM,direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return new StraightFdBoundaryProcessorGenerator2D
         < T,DESCRIPTOR,direction,orientation >  (x0, x1, y0, y1);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, T* uAv)
{
  return new StraightConvectionBoundaryProcessorGenerator2D
         < T,DESCRIPTOR,direction,orientation >  (x0, x1, y0, y1, uAv);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,DESCRIPTOR>*
InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,DESCRIPTOR>* InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,DESCRIPTOR>*
InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityCornerProcessor(int x, int y)
{
  return new OuterVelocityCornerProcessorGenerator2D<T,DESCRIPTOR, xNormal,yNormal> (x,y);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,DESCRIPTOR>*
InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM2D<T,DESCRIPTOR, xNormal,yNormal>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,DESCRIPTOR>* InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,DESCRIPTOR>*
InterpolationBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityCornerProcessor (int x, int y)
{
  return nullptr;
}

////////// Convenience wrappers for boundary functions ////////////////////////

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, T omega)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addVelocityBoundary(indicator, x0, x1, y0, y1, omega);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addVelocityBoundary(
  BlockIndicatorF2D<T>& indicator, T omega, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addVelocityBoundary(indicator,
                      margin, blockGeometryStructure.getNx()-1 -margin,
                      margin, blockGeometryStructure.getNy()-1 -margin,
                      omega);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega, bool includeOuterCells)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addVelocityBoundary(indicator, omega, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addSlipBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addSlipBoundary(indicator, x0, x1, y0, y1);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addSlipBoundary(
  BlockIndicatorF2D<T>& indicator, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addSlipBoundary(indicator,
                  margin, blockGeometryStructure.getNx()-1 -margin,
                  margin, blockGeometryStructure.getNy()-1 -margin);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addSlipBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, bool includeOuterCells)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addSlipBoundary(indicator, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPartialSlipBoundary(
  T tuner, BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addPartialSlipBoundary(tuner, indicator, x0, x1, y0, y1);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPartialSlipBoundary(
  T tuner, BlockIndicatorF2D<T>& indicator, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addPartialSlipBoundary(tuner, indicator,
                         margin, blockGeometryStructure.getNx()-1 -margin,
                         margin, blockGeometryStructure.getNy()-1 -margin);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPartialSlipBoundary(
  T tuner, BlockGeometryStructure2D<T>& blockGeometryStructure, int material, bool includeOuterCells)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addPartialSlipBoundary(tuner, indicator, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPressureBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, T omega)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addPressureBoundary(indicator, x0, x1, y0, y1, omega);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPressureBoundary(
  BlockIndicatorF2D<T>& indicator, T omega, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addPressureBoundary(indicator,
                      margin, blockGeometryStructure.getNx()-1 -margin,
                      margin, blockGeometryStructure.getNy()-1 -margin,
                      omega);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPressureBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega, bool includeOuterCells)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addPressureBoundary(indicator, omega, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addConvectionBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, T omega, T* uAv)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addConvectionBoundary(indicator, x0, x1, y0, y1, omega, uAv);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addConvectionBoundary(
  BlockIndicatorF2D<T>& indicator, T omega, T* uAv, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addConvectionBoundary(indicator,
                        margin, blockGeometryStructure.getNx()-1 -margin,
                        margin, blockGeometryStructure.getNy()-1 -margin,
                        omega, uAv);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addConvectionBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega, T* uAv, bool includeOuterCells)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addConvectionBoundary(indicator, omega, uAv, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, T addend, int latticeNumber)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addFreeEnergyWallBoundary(indicator, x0, x1, y0, y1, addend, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  BlockIndicatorF2D<T>& indicator, T addend, int latticeNumber, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addFreeEnergyWallBoundary(indicator,
                      margin, blockGeometryStructure.getNx()-1 -margin,
                      margin, blockGeometryStructure.getNy()-1 -margin,
                      addend, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T addend,
  int latticeNumber, bool includeOuterCells)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addFreeEnergyWallBoundary(indicator, addend, latticeNumber, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyInletBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, T omega, std::string type, int latticeNumber)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addFreeEnergyInletBoundary(indicator, x0, x1, y0, y1, omega, type, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyInletBoundary(
  BlockIndicatorF2D<T>& indicator, T omega, std::string type,
  int latticeNumber, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addFreeEnergyInletBoundary(indicator,
                      margin, blockGeometryStructure.getNx()-1 -margin,
                      margin, blockGeometryStructure.getNy()-1 -margin,
                      omega, type, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyInletBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  T omega, std::string type, int latticeNumber, bool includeOuterCells)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addFreeEnergyInletBoundary(indicator, omega, type, latticeNumber, includeOuterCells);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyOutletBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  int x0, int x1, int y0, int y1, T omega, std::string type, int latticeNumber)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addFreeEnergyOutletBoundary(indicator, x0, x1, y0, y1, omega, type, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyOutletBoundary(
  BlockIndicatorF2D<T>& indicator, T omega, std::string type,
  int latticeNumber, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  addFreeEnergyOutletBoundary(indicator,
                      margin, blockGeometryStructure.getNx()-1 -margin,
                      margin, blockGeometryStructure.getNy()-1 -margin,
                      omega, type, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void OnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyOutletBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
  T omega, std::string type, int latticeNumber, bool includeOuterCells)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addFreeEnergyOutletBoundary(indicator, omega, type, latticeNumber, includeOuterCells);
}

////////// Factory functions //////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
OnLatticeBoundaryCondition2D<T,DESCRIPTOR>* createLocalBoundaryCondition2D(BlockLatticeStructure2D<T,DESCRIPTOR>& block)
{
  return new BoundaryConditionInstantiator2D <
         T, DESCRIPTOR,
         RegularizedBoundaryManager2D<T,DESCRIPTOR, MixinDynamics> > (block);
}

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
OnLatticeBoundaryCondition2D<T,DESCRIPTOR>* createInterpBoundaryCondition2D(BlockLatticeStructure2D<T,DESCRIPTOR>& block)
{
  return new BoundaryConditionInstantiator2D <
         T, DESCRIPTOR,
         InterpolationBoundaryManager2D<T,DESCRIPTOR, MixinDynamics> > (block);
}

}  // namespace olb

#endif
