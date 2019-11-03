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
 * A helper for initialising 2D boundaries -- header file.
 */

#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_2D_H
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_2D_H

#include "momentaOnBoundaries2D.h"
#include "advectionDiffusionMomentaOnBoundaries.h"
#include "dynamics/dynamics.h"
#include "dynamics/advectionDiffusionDynamics.h"
#include "functors/lattice/indicator/blockIndicatorBaseF2D.h"

#include <vector>
#include <list>

namespace olb {

template<typename T, typename DESCRIPTOR>
class sOnLatticeBoundaryCondition2D;

template<typename T, typename DESCRIPTOR>
class OnLatticeAdvectionDiffusionBoundaryCondition2D {
public:
  virtual ~OnLatticeAdvectionDiffusionBoundaryCondition2D() { }

  virtual void addTemperatureBoundary0N(int x0, int x1, int y0, int y1,T omega) =0;
  virtual void addTemperatureBoundary0P(int x0, int x1, int y0, int y1,T omega) =0;
  virtual void addTemperatureBoundary1N(int x0, int x1, int y0, int y1,T omega) =0;
  virtual void addTemperatureBoundary1P(int x0, int x1, int y0, int y1,T omega) =0;

  virtual void addTemperatureCornerNN(int x, int y, T omega) =0;
  virtual void addTemperatureCornerNP(int x, int y, T omega) =0;
  virtual void addTemperatureCornerPN(int x, int y, T omega) =0;
  virtual void addTemperatureCornerPP(int x, int y, T omega) =0;

  BlockLatticeStructure2D<T,DESCRIPTOR>& getBlock();
  BlockLatticeStructure2D<T,DESCRIPTOR> const& getBlock() const;

  /// Add temperature boundary for indicated cells in range
  virtual void addTemperatureBoundary(BlockIndicatorF2D<T>& indicator,
                                      int x0, int x1, int y0, int y1,
                                      T omega) =0;

  /**
   * \name Convenience wrappers for temperature boundary functions
   * \{
   **/

  void addTemperatureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                              int x0, int x1, int y0, int y1,
                              T omega);
  void addTemperatureBoundary(BlockIndicatorF2D<T>& indicator,
                              T omega, bool includeOuterCells=false);
  void addTemperatureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                              T omega, bool includeOuterCells=false);

  ///\}

  /// adds a temperature boundary for one material or a range (x0-x1, y0-y1, z0-z1)
  virtual void addRegularizedTemperatureBoundary(BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addRegularizedTemperatureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addRegularizedTemperatureBoundary(BlockIndicatorF2D<T>& indicator, T omega) =0;
  virtual void addRegularizedTemperatureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega) =0;

  virtual void addRegularizedHeatFluxBoundary(BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1, T omega, T *heatFlux) =0;
  virtual void addRegularizedHeatFluxBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega, T *heatFlux) =0;
  virtual void addRegularizedHeatFluxBoundary(BlockIndicatorF2D<T>& indicator, T omega, T *heatFlux) =0;
  virtual void addRegularizedHeatFluxBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega, T *heatFlux) =0;
};

template<typename T, typename DESCRIPTOR, class MixinDynamics>
class AdvectionDiffusionBoundaryManager2D {
public:
  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getTemperatureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getTemperatureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int xNormal, int yNormal> static Momenta<T,DESCRIPTOR>*
  getTemperatureCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,DESCRIPTOR>*
  getTemperatureCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getTemperatureCornerProcessor(int x, int y);


  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getRegularizedTemperatureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getRegularizedTemperatureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getRegularizedTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int xNormal, int yNormal> static Momenta<T,DESCRIPTOR>*
  getRegularizedTemperatureCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,DESCRIPTOR>*
  getRegularizedTemperatureCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getRegularizedTemperatureCornerProcessor(int x, int y);

  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getRegularizedHeatFluxBoundaryMomenta(T *heatFlux);
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getRegularizedHeatFluxBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getRegularizedHeatFluxBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int xNormal, int yNormal> static Momenta<T,DESCRIPTOR>*
  getRegularizedHeatFluxCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,DESCRIPTOR>*
  getRegularizedHeatFluxCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getRegularizedHeatFluxCornerProcessor(int x, int y);
};

//////  Factory function for Regularized Thermal BC

/// blockLattice creator
template<typename T, typename DESCRIPTOR, typename MixinDynamics=AdvectionDiffusionRLBdynamics<T,DESCRIPTOR> >
OnLatticeAdvectionDiffusionBoundaryCondition2D<T,DESCRIPTOR>*
createAdvectionDiffusionBoundaryCondition2D(BlockLatticeStructure2D<T,DESCRIPTOR>& block);

/// superLattice creator, calls createAdvectionDiffusionBoundaryCondidtion3D from above.
template<typename T, typename DESCRIPTOR, typename MixinDynamics=AdvectionDiffusionRLBdynamics<T,DESCRIPTOR> >
void createAdvectionDiffusionBoundaryCondition2D(sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>& sBC);


} //namespace olb


#endif
