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
 * A helper for initialising 3D boundaries -- header file.
 */

#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_H
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_H


#include "dynamics/advectionDiffusionDynamics.h"
#include "geometry/blockGeometryStructure3D.h"

namespace olb {

template<typename T>
class BlockIndicatorF3D;

template<typename T, typename DESCRIPTOR>
class sOnLatticeBoundaryCondition3D;

/* Interface class whose childs are boundary instantiators.
*  Instantiators have a third template parameter BoundaryManager that finaly implements the interfaces from here.
*/
template<typename T, typename DESCRIPTOR>
class OnLatticeAdvectionDiffusionBoundaryCondition3D {
public:
  virtual ~OnLatticeAdvectionDiffusionBoundaryCondition3D() { }

  /**
   * \anchor AdvDiffBCimplI3D
   * \name Indicator-accepting boundary condition interfaces
   *
   * \param indicator Block indicator defining boundary cells
   * \param x0,x1,y0,y1,z0,z1 Range of cells to be traversed
   * \{
   **/

  virtual void addTemperatureBoundary(BlockIndicatorF3D<T>& indicator,
                                      int x0, int x1, int y0, int y1, int z0, int z1,
                                      T omega) =0;
  virtual void addConvectionBoundary(BlockIndicatorF3D<T>& indicator,
                                     int x0, int x1, int y0, int y1, int z0, int z1) =0;
  virtual void addExtFieldBoundary(BlockIndicatorF3D<T>& indicator,
                                   int offset,
                                   int x0, int x1, int y0, int y1, int z0, int z1) =0;
  /// Add BC that initializes zero distributions and computes the density that entered the boundary
  virtual void addZeroDistributionBoundary(BlockIndicatorF3D<T>& indicator,
      int x0, int x1, int y0, int y1, int z0, int z1) =0;

  ///\}

  /**
   * \name Convenience wrappers for boundary functions
   * In practice it is often preferable to define a boundary on a single material
   * number instead of instantiating an appropriate indicator by hand.
   *
   * These convenience functions are simple wrappers around the actual
   * \ref AdvDiffBCimplI3D "boundary condition interfaces".
   * \{
   **/

  void addTemperatureBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                              int x0, int x1, int y0, int y1, int z0, int z1,
                              T omega);
  void addTemperatureBoundary(BlockIndicatorF3D<T>& indicator,
                              T omega, bool includeOuterCells=false);
  void addTemperatureBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                              T omega, bool includeOuterCells=false);

  void addConvectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                             int x0, int x1, int y0, int y1, int z0, int z1);
  void addConvectionBoundary(BlockIndicatorF3D<T>& indicator, bool includeOuterCells=false);
  void addConvectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, bool includeOuterCells=false);

  void addExtFieldBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                           int offset,
                           int x0, int x1, int y0, int y1, int z0, int z1);
  void addExtFieldBoundary(BlockIndicatorF3D<T>& indicator,
                           int offset, bool includeOuterCells=false);
  void addExtFieldBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                           int offset, bool includeOuterCells=false);

  void addZeroDistributionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                                   int x0, int x1, int y0, int y1, int z0, int z1);
  void addZeroDistributionBoundary(BlockIndicatorF3D<T>& indicator, bool includeOuterCells=false);
  void addZeroDistributionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, bool includeOuterCells=false);

  ///\}
};

/// blockLattice creator
template<typename T, typename DESCRIPTOR, typename MixinDynamics=AdvectionDiffusionRLBdynamics<T,DESCRIPTOR> >
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR>*
createAdvectionDiffusionBoundaryCondition3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block);

/// superLattice creator, calls createAdvectionDiffusionBoundaryCondidtion3D from above.
template<typename T, typename DESCRIPTOR, typename MixinDynamics=AdvectionDiffusionRLBdynamics<T,DESCRIPTOR> >
void createAdvectionDiffusionBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& sBC);


} //namespace olb


#endif
