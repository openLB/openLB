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

#ifndef ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_3D_H
#define ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_3D_H

#include "advectionDiffusionBoundaryCondition3D.h"

namespace olb {

/* Implements methods by calling methods of template class BoundaryManager.
*  This class detects through a discreteNormal whether voxel is edge, corner or plate.
*  Additionally it detects orientation of plate, direction of corner, ...
*  Empose boundary model through template class BoundaryManager.
*/
template<typename T, typename DESCRIPTOR, class BoundaryManager>
class AdvectionDiffusionBoundaryConditionInstantiator3D : public OnLatticeAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR> {
public:
  AdvectionDiffusionBoundaryConditionInstantiator3D( BlockLatticeStructure3D<T,DESCRIPTOR>& block );
  ~AdvectionDiffusionBoundaryConditionInstantiator3D() override;

  void addTemperatureBoundary(BlockIndicatorF3D<T>& indicator,
                              int x0, int x1, int y0, int y1, int z0, int z1,
                              T omega) override;

  void addConvectionBoundary(BlockIndicatorF3D<T>& indicator,
                             int x0, int x1, int y0, int y1, int z0, int z1) override;

  void addExtFieldBoundary(BlockIndicatorF3D<T>& indicator,
                           int offset,
                           int x0, int x1, int y0, int y1, int z0, int z1) override;

  void addZeroDistributionBoundary(BlockIndicatorF3D<T>& indicator,
                                   int x0, int x1, int y0, int y1, int z0, int z1) override;

private:
  template<int direction, int orientation>
  void addTemperatureBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  template<int plane, int normal1, int normal2>
  void addTemperatureBoundaryEdge(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  template<int normalX, int normalY, int normalZ>
  void addTemperatureBoundaryCorner(int x, int y, int z, T omega);

private:
  BlockLatticeStructure3D<T,DESCRIPTOR>& _block;
  std::vector<Momenta<T,DESCRIPTOR>*>  momentaVector;
  std::vector<Dynamics<T,DESCRIPTOR>*> dynamicsVector;
};



} // namespace openlb


#endif
