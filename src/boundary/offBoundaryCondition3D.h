/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Kratzke, Mathias J. Krause
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

#ifndef OFF_BOUNDARY_CONDITION_3D_H
#define OFF_BOUNDARY_CONDITION_3D_H

#include <vector>
#include "core/blockLatticeStructure3D.h"
#include "core/blockLatticeStructure3D.h"
#include "offBoundaryCondition3D.h"
#include "dynamics/dynamics.h"
#include "geometry/blockGeometryStatistics3D.h"
#include "io/stlReader.h"
namespace olb {

/**
* This class provides a general off lattice boundary condition
*/

template<typename T, typename DESCRIPTOR>
class OffLatticeBoundaryCondition3D {

protected:
  T _epsFraction;

public:
  virtual ~OffLatticeBoundaryCondition3D() { }

  /// Using Bouzidi BC OnePoint corresponds to Bounce Back and TwoPoint to linear interpolation
  virtual void addOnePointZeroVelocityBoundary(int x, int y, int z, int iPop, T dist) =0;
  virtual void addTwoPointZeroVelocityBoundary(int x, int y, int z, int iPop, T dist) =0;
  virtual void addOnePointVelocityBoundary(int x, int y, int z, int iPop, T dist) =0;
  virtual void addTwoPointVelocityBoundary(int x, int y, int z, int iPop, T dist) =0;

  virtual void addOffDynamics(int x, int y, int z, T location[DESCRIPTOR::d]) =0;
  virtual void addOffDynamics(int x, int y, int z, T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q]) =0;

  virtual void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure,
                                       int x, int y, int z, int iPop, T dist) =0;
  virtual void defineU(int iX, int iY, int iZ, int iPop, const T u[DESCRIPTOR::d]) =0;

  virtual void addOffDynamics(BlockIndicatorF3D<T>& indicator) =0;
  virtual void addZeroVelocityBoundary(BlockIndicatorF3D<T>& boundaryIndicator,
                                       BlockIndicatorF3D<T>& bulkIndicator,
                                       IndicatorF3D<T>&      geometryIndicator) =0;
  virtual void addVelocityBoundary(BlockIndicatorF3D<T>& boundaryIndicator,
                                   BlockIndicatorF3D<T>& bulkIndicator,
                                   IndicatorF3D<T>&      geometryIndicator) =0;
  virtual void defineU(BlockIndicatorF3D<T>& indicator,
                       BlockIndicatorF3D<T>& bulkIndicator,
                       AnalyticalF3D<T,T>&   u) =0;

  /**
   * \name Convenience wrappers for boundary functions
   * \{
   **/

  void addOffDynamics(BlockGeometryStructure3D<T>& blockGeometryStructure, int material);
  void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                               IndicatorF3D<T>& geometryIndicator,
                               std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                           IndicatorF3D<T>& indicator,
                           std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void defineU(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
               AnalyticalF3D<T,T>& u,
               std::vector<int> bulkMaterials = std::vector<int>(1,1) );

  ///\}

  virtual BlockLatticeStructure3D<T,DESCRIPTOR>& getBlock() =0;
  virtual BlockLatticeStructure3D<T,DESCRIPTOR> const& getBlock() const =0;

  virtual void outputOn() =0;
  virtual void outputOff() =0;
};

////////// Factory functions //////////////////////////////////////////////////

/**
* Create specific off lattice boundary conditions
*/

template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR> >
OffLatticeBoundaryCondition3D<T,DESCRIPTOR>*
createBouzidiBoundaryCondition3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block);

}

#endif
