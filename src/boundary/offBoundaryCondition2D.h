/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2016 Jonas Kratzke, Mathias J. Krause
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

#ifndef OFF_BOUNDARY_CONDITION_2D_H
#define OFF_BOUNDARY_CONDITION_2D_H

#include <vector>
#include "core/blockLatticeStructure2D.h"
#include "core/blockLatticeStructure2D.h"
#include "offBoundaryCondition2D.h"
#include "dynamics/dynamics.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "io/stlReader.h"
namespace olb {

/**
* This class provides a general off lattice boundary condition
*/

template<typename T, typename DESCRIPTOR>
class OffLatticeBoundaryCondition2D {

protected:
  T _epsFraction;

public:
  virtual ~OffLatticeBoundaryCondition2D() { }

  /// Using Bouzidi BC OnePoint corresponds to Bounce Back and TwoPoint to linear interpolation
  virtual void addOnePointZeroVelocityBoundary(int iX, int iY, int iPop, T dist) =0;
  virtual void addTwoPointZeroVelocityBoundary(int iX, int iY, int iPop, T dist) =0;
  virtual void addOnePointVelocityBoundary(int iX, int iY, int iPop, T dist) =0;
  virtual void addTwoPointVelocityBoundary(int iX, int iY, int iPop, T dist) =0;

  virtual void addOffDynamics(int iX, int iY, T location[DESCRIPTOR::d]) =0;
  virtual void addOffDynamics(int iX, int iY, T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q]) =0;
  virtual void addOffDynamics(BlockIndicatorF2D<T>& indicator) =0;

  virtual void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, int iPop, T dist) =0;
  virtual void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator) =0;
  virtual void addZeroVelocityBoundary(BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator) =0;
  virtual void addZeroVelocityBoundary(BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator) =0;

  virtual void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int iPop, T dist) =0;
  virtual void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator) =0;
  virtual void addVelocityBoundary(BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator) =0;
  virtual void addVelocityBoundary(BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator) =0;

  virtual void addPressureBoundary(BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator) =0;
  virtual void addPressureBoundary(BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator) =0;

  virtual void defineU(int iX, int iY, int iPop, const T u[DESCRIPTOR::d]) =0;
  virtual void defineU(BlockIndicatorF2D<T>& indicator, BlockIndicatorF2D<T>& bulkIndicator, AnalyticalF2D<T,T>& u) =0;

  virtual void defineRho(int iX, int iY, int iPop, const T rho) =0;
  virtual void defineRho(BlockIndicatorF2D<T>& indicator, BlockIndicatorF2D<T>& bulkIndicator, AnalyticalF2D<T,T>& rho) =0;

  /**
   * \name Convenience wrappers for boundary functions
   * \{
   **/

  void addOffDynamics(BlockGeometryStructure2D<T>& blockGeometryStructure, int material);
  void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, IndicatorF2D<T>& geometryIndicator, std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& geometryIndicator, std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, IndicatorF2D<T>& geometryIndicator, std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& geometryIndicator, std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& geometryIndicator, std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void defineU(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, AnalyticalF2D<T,T>& u, std::vector<int> bulkMaterials = std::vector<int>(1,1) );
  void defineRho(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, AnalyticalF2D<T,T>& rho, std::vector<int> bulkMaterials = std::vector<int>(1,1) );

  ///\}

  virtual BlockLatticeStructure2D<T,DESCRIPTOR>& getBlock() =0;
  virtual BlockLatticeStructure2D<T,DESCRIPTOR> const& getBlock() const =0;

  virtual void outputOn() =0;
  virtual void outputOff() =0;
};

////////// Factory functions //////////////////////////////////////////////////

/**
* Create specific off lattice boundary conditions
*/

template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR> >
OffLatticeBoundaryCondition2D<T,DESCRIPTOR>*
createBouzidiBoundaryCondition2D(BlockLatticeStructure2D<T,DESCRIPTOR>& block);


template<typename T, typename DESCRIPTOR>
OffLatticeBoundaryCondition2D<T,DESCRIPTOR>*
createBounceBackBoundaryCondition2D(BlockLatticeStructure2D<T,DESCRIPTOR>& block)
{
  return createBounceBackBoundaryCondition2D<T,DESCRIPTOR>(block);
}

}

#endif
