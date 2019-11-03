/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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
 * Representation of the 3D block geometry view -- header file.
 */

#ifndef BLOCK_GEOMETRY_3D_H
#define BLOCK_GEOMETRY_3D_H

#include <vector>
#include <list>
#include "core/blockData3D.h"
#include "geometry/blockGeometryStatistics3D.h"
#include "geometry/blockGeometryStructure3D.h"
#include "geometry/cuboid3D.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

/// Representation of a block geometry
/** This class is derived from block geometry structure. It
 * holds the actual data with the materials. It stores pointers
 * to all dependent block geometry views.
 * It presents a volume of voxels where different types are
 * given by material numbers which is important e.g. to work
 * with different boundaries (like for inflow/output regions).
 *
 * This class is not intended to be derived from.
 */
template<typename T>
class BlockGeometry3D final : public BlockData3D<T,int>, public BlockGeometryStructure3D<T> {

private:
  /// Cuboid which charaterizes the block geometry
  Cuboid3D<T> _cuboid;
  /// List to all depending statistic status objects
  std::list<bool*> _statisticsUpdateNeeded;

public:
  /// Constructor
  BlockGeometry3D(T x0, T y0, T z0, T h, int nX, int nY, int nZ, int iCglob=-1);
  /// Constructor
  BlockGeometry3D(Cuboid3D<T>& cuboid, int iCglob=-1);
  /// Copy constructor
  BlockGeometry3D(BlockGeometry3D const& rhs);
  /// Copy assignment
  BlockGeometry3D& operator=(BlockGeometry3D const& rhs);

  /// Returns the underlying block structure
  BlockStructure3D& getBlockStructure() override;

  /// Write access to the associated block statistic
  BlockGeometryStatistics3D<T>& getStatistics(bool verbose=true) override;
  /// Read only access to the associated block statistic
  BlockGeometryStatistics3D<T> const& getStatistics(bool verbose=true) const override;

  /// Read only access to the origin position given in SI units (meter)
  Vector<T,3> getOrigin() const override;
  /// Read only access to the voxel size given in SI units (meter)
  const T getDeltaR() const override;
  /// Returns the extend (number of voxels) in X-direction
  int getNx() const override;
  /// Returns the extend (number of voxels) in Y-direction
  int getNy() const override;
  /// Returns the extend (number of voxels) in Z-direction
  int getNz() const override;
  /// Write access to a material number
  int& get(int iX, int iY, int iZ) override;  // override BlockGeometryStructure3D::get() by linking to BlockData3D<T,int>::get()
  /// Read only access to a material number
  int const& get(int iX, int iY, int iZ) const override;  // override BlockGeometryStructure3D::get() by linking to BlockData3D<T,int>::get()
  /// returns the (iX,iY,iZ) entry in the 3D scalar field
  int getMaterial(int iX, int iY, int iZ) const override; // TODO old

  /// Transforms lattice to physical coordinates (wrapped from cuboid geometry)
  void getPhysR(T physR[3], const int& iX, const int& iY, const int& iZ) const override;

  // returns the raw data field
  // olb::ScalarField3D<int>* getRawData();
  // void resize(int X, int Y, int Z, int nX, int nY, int nZ);
  // refine Mesh
  // void refineMesh(int level = 2);

  /// Adds a pointer to the list of dependent statistic classes
  void addToStatisticsList(bool* statisticStatus) override;
  /// Removes a pointer from the list of dependent statistic classes if existing
  void removeFromStatisticsList(bool* statisticStatus) override;

  /// Prints a chosen part of the block geometry
  void printLayer(int x0, int x1, int y0, int y1, int z0, int z1, bool linenumber = false);
  /// Prints a chosen part of the block geometry
  void printLayer(int direction, int layer, bool linenumber = false);
  /// Prints a chosen node and its neighbourhood
  void printNode(int x0, int y0, int z0);

private:

  /// Re-initialization of the block geometry
  /*void reInit(T originX, T originY, T originZ, T deltaR, int nX, int nY, int nZ, int iCglob = -1,
              olb::ScalarField3D<int>* geometryData = nullptr);*/
  /// Resets all depending statistic flags
  void resetStatistics();
};

} // namespace olb

#endif
