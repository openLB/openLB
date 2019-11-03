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

#ifndef BLOCK_GEOMETRY_VIEW_3D_H
#define BLOCK_GEOMETRY_VIEW_3D_H

#include <vector>

#include "core/blockStructure3D.h"
#include "geometry/blockGeometry3D.h"
#include "geometry/blockGeometryStatistics3D.h"
#include "geometry/blockGeometryStructure3D.h"

namespace olb {

/// Representation of a block geometry view
/** This class is derived from block structure and block geometry structure.
 * It holds a pointer to another block geometry structure. It operates
 * as a structure on a smaller intersection of the greater structure.
 * It presents a volume of voxels where different types are
 * given my material numbers which is imporant e.g. to work
 * with different boundaries (like for inflow/output regions).
 *
 * This class is not intended to be derived from.
 */

template<typename T> class BlockGeometryStatistics3D;
template<typename T> class BlockGeometryStructure3D;

/*
 * This class needs to be final as addToStatisticsList and removeFromStatisticsList in the constructor don't dispatch to derived classes.
 */
template<typename T>
class BlockGeometryView3D final : public BlockGeometryStructure3D<T>, public BlockStructure3D {

private:
  // Points to the structure where this view class is viewing at
  BlockGeometryStructure3D<T>* _originalBlockGeometry;
  // Offset of the data field with respect to the data of the original geometry
  int _x0, _y0, _z0;

public:
  /// Constructor
  BlockGeometryView3D(BlockGeometryStructure3D<T>& originalBlockGeometry, int x0, int x1, int y0, int y1, int z0, int z1);
  /// Copy constructor
  BlockGeometryView3D(BlockGeometryView3D const& rhs);
  /// Copy assignment
  BlockGeometryView3D& operator=(BlockGeometryView3D const& rhs);
  /// Destructor
  ~BlockGeometryView3D() override;

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
  int& get(int iX, int iY, int iZ) override;
  /// Read only access to a material number
  int const& get(int iX, int iY, int iZ) const override;
  /// returns the (iX,iY,iZ) entry in the 3D scalar field
  int getMaterial(int iX, int iY, int iZ) const override; // TODO old

  /// Transforms lattice to physical coordinates (wrapped from cuboid geometry)
  void getPhysR(T physR[3], const int& iX, const int& iY, const int& iZ) const override;

  /// Adds a pointer to the list of dependent statistic classes
  void addToStatisticsList(bool* statisticStatus) override;
  /// Removes a pointer from the list of dependent statistic classes if existing
  void removeFromStatisticsList(bool* statisticStatus) override;
};

} // namespace olb

#endif
