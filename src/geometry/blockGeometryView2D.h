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
 * Representation of the 2D block geometry view -- header file.
 */

#ifndef BLOCK_GEOMETRY_VIEW_2D_H
#define BLOCK_GEOMETRY_VIEW_2D_H

#include <vector>
#include "geometry/blockGeometry2D.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "geometry/blockGeometryStructure2D.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

/// Representation of a block geometry view
/** This class is derived from block geometry structure. It
 * holds a poniter to another block geometry structure. It operates
 * as a structure on a smaller intersection of teh greater structure.
 * It presents a volume of voxels where different types are
 * given my material numbers which is imporant e.g. to work
 * with different boundaries (like for inflow/output regions).
 *
 * This class is not intended to be derived from.
 */

template<typename T> class BlockGeometryStatistics2D;
template<typename T> class BlockGeometryStructure2D;

template<typename T>
class BlockGeometryView2D final : public BlockGeometryStructure2D<T>, public BlockStructure2D {

private:
  // Points to the structure where this view class is viewing at
  BlockGeometryStructure2D<T>* _originalBlockGeometry;
  // Offset of the data field with respect to the data of the original geometry
  int _x0, _y0;
  // Dimension of the view cuboid
  int _nx, _ny;

public:
  /// Constructor
  BlockGeometryView2D(BlockGeometryStructure2D<T>& originalBlockGeometry, int x0, int x1, int y0, int y1);
  /// Copy constructor
  BlockGeometryView2D(BlockGeometryView2D const& rhs);
  /// Copy assignment
  BlockGeometryView2D& operator=(BlockGeometryView2D const& rhs);
  /// Destructor
  ~BlockGeometryView2D() override;

  /// Returns the underlying block structure
  BlockStructure2D& getBlockStructure() override;

  /// Write access to the associated block statistic
  BlockGeometryStatistics2D<T>& getStatistics(bool verbose=true) override;
  /// Read only access to the associated block statistic
  BlockGeometryStatistics2D<T> const& getStatistics(bool verbose=true) const override;

  /// Read only access to the origin position given in SI units (meter)
  Vector<T,2> getOrigin() const override;
  /// Read only access to the voxel size given in SI units (meter)
  const T getDeltaR() const override;
  /// Returns the extend (number of voxels) in X-direction
  int getNx() const override;
  /// Returns the extend (number of voxels) in Y-direction
  int getNy() const override;

  /// Write access to a material number
  int& get(int iX, int iY) override;
  /// Read only access to a material number
  int const& get(int iX, int iY) const override;
  /// returns the (iX,iY,iZ) entry in the 2D scalar field
  int getMaterial(int iX, int iY) const override; // TODO old

  /// Transforms lattice to physical coordinates (wrapped from cuboid geometry)
  void getPhysR(T physR[2], const int& iX, const int& iY) const override;

  /// Adds a pointer to the list of dependent statistic classes
  void addToStatisticsList(bool* statisticStatus) override;
  /// Removes a pointer from the list of dependent statistic classes if existing
  void removeFromStatisticsList(bool* statisticStatus) override;
};

} // namespace olb

#endif
