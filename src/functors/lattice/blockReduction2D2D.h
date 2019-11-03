/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef BLOCK_REDUCTION_2D2D_H
#define BLOCK_REDUCTION_2D2D_H

#include "io/ostreamManager.h"
#include "core/blockData2D.h"
#include "core/vector.h"
#include "blockBaseF2D.h"
#include "superBaseF2D.h"
#include "superBaseF3D.h"
#include "utilities/functorPtr.h"
#include "utilities/blockDataSyncMode.h"
#include "utilities/hyperplaneLattice3D.h"

#include <tuple>

namespace olb {


class BlockStructure2D;

/// BlockReduction2D2D interpolates the data of a SuperF2D functor in a given resolution.
/**
 * This is primarily used for exporting GIF images via BlockGifWriter.
 **/
template<typename T>
class BlockReduction2D2D final : public BlockDataF2D<T,T> {
private:
  /// Functor to be reduced
  FunctorPtr<SuperF2D<T>> _f;
  /// Plane points scheduled for storage in _blockData
  /// i.e. Plane points whose physical location intersects the mother cuboid
  ///      and is nearest to a rank-local cuboid
  std::vector<std::tuple<int,int,int>> _rankLocalSubplane;

  Vector<T,2> _origin; /// origin of the cuboid
  T _h;                /// spacing according to given resolution
  int _nx;             /// horizontal size of the reduced cuboid
  int _ny;             /// vertical size of the reduced cuboid
  /// Data fields which hold the reduced data
  std::unique_ptr<BlockData2D<T,T>> _blockDataMemory;
  /// Synchronization mode, see BlockDataSyncMode enum for further information.
  /// This value only matters when PARALLEL_MODE_MPI is defined.
  const BlockDataSyncMode _syncMode;

  /// Updates _h, _nx, _ny such that the longest side is resolution voxels long
  void updateToWantedResolution(int resolution);

public:
  BlockReduction2D2D(FunctorPtr<SuperF2D<T>>&& f,
                     int resolution=600,
                     BlockDataSyncMode mode=BlockDataSyncMode::ReduceAndBcast);

  /// Transform lattice coordinates to their physical location
  Vector<T,2> getPhysR(const int& iX, const int& iY) const;
  /// Returns embedding of the discretized plane in 3D space
  /**
   * i.e. span vectors are X, Y unit vectors.
   *      Origin, spacing and resolution is exposed.
   *
   * Required for dimension agnostic implementation of GnuplotHeatMap<T>
   **/
  HyperplaneLattice3D<T> getPlaneDiscretizationIn3D() const;

  /// Initialize rank-local list of points to be stored in _blockData
  void initialize();
  /// Updates and writes the data to _blockData using _rankLocalSubplane
  void update();
  /// Overload of virtual function from class BlockF2D
  BlockStructure2D& getBlockStructure() override;

};


} // end namespace olb

#endif
