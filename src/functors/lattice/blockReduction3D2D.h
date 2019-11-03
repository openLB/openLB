/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef BLOCK_REDUCTION_3D2D_H
#define BLOCK_REDUCTION_3D2D_H

#include "core/blockData2D.h"
#include "core/vector.h"
#include "blockBaseF2D.h"
#include "superBaseF2D.h"
#include "superBaseF3D.h"
#include "utilities/hyperplane3D.h"
#include "utilities/hyperplaneLattice3D.h"
#include "utilities/functorPtr.h"
#include "utilities/blockDataSyncMode.h"
#include "utilities/blockDataReductionMode.h"

#include <tuple>

namespace olb {


/// BlockReduction3D2D reduces the data of a SuperF3D functor to the
/// intersection between a given hyperplane and the super geometry.
/**
 * This intersection is interpolated at a set of discrete points according to
 * the given resolution and exposed as a BlockDataF2D functor.
 *
 * The hyperplane is parametrized by a origin and two span vector u and v.
 * Definition of hyperplanes using e.g. origin and normal vectors is supported
 * via the Hyperplane3D interface.
 **/
template <typename T>
class BlockReduction3D2D final : public HyperplaneLattice3D<T>, public BlockDataF2D<T,T> {
private:
  /// Data fields to hold the reduced data
  std::unique_ptr<BlockData2D<T,T>> _blockDataMemory;
  /// Functor to be reduced
  FunctorPtr<SuperF3D<T>> _f;
  /// Plane points scheduled for storage in _blockData
  /// i.e. Plane points whose physical location intersects the mother cuboid
  ///      and is nearest to a rank-local cuboid
  std::vector<std::tuple<int,int,int>> _rankLocalSubplane;
  /// Synchronization mode, see BlockDataSyncMode enum for further information.
  /// This value only matters when PARALLEL_MODE_MPI is defined.
  const BlockDataSyncMode _syncMode;
  /// Reduction mode, see BlockDataReductionMode enum for further information.
  const BlockDataReductionMode _reductionMode;

  void updateBlockAnalytical(BlockData2D<T,T>& block);
  void updateBlockDiscrete(BlockData2D<T,T>& block);

public:
  /// Construction using functor and hyperplane lattice
  /**
   * \param f Functor to be reduced as a (non-)owning pointer or reference to SuperF3D<T>.
   * \param lattice Hyperplane lattice parametrization
   * \param syncMode
   *        Defines MPI synchronization strategy of the interpolated block data.
   * \param reductionMode
   *        Defines whether data is interpolated or read from discrete lattice locations.
   *        Note: BlockDataReductionMode::Analytical imposes restrictions on hyperplane
   *        definition and discretization.
   **/
  BlockReduction3D2D(FunctorPtr<SuperF3D<T>>&& f,
                     const HyperplaneLattice3D<T>& lattice,
                     BlockDataSyncMode syncMode=BlockDataSyncMode::ReduceAndBcast,
                     BlockDataReductionMode reductionMode=BlockDataReductionMode::Analytical);
  /// Construction using functor and hyperplane
  /**
   * \param f Functor to be reduced as a (non-)owning pointer or reference to SuperF3D<T>.
   * \param hyperplane Hyperplane parametrization
   * \param syncMode
   *        Defines MPI synchronization strategy of the interpolated block data.
   * \param reductionMode
   *        Defines whether data is interpolated or read from discrete lattice locations.
   **/
  BlockReduction3D2D(FunctorPtr<SuperF3D<T>>&& f,
                     const Hyperplane3D<T>& hyperplane,
                     BlockDataSyncMode syncMode=BlockDataSyncMode::ReduceAndBcast,
                     BlockDataReductionMode reductionMode=BlockDataReductionMode::Analytical);
  /// Construction using functor, hyperplane and resolution
  /**
   * \param f Functor to be reduced as a (non-)owning pointer or reference to SuperF3D<T>.
   * \param hyperplane Hyperplane parametrization
   * \param resolution Defines the number of voxel of the longest side.
   *                   If it equals zero, _h is set to the cuboid geometry's minDeltaR.
   * \param syncMode   Defines MPI synchronization strategy of the interpolated block data.
   **/
  BlockReduction3D2D(FunctorPtr<SuperF3D<T>>&& f,
                     const Hyperplane3D<T>& hyperplane,
                     int resolution=600,
                     BlockDataSyncMode syncMode=BlockDataSyncMode::ReduceAndBcast);

  /// Construction using functor, origin and span vectors as well as resolution
  /**
   * \param f Functor to be reduced as a (non-)owning pointer or reference to SuperF3D<T>.
   * \param origin Origin vector
   * \param u      Span vector
   * \param v      Span vector
   * \param resolution Defines the number of voxel of the longest side.
   *                   If it equals zero, _h is set to the cuboid geometry's minDeltaR.
   * \param syncMode Defines MPI synchronization strategy of the interpolated block data.
   **/
  BlockReduction3D2D(FunctorPtr<SuperF3D<T>>&& f,
                     const Vector<T,3>& origin, const Vector<T,3>& u, const Vector<T,3>& v,
                     int resolution=600,
                     BlockDataSyncMode syncMode=BlockDataSyncMode::ReduceAndBcast);
  /// Construction using functor, origin, normal and resolution
  /**
   * \param f Functor to be reduced as a (non-)owning pointer or reference to SuperF3D<T>.
   * \param origin Origin vector
   * \param normal Normal vector
   * \param resolution Defines the number of voxel of the longest side.
   *                   If it equals zero, _h is set to the cuboid geometry's minDeltaR.
   * \param syncMode Defines MPI synchronization strategy of the interpolated block data.
   **/
  BlockReduction3D2D(FunctorPtr<SuperF3D<T>>&& f,
                     const Vector<T,3>& origin, const Vector<T,3>& normal,
                     int resolution=600,
                     BlockDataSyncMode syncMode=BlockDataSyncMode::ReduceAndBcast);
  /// Construction using functor, normal vector and resolution
  /**
   * \param f Functor to be reduced as a (non-)owning pointer or reference to SuperF3D<T>.
   * \param normal     Normal vector
   * \param resolution Defines the number of voxel of the longest side.
   *                   If it equals zero, _h is set to the cuboid geometry's minDeltaR.
   * \param syncMode   Defines MPI synchronization strategy of the interpolated block data.
   **/
  BlockReduction3D2D(FunctorPtr<SuperF3D<T>>&& f,
                     const Vector<T,3>& normal,
                     int resolution=600,
                     BlockDataSyncMode syncMode=BlockDataSyncMode::ReduceAndBcast);

  /// Initialize rank-local list of plane points to be stored in _blockData
  void initialize();
  /// Updates and writes the data to _blockData using _rankLocalSubplane
  void update();
  /// Overload of virtual function from class BlockF2D
  BlockStructure2D& getBlockStructure() override;
  /// \return reference to the rank local list of discrete plane points, cuboid ids
  const std::vector<std::tuple<int,int,int>>& getRankLocalSubplane() const;

};


} // end namespace olb

#endif
