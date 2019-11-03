/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2011, 2014 Mathias J. Krause, Simon Zimny
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
 * Representation of a statistic for a 2D geometry -- header file.
 */

#ifndef BLOCK_GEOMETRY_STATISTICS_2D_H
#define BLOCK_GEOMETRY_STATISTICS_2D_H

#include <vector>
#include <map>
#include <string>

#include "io/ostreamManager.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

/// Representation of a statistic for a 2D geometry
/** A block geomety statistic computes different integral
 * values, like total number of different materials,
 * materials of any kind, min./max. physical position, of an
 * underlying block geoemtry structure.
 *
 * This class is not intended to be derived from.
 */

template<typename T>
class BlockGeometryStructure2D;

template<typename T>
class BlockGeometryStatistics2D {

private:

  /// Points to the underlying data from which the statistics is taken
  BlockGeometryStructure2D<T>* _blockGeometry;
  /// Specifies if an update is needed
  bool _statisticsUpdateNeeded;
  /// Number of voxels in each direction
  int _nX, _nY;
  /// Spacing
  T _h;

  /// Number of different material numbers
  int _nMaterials;
  /// Mapping a material number to the number of this kind found in the super geometry
  std::map<int, int> _material2n;
  /// Mapping a material number to the min. lattice position in each space direction
  std::map<int, std::vector<int> > _material2min;
  /// Mapping a material number to the max. lattice position in each space direction
  std::map<int, std::vector<int> > _material2max;

  /// class specific cout
  mutable OstreamManager clout;

public:

  /// Constructor
  BlockGeometryStatistics2D(BlockGeometryStructure2D<T>* blockGeometry);

  /// Read and write access to a flag, which indicates if an uptate is needed (=true)
  bool& getStatisticsStatus();
  /// Read only access to a flag, which indicates if an uptate is needed (=true)
  bool const & getStatisticsStatus() const;
  /// Returns the map with the numbers of voxels for each  material
  std::map<int, int> getMaterial2n();

  /// Updates the statistics if it is really needed
  void update(bool verbose=true);

  /// Returns the number of different materials
  int getNmaterials();
  /// Returns the number of voxels for a given material number
  int getNvoxel(int material);
  /// Returns the number of voxels with material!=0
  int getNvoxel();
  /// Returns the min. lattice position in each direction
  std::vector<int> getMinLatticeR(int material);
  /// Returns the max. lattice position in each direction
  std::vector<int> getMaxLatticeR(int material);
  /// Returns the min. phys position in each direction
  std::vector<T> getMinPhysR(int material);
  /// Returns the max. phys position in each direction
  std::vector<T> getMaxPhysR(int material);
  /// Returns the lattice extend as length in each direction
  std::vector<T> getLatticeExtend(int material);
  /// Returns the phys extend as length in each direction
  std::vector<T> getPhysExtend(int material);
  /// Returns the phys radius as length in each direction
  std::vector<T> getPhysRadius(int material);
  /// Returns the center position
  std::vector<T> getCenterPhysR(int material);
  /// Returns the boundary type which is characterized by a discrete normal (c.f. Zimny)
  std::vector<int> getType(int iX, int iY);

  /// Returns normal that points into the fluid for paraxial surfaces
  std::vector<int> computeNormal(int iX, int iY);
  /// Returns normal that points into the fluid for paraxial surfaces
  std::vector<T> computeNormal (int material);
  /// Returns discrete normal with norm maxNorm that points into the fluid for paraxial surfaces
  /// maxNorm=1.1 implies only normals parallel to the axises
  std::vector<int> computeDiscreteNormal (int material, T maxNorm = 1.1);

  // Returns true if at position (iX,iY,iZ) and in a neighbourhood of size (offsetX,offsetY,offsetZ) only voxels of the given material are found
  bool check(int material, int iX, int iY, unsigned offsetX, unsigned offsetY);
  // Returns true and a position (iX,iY,iZ) if there is a neighbourhood of size (offsetX,offsetY,offsetZ) around (iX,iY,iZ) with only voxels of the given material
  bool find(int material, unsigned offsetX, unsigned offsetY, int& iX, int& iY);

  /// Prints some statistic information, i.e. the number of voxels and min. max. physical position for each different material
  void print();

private:

  /// Helper function to simplify the implementation
  void takeStatistics(int iX, int iY);
};

} // namespace olb

#endif
