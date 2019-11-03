/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
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
 * The description of a vector of 2D cuboid -- header file.
 */


#ifndef CUBOID_GEOMETRY_2D_H
#define CUBOID_GEOMETRY_2D_H

#include <vector>
#include "geometry/cuboid2D.h"
#include "io/ostreamManager.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T> class IndicatorF2D;

/// A cuboid structure represents the grid of a considered domain
/** A cuboid structure is given by a number of cuboids. To represent
 * a connected domain it is required that the distance between two
 * neighbooring cuboids is less than the smallest delta of them.
 *
 * WARNING:
 * At the moment there are only cuboids with a constant delta possible
 * and the distance between two neighbooring cuboids must be delta
 * since an interpolation operator in time and space is missing in
 * cuboidNeigbourhood and superLattice.
 *
 * This class is not intended to be derived from.
 */
template<typename T>
class CuboidGeometry2D {

private:
  // TODO quick and dirty solution for 3d split to be removed after 2d clearence
  bool _oldApproach;

  /// Vector of the cuboids
  std::vector<Cuboid2D<T> > _cuboids;
  /// Cuboid which contains all other cuboids
  Cuboid2D<T> _motherCuboid;
  /// Periodicity flag
  std::vector<bool> _periodicityOn;
  /// class specific ostream
  mutable OstreamManager clout;

public:

  /// Constructs empty Geometry
  CuboidGeometry2D();
  /// Constructs a cuboid geometry with a cubic shape of size nX times nY with origin originR=(originX, originY) that consits of nC cuboids
  CuboidGeometry2D(T originX, T originY, T deltaR, int nX, int nY, int nC=1);
  /// Constructs a cuboid structure with a uniform spacing of voxelsize which consits of  nC cuboids, the cuboids not needed are removed and too big ones are shrinked
  CuboidGeometry2D(IndicatorF2D<T>& indicatorF, T voxelSize, int nC=1);

  /// Re init
  void reInit(T globPosX, T globPosY,
              T delta, int nX, int nY, int nC=1);
  /// Read and write access to the cuboids
  Cuboid2D<T>& get(int i);
  /// Read access to the cuboids
  Cuboid2D<T> const& get(int i) const;
  /// Set flag to enable/disable periodicity
  void setPeriodicity(bool periodicityX, bool periodicityY);

  /// Returns true and the cuboid number of the nearest lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  bool getC(std::vector<T> physR, int& iC) const;
  /// Returns true and the cuboid number of the nearest lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  bool getC(const Vector<T,2>& physR, int& iC) const;
  /// Returns true and the nearest lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  bool getLatticeR(std::vector<T> physR, std::vector<int>& latticeR) const;
  bool getLatticeR(int latticeR[], const T physR[]) const;
  bool getLatticeR(const Vector<T,2>& physR, Vector<int,3>& latticeR) const;
  /// Returns true and the floor lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  bool getFloorLatticeR(std::vector<T> physR, std::vector<int>& latticeR) const;
  /// Returns true and the floor lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  bool getFloorLatticeR(const Vector<T,2>& physR, Vector<int,3>& latticeR) const;
  /// Returns the physical position to the given lattice position respecting periodicity for the overlap nodes which are not in the mother cuboid for the case the flag periodicityOn[iDim]=true if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  std::vector<T> getPhysR(int iCglob, int iX, int iY) const;
  /// Returns the physical position to the given lattice position respecting periodicity for the overlap nodes which are not in the mother cuboid for the case the flag periodicityOn[iDim]=true
  std::vector<T> getPhysR(std::vector<int> latticeR) const;
  void getPhysR(T output[2], const int latticeR[3]) const;
  void getPhysR(T output[2], const int iCglob, const int iX, const int iY) const;

  /// Returns the number of cuboids in t < 2he structure
  int getNc() const;
  /// Returns the maximum/minimum of the ratio nX/NY in the structure
  T getMinRatio() const;
  T getMaxRatio() const;
  /// Returns the minimum coordinate in the structure
  std::vector<T> getMinPhysR() const;
  /// Returns the maximum coordinate in the structure
  std::vector<T> getMaxPhysR() const;
  /// Returns the maximum/minimum volume in the structure
  T getMinPhysVolume() const;
  T getMaxPhysVolume() const;
  /// Returns the maximum/minimum number of nodes in the structure
  size_t getMinLatticeVolume() const;
  size_t getMaxLatticeVolume() const;
  /// Returns the maximum/minimum delata in the structure
  T getMinDeltaR() const;
  T getMaxDeltaR() const;
  /// Returns the smallest cuboid that includes all cuboids of
  /// the structure
  Cuboid2D<T> getMotherCuboid() const;

  /// for a given point (globX/globY), returns the related cuboidID
  /// and _p if the point is not in any of the cuboid _childrenQ
  int get_iC(T globX, T globY, int offset = 0) const;
  /// This function checks if the points (globX/globY) and
  /// (globX + orientationX/delta /globY + orientationY/delta) is in
  /// a cuboid. It gives the related cuboidID and _p if the points are
  /// not in any of the cuboids.
  /// abs(orientationX) = abs(orientationY) = 1 must be satisfied
  int get_iC(T globX, T globY, int orientationX, int orientationY) const;

  /// Adds a cuboid
  void add(Cuboid2D<T> cuboid);
  /// Removes the cuboid iC
  void remove(int iC);
  /// Removes all cuboids where geometryData = 0
  //void remove(olb::ScalarField2D<int>* geometryData);
  /// Splits cuboid iC, removes it and add p cuboids
  void split(int iC, int p);
  /// Shrink all cuboids so that no empty planes are left
  void shrink(IndicatorF2D<T>& indicatorF);

  /// stores the neighbouring cuboids in array neighbours;
  void getNeighbourhood(int cuboid, std::vector<int> neighbours, int offset = 0);

  void refineArea(T x0, T x1, T y0, T y1, int coarse_level);

  /// Prints cuboid geometry details
  void print() const;
  /// Prints cuboid geometry details plus details of all cuboids
  void printExtended();
};

}  // namespace olb

#endif
