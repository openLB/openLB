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
 * The description of a single 3D cuboid -- header file.
 */

#ifndef CUBOID_3D_H
#define CUBOID_3D_H

#include <vector>
#include <math.h>

#include "core/serializer.h"
#include "io/ostreamManager.h"
#include "functors/analytical/indicator/indicatorF3D.h"


/// All OpenLB code is contained in this namespace.
namespace olb {


template<typename T> class IndicatorF3D;


/// A regular single 3D cuboid is the basic component of a 3D cuboid
/// structure which defines the grid.
/** A cuboid is given with its left lower corner, the number of nodes
 * in the direction x, y and z and the distance between two nodes.
 * Among after useful methods a cuboid can divide itself in a given
 * number of disjoint subcuboids. The number of nodes at the boundary
 * is minimized.
 *
 * This class is not intended to be derived from.
 */
template<typename T>
class Cuboid3D final : public Serializable {
private:
  /// Global position of the left lower corner of the cuboid
  T   _globPosX, _globPosY, _globPosZ;
  /// Distance to the next node
  T   _delta;
  /// Number of nodes in the direction x and y
  int _nX, _nY, _nZ;
  /// Number of full cells
  size_t _weight;
  /// refinement level, _delta = _delta0^_refinementLevel
  int _refinementLevel;

  /// Specific ostream for the classname in each line
  mutable OstreamManager clout;

public:
  /// Construction of an empty cuboid at position 0, 0, 0 with delta 0 and nX = nY = nZ = 0
  Cuboid3D();
  /// Construction of a cuboid
  Cuboid3D(T globPosX, T globPosY, T globPosZ, T delta, int nX, int nY, int nZ, int refinementLevel=0);
  /// Construction of a cuboid vector version
  Cuboid3D(std::vector<T> origin, T delta, std::vector<int> extend, int refinementLevel=0);
  /// Construction of a cuboid using indicator
  Cuboid3D(IndicatorF3D<T>& indicatorF, T voxelSize, int refinementLevel=0);
  /// Copy constructor
  Cuboid3D(Cuboid3D<T> const& rhs, int overlap=0);
  /// Copy assignment
  Cuboid3D& operator=(Cuboid3D const& rhs);
  /// Initializes the cuboid
  void init(T globPosX, T globPosY, T globPosZ, T delta, int nX, int nY, int nZ, int refinementLevel=0); //TODO: remove or private

  /// Read only access to left lower corner coordinates
  Vector<T,3> const getOrigin() const;
  /// Read only access to the distance of cuboid nodes
  T getDeltaR() const;
  /// Read access to cuboid width
  int getNx() const;
  /// Read access to cuboid height
  int getNy() const;
  /// Read access to cuboid depth
  int getNz() const;
  /// Read only access to the number of voxels in every dimension
  Vector<int,3> const getExtend() const;
  /// Returns the volume of the cuboid
  T getPhysVolume() const;
  /// Returns the actual value of weight (-1 for getLatticeVolume())
  int getWeightValue() const;
  /// Returns the number of full cells
  size_t getWeight() const;
  /// Sets the number of full cells
  void setWeight(size_t fullCells);
  /// Returns the refinementLevel
  int getRefinementLevel() const;
  /// Sets the refinementLevel
  void setRefinementLevel(int refLevel);
  /// Returns the number of Nodes in the volume
  size_t getLatticeVolume() const;
  /// Returns the perimeter of the cuboid
  T getPhysPerimeter() const;
  /// Returns the number of Nodes at the perimeter
  int getLatticePerimeter() const;

  /// equal operator
  bool operator==(const Cuboid3D<T>& rhs) const;

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer interface
  std::size_t getSerializableSize() const override;
  /// \return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

  /// Prints cuboid details
  void print() const;

  void getPhysR(T physR[3], const int latticeR[3]) const;
  void getPhysR(T physR[3], const int& iX, const int& iY, const int& iZ) const;

  void getLatticeR(int latticeR[3], const T physR[3]) const;
  void getLatticeR(int latticeR[3], const Vector<T,3>& physR) const;

  void getFloorLatticeR(const std::vector<T>& physR, std::vector<int>& latticeR) const;
  void getFloorLatticeR(int latticeR[3], const T physR[3]) const;

  /// Checks whether a point (globX/gloxY/globZ) is contained in the cuboid
  /// extended with an layer of size overlap*delta
  bool checkPoint(T globX, T globY, T globZ, int overlap = 0) const; //TODO globX-> x + additional interface: with (std::vector<T> physR, int overlap=0)
  /// Same for physical overlap
  bool physCheckPoint(T globX, T globY, T globZ, double overlap = 0) const;
  /// Checks whether a point (globX/gloxY/globZ) is contained and is a node
  /// in the cuboid extended with an layer of size overlap*delta and
  /// returns the local active node
  bool checkPoint(T globX, T globY, T globZ, int &locX, int &locY, int &locZ, int overlap = 0) const; //TODO globX-> x + additional interface: with (std::vector<T> physR, std::vector<int>& latticeR, int overlap=0)
  /// Checks whether there is an intersection with the cuboid extended
  /// with an layer of size overlap*delta
  bool checkInters(T globX0, T globX1, T globY0, T globY1, T globZ0, T globZ1, int overlap = 0) const;
  /// Checks whether a given point intersects the cuboid extended
  /// by a layer of size overlap*delta
  bool checkInters(T globX, T globY, T globZ, int overlap = 0) const;
  /// Checks whether there is an intersection and returns the local
  /// active node range which can be empty by means of locX0=1, locX1=0,
  /// locY0=1, locY1=0, locZ0=1; locZ1=0 of the cuboid extended with an
  /// layer of size overlap*delta
  bool checkInters(T globX0, T globX1, T globY0, T globY1, T globZ0, T globZ1,
                   int &locX0, int &locX1, int &locY0, int &locY1, int &locZ0, int &locZ1,
                   int overlap = 0) const;
  /// Divides the cuboid in p*q*r cuboids of equal volume and add them to the given vector
  void divide(int p, int q, int r, std::vector<Cuboid3D<T> > &childrenC) const;
  /// Divides the cuboid in p cuboids of equal volume and add them to the given vector
  void divide(int p, std::vector<Cuboid3D<T> > &childrenC) const;

  /// resize the cuboid to the passed size
  void resize(int X, int Y, int Z, int nX, int nY, int nZ);
};

}  // namespace olb

#endif
