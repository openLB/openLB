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
 * The description of a single 2D cuboid -- header file.
 */

#ifndef CUBOID_2D_H
#define CUBOID_2D_H

#include <vector>
#include <math.h>
#include "io/ostreamManager.h"
#include "core/vector.h"
#include "core/olbDebug.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

/// A regular single 2D cuboid is the basic component of a 2D cuboid
/// structure which defines the grid.
/** A cuboid is given with its left lower corner, the number of nodes
 * in the direction x and y and the distance between two nodes.
 * Among other useful methods, a cuboid can divide itself in a given
 * number of disjoint subcuboids. The number of nodes at the boundary
 * is minimized.
 *
 * This class is not intended to be derived from.
 */

template<typename T>
class Cuboid2D {
private:
  /// Global position of the left lower corner of the cuboid
  T   _globPosX, _globPosY;
  /// Distance to the next node
  T   _delta;
  /// Number of nodes in the direction x and y and the refinement Level
  int _nX, _nY;
  /// Number of full cells
  size_t _weight;
  /// refinement level, _delta = _delta0^_refinementLevel
  int _refinementLevel;

  /// Specific ostream for the classname in each line
  mutable OstreamManager clout;

public:
  /// Construction of a cuboid
  Cuboid2D(T globPosX, T globPosY, T delta, int nX, int nY, int refinementLevel=0);
  /// Construction of a cuboid vector version
  Cuboid2D(Vector<T,2> origin, T delta, Vector<int,2> extend, int refinementLevel=0);
  /// Copy constructor
  Cuboid2D(Cuboid2D<T> const& rhs, int overlap=0);
  /// Copy assignment
  Cuboid2D& operator=(Cuboid2D const& rhs);
  /// Initializes the cuboid
  void init(T globPosX, T globPosY, T delta, int nX, int nY, int refinementLevel=0);
  /// Read access to left lower corner coordinates
  T get_globPosX() const;
  T get_globPosY() const;
  /// Read only access to left lower corner coordinates
  Vector<T,2> const getOrigin() const;
  /// Read access to the distance of cuboid nodes
  T getDeltaR() const;
  /// Read access to cuboid width
  int getNx() const;
  /// Read access to cuboid height
  int getNy() const;
  /// Read only access to the number of voxels in every dimension
  Vector<int,2> const getExtend() const;

  /// Returns the volume of the cuboid
  T getPhysVolume() const;
  /// Returns the number of Nodes in the volume
  size_t getLatticeVolume() const;
  /// Returns the perimeter of the cuboid
  T getPhysPerimeter() const;
  /// Returns the number of Nodes at the perimeter
  int getLatticePerimeter() const;
  /// Prints cuboid details
  void print() const;

  void getPhysR(T physR[2], const int latticeR[2]) const;
  void getPhysR(T physR[2], const int& iX, const int& iY) const;

  void getLatticeR(int latticeR[2], const T physR[2]) const;
  void getLatticeR(int latticeR[2], const Vector<T,2>& physR) const;

  void getFloorLatticeR(const std::vector<T>& physR, std::vector<int>& latticeR) const
  {
    getFloorLatticeR(&latticeR[0], &physR[0]);
  }

  void getFloorLatticeR(int latticeR[2], const T physR[2]) const
  {
    latticeR[0] = (int)floor( (physR[0] - _globPosX)/_delta);
    latticeR[1] = (int)floor( (physR[1] - _globPosY)/_delta);
  }

  /// Returns the number of full cells
  size_t getWeight() const;
  /// Sets the number of full cells
  void setWeight(size_t fullCells);

  /// Checks whether a point (globX/globY) is contained in the cuboid
  /// extended with an layer of size overlap*delta
  bool checkPoint(T globX, T globY, int overlap = 0) const;
  /// Checks whether a point (globX/gloxY) is contained and is a node
  /// in the cuboid extended with an layer of size overlap*delta and
  /// returns the local active node
  bool checkPoint(T globX, T globY, int &locX, int &locY, int overlap = 0) const;
  /// Checks whether there is an intersection with the cuboid extended
  /// with an layer of size overlap*delta
  bool checkInters(T globX0, T globX1, T globY0, T globY1, int overlap = 0) const;
  /// Checks whether a given point intersects the cuboid extended
  /// by a layer of size overlap*delta
  bool checkInters(T globX, T globY, int overlap = 0) const;
  /// Checks whether there is an intersection and returns the local
  /// active node range which can be empty by means of locX0=1, locX1=0,
  /// locY0=1, locY1=0 of the cuboid extended with an layer of size
  /// overlap*delta
  bool checkInters(T globX0, T globX1, T globY0, T globY1,
                   int &locX0, int &locX1, int &locY0, int &locY1,
                   int overlap = 0) const;
  /// Divides the cuboid in p*q cuboids and adds them to the given vector
  void divide(int p, int q, std::vector<Cuboid2D<T> > &childrenC) const;
  /// Divides the cuboid in p cuboids and add them to the given vector
  void divide(int p, std::vector<Cuboid2D<T> > &childrenC) const;

  /// resize the cuboid to the passed size
  void resize(int X, int Y, int nX, int nY);

  /// Refines the cuboid with given refinement level
  void refineToLevel(unsigned int level);
  /// Refines one more level
  void refineIncrease();
  /// Refines one less level
  void refineDecrease();

  int get_refinementLevel() const;
};

}  // namespace olb

#endif
