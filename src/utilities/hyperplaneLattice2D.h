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

#ifndef HYPERPLANE_LATTICE_2D_H
#define HYPERPLANE_LATTICE_2D_H

#include "core/vector.h"
#include "geometry/cuboidGeometry2D.h"
#include "hyperplane2D.h"

namespace olb {

/// Parametrization of a hyperplane lattice (i.e. a line lattice).
/**
 * This class provides a common interface for describing how to discretize the intersection
 * of a hyperplane given by Hyperplane2D<T> and the mother cuboid of CuboidGeometry2D<T>.
 **/
template <typename T>
class HyperplaneLattice2D {
private:
  CuboidGeometry2D<T>&  _geometry;

  /// \return max possible distance
  int computeMaxLatticeDistance() const;
  /// Compute _hyperplane.origin, _n so that the cuboid is right inside the geometry
  void constructCuboid(int maxLatticeDistance);
  /// Update _h, _n _hyperplane.u so that the length matches the given resolution
  void setToResolution(int resolution);

protected:
  const Hyperplane2D<T> _hyperplane;

  /// Origin vector of the lattice
  /**
   * Note that this origin is set to a outermost point of the intersection between
   * cuboid geometry and hyperplane. Thus it is different from the Hyperplane2D<T>
   * origin vector in the general case.
   **/
  Vector<T,2> _origin;
  /// Direction vector of the lattice, normalized to grid width _h
  Vector<T,2> _u;

  /// Distance between discrete lattice points
  T _h;
  /// Number of lattice points in the direction of _u
  int _n;

public:
  /// Constructor for automatic discretization.
  /**
   * i.e. the grid width is set to CuboidGeometry2D<T>::getMinDeltaR.
   **/
  HyperplaneLattice2D(CuboidGeometry2D<T>& geometry,
                      Hyperplane2D<T>      hyperplane);
  /// Constructor for discretization of a given resolution.
  HyperplaneLattice2D(CuboidGeometry2D<T>& geometry,
                      Hyperplane2D<T>      hyperplane,
                      int                  resolution);
  /// Constructor for discretization of a given grid width.
  HyperplaneLattice2D(CuboidGeometry2D<T>& geometry,
                      Hyperplane2D<T>      hyperplane,
                      T                    h);

  HyperplaneLattice2D(const HyperplaneLattice2D&) = default;

  const Hyperplane2D<T>& getHyperplane() const;

  /// Transform 1d lattice coordinates to their physical 2d location
  Vector<T,2> getPhysR(const int& n) const;

  /// \return _n
  int getN() const;
  /// \return _h
  T getPhysSpacing() const;

  /// \return _origin
  Vector<T,2> getPhysOrigin() const;
  /// \return _u
  Vector<T,2> getVectorU() const;

};

}

#endif
