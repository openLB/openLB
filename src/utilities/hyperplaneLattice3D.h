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

#ifndef HYPERPLANE_LATTICE_3D_H
#define HYPERPLANE_LATTICE_3D_H

#include "core/vector.h"
#include "geometry/cuboidGeometry3D.h"
#include "hyperplane3D.h"

namespace olb {

/// Parametrization of a hyperplane lattice.
/**
 * i.e. the resolution / grid width of the discretization of a given hyperplane.
 *
 * This class provides a common interface for describing how to discretize the intersection
 * of a hyperplane given by Hyperplane3D<T> and the mother cuboid of CuboidGeometry3D<T>.
 **/
template <typename T>
class HyperplaneLattice3D {
private:
  /// \return max possible distance
  int computeMaxLatticeDistance(Cuboid3D<T>&& cuboid) const;
  /// Compute _hyperplane.origin, _nx, _ny so that the cuboid is right inside the geometry
  void constructCuboid(CuboidGeometry3D<T>& geometry, int maxLatticeDistance);
  /// Update _h, _nx, _ny, _hyperplane.(u,v) so that the longest side length matches the given resolution
  void setToResolution(int resolution);

protected:
  const Hyperplane3D<T> _hyperplane;

  /// Origin vector of the lattice
  /**
   * Note that this origin is set to a outermost point of the intersection between
   * cuboid geometry and hyperplane. Thus it is different from the Hyperplane3D<T>
   * origin vector in the general case.
   **/
  Vector<T,3> _origin;
  /// Span vector of the lattice, normalized to grid width _h
  Vector<T,3> _u;
  /// Span vector of the lattice, normalized to grid width _h
  Vector<T,3> _v;

  /// Distance between discrete lattice points
  T _h;
  /// Number of lattice points in the direction of _u
  int _nx;
  /// Number of lattice points in the direction of _v
  int _ny;

public:
  /// Constructor for automatic discretization.
  /**
   * i.e. the grid width is set to CuboidGeometry3D<T>::getMinDeltaR.
   **/
  HyperplaneLattice3D(CuboidGeometry3D<T>& geometry,
                      Hyperplane3D<T>      hyperplane);
  /// Constructor for discretization of a given resolution.
  HyperplaneLattice3D(CuboidGeometry3D<T>& geometry,
                      Hyperplane3D<T>      hyperplane,
                      int                  resolution);
  /// Constructor for discretization of a given grid width.
  HyperplaneLattice3D(CuboidGeometry3D<T>& geometry,
                      Hyperplane3D<T>      hyperplane,
                      T                    h);

  /// Constructor for manual discretization
  /**
   * \param hyperplane Hyperplane in 3D space
   * \param h          lattice point spacing
   * \param nx         X axis resolution
   * \param ny         Y axis resolution
   **/
  HyperplaneLattice3D(Hyperplane3D<T> hyperplane,
                      T h, int nx, int ny);

  HyperplaneLattice3D(const HyperplaneLattice3D&) = default;

  const Hyperplane3D<T>& getHyperplane() const;

  /// Transform 2d lattice coordinates to their physical 3d location
  Vector<T,3> getPhysR(const int& planeX, const int& planeY) const;

  /// \return _nx
  int getNx() const;
  /// \return _ny
  int getNy() const;
  /// \return _h
  T getPhysSpacing() const;

  /// \return _origin
  Vector<T,3> getPhysOrigin() const;
  /// \return _u
  Vector<T,3> getVectorU() const;
  /// \return _v
  Vector<T,3> getVectorV() const;

};

}

#endif
