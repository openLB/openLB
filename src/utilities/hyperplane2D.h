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

#ifndef HYPERPLANE_2D_H
#define HYPERPLANE_2D_H

#include "core/vector.h"
#include "geometry/cuboid2D.h"

namespace olb {

/// Definition of a analytical line embedded in 2D space
/**
 * Hyperplane2D defines a line using its origin and a direction vector.
 **/
template <typename T>
struct Hyperplane2D {
  Vector<T,2> origin;
  Vector<T,2> u;
  Vector<T,2> normal;

  Hyperplane2D() = default;

  /// Center the line at the given origin vector
  /// \return Hyperplane2D reference for further construction
  Hyperplane2D& originAt(const Vector<T,2>& origin);
  /// Center the line relative to the given cuboid
  /// \return Hyperplane2D reference for further construction
  Hyperplane2D& centeredIn(const Cuboid2D<T>& cuboid);
  /// Set the direction of the line parallel to a vector
  /// \return Hyperplane2D reference for further construction
  Hyperplane2D& parallelTo(const Vector<T,2>& direction);
  /// Calculate the direction vector of the line to be orthogonal to the given normal
  /// \return Hyperplane2D reference for further construction
  Hyperplane2D& normalTo(const Vector<T,2>& normal);

  /// \return true iff normal is orthogonal to X axis
  bool isParallelToX() const;
  /// \return true iff normal is orthogonal to Y axis
  bool isParallelToY() const;
};

}

#endif
