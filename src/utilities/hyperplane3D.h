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

#ifndef HYPERPLANE_3D_H
#define HYPERPLANE_3D_H

#include "core/vector.h"
#include "geometry/cuboid3D.h"

namespace olb {

/// Definition of a analytical 2D plane embedded in 3D space
/**
 * Hyperplane3D defines a hyperplane using its origin and two span vectors.
 *
 * In practice it might be preferable to define a hyperplane using a normal
 * vector or to automatically center the origin in a given cuboid.
 * For this purpose a fluent construction interface is offered:
 * \code{.cpp}
 * // construct a hyperplane positioned at (1,1,1) and normal to (1,0,0)
 * auto plane = Hyperplane3D<T>().originAt({1,1,1})
 *                               .normalTo({1,0,0});
 * \endcode
 *
 * The primary reason for this development was the increasing constructor
 * clutter in BlockReduction3D2D: Instead of using the correct Vector<T,3> types
 * for passing span and origin vectors they were passed as a mix between raw values
 * and array types to prevent the constructor from becoming ambiguous. e.g. on the
 * type level passing two span vectors is indistinguishable from passing normal and
 * origin vectors.
 **/
template <typename T>
struct Hyperplane3D {
  Vector<T,3> origin;
  Vector<T,3> u;
  Vector<T,3> v;
  Vector<T,3> normal;

  Hyperplane3D() = default;

  /// Center the hyperplane at the given origin vector
  /// \return Hyperplane3D reference for further construction
  Hyperplane3D& originAt(const Vector<T,3>& origin);
  /// Center the hyperplane relative to the given cuboid
  /// \return Hyperplane3D reference for further construction
  Hyperplane3D& centeredIn(const Cuboid3D<T>& cuboid);
  /// Span the hyperplane using two span vectors
  /// \return Hyperplane3D reference for further construction
  Hyperplane3D& spannedBy(const Vector<T,3>& u, const Vector<T,3>& v);
  /// Calculate the spanning vectors of the hyperplane to be orthogonal to the given normal
  /// \return Hyperplane3D reference for further construction
  Hyperplane3D& normalTo(const Vector<T,3>& normal);

  /// Apply a matrix given by its row vectors to both span vectors
  /// \return Hyperplane3D reference for further construction
  Hyperplane3D& applyMatrixToSpan(const Vector<T,3>& row0,
                                  const Vector<T,3>& row1,
                                  const Vector<T,3>& row2);
  /// Rotate the spanning vectors around the X axis
  /// \return Hyperplane3D reference for further construction
  Hyperplane3D& rotateSpanAroundX(T r);
  /// Rotate the spanning vectors around the Y axis
  /// \return Hyperplane3D reference for further construction
  Hyperplane3D& rotateSpanAroundY(T r);
  /// Rotate the spanning vectors around the Z axis
  /// \return Hyperplane3D reference for further construction
  Hyperplane3D& rotateSpanAroundZ(T r);

  /// \return true iff normal is orthogonal to X, Y axis
  bool isXYPlane() const;
  /// \return true iff normal is orthogonal to X, Z axis
  bool isXZPlane() const;
  /// \return true iff normal is orthogonal to Y, Z axis
  bool isYZPlane() const;
};

}

#endif
