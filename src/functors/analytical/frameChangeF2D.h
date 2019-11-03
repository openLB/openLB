/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2015 Mathias J. Krause, Marie-Luise Maier
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

#ifndef ANALYTICAL_FRAME_CHANGE_F_2D_H
#define ANALYTICAL_FRAME_CHANGE_F_2D_H

#include<vector>
#include<cmath>
#include<string>
#include"math.h"

#include "functors/genericF.h"
#include "analyticalF.h"
#include "core/superLattice2D.h"

/** \file
  This file contains two different classes of functors, in the FIRST part
  - for simulations in a rotating frame
  - different functors for
      velocity (3d, RotatingLinear3D),
      pressure (1d, RotatingQuadratic1D) and
      force    (3d, RotatingForceField3D)
  The functors return the displacement of a point x in a fixed amount of time.

  The ones in the SECOND part are useful to set Poiseuille velocity profiles on
  - pipes with round cross-section and
  - pipes with square-shaped cross-section.
*/

/** To enable simulations in a rotating frame, the axis is set in the
  * constructor with axisPoint and axisDirection. The axisPoint can be the
  * coordinate of any point on the axis. The axisDirection has to be a normed to
  * 1. The pulse w is in rad/s. It determines the pulse speed by its norm while
  * the trigonometric or clockwise direction is determined by its sign: When the
  * axisDirection is pointing "towards you", a positive pulse makes it turn in
  * the trigonometric way. It has to be noticed that putting both axisDirection
  * into -axisDirection and w into -w yields an exactly identical situation.
  */


namespace olb {

template <typename T>
class PowerLaw2D : public AnalyticalF2D<T,T> {
protected:
  std::vector<T> _axisPoint;
  std::vector<T> _axisDirection;
  T _maxVelocity;
  T _radius;
  T _exponent;

public:
  PowerLaw2D(std::vector<T> axisPoint, std::vector<T> axisDirection,  T maxVelocity, T radius, T exponent);
  /// construct from material number, note: untested
  PowerLaw2D(SuperGeometry2D<T>& superGeometry, int material, T maxVelocity, T distance2Wall, T exponent);
  bool operator()(T output[], const T input[]) override;
};

template <typename T>
class Poiseuille2D : public PowerLaw2D<T> {

public:
  Poiseuille2D(std::vector<T> axisPoint, std::vector<T> axisDirection,  T maxVelocity, T radius);
  /// construct from material number, note: untested
  Poiseuille2D(SuperGeometry2D<T>& superGeometry, int material, T maxVelocity, T distance2Wall);
  //bool operator()(T output[], const T x[]);
};


template <typename T, typename S, typename DESCRIPTOR>
class PoiseuilleStrainRate2D : public AnalyticalF2D<T,S> {
protected:
  T _lengthY;
  T _maxVelocity;

public:
  PoiseuilleStrainRate2D(UnitConverter<T,DESCRIPTOR> const& converter, T ly);
  bool operator()(T output[], const S input[]) override;
};

/// Analytical solution of porous media channel flow with low Reynolds number
/// See Spaid and Phelan (doi:10.1063/1.869392)
template <typename T>
class AnalyticalPorousVelocity2D : public AnalyticalF2D<T,T> {
protected:
  std::vector<T> axisDirection;
  T K, mu, gradP, radius;
  T eps;
public:
  AnalyticalPorousVelocity2D(std::vector<T> axisDirection_, T K_, T mu_, T gradP_, T radius_, T eps_=T(1));
  T getPeakVelocity();
  bool operator()(T output[], const T input[]) override;
};


////////// Polar & Cartesian ///////////////


/// This class converts polar coordinates of point x (x[0] = radius, x[1] = phi)
/// to Cartesian coordinates (wrote into output field).
/// Initial situation for the Cartesian coordinate system is that angle phi
/// lies in the x-y-plane and turns round the _polarOrigin in math. pos sense.
template <typename T, typename S>
class PolarToCartesian2D : public AnalyticalF2D<T,S> {
protected:
  /// origin of the polar coordinate system to which point x is related
  std::vector<T> _polarOrigin;
  //  std::vector<T> _normalToPlane;
public:
  /// constructor to obtain Cartesian coordinates of polar coordinates
  PolarToCartesian2D(std::vector<T> polarOrigin);
  /// operator writes Cartesian coordinates of polar coordinates
  /// x[0] = radius >= 0, x[1] = phi in [0, 2*Pi) into output field
  bool operator()(T output[], const S x[]) override;
};


/// This class converts Cartesian coordinates of point x to polar coordinates
/// wrote into output field
/// (output[0] = radius>= 0, output[1] = phi in [0, 2Pi).
/// Initial situation for the polar coordinate system is that angle phi
/// lies in plane perpendicular to the _axisDirection and turns around
/// the _cartesianOrigin. The radius is the distance of point x to the
/// _axisDirection.
template <typename T, typename S>
class CartesianToPolar2D : public AnalyticalF2D<T,S> {
protected:
  /// origin of the Cartesian coordinate system
  std::vector<T> _cartesianOrigin;
  /// direction of the axis along which the polar coordinates are calculated
  std::vector<T> _axisDirection;
  /// direction to know orientation for math positive to obtain angle phi
  /// of Cartesian point x
  std::vector<T> _orientation;
public:
  CartesianToPolar2D(std::vector<T> cartesianOrigin,
                     std::vector<T> axisDirection,
                     std::vector<T> orientation = {T(1),T(),T()});
  CartesianToPolar2D(T cartesianOriginX, T cartesianOriginY,
                     T cartesianOriginZ,
                     T axisDirectionX, T axisDirectionY, T axisDirectionZ,
                     T orientationX = T(1), T orientationY = T(),
                     T orientationZ = T());
  /// operator writes polar coordinates of Cartesian point x into output
  /// field,
  /// returns output[0] = radius ( >= 0 ), output[1] = phi in [0, 2Pi)
  bool operator()(T output[], const S x[]) override;
};

} // end namespace olb
#endif
