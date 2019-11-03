/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2015 Gilles Zahnd, Mathias J. Krause
 *  Marie-Luise Maier
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

#ifndef FRAME_CHANGE_F_3D_H
#define FRAME_CHANGE_F_3D_H

#include<vector>
#include<string>

#include "analyticalF.h"

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


template<typename T> class SuperGeometry3D;

// PART 1: /////////////////////////////////////////////////////////////////////
// Functors for rotating the coordinate system (velocity, pressure, force,...)

/**
  * This functor gives a linar profile for a given point x as it computes
  * the distance between x and the axis.
  *
  * The field in outcome is the velocity field of q rotating solid
  */

/// Functor with a linear profile e.g. for rotating velocity fields.
template <typename T>
class RotatingLinear3D final : public AnalyticalF3D<T,T> {
protected:
  std::vector<T> axisPoint;
  std::vector<T> axisDirection;
  T w;
  T scale;
public:
  RotatingLinear3D(std::vector<T> axisPoint_, std::vector<T> axisDirection_, T w_, T scale_=1);
  bool operator()(T output[], const T x[]) override;
};

/**
  * This functor gives a linar profile in an annulus for a given point x between the inner and outer radius as it computes
  * the distance between x and the inner and outer radius.
  *
  * The field in outcome is the velocity field of q rotating solid in an annulus
  */

/// Functor with a linear profile e.g. for rotating velocity fields.
template <typename T>
class RotatingLinearAnnulus3D final : public AnalyticalF3D<T,T> {
protected:
  std::vector<T> axisPoint;
  std::vector<T> axisDirection;
  T w;
  T ri;
  T ro;
  T scale;
public:
  RotatingLinearAnnulus3D(std::vector<T> axisPoint_, std::vector<T> axisDirection_, T w_, T ri_, T ro_, T scale_=1);
  bool operator()(T output[], const T x[]);
};


/**
  * This functor gives a parabolic profile for a given point x as it computes
  * the distance between x and the axis.
  *
  * This field is a scalar field, a vector with one component will be used
  */

/// Functor with a parabolic profile e.g. for rotating pressure fields.
template <typename T>
class RotatingQuadratic1D final : public AnalyticalF3D<T,T> {
protected:
  std::vector<T> axisPoint;
  std::vector<T> axisDirection;
  T w;
  T scale;
  T additive;
public:
  RotatingQuadratic1D(std::vector<T> axisPoint_, std::vector<T> axisDirection_,
                      T w_, T scale_=1, T additive_=0);
  bool operator()(T output[], const T x[]) override;
};


// PART 2: /////////////////////////////////////////////////////////////////////
// Functors for setting velocities on a velocity boundary of a pipe

/**
  * This functor returns a quadratic Poiseuille profile for use with a pipe with
  * round cross-section. It uses cylinder coordinates and is valid for the
  * entire length of the pipe.
  *
  * This functor gives a parabolic velocity profile for a given point x as it
  * computes the distance between x and the axis.
  *
  * The axis is set in the input with axisPoint and axisDirection. The axisPoint
  * can be the coordinate of any point where the axis passes.
  * axisDirection has to be normed to 1.
  * Once the axis is set in the middle of the pipe, the radius of the
  * pipe "radius" and the velocity in the middle of the pipe "maxVelocity"
  * determine the Poisseuille profile entierly.
  */

/// Velocity profile for round pipes and power law fluids: u(r)=u_max*(1-(r/R)^((n+1)/n)). The exponent n characterizes the fluid behavior.
/// n<1: Pseudoplastic, n=1: Newtonian fluid, n>1: Dilatant
template <typename T>
class CirclePowerLaw3D : public AnalyticalF3D<T,T> {
protected:
  std::vector<T> _center;
  std::vector<T> _normal;
  T _radius;
  T _maxVelocity;
  T _n;
  T _scale;

public:
  CirclePowerLaw3D(std::vector<T> axisPoint, std::vector<T> axisDirection,  T maxVelocity, T radius, T n, T scale = T(1));
  CirclePowerLaw3D(T center0, T center1, T center2, T normal0, T normal1, T normal2, T radius, T maxVelocity, T n, T scale = T(1));
  CirclePowerLaw3D(SuperGeometry3D<T>& superGeometry, int material, T maxVelocity, T n, T scale = T(1));

  CirclePowerLaw3D(bool useMeanVelocity, std::vector<T> axisPoint, std::vector<T> axisDirection,  T Velocity, T radius, T n, T scale = T(1));
  CirclePowerLaw3D(bool useMeanVelocity, T center0, T center1, T center2, T normal0, T normal1, T normal2, T radius, T Velocity, T n, T scale = T(1));
  CirclePowerLaw3D(bool useMeanVelocity, SuperGeometry3D<T>& superGeometry, int material, T Velocity, T n, T scale = T(1));

  /// Returns centerpoint vector
  std::vector<T> getCenter()
  {
    return _center;
  };
  /// Returns normal vector
  std::vector<T> getNormal()
  {
    return _normal;
  };
  /// Returns radi
  T getRadius()
  {
    return _radius;
  };

  bool operator()(T output[], const T x[]) override;
};

/// Velocity profile for round pipes and turbulent flows: u(r)=u_max*(1-r/R)^(1/n) The exponent n can be calculated by n = 1.03 * ln(Re) - 3.6
/// n=7 is used for many flow applications
template <typename T>
class CirclePowerLawTurbulent3D : public CirclePowerLaw3D<T> {
public:
  CirclePowerLawTurbulent3D(std::vector<T> axisPoint_, std::vector<T> axisDirection,  T maxVelocity, T radius, T n, T scale = T(1));
  CirclePowerLawTurbulent3D(T center0, T center1, T center2, T normal0, T normal1, T normal2, T radius, T maxVelocity, T n, T scale = T(1));
  CirclePowerLawTurbulent3D(SuperGeometry3D<T>& superGeometry, int material, T maxVelocity, T n, T scale = T(1));

  CirclePowerLawTurbulent3D(bool useMeanVelocity, std::vector<T> axisPoint, std::vector<T> axisDirection,  T Velocity, T radius, T n, T scale = T(1));
  CirclePowerLawTurbulent3D(bool useMeanVelocity, T center0, T center1, T center2, T normal0, T normal1, T normal2, T radius, T Velocity, T n, T scale = T(1));
  CirclePowerLawTurbulent3D(bool useMeanVelocity, SuperGeometry3D<T>& superGeometry, int material, T Velocity, T n, T scale = T(1));

  bool operator()(T output[], const T x[]) override;
};

/// Velocity profile for round pipes and a laminar flow of a Newtonian fluid: u(r)=u_max*(1-(r/R)^2)
template <typename T>
class CirclePoiseuille3D final : public CirclePowerLaw3D<T> {

public:
  CirclePoiseuille3D(std::vector<T> axisPoint, std::vector<T> axisDirection,  T maxVelocity, T radius, T scale = T(1));
  CirclePoiseuille3D(T center0, T center1, T center2, T normal0, T normal1, T normal2, T radius, T maxVelocity, T scale = T(1));
  CirclePoiseuille3D(SuperGeometry3D<T>& superGeometry, int material, T maxVelocity, T scale = T(1));

  CirclePoiseuille3D(bool useMeanVelocity, std::vector<T> axisPoint, std::vector<T> axisDirection,  T Velocity, T radius, T scale = T(1));
  CirclePoiseuille3D(bool useMeanVelocity, T center0, T center1, T center2, T normal0, T normal1, T normal2, T radius, T Velocity, T scale = T(1));
  CirclePoiseuille3D(bool useMeanVelocity, SuperGeometry3D<T>& superGeometry, int material, T Velocity, T scale = T(1));
};

/// Strain rate for round pipes and laminar flow of a Newtonian fluid
template < typename T,typename DESCRIPTOR>
class CirclePoiseuilleStrainRate3D : public AnalyticalF3D<T,T> {
protected:
 T lengthY;
 T lengthZ;
 T maxVelocity;

public:
 CirclePoiseuilleStrainRate3D(UnitConverter<T, DESCRIPTOR> const& converter, T ly);
 bool operator()(T output[], const T input[]);
};

/**
  This functor returns a Poiseuille profile for use with a pipe with square shaped cross-section.
*/
/// velocity profile for pipes with rectangular cross-section.
template <typename T>
class RectanglePoiseuille3D final : public AnalyticalF3D<T,T> {
protected:
  OstreamManager clout;
  std::vector<T> x0;
  std::vector<T> x1;
  std::vector<T> x2;
  std::vector<T> maxVelocity;

public:
  RectanglePoiseuille3D(std::vector<T>& x0_, std::vector<T>& x1_, std::vector<T>& x2_, std::vector<T>& maxVelocity_);
  /// constructor from material numbers
  /** offsetX,Y,Z is a positive number denoting the distance from border
    * voxels of material_ to the zerovelocity boundary */
  RectanglePoiseuille3D(SuperGeometry3D<T>& superGeometry_, int material_, std::vector<T>& maxVelocity_, T offsetX, T offsetY, T offsetZ);
  bool operator()(T output[], const T x[]) override;
};


/**
  * This functor returns a quadratic Poiseuille profile for use with a pipe with
  * elliptic cross-section.
  *
  * Ellpise is orthogonal to z-direction.
  * a is in x-direction
  * b is in y-direction
  *
  */
/// Velocity profile for round pipes.
template <typename T>
class EllipticPoiseuille3D final : public AnalyticalF3D<T,T> {
protected:
  OstreamManager clout;
  std::vector<T> _center;
  T _a2, _b2;
  T _maxVel;

public:
  EllipticPoiseuille3D(std::vector<T> center, T a, T b, T maxVel);
  bool operator()(T output[], const T x[]) override;
};


/// Analytical solution of porous media channel flow with low Reynolds number
/// See Spaid and Phelan (doi:10.1063/1.869392)
template <typename T>
class AnalyticalPorousVelocity3D : public AnalyticalF3D<T,T> {
protected:
  std::vector<T> center;
  std::vector<T> normal;
  T K, mu, gradP, radius;
  T eps;
public:
  AnalyticalPorousVelocity3D(SuperGeometry3D<T>& superGeometry, int material, T K_, T mu_, T gradP_, T radius_, T eps_=T(1));
  T getPeakVelocity();
  bool operator()(T output[], const T input[]) override;
};


////////////////Coordinate Transformation//////////////


/// This class calculates the angle alpha between vector _referenceVector and
/// any other vector x.
/// Vector x has to be turned by alpha in mathematical positive sense depending
/// to _orientation to lie over vector _referenceVector
template <typename T, typename S>
class AngleBetweenVectors3D final : public AnalyticalF3D<T,S> {
protected:
  /// between _referenceVector and vector x, angle is calculated
  std::vector<T> _referenceVector;
  /// direction around which x has to be turned with angle is in math pos sense
  /// to lie over _referenceVector
  std::vector<T> _orientation;
public:
  /// constructor defines _referenceVector and _orientation
  AngleBetweenVectors3D(std::vector<T> referenceVector,
                        std::vector<T> orientation);
  /// operator writes angle between x and _referenceVector inro output field
  bool operator()(T output[], const S x[]) override;
};


/// This class saves coordinates of the resulting point after rotation in the
/// output field. The resulting point is achieved after rotation of x with
/// angle _alpha in math positive sense relative to _rotAxisDirection with a
/// given _origin.
template <typename T, typename S>
class RotationRoundAxis3D final : public AnalyticalF3D<T,S> {
protected:
  /// origin, around which x is turned
  std::vector<T> _origin;
  /// direction, around which x is turned
  std::vector<T> _rotAxisDirection;
  /// angle, by which vector x is rotated around _rotAxisDirection
  T _alpha;
public:
  /// constructor defines _rotAxisDirection, _alpha and _origin
  RotationRoundAxis3D(std::vector<T> origin, std::vector<T> rotAxisDirection,
                      T alpha);
  /// operator writes coordinates of the resulting point after rotation
  /// of x by angle _alpha in math positive sense to _rotAxisDirection
  /// with _origin into output field
  bool operator()(T output[], const S x[]) override;
};


////////// Cylinder & Cartesian ///////////////


/// This class converts cylindrical of point x (x[0] = radius, x[1] = phi,
/// x[2] = z) to Cartesian coordinates (wrote into output field).
/// Initial situation for the Cartesian coordinate system is that angle phi
/// lies in the x-y-plane and turns round the _cylinderOrigin, the z-axis
/// direction equals the axis direction of the cylinder.
template <typename T, typename S>
class CylinderToCartesian3D final : public AnalyticalF3D<T,S> {
protected:
  std::vector<T> _cylinderOrigin;
public:
  /// constructor defines _cylinderOrigin
  CylinderToCartesian3D(std::vector<T> cylinderOrigin);
  /// operator writes Cartesian coordinates of cylinder coordinates
  /// x[0] = radius >= 0, x[1] = phi in [0, 2*Pi), x[2] = z into output field
  bool operator()(T output[], const S x[]) override;
};


/// This class converts Cartesian coordinates of point x to cylindrical
/// coordinates wrote into output field
/// (output[0] = radius, output[1] = phi, output[2] = z).
/// Initial situation for the cylindrical coordinate system is that angle phi
/// lies in plane perpendicular to the _axisDirection and turns around the
/// _cartesianOrigin. The radius is the distance of point x to the
/// _axisDirection and z the distance to _cartesianOrigin along _axisDirection.
template <typename T, typename S>
class CartesianToCylinder3D final : public AnalyticalF3D<T,S> {
protected:
  /// origin of the Cartesian coordinate system
  std::vector<T> _cartesianOrigin;
  /// direction of the axis along which the cylindrical coordinates are calculated
  std::vector<T> _axisDirection;
  /// direction to know orientation for math positive to obtain angle phi
  /// of Cartesian point x
  std::vector<T> _orientation;
public:
  CartesianToCylinder3D(std::vector<T> cartesianOrigin, T& angle,
                        std::vector<T> orientation = {T(1),T(),T()});
  CartesianToCylinder3D(std::vector<T> cartesianOrigin,
                        std::vector<T> axisDirection,
                        std::vector<T> orientation = {T(1),T(),T()});
  CartesianToCylinder3D(T cartesianOriginX, T cartesianOriginY,
                        T cartesianOriginZ,
                        T axisDirectionX, T axisDirectionY,
                        T axisDirectionZ,
                        T orientationX = T(1), T orientationY = T(),
                        T orientationZ = T());
  /// operator writes cylindrical coordinates of Cartesian point x into output
  /// field,
  /// output[0] = radius ( >= 0), output[1] = phi (in [0, 2Pi)), output[2] = z
  bool operator()(T output[], const S x[]) override;
  /// Read only access to _axisDirection
  std::vector<T> const& getAxisDirection() const;
  /// Read and write access to _axisDirection
  std::vector<T>& getAxisDirection();
  /// Read only access to _catesianOrigin
  std::vector<T> const& getCartesianOrigin() const;
  /// Read and write access to _catesianOrigin
  std::vector<T>& getCartesianOrigin();
  /// Read and write access to _axisDirection using angle
  void setAngle(T angle);
};


////////// Spherical & Cartesian ///////////////


/// This class converts spherical coordinates of point x
/// (x[0] = radius, x[1] = phi, x[2] = theta) to Cartesian coordinates wrote
/// into output field.
/// Initial situation for the Cartesian coordinate system is that angle phi
/// (phi in [0, 2Pi)) lies in x-y-plane and turns around the Cartesian origin
/// (0,0,0). Angle theta (theta in [0, Pi]) lies in y-z-plane. z axis equals
/// the _orientation of the spherical coordinate system.
/// The radius is distance of point x to the Cartesian _origin
template <typename T, typename S>
class SphericalToCartesian3D final : public AnalyticalF3D<T,S> {
  //protected:
public:
  SphericalToCartesian3D();
  /// operator writes spherical coordinates of Cartesian coordinates of x
  /// (x[0] = radius, x[1] = phi, x[2] = theta) into output field.
  /// phi is in x-y-plane, theta in x-z-plane, z axis as orientation
  bool operator()(T output[], const S x[]) override;
};


/// This class converts Cartesian coordinates of point x to spherical
/// coordinates wrote into output field
/// (output[0] = radius, output[1] = phi, output[2] = theta).
/// Initial situation for the spherical coordinate system is that angle phi
/// lies in x-y-plane and theta in x-z-plane.
/// The z axis is the orientation for the mathematical positive sense of phi.
template <typename T, typename S>
class CartesianToSpherical3D final : public AnalyticalF3D<T,S> {
protected:
  /// origin of the Cartesian coordinate system
  std::vector<T> _cartesianOrigin;
  /// direction of the axis along which the spherical coordinates are calculated
  std::vector<T> _axisDirection;
public:
  CartesianToSpherical3D(std::vector<T> cartesianOrigin,
                         std::vector<T> axisDirection);//, std::vector<T> orientation);
  /// operator writes shperical coordinates of Cartesian point x into output
  /// field,
  /// output[0] = radius ( >= 0), output[1] = phi (in [0, 2Pi)),
  /// output[2] = theta (in [0, Pi])
  bool operator()(T output[], const S x[]) override;
};


////////// Magnetic Field ///////////////


/// Magnetic field that creates magnetization in wire has to be orthogonal to
/// the wire. To calculate the magnetic force on particles around a cylinder
/// (J. Lindner et al. / Computers and Chemical Engineering 54 (2013) 111-121)
///  Fm = mu0*4/3.*PI*radParticle^3*_Mp*radWire^2/r^3 *
///       [radWire^2/r^2+cos(2*theta), sin(2*theta), 0]
template <typename T, typename S>
class MagneticForceFromCylinder3D final : public AnalyticalF3D<T,S> {
protected:
  CartesianToCylinder3D<T,S>& _car2cyl;
  /// length of the wire, from origin to _car2cyl.axisDirection
  T _length;
  /// wire radius
  T _radWire;
  /// maximal distance from wire cutoff/threshold
  T _cutoff;
  /// saturation magnetization wire, linear scaling factor
  T _Mw;
  /// magnetization particle, linear scaling factor
  T _Mp;
  /// particle radius, cubic scaling factor
  T _radParticle;
  /// permeability constant, linear scaling factor
  T _mu0;
  /// factor = mu0*4/3.*PI*radParticle^3*_Mp*radWire^2/r^3
  T _factor;
public:
  MagneticForceFromCylinder3D(CartesianToCylinder3D<T,S>& car2cyl, T length,
                              T radWire, T cutoff, T Mw, T Mp = T(1),
                              T radParticle = T(1), T mu0 = 4*3.14159265e-7);
  /// operator writes the magnetic force in a point x round a cylindrical wire into output field
  bool operator()(T output[], const S x[]) override;
};

template <typename T, typename S>
class MagneticFieldFromCylinder3D final : public AnalyticalF3D<T, S> {
protected:
  CartesianToCylinder3D<T, S>& _car2cyl;
  /// length of the wire, from origin to _car2cyl.axisDirection
  T _length;
  /// wire radius
  T _radWire;
  /// maximal distance from wire cutoff/threshold
  T _cutoff;
  /// saturation magnetization wire, linear scaling factor
  T _Mw;
  /// factor = mu0*4/3.*PI*radParticle^3*_Mp*radWire^2/r^3
  T _factor;
public:
  MagneticFieldFromCylinder3D(CartesianToCylinder3D<T, S>& car2cyl,
                              T length,
                              T radWire,
                              T cutoff,
                              T Mw
                             );
  /// operator writes the magnetic force in a point x round a cylindrical wire into output field
  bool operator()(T output[], const S x[]);
};

} // end namespace olb

#endif
