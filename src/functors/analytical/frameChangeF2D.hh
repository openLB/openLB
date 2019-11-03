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

#ifndef ANALYTICAL_FRAME_CHANGE_F_2D_HH
#define ANALYTICAL_FRAME_CHANGE_F_2D_HH

#include<vector>
#include<cmath>
#include<string>

#include "frameChangeF2D.h"
#include "frameChangeF3D.h"
#include "functors/genericF.h"
#include "analyticalF.h"
#include "functors/lattice/superBaseF2D.h"
#include "geometry/superGeometry2D.h"
#include "functors/lattice/superLatticeLocalF2D.h"
#include "core/superLattice2D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "utilities/vectorHelpers.h"  // for normalize

using namespace olb::util;

namespace olb {

template <typename T>
PowerLaw2D<T>::PowerLaw2D(std::vector<T> axisPoint, std::vector<T> axisDirection, T maxVelocity, T radius, T exponent) : AnalyticalF2D<T,T>(2)
{
  this->getName() = "PowerLaw2D";
  _axisPoint.resize(2);
  _axisDirection.resize(2);
  for (int i = 0; i < 2; ++i) {
    _axisDirection[i] = axisDirection[i];
    _axisPoint[i] = axisPoint[i];
  }
  _maxVelocity = maxVelocity;
  _radius = radius;
  _exponent = exponent;
}

template <typename T>
PowerLaw2D<T>::PowerLaw2D(SuperGeometry2D<T>& superGeometry, int material, T maxVelocity, T distance2Wall, T exponent) : AnalyticalF2D<T,T>(2)
{
  this->getName() = "PowerLaw2D";
  _axisPoint.resize(2);
  _axisDirection.resize(2);
  _axisPoint = superGeometry.getStatistics().getCenterPhysR(material);
  std::vector<int> discreteNormal = superGeometry.getStatistics().computeDiscreteNormal(material);
  for (int i = 0; i < 2; ++i) {
    _axisDirection[i] = discreteNormal[i];
  }

  _radius = T(distance2Wall);
  for (int iD = 0; iD < 2; iD++) {
    _radius += (superGeometry.getStatistics().getPhysRadius(material)[iD]);
  }
  _maxVelocity = maxVelocity;
  _exponent = exponent;
}


template <typename T>
bool PowerLaw2D<T>::operator()(T output[], const T x[])
{
  T d = fabs(_axisDirection[1]*(x[0] - _axisPoint[0]) - _axisDirection[0]*(x[1] - _axisPoint[1]));
  output[0] = _maxVelocity*_axisDirection[0]*(1. - pow(d/_radius,_exponent));
  output[1] = _maxVelocity*_axisDirection[1]*(1. - pow(d/_radius,_exponent));
  if ( 1. - pow(d/_radius,_exponent)  < 0.) {
    output[0] = T();
    output[1] = T();
  }
  return true;
}


template <typename T>
Poiseuille2D<T>::Poiseuille2D(std::vector<T> axisPoint, std::vector<T> axisDirection, T maxVelocity, T radius) : PowerLaw2D<T>(axisPoint, axisDirection, maxVelocity, radius, 2)
{
  this->getName() = "Poiseuille2D";
}


template <typename T>
Poiseuille2D<T>::Poiseuille2D(SuperGeometry2D<T>& superGeometry, int material, T maxVelocity, T distance2Wall) : PowerLaw2D<T>(superGeometry, material, maxVelocity, distance2Wall, 2)
{
  this->getName() = "Poiseuille2D";
}


template <typename T, typename S, typename DESCRIPTOR>
PoiseuilleStrainRate2D<T,S,DESCRIPTOR>::PoiseuilleStrainRate2D(UnitConverter<T,DESCRIPTOR> const& converter, T ly)
  : AnalyticalF2D<T,S>(4), _lengthY(ly), _maxVelocity(converter.getCharPhysVelocity())
{
  this->getName() = "PoiseuilleStrainRate2D";
}


template <typename T, typename S, typename DESCRIPTOR>
bool PoiseuilleStrainRate2D<T,S,DESCRIPTOR>::operator()(T output[], const S input[])
{
  T y = input[1];

  T DuDx = T();
  T DuDy = (T) _maxVelocity*(-2.*(y-(_lengthY/2.))/((_lengthY/2.)*(_lengthY/2.)));
  T DvDx = T();
  T DvDy = T();

  output[0] = (DuDx + DuDx)/2.;
  output[1] = (DuDy + DvDx)/2.;
  output[2] = (DvDx + DuDy)/2.;
  output[3] = (DvDy + DvDy)/2.;

  return true;
};

template <typename T>
AnalyticalPorousVelocity2D<T>::AnalyticalPorousVelocity2D(std::vector<T> axisDirection_, T K_, T mu_, T gradP_, T radius_, T eps_)
  :  AnalyticalF2D<T,T>(2), axisDirection(axisDirection_), K(K_), mu(mu_), gradP(gradP_),radius(radius_), eps(eps_)
{
  this->getName() = "AnalyticalPorousVelocity2D";
};


template <typename T>
T AnalyticalPorousVelocity2D<T>::getPeakVelocity()
{
  T uMax = K / mu*gradP*(1. - 1./(cosh((sqrt(1./K))*radius)));

  return uMax/eps;
};


template <typename T>
bool AnalyticalPorousVelocity2D<T>::operator()(T output[], const T input[])
{
  output[0] = K / mu*gradP*(1. - (cosh((sqrt(1./K))*(input[1] - radius)))/(cosh((sqrt(1./K))*radius)));
  output[1] = K / mu*gradP*(1. - (cosh((sqrt(1./K))*(input[0] - radius)))/(cosh((sqrt(1./K))*radius)));

  output[0] *= axisDirection[0]/eps;
  output[1] *= axisDirection[1]/eps;

  return true;
};


////////// Polar & Cartesian ///////////////


// constructor to obtain Cartesian coordinates of polar coordinates,
// with _polarOrigin, in x-y-plane
template<typename T, typename S>
PolarToCartesian2D<T, S>::PolarToCartesian2D(std::vector<T> polarOrigin)
  : AnalyticalF2D<T, S>(2),
    _polarOrigin(polarOrigin)
{
  this->getName() = "const";
}

// operator returns Cartesian coordinates of polar coordinates,
// input is x[0] = radius , x[1] = phi
template<typename T, typename S>
bool PolarToCartesian2D<T, S>::operator()(T output[], const S x[])
{
  output[0] = x[0]*cos(x[1]) + _polarOrigin[0];
  output[1] = x[0]*sin(x[1]) + _polarOrigin[1];
  return true;
}


// constructor to obtain polar coordinates of Cartesian coordinates
template<typename T, typename S>
CartesianToPolar2D<T, S>::CartesianToPolar2D(T cartesianOriginX,
    T cartesianOriginY, T cartesianOriginZ,
    T axisDirectionX, T axisDirectionY, T axisDirectionZ,
    T orientationX, T orientationY, T orientationZ)
  : AnalyticalF2D<T, S>(2)
{
  _cartesianOrigin.push_back(cartesianOriginX);
  _cartesianOrigin.push_back(cartesianOriginY);
  _cartesianOrigin.push_back(cartesianOriginZ);

  _axisDirection.push_back(axisDirectionX);
  _axisDirection.push_back(axisDirectionY);
  _axisDirection.push_back(axisDirectionZ);

  _orientation.push_back(orientationX);
  _orientation.push_back(orientationY);
  _orientation.push_back(orientationZ);
}

template<typename T, typename S>
CartesianToPolar2D<T, S>::CartesianToPolar2D(std::vector<T> cartesianOrigin,
    std::vector<T> axisDirection,
    std::vector<T> orientation)
  : AnalyticalF2D<T, S>(2),
    _cartesianOrigin(cartesianOrigin),
    _axisDirection(axisDirection),
    _orientation(orientation)
{
  this->getName() = "const";
}

// operator returns polar coordinates, output[0] = radius(>= 0),
// output[1] = phi in [0, 2Pi)
template<typename T, typename S>
bool CartesianToPolar2D<T, S>::operator()(T output[], const S x[])
{
  //Vector<T, 3> x_rot(x[this->getSourceDim()]); //// CAUSES ERROR
  int dim = this->getSourceDim();
  //std::cout<<"Source dim return: "<<dim<<std::endl; //sourceDim = 2
  Vector<T, 3> x_rot;
  for (int i = 0; i < dim; ++i) {
    x_rot[i] = x[i];
  }
  Vector<T, 3> e3(T(), T(), 1);
  Vector<T, 3> axisDirection(_axisDirection);

  Vector<T, 3> normal = crossProduct3D(axisDirection, e3);
  Vector<T, 3> normalAxisDir(axisDirection);
  normalAxisDir.normalize();

  // if axis has to be rotated
  if (!( util::nearZero(normalAxisDir[0]) && util::nearZero(normalAxisDir[1]) && util::nearZero(normalAxisDir[2]-1) ) ) {

    if ( !util::nearZero(util::norm(_orientation)) ) {
      normal = _orientation;
    }

    // normal is orientation direction
    AngleBetweenVectors3D<T, S> angle(fromVector3(e3), fromVector3(normal));
    T tmp[3] = {_axisDirection[0],_axisDirection[1],_axisDirection[2]};
    T alpha[1] = {};
    angle(alpha, tmp);

    // cross is rotation axis
    //Vector<T, 3> cross = crossProduct3D(e3, axisDirection);
    // rotation with angle alpha to rotAxisDir
    RotationRoundAxis3D<T, S> rotRAxis(_cartesianOrigin, fromVector3(normal), -alpha[0]);
    T x_tmp[3] = {};
    x_tmp[0] = x[0];
    x_tmp[1] = x[1];
    x_tmp[2] = x[2];
    T tmp2[3] = {};
    rotRAxis(tmp2,x_tmp);
    x_rot[0] = tmp2[0];
    x_rot[1] = tmp2[1];
    x_rot[2] = tmp2[2];
  }

  // calculates r, phi related to cartesianAxisDirection and their origin
  Vector<T, 2> difference;
  difference[0] = x_rot[0] - _cartesianOrigin[0];
  difference[1] = x_rot[1] - _cartesianOrigin[1];
  T distance = T();

  for (int i = 0; i < 2; ++i) {
    distance += difference[i]*difference[i];
  }

  distance = sqrt(distance);
  T phi[1] = {};

  if (distance > T()) {
    Vector<T, 3> e1(T(1), T(), T());
    AngleBetweenVectors3D<T, S> angle(fromVector3(e1), fromVector3(e3));

    T x_help[3] = {difference[0],difference[1],0.};
    angle(phi, x_help);
  }

  output[0] = distance;
  output[1] = phi[0];
  return true;
}

} // end namespace olb
#endif
