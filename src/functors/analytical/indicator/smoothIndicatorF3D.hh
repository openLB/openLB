/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Mathias J. Krause, Albert Mink
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

#ifndef SMOOTH_INDICATOR_F_3D_HH
#define SMOOTH_INDICATOR_F_3D_HH

#include <vector>
#include <cmath>
#include <sstream>

#include "smoothIndicatorF3D.h"
#include "smoothIndicatorBaseF3D.h"
#include "smoothIndicatorCalcF3D.h"
#include "utilities/vectorHelpers.h"
#include "functors/analytical/interpolationF3D.h"
#include "functors/lattice/reductionF3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI2
#define M_PI2 1.57079632679489661923
#endif

namespace olb {

//Constructor: SmoothIndicatorCuboid3D
template <typename T, typename S, bool HLBM>
SmoothIndicatorCuboid3D<T,S,HLBM>::SmoothIndicatorCuboid3D(Vector<S,3> center, S xLength, S yLength, S zLength, S epsilon, Vector<S,3> theta, S density, Vector<S,3> vel)
  : _xLength(xLength),_yLength(yLength),_zLength(zLength)
{
  this->_pos = center;
  this->_circumRadius = .5*(sqrt(xLength*xLength+yLength*yLength+zLength*zLength))+0.5*epsilon;
  this->_myMin = {
    center[0] - this->getCircumRadius(), 
    center[1] - this->getCircumRadius(),
    center[2] - this->getCircumRadius()
  };
  this->_myMax = {
    center[0] + this->getCircumRadius(),
    center[1] + this->getCircumRadius(),
    center[2] + this->getCircumRadius()
  };
  this->_epsilon = epsilon;
  this->_theta = {
    theta[0] * (M_PI/180.),
    theta[1] * (M_PI/180.),
    theta[2] * (M_PI/180.)
  };
  T mass = xLength*yLength*zLength*density;
  T xLength2 = xLength*xLength;
  T yLength2 = yLength*yLength;
  T zLength2 = zLength*zLength;
  Vector<S,3> mofi;
  mofi[0] = 1./12.*mass*(yLength2+zLength2);
  mofi[1] = 1./12.*mass*(xLength2+zLength2);
  mofi[2] = 1./12.*mass*(yLength2+xLength2);
  this->init(this->_theta, vel, mass, mofi);
}

template <typename T, typename S, bool HLBM>
bool SmoothIndicatorCuboid3D<T,S,HLBM>::operator()(T output[], const S input[])
{
  //1.Calculate distance between point and center of unrotated indicator
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];
  T zDist = input[2] - this->getPos()[2];

  //2.Calculate point projected to rotated indicator
  // counter-clockwise rotation by _theta=-theta around center
  T x = this->getPos()[0] + this->getRotationMatrix()[0]*xDist + this->getRotationMatrix()[3]*yDist + this->getRotationMatrix()[6]*zDist;
  T y = this->getPos()[1] + this->getRotationMatrix()[1]*xDist + this->getRotationMatrix()[4]*yDist + this->getRotationMatrix()[7]*zDist;
  T z = this->getPos()[2] + this->getRotationMatrix()[2]*xDist + this->getRotationMatrix()[5]*yDist + this->getRotationMatrix()[8]*zDist;

  //3.Calculate distance between projected point and rotated indicator bounds
  xDist = fabs(x-this->getPos()[0]) - 0.5*(_xLength-this->getEpsilon());
  yDist = fabs(y-this->getPos()[1]) - 0.5*(_yLength-this->getEpsilon());
  zDist = fabs(z-this->getPos()[2]) - 0.5*(_zLength-this->getEpsilon());

  //4.Evaluate distance
  if (xDist <= 0 && yDist <= 0 && zDist <= 0) {
    output[0] = 1.;
    return true;
  }
  if (xDist >= this->getEpsilon() || yDist >= this->getEpsilon() || zDist >= this->getEpsilon()) {
    output[0] = 0.;
    return false;
  }
  //Evaluate epsilon on edges and borders
  T dist2 = 0.;
  T epsilon2 = this->getEpsilon()*this->getEpsilon();
  if (xDist < this->getEpsilon() && xDist > 0) {
    dist2 += xDist*xDist;
  }
  if (yDist < this->getEpsilon() && yDist > 0) {
    dist2 += yDist*yDist;
  }
  if (zDist < this->getEpsilon() && zDist > 0) {
    dist2 += zDist*zDist;
  }
  if (dist2 > 0 && dist2 <= epsilon2) {
    output[0] = T(cos(M_PI2*sqrt(dist2)/this->getEpsilon())*cos(M_PI2*sqrt(dist2)/this->getEpsilon()));
    return true;
  }
  output[0] = 0.;
  return false;
}


//Constructor: SmoothIndicatorSphere3D
template <typename T, typename S, bool HLBM>
SmoothIndicatorSphere3D<T,S,HLBM>::SmoothIndicatorSphere3D(Vector<S,3> center,
    S radius, S epsilon, S density, Vector<S,3> vel)
  : _radius(radius)
{
  this->_pos = center;
  this->_circumRadius = radius + 0.5*epsilon;
  this->_myMin = {
    center[0] - this->getCircumRadius(), 
    center[1] - this->getCircumRadius(),
    center[2] - this->getCircumRadius()
  };
  this->_myMax = {
    center[0] + this->getCircumRadius(),
    center[1] + this->getCircumRadius(),
    center[2] + this->getCircumRadius()
  };
  this->_epsilon = epsilon;

  T mass = 4./3.*M_PI*std::pow(radius, 3.)*density;
  T radius2 = radius * radius;
  Vector<S,3> mofi;
  mofi[0] = 2./5.*mass*radius2;
  mofi[1] = 2./5.*mass*radius2;
  mofi[2] = 2./5.*mass*radius2;
  Vector<S,3> theta(0.,0.,0.);
  this->init(theta, vel, mass, mofi);
}

template <typename T, typename S, bool HLBM>
bool SmoothIndicatorSphere3D<T, S, HLBM>::operator()(T output[], const S input[])
{
  //1.Calculate distance between point and center of indicator
  T distToCenter2 = (this->getPos()[0]-input[0])*(this->getPos()[0]-input[0]) +
                         (this->getPos()[1]-input[1])*(this->getPos()[1]-input[1]) +
                         (this->getPos()[2]-input[2])*(this->getPos()[2]-input[2]); 
  
  //3.Calculate distance between point and indicator bounds
  T rDist = sqrt(distToCenter2) - (_radius-this->getEpsilon()*0.5);

  //4. Evaluate distance
  if (rDist <= 0) {
    output[0] = 1.;
    return true;
  } else if (rDist > 0 && rDist < this->getEpsilon()) {
    output[0] = T(cos(M_PI2*rDist/this->getEpsilon())*cos(M_PI2*rDist/this->getEpsilon()));
    return true;
  }
  output[0] = 0.;
  return false;
}


template <typename T, typename S, bool HLBM>
void SmoothIndicatorCylinder3D<T,S,HLBM>::initIndicatorCylinder3D(Vector<S,3> normal, Vector<S,3> theta, S density, Vector<S,3> vel)
{ 
  T normLength = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
  T dx = normal[0]*(_length/normLength); 
  T dy = normal[1]*(_length/normLength); 
  T dz = normal[2]*(_length/normLength); 

  //Rotate according to normal orientation 
  theta[0] -= asin(dy/_length);
  if (dz >= 0){
    theta[1] += asin(dx/_length);
  } else {  
    theta[1] += M_PI;
    theta[1] -= asin(dx/_length);
  }
  
  this->_circumRadius = std::sqrt(_radius*_radius+(0.5*_length)*(0.5*_length))+0.5*this->getEpsilon();
  this->_myMin = {
    this->_pos[0] - this->getCircumRadius(),
    this->_pos[1] - this->getCircumRadius(),
    this->_pos[2] - this->getCircumRadius()
  };
  this->_myMax = {
    this->_pos[0] + this->getCircumRadius(),
    this->_pos[1] + this->getCircumRadius(),
    this->_pos[2] + this->getCircumRadius()
  };
  
  T radius2 = _radius * _radius;
  T mass = M_PI*radius2*_length*density; 
  Vector<S,3> mofi;
  mofi[0] = 0.5*mass*radius2;
  mofi[1] = 1/12.*mass*(_length*_length+3.*radius2);
  mofi[2] = 1/12.*mass*(_length*_length+3.*radius2);
  
  this->init(theta, vel, mass, mofi);
}

//Constructor 1: SmoothIndicatorCylinder3D
template <typename T, typename S, bool HLBM>
SmoothIndicatorCylinder3D<T,S,HLBM>::SmoothIndicatorCylinder3D(Vector<S,3> pointA, Vector<S,3> pointB, S radius, S epsilon, Vector<S,3> theta, S density, Vector<S,3> vel)
 : _radius(radius) 
{
  this->_epsilon = epsilon;
  this->_pos[0] = 0.5*(pointA[0]+pointB[0]);
  this->_pos[1] = 0.5*(pointA[1]+pointB[1]);
  this->_pos[2] = 0.5*(pointA[2]+pointB[2]);
  T dx = pointB[0]-pointA[0];
  T dy = pointB[1]-pointA[1];
  T dz = pointB[2]-pointA[2];
  theta[0] *= M_PI/180.;
  theta[1] *= M_PI/180.;
  theta[2] *= M_PI/180.;
  this->_theta = theta;
  _length = sqrt( dx*dx + dy*dy + dz*dz );
  Vector<S,3> normal (dx/_length, dy/_length, dz/_length);
  initIndicatorCylinder3D(normal, theta, density, vel);
}

//Constructor 2: SmoothIndicatorCylinder3D
template <typename T, typename S, bool HLBM>
SmoothIndicatorCylinder3D<T,S,HLBM>::SmoothIndicatorCylinder3D(Vector<S,3> center, Vector<S,3> normal, S radius, S length, S epsilon, Vector<S,3> theta, S density, Vector<S,3> vel)
: _radius(radius), _length(length) 
{ 
  this->_pos = center;
  this->_epsilon = epsilon;
  theta[0] *= M_PI/180.;
  theta[1] *= M_PI/180.;
  theta[2] *= M_PI/180.;
  this->_theta = theta;
  initIndicatorCylinder3D(normal, theta, density, vel);
}


template <typename T, typename S, bool HLBM>
bool SmoothIndicatorCylinder3D<T,S,HLBM>::operator()(T output[],const S input[])
{
  //1.Calculate distance between point and center of unrotated (z aligned) indicator
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];
  T zDist = input[2] - this->getPos()[2];

  //2.Calculate point projected to rotated indicator
  // counter-clockwise rotation by _theta=-theta around center
  T x= this->getPos()[0] + this->getRotationMatrix()[0]*xDist + this->getRotationMatrix()[3]*yDist + this->getRotationMatrix()[6]*zDist;
  T y= this->getPos()[1] + this->getRotationMatrix()[1]*xDist + this->getRotationMatrix()[4]*yDist + this->getRotationMatrix()[7]*zDist;
  T z= this->getPos()[2] + this->getRotationMatrix()[2]*xDist + this->getRotationMatrix()[5]*yDist + this->getRotationMatrix()[8]*zDist;

  //3.Calculate distance between projected point and indicator bounds
  T xyDistToCenter2 = (this->getPos()[0]-x)*(this->getPos()[0]-x)
                    + (this->getPos()[1]-y)*(this->getPos()[1]-y);
  T rDist = sqrt(xyDistToCenter2) - (_radius-this->getEpsilon()*0.5);
  zDist = fabs(z -this-> getPos()[2]) - 0.5*(_length-this->getEpsilon());

  //4.Evaluate distance
  if ( zDist <= 0 && rDist <= 0) {
    output[0] = 1.;
    return true;
  }  
  if (zDist >= this->getEpsilon() || rDist >= this->getEpsilon()) {
    output[0] = 0.;
    return false;
  }
  //Evaluate epsilon on edges and borders
  T dist2 = 0.;
  T epsilon2 = this->getEpsilon()*this->getEpsilon();
  if (zDist < this->getEpsilon() && zDist > 0) {
    dist2 += zDist*zDist;
  }
  if (rDist < this->getEpsilon() && rDist > 0) {
    dist2 += rDist*rDist;
  }
  if (dist2 > 0 && dist2 < epsilon2) {
    output[0] = T(cos(M_PI2*sqrt(dist2)/this->getEpsilon())*cos(M_PI2*sqrt(dist2)/this->getEpsilon()));
    return true;
  }
  output[0] = 0.;
  return false;
}


template <typename T, typename S, bool HLBM>
void SmoothIndicatorCone3D<T,S,HLBM>::initIndicatorCone3D(Vector<S,3> normal, Vector<S,3> theta, S density, Vector<S,3> vel)
{ 
  T normLength = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
  T dx = normal[0]*(_length/normLength); 
  T dy = normal[1]*(_length/normLength); 
  T dz = normal[2]*(_length/normLength); 
  
  //Rotate according to normal orientation 
  theta[0] -= asin(dy/_length);
  if (dz >= 0){
    theta[1] += asin(dx/_length);
  } else {  
    theta[1] += M_PI;
    theta[1] -= asin(dx/_length);
  }
   
  T radiusA2 = _radiusA *_radiusA;
  T radiusB2 = _radiusB *_radiusB;
  T mass = (1/3.)*M_PI*(radiusA2+_radiusA*_radiusB+radiusB2)*_length*density; 
  Vector<S,3> mofi;

  //FOR NOW ONLY VALID FOR SIMPLE CONE. NOT FOR TRUNCATED ONE 
  mofi[0] = 3/10.*mass*radiusA2;
  mofi[1] = 3/20.*mass*(4*_length*_length+radiusA2);
  mofi[2] = 3/20.*mass*(4*_length*_length+radiusA2);

  if (_radiusA >= _radiusB){
    this->_circumRadius = std::sqrt(_radiusA*_radiusA+(0.5*_length)*(0.5*_length))+0.5*this->getEpsilon();
  } else {
    this->_circumRadius = std::sqrt(_radiusB*_radiusB+(0.5*_length)*(0.5*_length))+0.5*this->getEpsilon();
  }

  this->_myMin = {
    this->_pos[0] - this->getCircumRadius(),
    this->_pos[1] - this->getCircumRadius(),
    this->_pos[2] - this->getCircumRadius()
  };
  this->_myMax = {
    this->_pos[0] + this->getCircumRadius(),
    this->_pos[1] + this->getCircumRadius(),
    this->_pos[2] + this->getCircumRadius()
  };
  
  this->init(theta, vel, mass, mofi);
}

//Constructor 1: SmoothIndicatorCone3D
template <typename T, typename S, bool HLBM>
SmoothIndicatorCone3D<T,S,HLBM>::SmoothIndicatorCone3D(Vector<S,3> pointA, Vector<S,3> pointB, S radiusA, S radiusB, S epsilon, Vector<S,3> theta, S density, Vector<S,3> vel)
 : _radiusA(radiusA), _radiusB(radiusB)
{
  this->_epsilon = epsilon;
  this->_pos[0] = 0.5*(pointA[0]+pointB[0]);
  this->_pos[1] = 0.5*(pointA[1]+pointB[1]);
  this->_pos[2] = 0.5*(pointA[2]+pointB[2]);
  T dx = pointB[0]-pointA[0];
  T dy = pointB[1]-pointA[1];
  T dz = pointB[2]-pointA[2];
  theta[0] *= M_PI/180.;
  theta[1] *= M_PI/180.;
  theta[2] *= M_PI/180.;
  this->_theta = theta;
  _length = sqrt( dx*dx + dy*dy + dz*dz );
  Vector<S,3> normal (dx/_length, dy/_length, dz/_length);
  initIndicatorCone3D(normal, theta, density, vel);
}

//Constructor 2: SmoothIndicatorCone3D
template <typename T, typename S, bool HLBM>
SmoothIndicatorCone3D<T,S,HLBM>::SmoothIndicatorCone3D(Vector<S,3> center, Vector<S,3> normal, S radiusA, S radiusB, S length, S epsilon, Vector<S,3> theta, S density, Vector<S,3> vel)
: _radiusA(radiusA), _radiusB(radiusB), _length(length)
{ 
  this->_pos = center;
  this->_epsilon = epsilon;
  theta[0] *= M_PI/180.;
  theta[1] *= M_PI/180.;
  theta[2] *= M_PI/180.;
  this->_theta = theta;
  initIndicatorCone3D(normal, theta, density, vel);
}

template <typename T, typename S, bool HLBM>
bool SmoothIndicatorCone3D<T,S,HLBM>::operator()(T output[],const S input[])
{
  //1.Calculate distance between point and center of unrotated (z aligned) indicator
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];
  T zDist = input[2] - this->getPos()[2];

  //2.Calculate point projected to rotated indicator
  // counter-clockwise rotation by _theta=-theta around center
  T x= this->getPos()[0] + this->getRotationMatrix()[0]*xDist + this->getRotationMatrix()[3]*yDist + this->getRotationMatrix()[6]*zDist;
  T y= this->getPos()[1] + this->getRotationMatrix()[1]*xDist + this->getRotationMatrix()[4]*yDist + this->getRotationMatrix()[7]*zDist;
  T z= this->getPos()[2] + this->getRotationMatrix()[2]*xDist + this->getRotationMatrix()[5]*yDist + this->getRotationMatrix()[8]*zDist;

  //3.Calculate distance between projected point and indicator bounds
  zDist = fabs(z -this-> getPos()[2]) - 0.5*(_length-this->getEpsilon()); 
  T axDist = z -this-> getPos()[2];
  T radiusAtZ = _radiusA - (axDist/_length+0.5)*(_radiusA-_radiusB);
  T xyDistToCenter2 = (this->getPos()[0]-x)*(this->getPos()[0]-x)
                     + (this->getPos()[1]-y)*(this->getPos()[1]-y);
  T rDist = sqrt(xyDistToCenter2) - (radiusAtZ-this->getEpsilon()*0.5);


   //4.Evaluate distance
  if ( zDist <= 0 && rDist <= 0) {
    output[0] = 1.;
    return true;
  }
  if (zDist >= this->getEpsilon() || rDist >= this->getEpsilon()) {
    output[0] = 0.;
    return false;
  }
  //Evaluate epsilon on edges and borders
  T dist2 = 0.;
  T epsilon2 = this->getEpsilon()*this->getEpsilon();
  if (zDist < this->getEpsilon() && zDist > 0) {
    dist2 += zDist*zDist;
  }
  if (rDist < this->getEpsilon() && rDist > 0) {
    dist2 += rDist*rDist;
  }
  if (dist2 > 0 && dist2 < epsilon2) {
    output[0] = T(cos(M_PI2*sqrt(dist2)/this->getEpsilon())*cos(M_PI2*sqrt(dist2)/this->    getEpsilon()));
    return true;
  }
  output[0] = 0.;
  return false;
}

} // namespace olb

#endif
