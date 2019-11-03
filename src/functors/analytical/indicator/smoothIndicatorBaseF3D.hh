/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Albert Mink, Mathias J. Krause, Benjamin Förster
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

#ifndef SMOOTH_INDICATOR_BASE_F_3D_HH
#define SMOOTH_INDICATOR_BASE_F_3D_HH

#include <cmath>

#include "smoothIndicatorBaseF3D.h"
#include "utilities/vectorHelpers.h"

namespace olb {


template <typename T, typename S>
SmoothIndicatorF3D<T, S, false>::SmoothIndicatorF3D()
  : AnalyticalF3D<T,S>(1),
    _myMin(S()), _myMax(S()), _pos(S()), _rotMat(S()), _circumRadius(S()),
    _theta(S()), _epsilon(S())
{ }

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,false>::init(Vector<S,3> theta, Vector<S,3> vel, T mass, Vector<S,3> mofi) {
  this->_rotMat[0] = std::cos(theta[1])*std::cos(theta[2]);
  this->_rotMat[1] = std::sin(theta[0])*std::sin(theta[1])*std::cos(theta[2]) - std::cos(theta[0])*std::sin(theta[2]);
  this->_rotMat[2] = std::cos(theta[0])*std::sin(theta[1])*std::cos(theta[2]) + std::sin(theta[0])*std::sin(theta[2]);
  this->_rotMat[3] = std::cos(theta[1])*std::sin(theta[2]);
  this->_rotMat[4] = std::sin(theta[0])*std::sin(theta[1])*std::sin(theta[2]) + std::cos(theta[0])*std::cos(theta[2]);
  this->_rotMat[5] = std::cos(theta[0])*std::sin(theta[1])*std::sin(theta[2]) - std::sin(theta[0])*std::cos(theta[2]);
  this->_rotMat[6] = -std::sin(theta[1]);
  this->_rotMat[7] = std::sin(theta[0])*std::cos(theta[1]);
  this->_rotMat[8] = std::cos(theta[0])*std::cos(theta[1]);
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T, S, false>::getMin() const
{
  return _myMin;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T, S, false>::getMax() const
{
  return _myMax;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T,S,false>::getPos() const {
  return _pos;
}

template <typename T, typename S>
const S& SmoothIndicatorF3D<T,S,false>::getCircumRadius() const 
{
  return _circumRadius;
}

template <typename T, typename S>
const Vector<S,9>& SmoothIndicatorF3D<T,S,false>::getRotationMatrix() const {
  return _rotMat;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T,S,false>::getTheta() const {
  return _theta;
}

template <typename T, typename S>
const S& SmoothIndicatorF3D<T,S,false>::getEpsilon() const {
  return _epsilon;
}

template <typename T, typename S>
std::string SmoothIndicatorF3D<T,S,false>::name() {
  return _name;
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,false>::setTheta(Vector<S,3> theta) {
  _theta[0] = theta[0];
  _theta[1] = theta[1];
  _theta[2] = theta[2];
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,false>::setEpsilon(S epsilon) {
  _epsilon = epsilon;
}

// identity to "store results"
template <typename T, typename S>
SmoothIndicatorIdentity3D<T,S>::SmoothIndicatorIdentity3D(SmoothIndicatorF3D<T,S,false>& f)
  : _f(f) {
  this->_myMin = _f.getMin();
  this->_myMax = _f.getMax();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool SmoothIndicatorIdentity3D<T,S>::operator() (T output[], const S input[]) {
  _f(output, input);
  return true;
}

///////////////////////// for template specialisation HLBM=true

template <typename T, typename S>
SmoothIndicatorF3D<T,S,true>::SmoothIndicatorF3D()
  : AnalyticalF3D<T,S>(1),
    _myMin(S()), _myMax(S()), _pos(S()), _vel(S()), _acc(S()), _acc2(S()), _force(S()), 
    _rotMat(S()), _circumRadius(S()), _theta(S()), _omega(S()), _alpha(S()), 
    _alpha2(S()), _mass(S()), _mofi(S()), _epsilon(S())
{ }

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::init(Vector<S,3> theta, Vector<S,3> vel, T mass, Vector<S,3> mofi) {
  this->_rotMat[0] = std::cos(theta[1])*std::cos(theta[2]);
  this->_rotMat[1] = std::sin(theta[0])*std::sin(theta[1])*std::cos(theta[2]) - std::cos(theta[0])*std::sin(theta[2]);
  this->_rotMat[2] = std::cos(theta[0])*std::sin(theta[1])*std::cos(theta[2]) + std::sin(theta[0])*std::sin(theta[2]);
  this->_rotMat[3] = std::cos(theta[1])*std::sin(theta[2]);
  this->_rotMat[4] = std::sin(theta[0])*std::sin(theta[1])*std::sin(theta[2]) + std::cos(theta[0])*std::cos(theta[2]);
  this->_rotMat[5] = std::cos(theta[0])*std::sin(theta[1])*std::sin(theta[2]) - std::sin(theta[0])*std::cos(theta[2]);
  this->_rotMat[6] = -std::sin(theta[1]);
  this->_rotMat[7] = std::sin(theta[0])*std::cos(theta[1]);
  this->_rotMat[8] = std::cos(theta[0])*std::cos(theta[1]);
  _vel = vel;
  _mass = mass;
  _mofi = mofi;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T,S,true>::getMin() const {
  return _myMin;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T,S,true>::getMax() const {
  return _myMax;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T,S,true>::getPos() const {
  return _pos;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T,S,true>::getVel() const {
   return _vel;
 }

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T, S, true>::getAcc() const
{
  return _acc;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T, S, true>::getAcc2() const
{
  return _acc2;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T,S,true>::getForce() const {
  return _force;
}

template <typename T, typename S>
const Vector<S,9>& SmoothIndicatorF3D<T,S,true>::getRotationMatrix() const {
  return _rotMat;
}

template <typename T, typename S>
const S& SmoothIndicatorF3D<T,S,true>::getCircumRadius() const {
  return _circumRadius;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T, S, true>::getTheta() const
{
  return _theta;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T, S, true>::getOmega() const
{
  return _omega;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T, S, true>::getAlpha() const
{
  return _alpha;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T,S,true>::getAlpha2() const {
  return _alpha2;
}

template <typename T, typename S>
const S& SmoothIndicatorF3D<T, S, true>::getMass() const
{
  return _mass;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T, S, true>::getMofi() const
{
  return _mofi;
}

template <typename T, typename S>
const S& SmoothIndicatorF3D<T, S, true>::getEpsilon() const
{
  return _epsilon;
} 

template <typename T, typename S>
std::string SmoothIndicatorF3D<T,S,true>::name() {
  return _name;
}

// identity to "store results"
template <typename T, typename S>
void SmoothIndicatorF3D<T, S, true>::setPos(Vector<S, 3> pos) {
  _pos[0] = pos[0];
  _pos[1] = pos[1];
  _pos[2] = pos[2];
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setVel(Vector<S,3> vel) {
  _vel[0] = vel[0];
  _vel[1] = vel[1];
  _vel[2] = vel[2];
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setAcc(Vector<S,3> acc) {
  _acc[0] = acc[0];
  _acc[1] = acc[1];
  _acc[2] = acc[2];
} 

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setAcc2(Vector<S,3> acc2) {
  _acc2[0] = acc2[0];
  _acc2[1] = acc2[1];
  _acc2[2] = acc2[2];
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setForce(Vector<S,3> force) {
  _force[0] = force[0];
  _force[1] = force[1];
  _force[2] = force[2];
}


template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setRotationMatrix(Vector<S,9> rotMat) {
  _rotMat[0] = rotMat[0];
  _rotMat[1] = rotMat[1];
  _rotMat[2] = rotMat[2];
  _rotMat[3] = rotMat[3];
  _rotMat[4] = rotMat[4];
  _rotMat[5] = rotMat[5];
  _rotMat[6] = rotMat[6];
  _rotMat[7] = rotMat[7];
  _rotMat[8] = rotMat[8];
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setTheta(Vector<S,3> theta) {
  _theta[0] = theta[0];
  _theta[1] = theta[1];
  _theta[2] = theta[2];
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setOmega(Vector<S,3> omega) {
  _omega[0] = omega[0];
  _omega[1] = omega[1];
  _omega[2] = omega[2];
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setAlpha(Vector<S,3> alpha) {
  _alpha[0] = alpha[0];
  _alpha[1] = alpha[1];
  _alpha[2] = alpha[2];
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setAlpha2(Vector<S,3> alpha2) {
  _alpha2[0] = alpha2[0];
  _alpha2[1] = alpha2[1];
  _alpha2[2] = alpha2[2];
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setMass(S mass) {
  _mass = mass;
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setMofi(Vector<S, 3> mofi) {
  _mofi[0] = mofi[0];
  _mofi[1] = mofi[1];
  _mofi[2] = mofi[2];
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setEpsilon(S epsilon) {
  _epsilon = epsilon;
}

} // namespace olb


#endif
