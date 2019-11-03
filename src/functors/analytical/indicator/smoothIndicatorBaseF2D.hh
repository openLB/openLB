/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Mathias J. Krause, Benjamin Förster
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

#ifndef SMOOTH_INDICATOR_BASE_F_2D_HH
#define SMOOTH_INDICATOR_BASE_F_2D_HH

#include <cmath>

#include "smoothIndicatorBaseF2D.h"

namespace olb {


template <typename T, typename S>
SmoothIndicatorF2D<T,S,false>::SmoothIndicatorF2D()
  : AnalyticalF2D<T,S>(1),
    _myMin(S()), _myMax(S()), _pos(S()), _rotMat(S()), _circumRadius(S()), _theta(S()), _epsilon(S())
{ }

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,false>::init(T theta, Vector<S,2> vel, T mass, T mofi) {
  _rotMat[0] = std::cos(theta);
  _rotMat[1] = std::sin(theta);
  _rotMat[2] = -std::sin(theta);
  _rotMat[3] = std::cos(theta);
}

template <typename T, typename S>
const Vector<S,2>& SmoothIndicatorF2D<T,S,false>::getMin() const {
  return _myMin;
}

template <typename T, typename S>
const Vector<S,2>& SmoothIndicatorF2D<T,S,false>::getMax() const {
  return _myMax;
}

template <typename T, typename S>
const Vector<S,2>& SmoothIndicatorF2D<T,S,false>::getPos() const {
  return _pos;
}

template <typename T, typename S>
const Vector<S,4>& SmoothIndicatorF2D<T,S,false>::getRotationMatrix() const {
  return _rotMat;
}

template <typename T, typename S>
const S& SmoothIndicatorF2D<T,S,false>::getCircumRadius() const {
  return _circumRadius;
}

template <typename T, typename S>
const S& SmoothIndicatorF2D<T,S,false>::getTheta() const {
  return _theta;
}

template <typename T, typename S>
const S& SmoothIndicatorF2D<T,S,false>::getEpsilon() const {
  return _epsilon;
}

template <typename T, typename S>
std::string SmoothIndicatorF2D<T,S,false>::name() {
  return _name;
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,false>::setRotationMatrix(Vector<S,4> rotMat) {
  _rotMat[0] = rotMat[0];
  _rotMat[1] = rotMat[1];
  _rotMat[2] = rotMat[2];
  _rotMat[3] = rotMat[3];
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,false>::setTheta(S theta) {
  _theta = theta;
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,false>::setEpsilon(S epsilon) {
  _epsilon = epsilon;
}

// identity to "store results"
template <typename T, typename S>
SmoothIndicatorIdentity2D<T,S>::SmoothIndicatorIdentity2D(SmoothIndicatorF2D<T,S,false>& f)
  : _f(f) {
  this->_myMin = _f.getMin();
  this->_myMax = _f.getMax();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool SmoothIndicatorIdentity2D<T,S>::operator() (T output[], const S input[]) {
  _f(output, input);
  return true;
}

///////////////////////// for template specialisation HLBM=true

template <typename T, typename S>
SmoothIndicatorF2D<T,S,true>::SmoothIndicatorF2D()
  : AnalyticalF2D<T,S>(1),
    _myMin(S()), _myMax(S()), _pos(S()), _vel(S()), _acc(S()), _acc2(S()), _force(S()),
    _rotMat(S()), _circumRadius(S()), _theta(S()), _omega(S()), _alpha(S()), 
    _alpha2(S()), _mass(S()), _mofi(S()), _epsilon(S())
{ }

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::init(T theta, Vector<S,2> vel, T mass, T mofi) {
  _rotMat[0] = std::cos(theta);
  _rotMat[1] = std::sin(theta);
  _rotMat[2] = -std::sin(theta);
  _rotMat[3] = std::cos(theta);
  _vel = vel;
  _mass = mass;
  _mofi = mofi;
}

template <typename T, typename S>
const Vector<S,2>& SmoothIndicatorF2D<T,S,true>::getMin() const {
  return _myMin;
}

template <typename T, typename S>
const Vector<S,2>& SmoothIndicatorF2D<T,S,true>::getMax() const {
  return _myMax;
}

template <typename T, typename S>
const Vector<S,2>& SmoothIndicatorF2D<T,S,true>::getPos() const {
  return _pos;
}

template <typename T, typename S>
const Vector<S,2>& SmoothIndicatorF2D<T,S,true>::getVel() const {
  return _vel;
}

template <typename T, typename S>
const Vector<S,2>& SmoothIndicatorF2D<T,S,true>::getAcc() const {
  return _acc;
}

template <typename T, typename S>
const Vector<S,2>& SmoothIndicatorF2D<T,S,true>::getAcc2() const {
  return _acc2;
}

template <typename T, typename S>
const Vector<S,2>& SmoothIndicatorF2D<T,S,true>::getForce() const {
  return _force;
}

template <typename T, typename S>
const Vector<S,4>& SmoothIndicatorF2D<T,S,true>::getRotationMatrix() const {
  return _rotMat;
}

template <typename T, typename S>
const S& SmoothIndicatorF2D<T,S,true>::getCircumRadius() const {
  return _circumRadius;
}

template <typename T, typename S>
const S& SmoothIndicatorF2D<T,S,true>::getTheta() const {
  return _theta;
}

template <typename T, typename S>
const S& SmoothIndicatorF2D<T,S,true>::getOmega() const {
  return _omega;
}

template <typename T, typename S>
const S& SmoothIndicatorF2D<T,S,true>::getAlpha() const {
  return _alpha;
}

template <typename T, typename S>
const S& SmoothIndicatorF2D<T,S,true>::getAlpha2() const {
  return _alpha2;
}

template <typename T, typename S>
const S& SmoothIndicatorF2D<T,S,true>::getMass() const {
  return _mass;
}

template <typename T, typename S>
const S& SmoothIndicatorF2D<T,S,true>::getMofi() const {
  return _mofi;
}

template <typename T, typename S>
const S& SmoothIndicatorF2D<T,S,true>::getEpsilon() const {
  return _epsilon;
}

template <typename T, typename S>
std::string SmoothIndicatorF2D<T,S,true>::name() {
  return _name;
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setPos(Vector<S,2> pos) {
  _pos[0] = pos[0];
  _pos[1] = pos[1];
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setVel(Vector<S,2> vel) {
  _vel[0] = vel[0];
  _vel[1] = vel[1];
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setAcc(Vector<S,2> acc) {
  _acc[0] = acc[0];
  _acc[1] = acc[1];
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setAcc2(Vector<S,2> acc2) {
  _acc2[0] = acc2[0];
  _acc2[1] = acc2[1];
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setForce(Vector<S,2> force) {
  _force[0] = force[0];
  _force[1] = force[1];
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setRotationMatrix(Vector<S,4> rotMat) {
  _rotMat[0] = rotMat[0];
  _rotMat[1] = rotMat[1];
  _rotMat[2] = rotMat[2];
  _rotMat[3] = rotMat[3];
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setTheta(S theta) {
  _theta = theta;
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setOmega(S omega) {
  _omega = omega;
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setAlpha(S alpha) {
  _alpha = alpha;
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setAlpha2(S alpha2) {
  _alpha2 = alpha2;
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setMass(S mass) {
  _mass = mass;
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setMofi(S mofi) {
  _mofi = mofi;
}

template <typename T, typename S>
void SmoothIndicatorF2D<T,S,true>::setEpsilon(S epsilon) {
  _epsilon = epsilon;
}

} // namespace olb

#endif
