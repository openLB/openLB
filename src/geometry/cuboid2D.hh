/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
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

/** \file
 * The description of a single 2D cuboid -- generic implementation.
 */


#ifndef CUBOID_2D_HH
#define CUBOID_2D_HH

#include <iostream>
#include <math.h>
#include <vector>
#include "cuboid2D.h"
#include "dynamics/lbHelpers.h"
#include "io/ostreamManager.h"


namespace olb {

////////////////////// Class Cuboid2D /////////////////////////

template<typename T>
Cuboid2D<T>::Cuboid2D(T globPosX, T globPosY, T delta ,int nX, int nY, int refinementLevel)
  : _weight(std::numeric_limits<size_t>::max()), clout(std::cout,"Cuboid2D")
{
  init(globPosX, globPosY, delta, nX, nY, refinementLevel);
}

template<typename T>
Cuboid2D<T>::Cuboid2D(Vector<T,2> origin, T delta, Vector<int,2> extend, int refinementLevel)
  : _weight(std::numeric_limits<size_t>::max()), clout(std::cout,"Cuboid2D")
{
  this->init(origin[0], origin[1], delta, extend[0], extend[1], refinementLevel);
}

// copy constructor
template<typename T>
Cuboid2D<T>::Cuboid2D(Cuboid2D<T> const& rhs, int overlap) : clout(std::cout,"Cuboid2D")
{
  this->init(rhs._globPosX - rhs._delta*overlap, rhs._globPosY - rhs._delta*overlap,
             rhs._delta, rhs._nX + 2*overlap, rhs._nY + 2*overlap);
  _weight = rhs._weight;
  _refinementLevel = rhs._refinementLevel;
}

// copy assignment
template<typename T>
Cuboid2D<T>& Cuboid2D<T>::operator=(Cuboid2D<T> const& rhs)
{
  this->init(rhs._globPosX, rhs._globPosY, rhs._delta, rhs._nX, rhs._nY);
  _weight = rhs._weight;
  _refinementLevel = rhs._refinementLevel;
  return *this;
}

template<typename T>
void Cuboid2D<T>::init(T globPosX, T globPosY, T delta, int nX, int nY, int refinementLevel)
{
  _globPosX = globPosX;
  _globPosY = globPosY;
  _delta    = delta;
  _nX       = nX;
  _nY       = nY;
  _refinementLevel = refinementLevel;
}


template<typename T>
T Cuboid2D<T>::get_globPosX() const
{
  return _globPosX;
}

template<typename T>
T Cuboid2D<T>::get_globPosY() const
{
  return _globPosY;
}

template<typename T>
Vector<T,2> const Cuboid2D<T>::getOrigin() const
{
  return Vector<T,2> (_globPosX,_globPosY);
}

template<typename T>
T Cuboid2D<T>::getDeltaR() const
{
  return _delta;
}

template<typename T>
int Cuboid2D<T>::getNx() const
{
  return _nX;
}

template<typename T>
int Cuboid2D<T>::getNy() const
{
  return _nY;
}

template<typename T>
Vector<int,2> const Cuboid2D<T>::getExtend() const
{
  return Vector<int,2> (_nX, _nY);
}

template<typename T>
T Cuboid2D<T>::getPhysVolume() const
{
  return _nY*_nX*_delta*_delta;
}

template<typename T>
size_t Cuboid2D<T>::getLatticeVolume() const
{
  return static_cast<size_t>(_nY)*static_cast<size_t>(_nX);
}

template<typename T>
T Cuboid2D<T>::getPhysPerimeter() const
{
  return 2*_nY*_delta + 2*_nX*_delta;
}

template<typename T>
int Cuboid2D<T>::getLatticePerimeter() const
{
  return 2*_nY + 2*_nX -4;
}

template<typename T>
void Cuboid2D<T>::print() const
{
  clout << "--------Cuboid Details----------" << std::endl;
  clout << " Left Corner (x/y): " << "\t" << "(" << this->get_globPosX() << "/" << this->get_globPosY() << ")" << std::endl;
  clout << " Delta: " << "\t" << "\t" << this->getDeltaR() << std::endl;
  clout << " Perimeter: " << "\t" << "\t" << this->getPhysPerimeter() << std::endl;
  clout << " Volume: " << "\t" << "\t" << this->getPhysVolume() << std::endl;
  clout << " Extent (x/y): " << "\t" << "\t" << "(" << this->getNx() << "/" << this->getNy() << ")" << std::endl;
  clout << " Nodes at Perimeter: " << "\t" << this->getLatticePerimeter() << std::endl;
  clout << " Nodes in Volume: " << "\t" << this->getLatticeVolume() << std::endl;
  clout << "--------------------------------" << std::endl;
}

template<typename T>
void Cuboid2D<T>::getPhysR(T physR[2], const int latticeR[2]) const
{
  physR[0] = _globPosX + latticeR[0]*_delta;
  physR[1] = _globPosY + latticeR[1]*_delta;
}

template<typename T>
void Cuboid2D<T>::getPhysR(T physR[2], const int& iX, const int& iY) const
{
  physR[0] = _globPosX + iX*_delta;
  physR[1] = _globPosY + iY*_delta;
}

template<typename T>
void Cuboid2D<T>::getLatticeR(int latticeR[2], const T physR[2]) const
{
  latticeR[0] = (int)floor( (physR[0] - _globPosX )/_delta +.5);
  latticeR[1] = (int)floor( (physR[1] - _globPosY )/_delta +.5);
}

template<typename T>
void Cuboid2D<T>::getLatticeR(int latticeR[2], const Vector<T,2>& physR) const
{
  latticeR[0] = (int)floor( (physR[0] - _globPosX )/_delta +.5);
  latticeR[1] = (int)floor( (physR[1] - _globPosY )/_delta +.5);
}

template<typename T>
bool Cuboid2D<T>::checkPoint(T globX, T globY, int overlap) const
{

  if (_globPosX - T(0.5 + overlap)*_delta <= globX &&
      _globPosX + T(_nX-0.5+overlap)*_delta  > globX &&
      _globPosY - T(0.5 + overlap)*_delta <= globY &&
      _globPosY + T(_nY-0.5+overlap)*_delta  > globY ) {
    return true;
  } else {
    return false;
  }
}

template<typename T>
bool Cuboid2D<T>::checkPoint(T globX, T globY, int &locX, int &locY, int overlap) const
{
  if (overlap!=0) {
    Cuboid2D tmp(_globPosX - overlap*_delta, _globPosY - overlap*_delta, _delta , _nX + overlap*2, _nY + overlap*2);
    return tmp.checkPoint(globX, globY, locX, locY);
  } else if (!checkPoint(globX, globY)) {
    return false;
  } else {
    locX = (int)floor((globX - (T)_globPosX)/_delta + .5);
    locY = (int)floor((globY - (T)_globPosY)/_delta + .5);
    return true;
  }
}

template<typename T>
bool Cuboid2D<T>::checkInters(T globX0, T globX1, T globY0, T globY1, int overlap) const
{

  T locX0d = std::max(_globPosX-overlap*_delta,globX0);
  T locY0d = std::max(_globPosY-overlap*_delta,globY0);
  T locX1d = std::min(_globPosX+(_nX+overlap-1)*_delta,globX1);
  T locY1d = std::min(_globPosY+(_nY+overlap-1)*_delta,globY1);

  if (!(locX1d>=locX0d && locY1d>=locY0d)) {
    return false;
  }
  return true;
}

template<typename T>
bool Cuboid2D<T>::checkInters(T globX, T globY, int overlap) const
{
  return checkInters(globX, globX, globY, globY, overlap);
}

template<typename T>
bool Cuboid2D<T>::checkInters(T globX0, T globX1, T globY0, T globY1,
                              int &locX0, int &locX1, int &locY0, int &locY1, int overlap) const
{
  if (overlap!=0) {
    Cuboid2D tmp(_globPosX - overlap*_delta, _globPosY - overlap*_delta, _delta, _nX + overlap*2, _nY + overlap*2);
    return tmp.checkInters(globX0, globX1, globY0, globY1, locX0, locX1, locY0, locY1);
  } else if (!checkInters(globX0, globX1, globY0, globY1)) {
    locX0 = 1;
    locX1 = 0;
    locY0 = 1;
    locY1 = 0;
    return false;
  } else {
    locX0 = 0;
    for (int i=0; _globPosX + i*_delta < globX0; i++) {
      locX0 = i+1;
    }
    locX1 = _nX-1;
    for (int i=_nX-1; _globPosX + i*_delta > globX1; i--) {
      locX1 = i-1;
    }
    locY0 = 0;
    for (int i=0; _globPosY + i*_delta < globY0; i++) {
      locY0 = i+1;
    }
    locY1 = _nY-1;
    for (int i=_nY-1; _globPosY + i*_delta > globY1; i--) {
      locY1 = i-1;
    }
    return true;
  }
}

template<typename T>
void Cuboid2D<T>::divide(int nX, int nY, std::vector<Cuboid2D<T> > &childrenC) const
{
  T globPosX_child, globPosY_child;
  int xN_child = 0;
  int yN_child = 0;

  globPosX_child = _globPosX;
  globPosY_child = _globPosY;

  for (int iX=0; iX<nX; iX++) {
    for (int iY=0; iY<nY; iY++) {
      xN_child       = (_nX+nX-iX-1)/nX;
      yN_child       = (_nY+nY-iY-1)/nY;
      Cuboid2D<T> child(globPosX_child, globPosY_child, _delta, xN_child, yN_child, _refinementLevel);
      childrenC.push_back(child);
      globPosY_child += yN_child*_delta;
    }
    globPosY_child = _globPosY;
    globPosX_child += xN_child*_delta;
  }
}

template<typename T>
void Cuboid2D<T>::divide(int p, std::vector<Cuboid2D<T> > &childrenC) const
{

  OLB_PRECONDITION(p>0);
  int nX = 0;
  int nY = 0;
  T ratio;
  T bestRatio = (T)_nX/(T)_nY;
  T difRatio = fabs(bestRatio - 1) + 1;
  for (int i=1; i<= p; i++) {
    int j = p / i;
    if (i*j<=p) {
      if ( fabs(bestRatio - (T)i/(T)j) <= difRatio) {
        difRatio = fabs(bestRatio - (T)i/(T)j);
        nX = i;
        nY = j;
      }
    }
  }

  ratio = T(nX)/(T)nY;
  int rest = p - nX*nY;

  if (rest==0) {
    divide(nX,nY,childrenC);
    return;
  }

  if (ratio < bestRatio && (nY-rest) >= 0) {
    int n_QNoInsertions = nX*(nY-rest);
    T bestVolume_QNoInsertions = (T)_nX*_nY * n_QNoInsertions/(T)p;
    int yN_QNoInsertions = (int)(bestVolume_QNoInsertions / (T)_nX);
    int xN_QNoInsertions = _nX;
    int yN_QInsertions = _nY-yN_QNoInsertions;
    int xN_QInsertions = _nX;
    Cuboid2D<T> firstChildQ(_globPosX, _globPosY, _delta, xN_QNoInsertions, yN_QNoInsertions, _refinementLevel);
    Cuboid2D<T> secondChildQ(_globPosX, _globPosY+yN_QNoInsertions*_delta, _delta, xN_QInsertions, yN_QInsertions, _refinementLevel);
    firstChildQ.divide(nX, nY-rest, childrenC);
    secondChildQ.divide(nX+1,rest, childrenC);
  } else {
    int n_QNoInsertions = nY*(nX-rest);
    T bestVolume_QNoInsertions = (T)_nX*_nY * n_QNoInsertions/(T)p;
    int xN_QNoInsertions = (int)(bestVolume_QNoInsertions / (T)_nY + 0.9999);
    int yN_QNoInsertions = _nY;
    int xN_QInsertions = _nX-xN_QNoInsertions;
    int yN_QInsertions = _nY;
    Cuboid2D<T> firstChildQ(_globPosX, _globPosY, _delta, xN_QNoInsertions, yN_QNoInsertions, _refinementLevel);
    Cuboid2D<T> secondChildQ(_globPosX+xN_QNoInsertions*_delta, _globPosY, _delta, xN_QInsertions, yN_QInsertions, _refinementLevel);
    firstChildQ.divide(nX-rest, nY, childrenC);
    secondChildQ.divide(rest,nY+1, childrenC);
  }
}


template<typename T>
void Cuboid2D<T>::resize(int iX, int iY, int nX, int nY)
{
  _globPosX = _globPosX+iX*_delta;
  _globPosY = _globPosY+iY*_delta;
  _nX = nX;
  _nY = nY;
}

template<typename T>
void Cuboid2D<T>::refineToLevel(unsigned int level)
{
  int leveldiffabs = int (std::pow(2, level - _refinementLevel ));
  _delta /= (T)(leveldiffabs);
  _nX *= leveldiffabs;
  _nY *= leveldiffabs;
  _refinementLevel = level;
}

template<typename T>
void Cuboid2D<T>::refineIncrease()
{
  _delta /= (T)2;
  _nX = _nX*2;
  _nY = _nY*2;
  _refinementLevel++;
}

template<typename T>
void Cuboid2D<T>::refineDecrease()
{
  if (_refinementLevel == 0) {
    return;
  }
  _delta *= (T)2;
  _nX = (_nX)/2;
  _nY = (_nY)/2;
  _refinementLevel--;
}

template<typename T>
int Cuboid2D<T>::get_refinementLevel() const
{
  return _refinementLevel;
}

template<typename T>
size_t Cuboid2D<T>::getWeight() const
{
  if (_weight == std::numeric_limits<size_t>::max()) {
    return getLatticeVolume();
  } else {
    return _weight;
  }
}

template<typename T>
void Cuboid2D<T>::setWeight(size_t fullCells)
{
  _weight = fullCells;
}


}  // namespace olb

#endif
