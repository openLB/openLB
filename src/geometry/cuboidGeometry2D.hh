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
 * The description of a vector of 2D cuboid  -- generic implementation.
 */


#ifndef CUBOID_GEOMETRY_2D_HH
#define CUBOID_GEOMETRY_2D_HH


#include <iostream>
#include <math.h>
#include "geometry/cuboidGeometry2D.h"
#include "functors/analytical/indicator/indicatorF2D.h"

namespace olb {


////////////////////// Class CuboidGeometry2D /////////////////////////

template<typename T>
CuboidGeometry2D<T>::CuboidGeometry2D()
  : _motherCuboid(0,0,0,0,0), _periodicityOn(3, bool(false)), clout(std::cout, "CuboidGeometry3D")
{
  add(_motherCuboid);
  split(0, 1);
}

template<typename T>
CuboidGeometry2D<T>::CuboidGeometry2D(T originX, T originY, T deltaR, int nX, int nY, int nC)
  : _motherCuboid(originX, originY, deltaR, nX, nY), _periodicityOn(2, bool(false)),
    clout(std::cout, "CuboidGeometry2D")
{
  add(_motherCuboid);
  split(0, nC);
}

template<typename T>
CuboidGeometry2D<T>::CuboidGeometry2D(IndicatorF2D<T>& indicatorF, T voxelSize, int nC)
  : _motherCuboid(indicatorF.getMin()[0],  indicatorF.getMin()[1], voxelSize,
                  (int)((indicatorF.getMax()[0] - indicatorF.getMin()[0]) / voxelSize + 1.5),
                  (int)((indicatorF.getMax()[1] - indicatorF.getMin()[1]) / voxelSize + 1.5)),
    _periodicityOn(2, bool(false)), clout(std::cout, "CuboidGeometry2D")
{

  add(_motherCuboid);
  split(0, nC);
  shrink(indicatorF);
}


template<typename T>
void CuboidGeometry2D<T>::reInit(T globPosX, T globPosY, T delta, int nX, int nY, int nC)
{
  _cuboids.clear();
  _motherCuboid = Cuboid2D<T>(globPosX, globPosY, delta, nX, nY);
  Cuboid2D<T> cuboid(globPosX, globPosY, delta, nX, nY);
  if (_oldApproach) {
    cuboid.init(0, 0, 1, nX, nY);
  }

  add(cuboid);
  split(0, nC);
}

template<typename T>
Cuboid2D<T>& CuboidGeometry2D<T>::get(int i)
{
  return _cuboids[i];
}

template<typename T>
Cuboid2D<T> const& CuboidGeometry2D<T>::get(int i) const
{
  return _cuboids[i];
}

template<typename T>
void CuboidGeometry2D<T>::setPeriodicity(bool periodicityX, bool periodicityY)
{
  _periodicityOn.resize(2);
  _periodicityOn[0] = periodicityX;
  _periodicityOn[1] = periodicityY;
}


template<typename T>
bool CuboidGeometry2D<T>::getC(std::vector<T> physR, int& iC) const
{
  int iCtmp = get_iC(physR[0], physR[1]);
  if (iCtmp < getNc()) {
    iC = iCtmp;
    return true;
  }
  else {
    return false;
  }
}

template<typename T>
bool CuboidGeometry2D<T>::getC(const Vector<T,2>& physR, int& iC) const
{
  int iCtmp = get_iC(physR[0], physR[1]);
  if (iCtmp < getNc()) {
    iC = iCtmp;
    return true;
  }
  else {
    return false;
  }
}

template<typename T>
bool CuboidGeometry2D<T>::getLatticeR(std::vector<T> physR, std::vector<int>& latticeR) const
{
  return getLatticeR(&latticeR[0], &physR[0]);
  /*  int iCtmp = get_iC(physR[0], physR[1]);
    if (iCtmp < getNc()) {
      latticeR[0] = iCtmp;
      latticeR[1] = (int)floor( (physR[0] - _cuboids[latticeR[0]].getOrigin()[0] ) / _cuboids[latticeR[0]].getDeltaR() + .5);
      latticeR[2] = (int)floor( (physR[1] - _cuboids[latticeR[0]].getOrigin()[1] ) / _cuboids[latticeR[0]].getDeltaR() + .5);
      return true;
    } else {
      return false;
    }*/
}

template<typename T>
bool CuboidGeometry2D<T>::getLatticeR(int latticeR[], const T physR[]) const
{
  int iCtmp = get_iC(physR[0], physR[1]);
  if (iCtmp < getNc()) {
    latticeR[0] = iCtmp;
    latticeR[1] = (int)floor( (physR[0] - _cuboids[latticeR[0]].getOrigin()[0] ) / _cuboids[latticeR[0]].getDeltaR() + .5);
    latticeR[2] = (int)floor( (physR[1] - _cuboids[latticeR[0]].getOrigin()[1] ) / _cuboids[latticeR[0]].getDeltaR() + .5);
    return true;
  }
  else {
    return false;
  }
}

template<typename T>
bool CuboidGeometry2D<T>::getLatticeR(
  const Vector<T,2>& physR, Vector<int,3>& latticeR) const
{
  int iCtmp = get_iC(physR[0], physR[1]);
  if (iCtmp < getNc()) {
    latticeR[0] = iCtmp;
    latticeR[1] = (int)floor( (physR[0] - _cuboids[latticeR[0]].getOrigin()[0] ) / _cuboids[latticeR[0]].getDeltaR() + .5);
    latticeR[2] = (int)floor( (physR[1] - _cuboids[latticeR[0]].getOrigin()[1] ) / _cuboids[latticeR[0]].getDeltaR() + .5);
    return true;
  }
  else {
    return false;
  }
}

template<typename T>
bool CuboidGeometry2D<T>::getFloorLatticeR(std::vector<T> physR, std::vector<int>& latticeR) const
{
  int iCtmp = get_iC(physR[0], physR[1]);
  if (iCtmp < getNc()) {
    latticeR[0] = iCtmp;
    latticeR[1] = (int)floor( (physR[0] - _cuboids[latticeR[0]].getOrigin()[0] ) / _cuboids[latticeR[0]].getDeltaR() );
    latticeR[2] = (int)floor( (physR[1] - _cuboids[latticeR[0]].getOrigin()[1] ) / _cuboids[latticeR[0]].getDeltaR() );
    return true;
  }
  else {
    return false;
  }
}

template<typename T>
bool CuboidGeometry2D<T>::getFloorLatticeR(
  const Vector<T,2>& physR, Vector<int,3>& latticeR) const
{
  int iCtmp = get_iC(physR[0], physR[1]);
  if (iCtmp < getNc()) {
    latticeR[0] = iCtmp;
    latticeR[1] = (int)floor( (physR[0] - _cuboids[latticeR[0]].getOrigin()[0] ) / _cuboids[latticeR[0]].getDeltaR() );
    latticeR[2] = (int)floor( (physR[1] - _cuboids[latticeR[0]].getOrigin()[1] ) / _cuboids[latticeR[0]].getDeltaR() );
    return true;
  }
  else {
    return false;
  }
}

template<typename T>
std::vector<T> CuboidGeometry2D<T>::getPhysR(int iCglob, int iX, int iY) const
{
  std::vector<T> physR(2,T());
  _cuboids[iCglob].getPhysR(&(physR[0]), iX, iY);
  for (int iDim = 0; iDim < 2; iDim++) {
    if (_periodicityOn[iDim]) {
      //std::cout << iDim << _periodicityOn[iDim] <<":"<< _motherCuboid.getDeltaR()*(_motherCuboid.getExtend()[iDim]) << std::endl;
      physR[iDim] = remainder( physR[iDim] - _motherCuboid.getOrigin()[iDim]
                               + _motherCuboid.getDeltaR() * (_motherCuboid.getExtend()[iDim]),
                               _motherCuboid.getDeltaR() * (_motherCuboid.getExtend()[iDim]));
      // solving the rounding error problem for double
      if ( physR[iDim]*physR[iDim] < 0.001 * _motherCuboid.getDeltaR()*_motherCuboid.getDeltaR() ) {
        physR[iDim] = T();
      }
      // make it to mod instead remainer
      if ( physR[iDim] < 0 ) {
        physR[iDim] += _motherCuboid.getDeltaR() * _motherCuboid.getExtend()[iDim];
      }
      // add origin
      physR[iDim] += _motherCuboid.getOrigin()[iDim];
    }
  }
  return physR;
}

template<typename T>
std::vector<T> CuboidGeometry2D<T>::getPhysR(std::vector<int> latticeR) const
{
  return getPhysR(latticeR[0], latticeR[1], latticeR[2]);
}

template<typename T>
void CuboidGeometry2D<T>::getPhysR(T output[2], const int latticeR[3]) const
{
  getPhysR(output, latticeR[0], latticeR[1], latticeR[2]);
}

template<typename T>
void CuboidGeometry2D<T>::getPhysR(T output[2], const int iCglob, const int iX, const int iY) const
{
  _cuboids[iCglob].getPhysR(output, iX, iY);
  for (int iDim = 0; iDim < 2; iDim++) {
    if (_periodicityOn[iDim]) {
      //std::cout << iDim << _periodicityOn[iDim] <<":"<< _motherCuboid.getDeltaR()*(_motherCuboid.getExtend()[iDim]) << std::endl;
      output[iDim] = remainder( output[iDim] - _motherCuboid.getOrigin()[iDim]
                                + _motherCuboid.getDeltaR() * (_motherCuboid.getExtend()[iDim]),
                                _motherCuboid.getDeltaR() * (_motherCuboid.getExtend()[iDim]));
      // solving the rounding error problem for double
      if ( output[iDim]*output[iDim] < 0.001 * _motherCuboid.getDeltaR()*_motherCuboid.getDeltaR() ) {
        if ( output[iDim] > 0 ) {
          output[iDim] = _motherCuboid.getDeltaR() * _motherCuboid.getExtend()[iDim];
        }
        else {
          output[iDim] = T();
        }
      }
      // make it to mod instead remainer
      if ( output[iDim] < 0 ) {
        output[iDim] += _motherCuboid.getDeltaR() * _motherCuboid.getExtend()[iDim];
      }
      // add origin
      output[iDim] += _motherCuboid.getOrigin()[iDim];
    }
  }
}

template<typename T>
int CuboidGeometry2D<T>::getNc() const
{
  return _cuboids.size();
}

template<typename T>
T CuboidGeometry2D<T>::getMinRatio() const
{
  T minRatio = 1.;
  for ( auto& cuboid : _cuboids ) {
    if ((T)cuboid.getNx() / (T)cuboid.getNy() < minRatio) {
      minRatio = (T)cuboid.getNx() / (T)cuboid.getNy();
    }
  }
  return minRatio;
}

template<typename T>
T CuboidGeometry2D<T>::getMaxRatio() const
{
  T maxRatio = 1.;
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if ((T)_cuboids[i].getNx() / (T)_cuboids[i].getNy() > maxRatio) {
      maxRatio = (T)_cuboids[i].getNx() / (T)_cuboids[i].getNy();
    }
  }
  return maxRatio;
}

template<typename T>
std::vector<T> CuboidGeometry2D<T>::getMinPhysR() const
{
  Vector<T,2> output(_cuboids[0].getOrigin());
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getOrigin()[0] < output[0]) {
      output[0] = _cuboids[i].getOrigin()[0];
    }
    if (_cuboids[i].getOrigin()[1] < output[1]) {
      output[1] = _cuboids[i].getOrigin()[1];
    }
  }
  return std::vector<T> {output[0],output[1]};
}

template<typename T>
std::vector<T> CuboidGeometry2D<T>::getMaxPhysR() const
{
  Vector<T,2> output(_cuboids[0].getOrigin());
  output[0] += _cuboids[0].getNx()*_cuboids[0].getDeltaR();
  output[1] += _cuboids[0].getNy()*_cuboids[0].getDeltaR();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getOrigin()[0] + _cuboids[i].getNx()*_cuboids[i].getDeltaR() > output[0]) {
      output[0] = _cuboids[i].getOrigin()[0] + _cuboids[i].getNx()*_cuboids[i].getDeltaR();
    }
    if (_cuboids[i].getOrigin()[1] + _cuboids[i].getNy()*_cuboids[i].getDeltaR() > output[1]) {
      output[1] = _cuboids[i].getOrigin()[1] + _cuboids[i].getNy()*_cuboids[i].getDeltaR();
    }
  }
  return std::vector<T> {output[0],output[1]};
}

template<typename T>
T CuboidGeometry2D<T>::getMinPhysVolume() const
{
  T minVolume = _cuboids[0].getPhysVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getPhysVolume() < minVolume) {
      minVolume = _cuboids[i].getPhysVolume();
    }
  }
  return minVolume;
}

template<typename T>
T CuboidGeometry2D<T>::getMaxPhysVolume() const
{
  T maxVolume = _cuboids[0].getPhysVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getPhysVolume() > maxVolume) {
      maxVolume = _cuboids[i].getPhysVolume();
    }
  }
  return maxVolume;
}

template<typename T>
size_t CuboidGeometry2D<T>::getMinLatticeVolume() const
{
  size_t minNodes = _cuboids[0].getLatticeVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getLatticeVolume() < minNodes) {
      minNodes = _cuboids[i].getLatticeVolume();
    }
  }
  return minNodes;
}

template<typename T>
size_t CuboidGeometry2D<T>::getMaxLatticeVolume() const
{
  size_t maxNodes = _cuboids[0].getLatticeVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getLatticeVolume() > maxNodes) {
      maxNodes = _cuboids[i].getLatticeVolume();
    }
  }
  return maxNodes;
}

template<typename T>
T CuboidGeometry2D<T>::getMinDeltaR() const
{
  T minDelta = _cuboids[0].getDeltaR();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getDeltaR() < minDelta) {
      minDelta = _cuboids[i].getDeltaR();
    }
  }
  return minDelta;
}

template<typename T>
T CuboidGeometry2D<T>::getMaxDeltaR() const
{
  T maxDelta = _cuboids[0].getDeltaR();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getDeltaR() > maxDelta) {
      maxDelta = _cuboids[i].getDeltaR();
    }
  }
  return maxDelta;
}

template<typename T>
Cuboid2D<T> CuboidGeometry2D<T>::getMotherCuboid() const
{

  /*Cuboid2D<T> found;
  if(_cuboids.size()==0) {
      found.init(0, 0, 0, 0, 0);
      return found;
  }

  T delta = _cuboids[0].getDeltaR();
  T globPosXmin = _cuboids[0].get_globPosX();
  T globPosYmin = _cuboids[0].get_globPosY();
  T globPosXmax = _cuboids[0].get_globPosX() + delta*(_cuboids[0].getNx()-1);
  T globPosYmax = _cuboids[0].get_globPosY() + delta*(_cuboids[0].getNy()-1);

  for (unsigned i=1; i<_cuboids.size(); i++) {
      if(delta > _cuboids[i].getDeltaR() ) {
          delta = _cuboids[i].getDeltaR();
      }
      if(globPosXmin > _cuboids[i].get_globPosX() ) {
          globPosXmin = _cuboids[i].get_globPosX();
      }
      if(globPosYmin > _cuboids[i].get_globPosY() ) {
          globPosYmin = _cuboids[i].get_globPosY();
      }
      if(globPosXmax < _cuboids[i].get_globPosX()
                          + delta*(_cuboids[i].getNx()-1)) {
          globPosXmax = _cuboids[i].get_globPosX()
                          + delta*(_cuboids[i].getNx()-1);
      }
      if(globPosYmax < _cuboids[i].get_globPosY()
                          + delta*(_cuboids[i].getNy()-1)) {
          globPosYmax = _cuboids[i].get_globPosY()
                          + delta*(_cuboids[i].getNy()-1);
      }
  }
  int nX = int(ceil((globPosXmax - globPosXmin)/delta))+1;
  int nY = int(ceil((globPosYmax - globPosYmin)/delta))+1;

  found.init(globPosXmin, globPosYmin, delta, nX, nY);

  return found;*/
  return _motherCuboid;
}

template<typename T>
int CuboidGeometry2D<T>::get_iC(T globX, T globY, int offset) const
{
  unsigned i;
  for (i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].checkPoint(globX, globY, offset)) {
      return (int)i;
    }
  }
  return (int)i;
}

template<typename T>
int CuboidGeometry2D<T>::get_iC(T globX, T globY, int orientationX, int orientationY) const
{
  unsigned i;
  for (i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].checkPoint(globX, globY) &&
        _cuboids[i].checkPoint(globX + orientationX / _cuboids[i].getDeltaR(),
                               globY + orientationY / _cuboids[i].getDeltaR())) {
      return (int)i;
    }
  }
  return (int)i;
}


template<typename T>
void CuboidGeometry2D<T>::add(Cuboid2D<T> cuboid)
{

  _cuboids.push_back(cuboid);
}

template<typename T>
void CuboidGeometry2D<T>::remove(int iC)
{

  _cuboids.erase(_cuboids.begin() + iC);
}
/*
template<typename T>
void CuboidGeometry2D<T>::remove(olb::ScalarField2D<int>* geometryData)
{

  std::vector<Cuboid2D<T> > cuboids;
  unsigned size = _cuboids.size();

  std::vector<bool> allZero;
  for (unsigned i = 0; i < size; i++) {
    allZero.push_back(1);
    for (int iX = 0; iX < _cuboids[i].getNx(); iX++) {
      for (int iY = 0; iY < _cuboids[i].getNy(); iY++) {
        if (geometryData->get(_cuboids[i].get_globPosX() + iX,
                              _cuboids[i].get_globPosY() + iY) != 0 ) {
          allZero[i] = 0;
        }
      }
    }
  }
  for (unsigned i = 0; i < size; i++) {
    if (!allZero[i] ) {
      cuboids.push_back(_cuboids[i]);
    }
  }
  _cuboids.clear();
  for (unsigned i = 0; i < cuboids.size(); i++) {
    _cuboids.push_back(cuboids[i]);
  }
}*/

template<typename T>
void CuboidGeometry2D<T>::shrink(IndicatorF2D<T>& indicatorF)
{
  //IndicatorIdentity3D<T> tmpIndicatorF(indicatorF);
  int newX, newY, maxX, maxY;
  int nC = getNc();
  std::vector<int> latticeR(3, 0);
  std::vector<T> physR(2, T());
  bool inside[1];
  for (int iC = nC - 1; iC >= 0; iC--) {
    latticeR[0] = iC;
    int fullCells = 0;
    int xN = get(iC).getNx();
    int yN = get(iC).getNy();
    maxX = 0;
    maxY = 0;
    newX = xN - 1;
    newY = yN - 1;
    for (int iX = 0; iX < xN; iX++) {
      for (int iY = 0; iY < yN; iY++) {
        latticeR[1] = iX;
        latticeR[2] = iY;
        physR = getPhysR(latticeR);
        indicatorF(inside,&physR[0]);
        if (inside[0]) {
          fullCells++;
          maxX = std::max(maxX, iX);
          maxY = std::max(maxY, iY);
          newX = std::min(newX, iX);
          newY = std::min(newY, iY);
        }
      }
    }
    if (fullCells > 0) {
      get(iC).setWeight(fullCells);
      _cuboids[iC].resize(newX, newY, maxX - newX + 1, maxY - newY + 1);
    }
    else {
      remove(iC);
    }
  }
  // shrink mother cuboid
  std::vector<T> minPhysR = getMinPhysR();
  std::vector<T> maxPhysR = getMaxPhysR();
  T minDelataR = getMinDeltaR();
  _motherCuboid = Cuboid2D<T>(minPhysR[0], minPhysR[1], minDelataR, (int)((maxPhysR[0]-minPhysR[0])/minDelataR + 0.5), (int)((maxPhysR[1]-minPhysR[1])/minDelataR + 0.5));
}

template<typename T>
void CuboidGeometry2D<T>::split(int iC, int p)
{

  Cuboid2D<T> temp(_cuboids[iC].get_globPosX(), _cuboids[iC].get_globPosY(),
                   _cuboids[iC].getDeltaR(), _cuboids[iC].getNx(), _cuboids[iC].getNy());
  temp.divide(p, _cuboids);
  remove(iC);
}

template<typename T>
void CuboidGeometry2D<T>::getNeighbourhood(int cuboid, std::vector<int> neighbours, int offset)
{
  for (int iC = 0; iC < getNc(); iC++) {
    if (cuboid == iC) {
      continue;
    }
    T globX = get(iC).get_globPosX();
    T globY = get(iC).get_globPosY();
    T nX = get(iC).getNx();
    T nY = get(iC).getNy();
    if (get(cuboid).checkInters(globX, globX + nX, globY, globY + nY, offset)) {
      neighbours.push_back(iC);
    }
  }
}

template<typename T>
void CuboidGeometry2D<T>::refineArea(T x0, T x1, T y0, T y1, int coarse_level)
{

  for (int iC = 0; iC < getNc(); iC++) {
    if (get(iC).get_refinementLevel() != coarse_level) {
      continue;
    }
    int locX0, locX1, locY0, locY1;
    bool inter = get(iC).checkInters(x0, y1, y0, y1, locX0, locX1, locY0, locY1, 0);
    if (!inter) {
      continue;
    }

    T globX = get(iC).get_globPosX();
    T globY = get(iC).get_globPosY();
    T delta = get(iC).getDeltaR();
    int nx = get(iC).getNx();
    int ny = get(iC).getNy();

    if (locX0 != 0) {
      Cuboid2D<T> right_side(globX, globY, delta, locX0, ny, coarse_level);
      add(right_side);
    }

    if (locY0 != 0) {
      Cuboid2D<T> down_side(globX + locX0 * delta, globY, delta, nx - locX0, locY0, coarse_level);
      add(down_side);
    }

    if (locX1 != get(iC).getNx() - 1) {
      Cuboid2D<T> left_side(globX + (locX1 + 1)*delta, globY + locY0 * delta, delta, nx - locX1 - 1, ny - locY0, coarse_level);
      add(left_side);
    }

    if (locY1 != get(iC).getNy() - 1) {
      Cuboid2D<T> top_side(globX + locX0 * delta, globY + (locY1 + 1)*delta, delta, locX1 - locX0 + 1, ny - locY1 - 1, coarse_level);
      add(top_side);
    }
    get(iC).init(globX + locX0 * delta, globY + locY0 * delta, delta, locX1 - locX0 + 1, locY1 - locY0 + 1, coarse_level);
    get(iC).refineIncrease();
  }
}

template<typename T>
void CuboidGeometry2D<T>::print() const
{
  clout << "---Cuboid Stucture Statistics---" << std::endl;
  clout << " Number of Cuboids: " << "\t" << getNc() << std::endl;
  clout << " Delta (min): " << "\t" << "\t" << getMinDeltaR() << std::endl;
  clout << "       (max): " << "\t" << "\t" << getMaxDeltaR() << std::endl;
  clout << " Ratio (min): " << "\t" << "\t" << getMinRatio() << std::endl;
  clout << "       (max): " << "\t" << "\t" << getMaxRatio() << std::endl;
  clout << " Nodes (min): " << "\t" << "\t" << getMinLatticeVolume() << std::endl;
  clout << "       (max): " << "\t" << "\t" << getMaxLatticeVolume() << std::endl;
  clout << "--------------------------------" << std::endl;
}

template<typename T>
void CuboidGeometry2D<T>::printExtended()
{
  clout << "Mothercuboid :" << std::endl;
  getMotherCuboid().print();

  for (int iC = 0; iC < getNc(); iC++) {
    clout << "Cuboid #" << iC << ": " << std::endl;
    get(iC).print();
  }
}

}  // namespace olb

#endif
