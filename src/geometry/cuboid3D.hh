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
 * The description of a single 3D cuboid -- generic implementation.
 */


#ifndef CUBOID_3D_HH
#define CUBOID_3D_HH


#include <vector>
#include "geometry/cuboidGeometry2D.h"
#include "cuboid3D.h"
#include "dynamics/lbHelpers.h"
#include "utilities/vectorHelpers.h"


namespace olb {

////////////////////// Class Cuboid3D /////////////////////////

/*template<typename T>
Cuboid3D<T>::Cuboid3D() : clout(std::cout,"Cuboid3D") {
  _weight = -1;
}*/

template<typename T>
Cuboid3D<T>::Cuboid3D() : Cuboid3D<T>(0, 0, 0, 0, 0, 0, 0, 0)
{
}

template<typename T>
Cuboid3D<T>::Cuboid3D(T globPosX, T globPosY, T globPosZ, T delta , int nX, int nY,
                      int nZ, int refinementLevel)
  : _weight(std::numeric_limits<size_t>::max()), clout(std::cout,"Cuboid3D")
{
  init(globPosX, globPosY, globPosZ, delta, nX, nY, nZ, refinementLevel);
}

template<typename T>
Cuboid3D<T>::Cuboid3D(std::vector<T> origin, T delta , std::vector<int> extend,
                      int refinementLevel)
  : clout(std::cout,"Cuboid3D")
{
  _weight = std::numeric_limits<size_t>::max();
  init(origin[0], origin[1], origin[2], delta, extend[0], extend[1], extend[2], refinementLevel);
}

template<typename T>
Cuboid3D<T>::Cuboid3D(IndicatorF3D<T>& indicatorF, T voxelSize, int refinementLevel)
  : clout(std::cout,"Cuboid3D")
{
  _weight = std::numeric_limits<size_t>::max();
  init(indicatorF.getMin()[0],  indicatorF.getMin()[1],  indicatorF.getMin()[2], voxelSize,
       (int)((indicatorF.getMax()[0]-indicatorF.getMin()[0])/voxelSize+1.5),
       (int)((indicatorF.getMax()[1]-indicatorF.getMin()[1])/voxelSize+1.5),
       (int)((indicatorF.getMax()[2]-indicatorF.getMin()[2])/voxelSize+1.5), refinementLevel);
}

// copy constructor
template<typename T>
Cuboid3D<T>::Cuboid3D(Cuboid3D<T> const& rhs, int overlap) : clout(std::cout,"Cuboid3D")
{
  init( rhs._globPosX-rhs._delta*overlap, rhs._globPosY-rhs._delta*overlap,
        rhs._globPosZ-rhs._delta*overlap, rhs._delta, rhs._nX+2*overlap,
        rhs._nY+2*overlap, rhs._nZ+2*overlap);
  _weight = rhs._weight;
  _refinementLevel = rhs._refinementLevel;
}


// copy assignment
template<typename T>
Cuboid3D<T>& Cuboid3D<T>::operator=(Cuboid3D<T> const& rhs)
{
  init( rhs._globPosX, rhs._globPosY, rhs._globPosZ, rhs._delta, rhs._nX,
        rhs._nY, rhs._nZ);
  _weight = rhs._weight;
  _refinementLevel = rhs._refinementLevel;
  return *this;
}

template<typename T>
void Cuboid3D<T>::init(T globPosX, T globPosY, T globPosZ, T delta, int nX, int nY,
                       int nZ, int refinementLevel)
{
  _globPosX = globPosX;
  _globPosY = globPosY;
  _globPosZ = globPosZ;
  _delta    = delta;
  _nX       = nX;
  _nY       = nY;
  _nZ       = nZ;
  _refinementLevel = refinementLevel;
}

template<typename T>
std::size_t Cuboid3D<T>::getNblock() const
{
  return 9;
}

template<typename T>
std::size_t Cuboid3D<T>::getSerializableSize() const
{
  return ( 4 * sizeof(T) )  + ( 5 * sizeof(int) );
}

template<typename T>
Vector<T,3> const Cuboid3D<T>::getOrigin() const
{
  return Vector<T,3> (_globPosX, _globPosY, _globPosZ);
}

template<typename T>
T Cuboid3D<T>::getDeltaR() const
{
  return _delta;
}

template<typename T>
int Cuboid3D<T>::getNx() const
{
  return _nX;
}

template<typename T>
int Cuboid3D<T>::getNy() const
{
  return _nY;
}

template<typename T>
int Cuboid3D<T>::getNz() const
{
  return _nZ;
}

template<typename T>
Vector<int,3> const Cuboid3D<T>::getExtend() const
{
  return Vector<int,3> (_nX, _nY, _nZ);
}

template<typename T>
T Cuboid3D<T>::getPhysVolume() const
{
  return _nY*_nX*_nZ*_delta*_delta*_delta;
}

template<typename T>
int Cuboid3D<T>::getWeightValue() const
{
  return _weight;
}

template<typename T>
size_t Cuboid3D<T>::getWeight() const
{
  if (_weight == std::numeric_limits<size_t>::max()) {
    return getLatticeVolume();
  } else {
    return _weight;
  }
}

template<typename T>
void Cuboid3D<T>::setWeight(size_t fullCells)
{
  _weight = fullCells;
}

template<typename T>
int Cuboid3D<T>::getRefinementLevel() const
{
  return _refinementLevel;
}

template<typename T>
void Cuboid3D<T>::setRefinementLevel(int refLevel)
{
  _refinementLevel = refLevel;
}

template<typename T>
size_t Cuboid3D<T>::getLatticeVolume() const
{
  return static_cast<size_t>(_nY)*static_cast<size_t>(_nX)*static_cast<size_t>(_nZ);
}

template<typename T>
T Cuboid3D<T>::getPhysPerimeter() const
{
  return 2*_delta*_delta*(_nX*_nY + _nY*_nZ + _nZ*_nX);
}

template<typename T>
int Cuboid3D<T>::getLatticePerimeter() const
{
  return 2*((_nX-1)*(_nY-1) + (_nY-1)*(_nZ-1) + (_nZ-1)*(_nX-1));
}


template<typename T>
bool Cuboid3D<T>::operator==(const Cuboid3D<T>& rhs) const
{
  return util::nearZero<T>(_globPosX - rhs._globPosX) &&
         util::nearZero<T>(_globPosY - rhs._globPosY) &&
         util::nearZero<T>(_globPosZ - rhs._globPosZ) &&
         util::nearZero<T>(_delta    - rhs._delta) &&
         _nX              == rhs._nX &&
         _nY              == rhs._nY &&
         _nZ              == rhs._nZ &&
         _weight          == rhs._weight &&
         _refinementLevel == rhs._refinementLevel;
}


template<typename T>
bool* Cuboid3D<T>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _globPosX);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _globPosY);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _globPosZ);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _delta);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _nX);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _nY);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _nZ);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _weight);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _refinementLevel);

  return dataPtr;
}

template<typename T>
void Cuboid3D<T>::print() const
{
  clout << "--------Cuboid Details----------" << std::endl;
  clout << " Corner (x/y/z): " << "\t" << "("
        << _globPosX-_delta/2. << "/" << _globPosY - _delta/2.
        << "/" << _globPosZ - _delta/2. << ")" << std::endl;
  clout << " Delta: " << "\t" << "\t" << getDeltaR() << std::endl;
  clout << " Perimeter: " << "\t" << "\t" << getPhysPerimeter() << std::endl;
  clout << " Volume: " << "\t" << "\t" << getPhysVolume() << std::endl;
  clout << " Extent (x/y/z): " << "\t" << "("
        << getNx() << "/" << getNy() << "/"
        << getNz() << ")" << std::endl;
  clout << " Nodes at Perimeter: " << "\t" << getLatticePerimeter() << std::endl;
  clout << " Nodes in Volume: " << "\t" << getLatticeVolume() << std::endl;
  clout << " Nodes in Indicator: " << "\t" << getWeight() << std::endl;
  clout << " Other Corner: "  << "\t"<<  "(" << _globPosX + T(_nX-0.5)*_delta << "/"
        << _globPosY + T(_nY-0.5)*_delta << "/"
        << _globPosZ + T(_nZ-0.5)*_delta <<  ")" << std::endl;
  clout << "--------------------------------" << std::endl;

}

template<typename T>
void Cuboid3D<T>::getPhysR(T physR[3], const int latticeR[3]) const
{
  physR[0] = _globPosX + latticeR[0]*_delta;
  physR[1] = _globPosY + latticeR[1]*_delta;
  physR[2] = _globPosZ + latticeR[2]*_delta;
}

template<typename T>
void Cuboid3D<T>::getPhysR(T physR[3], const int& iX, const int& iY, const int& iZ) const
{
  physR[0] = _globPosX + iX*_delta;
  physR[1] = _globPosY + iY*_delta;
  physR[2] = _globPosZ + iZ*_delta;
}

template<typename T>
void Cuboid3D<T>::getLatticeR(int latticeR[3], const T physR[3]) const
{
  latticeR[0] = (int)floor( (physR[0] - _globPosX )/_delta +.5);
  latticeR[1] = (int)floor( (physR[1] - _globPosY )/_delta +.5);
  latticeR[2] = (int)floor( (physR[2] - _globPosZ )/_delta +.5);
}

template<typename T>
void Cuboid3D<T>::getLatticeR(int latticeR[3], const Vector<T,3>& physR) const
{
  latticeR[0] = (int)floor( (physR[0] - _globPosX )/_delta +.5);
  latticeR[1] = (int)floor( (physR[1] - _globPosY )/_delta +.5);
  latticeR[2] = (int)floor( (physR[2] - _globPosZ )/_delta +.5);
}

template<typename T>
void Cuboid3D<T>::getFloorLatticeR(const std::vector<T>& physR, std::vector<int>& latticeR) const
{
  getFloorLatticeR(&latticeR[0], &physR[0]);
}

template<typename T>
void Cuboid3D<T>::getFloorLatticeR(int latticeR[3], const T physR[3]) const
{
  latticeR[0] = (int)floor( (physR[0] - _globPosX)/_delta);
  latticeR[1] = (int)floor( (physR[1] - _globPosY)/_delta);
  latticeR[2] = (int)floor( (physR[2] - _globPosZ)/_delta);
}

template<typename T>
bool Cuboid3D<T>::checkPoint(T globX, T globY, T globZ, int overlap) const
{
  return _globPosX - _delta / 2. <= globX + overlap * _delta &&
         _globPosX + T(_nX-0.5+overlap) * _delta  > globX &&
         _globPosY - _delta/2. <= globY + overlap * _delta &&
         _globPosY + T(_nY-0.5+overlap) * _delta  > globY &&
         _globPosZ - _delta/2. <= globZ + overlap * _delta &&
         _globPosZ + T(_nZ-0.5+overlap) * _delta  > globZ;
}

template<typename T>
bool Cuboid3D<T>::physCheckPoint(T globX, T globY, T globZ, double overlap) const
{
  return _globPosX - T(0.5 + overlap) * _delta <= globX &&
         _globPosX + T(_nX-0.5+overlap)*_delta  > globX &&
         _globPosY - T(0.5 + overlap)*_delta <= globY &&
         _globPosY + T(_nY-0.5+overlap)*_delta  > globY &&
         _globPosZ - T(0.5 + overlap)*_delta <= globZ &&
         _globPosZ + T(_nZ-0.5+overlap)*_delta  > globZ;
}

template<typename T>
bool Cuboid3D<T>::checkPoint(T globX, T globY, T globZ, int &locX, int &locY,
                             int &locZ, int overlap) const
{
  if (overlap!=0) {
    Cuboid3D tmp(_globPosX - overlap*_delta, _globPosY - overlap*_delta, _globPosZ - overlap*_delta,
                 _delta , _nX + overlap*2, _nY + overlap*2, _nZ + overlap*2);
    return tmp.checkPoint(globX, globY, globZ, locX, locY, locZ);
  } else if (!checkPoint(globX, globY, globZ)) {
    return false;
  } else {
    locX = (int)floor((globX - (T)_globPosX)/_delta + .5);
    locY = (int)floor((globY - (T)_globPosY)/_delta + .5);
    locZ = (int)floor((globZ - (T)_globPosZ)/_delta + .5);
    return true;
  }
}

template<typename T>
bool Cuboid3D<T>::checkInters(T globX0, T globX1, T globY0, T globY1, T globZ0,
                              T globZ1, int overlap) const
{
  T locX0d = std::max(_globPosX-overlap*_delta,globX0);
  T locY0d = std::max(_globPosY-overlap*_delta,globY0);
  T locZ0d = std::max(_globPosZ-overlap*_delta,globZ0);
  T locX1d = std::min(_globPosX+(_nX+overlap-1)*_delta,globX1);
  T locY1d = std::min(_globPosY+(_nY+overlap-1)*_delta,globY1);
  T locZ1d = std::min(_globPosZ+(_nZ+overlap-1)*_delta,globZ1);

  return locX1d >= locX0d
      && locY1d >= locY0d
      && locZ1d >= locZ0d;
}

template<typename T>
bool Cuboid3D<T>::checkInters(T globX, T globY, T globZ, int overlap) const
{
  return checkInters(globX, globX, globY, globY, globZ, globZ, overlap);
}

template<typename T>
bool Cuboid3D<T>::checkInters(T globX0, T globX1, T globY0, T globY1, T globZ0, T globZ1,
                              int &locX0, int &locX1, int &locY0, int &locY1, int &locZ0, int &locZ1, int overlap) const
{
  if (overlap!=0) {
    Cuboid3D tmp(_globPosX - overlap*_delta, _globPosY - overlap*_delta, _globPosZ - overlap*_delta,
                 _delta , _nX + overlap*2, _nY + overlap*2, _nZ + overlap*2);
    return tmp.checkInters(globX0, globX1, globY0, globY1, globZ0, globZ1,
                           locX0, locX1, locY0, locY1, locZ0, locZ1);
  } else if (!checkInters(globX0, globX1, globY0, globY1, globZ0, globZ1)) {
    locX0 = 1;
    locX1 = 0;
    locY0 = 1;
    locY1 = 0;
    locZ0 = 1;
    locZ1 = 0;
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
    locZ0 = 0;
    for (int i=0; _globPosZ + i*_delta < globZ0; i++) {
      locZ0 = i+1;
    }
    locZ1 = _nZ-1;
    for (int i=_nZ-1; _globPosZ + i*_delta > globZ1; i--) {
      locZ1 = i-1;
    }
    return true;
  }
}

template<typename T>
void Cuboid3D<T>::divide(int nX, int nY, int nZ, std::vector<Cuboid3D<T> > &childrenC) const
{
  int xN_child = 0;
  int yN_child = 0;
  int zN_child = 0;

  T globPosX_child = _globPosX;
  T globPosY_child = _globPosY;
  T globPosZ_child = _globPosZ;

  for (int iX=0; iX<nX; iX++) {
    for (int iY=0; iY<nY; iY++) {
      for (int iZ=0; iZ<nZ; iZ++) {
        xN_child       = (_nX+nX-iX-1)/nX;
        yN_child       = (_nY+nY-iY-1)/nY;
        zN_child       = (_nZ+nZ-iZ-1)/nZ;
        Cuboid3D<T> child(globPosX_child, globPosY_child, globPosZ_child,
                          _delta, xN_child, yN_child, zN_child);
        childrenC.push_back(child);
        globPosZ_child += zN_child*_delta;
      }
      globPosZ_child = _globPosZ;
      globPosY_child += yN_child*_delta;
    }
    globPosY_child = _globPosY;
    globPosX_child += xN_child*_delta;
  }
}


template<typename T>
void Cuboid3D<T>::resize(int iX, int iY, int iZ, int nX, int nY, int nZ)
{
  _globPosX = _globPosX+iX*_delta;
  _globPosY = _globPosY+iY*_delta;
  _globPosZ = _globPosZ+iZ*_delta;

  _nX = nX;
  _nY = nY;
  _nZ = nZ;
}


template<typename T>
void Cuboid3D<T>::divide(int p, std::vector<Cuboid3D<T> > &childrenC) const
{

  OLB_PRECONDITION(p>0);

  int iXX = 1;
  int iYY = 1;
  int iZZ = p;
  int nX = _nX/iXX;
  int bestIx = iXX;
  int nY = _nY/iYY;
  int bestIy = iYY;
  int nZ = _nZ/iZZ;
  int bestIz = iZZ;
  T bestRatio = ((T)(_nX/iXX)/(T)(_nY/iYY)-1)*((T)(_nX/iXX)/(T)(_nY/iYY)-1)
                + ((T)(_nY/iYY)/(T)(_nZ/iZZ)-1)*((T)(_nY/iYY)/(T)(_nZ/iZZ)-1)
                + ((T)(_nZ/iZZ)/(T)(_nX/iXX)-1)*((T)(_nZ/iZZ)/(T)(_nX/iXX)-1);

  for (int iX=1; iX<=p; iX++) {
    for (int iY=1; iY*iX<=p; iY++) {
      for (int iZ=p/(iX*iY); iZ*iY*iX<=p; iZ++) {
        if ((iX+1)*iY*iZ>p && iX*(iY+1)*iZ>p ) {
          T ratio = ((T)(_nX/iX)/(T)(_nY/iY)-1)*((T)(_nX/iX)/(T)(_nY/iY)-1)
                    + ((T)(_nY/iY)/(T)(_nZ/iZ)-1)*((T)(_nY/iY)/(T)(_nZ/iZ)-1)
                    + ((T)(_nZ/iZ)/(T)(_nX/iX)-1)*((T)(_nZ/iZ)/(T)(_nX/iX)-1);
          if (ratio<bestRatio) {
            bestRatio = ratio;
            bestIx = iX;
            bestIy = iY;
            bestIz = iZ;
            nX = _nX/iX;
            nY = _nY/iY;
            nZ = _nZ/iZ;
          }
        }
      }
    }
  }

  int rest = p - bestIx*bestIy*bestIz;

  // split in one cuboid
  if (rest==0) {
    divide(bestIx, bestIy, bestIz, childrenC);
    return;
  } else {

    // add in z than in y direction
    if (nZ>nY && nZ>nX) {

      int restY = rest%bestIy;
      // split in two cuboid
      if (restY==0) {
        int restX = rest/bestIy;
        CuboidGeometry2D<T> helpG(_globPosX, _globPosZ, _delta, _nX, _nZ, bestIx*bestIz+restX);

        int yN_child = 0;
        T globPosY_child = _globPosY;

        for (int iY=0; iY<bestIy; iY++) {
          yN_child         = (_nY+bestIy-iY-1)/bestIy;
          for (int iC=0; iC<helpG.getNc(); iC++) {
            int xN_child     = helpG.get(iC).getNx();
            int zN_child     = helpG.get(iC).getNy();
            T globPosX_child = helpG.get(iC).get_globPosX();
            T globPosZ_child = helpG.get(iC).get_globPosY();

            Cuboid3D<T> child(globPosX_child, globPosY_child, globPosZ_child,
                              _delta, xN_child, yN_child, zN_child);
            childrenC.push_back(child);

          }
          globPosY_child += yN_child*_delta;
        }
        return;
      }

      // split in four cuboid

      int restX = rest/bestIy+1;
      int yN_child = 0;
      T globPosY_child = _globPosY;
      int splited_nY = (int) (_nY * (T)((bestIx*bestIz+restX)*restY)/(T)p);
      CuboidGeometry2D<T> helpG0(_globPosX, _globPosZ, _delta, _nX, _nZ, bestIx*bestIz+restX);

      for (int iY=0; iY<restY; iY++) {
        yN_child         = (splited_nY+restY-iY-1)/restY;
        for (int iC=0; iC<helpG0.getNc(); iC++) {
          int xN_child     = helpG0.get(iC).getNx();
          int zN_child     = helpG0.get(iC).getNy();
          T globPosX_child = helpG0.get(iC).get_globPosX();
          T globPosZ_child = helpG0.get(iC).get_globPosY();

          Cuboid3D<T> child(globPosX_child, globPosY_child, globPosZ_child,
                            _delta, xN_child, yN_child, zN_child);
          childrenC.push_back(child);
        }
        globPosY_child += yN_child*_delta;
      }

      splited_nY = _nY - splited_nY;
      restX = rest/bestIy;
      CuboidGeometry2D<T> helpG1(_globPosX, _globPosZ, _delta, _nX, _nZ, bestIx*bestIz+restX);
      yN_child = 0;

      for (int iY=0; iY<bestIy-restY; iY++) {
        yN_child         = (splited_nY+bestIy-restY-iY-1)/(bestIy-restY);
        for (int iC=0; iC<helpG1.getNc(); iC++) {
          int xN_child     = helpG1.get(iC).getNx();
          int zN_child     = helpG1.get(iC).getNy();
          T globPosX_child = helpG1.get(iC).get_globPosX();
          T globPosZ_child = helpG1.get(iC).get_globPosY();

          Cuboid3D<T> child(globPosX_child, globPosY_child, globPosZ_child,
                            _delta, xN_child, yN_child, zN_child);
          childrenC.push_back(child);
        }
        globPosY_child += yN_child*_delta;
      }
      return;
    }

    // add in x than in y direction
    else if (nX>nY && nX>nZ) {
      int restY = rest%bestIy;
      // split in two cuboid
      if (restY==0) {
        int restZ = rest/bestIy;
        CuboidGeometry2D<T> helpG(_globPosX, _globPosZ, _delta, _nX, _nZ, bestIx*bestIz+restZ);

        int yN_child = 0;
        T globPosY_child = _globPosY;

        for (int iY=0; iY<bestIy; iY++) {
          yN_child         = (_nY+bestIy-iY-1)/bestIy;
          for (int iC=0; iC<helpG.getNc(); iC++) {
            int xN_child     = helpG.get(iC).getNx();
            int zN_child     = helpG.get(iC).getNy();
            T globPosX_child = helpG.get(iC).get_globPosX();
            T globPosZ_child = helpG.get(iC).get_globPosY();

            Cuboid3D<T> child(globPosX_child, globPosY_child, globPosZ_child,
                              _delta, xN_child, yN_child, zN_child);
            childrenC.push_back(child);

          }
          globPosY_child += yN_child*_delta;
        }
        return;
      }

      // split in four cuboid

      int restZ = rest/bestIy+1;

      int yN_child = 0;
      T globPosY_child = _globPosY;
      int splited_nY = (int) (_nY * (T)((bestIx*bestIz+restZ)*restY)/(T)p);
      CuboidGeometry2D<T> helpG0(_globPosX, _globPosZ, _delta, _nX, _nZ, bestIx*bestIz+restZ);

      for (int iY=0; iY<restY; iY++) {
        yN_child         = (splited_nY+restY-iY-1)/restY;
        for (int iC=0; iC<helpG0.getNc(); iC++) {
          int xN_child     = helpG0.get(iC).getNx();
          int zN_child     = helpG0.get(iC).getNy();
          T globPosX_child = helpG0.get(iC).get_globPosX();
          T globPosZ_child = helpG0.get(iC).get_globPosY();

          Cuboid3D<T> child(globPosX_child, globPosY_child, globPosZ_child,
                            _delta, xN_child, yN_child, zN_child);
          childrenC.push_back(child);
        }
        globPosY_child += yN_child*_delta;
      }

      splited_nY = _nY - splited_nY;
      restZ = rest/bestIy;

      CuboidGeometry2D<T> helpG1(_globPosX, _globPosZ, _delta, _nX, _nZ, bestIx*bestIz+restZ);
      yN_child = 0;

      for (int iY=0; iY<bestIy-restY; iY++) {
        yN_child         = (splited_nY+bestIy-restY-iY-1)/(bestIy-restY);
        for (int iC=0; iC<helpG1.getNc(); iC++) {
          int xN_child     = helpG1.get(iC).getNx();
          int zN_child     = helpG1.get(iC).getNy();
          T globPosX_child = helpG1.get(iC).get_globPosX();
          T globPosZ_child = helpG1.get(iC).get_globPosY();

          Cuboid3D<T> child(globPosX_child, globPosY_child, globPosZ_child,
                            _delta, xN_child, yN_child, zN_child);
          childrenC.push_back(child);
        }
        globPosY_child += yN_child*_delta;
      }
      return;
    }

    // add in y than in x direction
    else {
      int restX = rest%bestIx;
      // split in two cuboid
      if (restX==0) {
        int restZ = rest/bestIx;
        CuboidGeometry2D<T> helpG(_globPosZ, _globPosY, _delta, _nZ, _nY, bestIz*bestIy+restZ);


        int xN_child = 0;
        T globPosX_child = _globPosX;

        for (int iX=0; iX<bestIx; iX++) {
          xN_child         = (_nX+bestIx-iX-1)/bestIx;
          for (int iC=0; iC<helpG.getNc(); iC++) {
            int zN_child     = helpG.get(iC).getNx();
            int yN_child     = helpG.get(iC).getNy();
            T globPosZ_child = helpG.get(iC).get_globPosX();
            T globPosY_child = helpG.get(iC).get_globPosY();

            Cuboid3D<T> child(globPosX_child, globPosY_child, globPosZ_child,
                              _delta, xN_child, yN_child, zN_child);
            childrenC.push_back(child);
          }
          globPosX_child += xN_child*_delta;
        }
        return;
      }

      // split in four cuboid

      int restZ = rest/bestIx+1;
      int xN_child = 0;
      T globPosX_child = _globPosX;
      int splited_nX = (int) (_nX * (T)((bestIz*bestIy+restZ)*restX)/(T)p);
      CuboidGeometry2D<T> helpG0(_globPosZ, _globPosY, _delta, _nZ, _nY, bestIz*bestIy+restZ);

      for (int iX=0; iX<restX; iX++) {
        xN_child         = (splited_nX+restX-iX-1)/restX;
        for (int iC=0; iC<helpG0.getNc(); iC++) {
          int zN_child     = helpG0.get(iC).getNx();
          int yN_child     = helpG0.get(iC).getNy();
          T globPosZ_child = helpG0.get(iC).get_globPosX();
          T globPosY_child = helpG0.get(iC).get_globPosY();

          Cuboid3D<T> child(globPosX_child, globPosY_child, globPosZ_child,
                            _delta, xN_child, yN_child, zN_child);
          childrenC.push_back(child);
        }
        globPosX_child += xN_child*_delta;
      }

      splited_nX = _nX - splited_nX;
      restZ = rest/bestIx;
      CuboidGeometry2D<T> helpG1(_globPosZ, _globPosY, _delta, _nZ, _nY, bestIz*bestIy+restZ);
      xN_child = 0;

      for (int iX=0; iX<bestIx-restX; iX++) {
        xN_child         = (splited_nX+bestIx-restX-iX-1)/(bestIx-restX);
        for (int iC=0; iC<helpG1.getNc(); iC++) {
          int zN_child     = helpG1.get(iC).getNx();
          int yN_child     = helpG1.get(iC).getNy();
          T globPosZ_child = helpG1.get(iC).get_globPosX();
          T globPosY_child = helpG1.get(iC).get_globPosY();

          Cuboid3D<T> child(globPosX_child, globPosY_child, globPosZ_child,
                            _delta, xN_child, yN_child, zN_child);
          childrenC.push_back(child);
        }
        globPosX_child += xN_child*_delta;
      }
      return;
    }
  }
}

}  // namespace olb

#endif
