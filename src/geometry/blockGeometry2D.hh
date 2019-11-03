/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2011, 2014 Mathias J. Krause, Simon Zimny
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
 * Representation of the 2D block geometry -- generic implementation.
 */

#ifndef BLOCK_GEOMETRY_2D_HH
#define BLOCK_GEOMETRY_2D_HH


#include "geometry/blockGeometry2D.h"

namespace olb {

template<typename T>
BlockGeometry2D<T>::BlockGeometry2D(T x0, T y0, T h, int nX, int nY, int iCglob)
  : BlockData2D<T,int>(nX,nY), BlockGeometryStructure2D<T>(iCglob), _cuboid(x0, y0, h, nX, nY)
{
  this->_statistics = BlockGeometryStatistics2D<T>(this);
  addToStatisticsList( &(this->_statistics.getStatisticsStatus()) );
}

template<typename T>
BlockGeometry2D<T>::BlockGeometry2D(Cuboid2D<T>& cuboid, int iCglob)
  : BlockData2D<T,int>(cuboid.getNx(),cuboid.getNy()), BlockGeometryStructure2D<T>(iCglob),
    _cuboid(cuboid)
{
  this->_statistics = BlockGeometryStatistics2D<T>(this);
  addToStatisticsList( &(this->_statistics.getStatisticsStatus()) );
}

template<typename T>
BlockGeometry2D<T>::BlockGeometry2D(BlockGeometry2D const& rhs)
  : BlockData2D<T,int>(rhs.getNx(), rhs.getNy()), BlockGeometryStructure2D<T>(rhs._iCglob),
    _cuboid(rhs._cuboid)
{
  this->_statistics = BlockGeometryStatistics2D<T>(this);
  addToStatisticsList( &(this->_statistics.getStatisticsStatus()) );
}

template<typename T>
BlockGeometry2D<T>& BlockGeometry2D<T>::operator=(BlockGeometry2D const& rhs)
{
  this->_nx = rhs._nx;
  this->_ny = rhs._ny;
  this->_iCglob = rhs._iCglob;
  _cuboid = rhs._cuboid;
  this->_statistics = BlockGeometryStatistics2D<T>(this);
  addToStatisticsList( &(this->_statistics.getStatisticsStatus()) );
  return *this;
}

template<typename T>
BlockStructure2D& BlockGeometry2D<T>::getBlockStructure()
{
  return *this;
}

template<typename T>
BlockGeometryStatistics2D<T>& BlockGeometry2D<T>::getStatistics(bool verbose)
{
  return this->_statistics;
}

template<typename T>
BlockGeometryStatistics2D<T> const& BlockGeometry2D<T>::getStatistics(bool verbose) const
{
  return this->_statistics;
}

template<typename T>
Vector<T,2> BlockGeometry2D<T>::getOrigin() const
{
  return _cuboid.getOrigin();
}

template<typename T>
const T BlockGeometry2D<T>::getDeltaR() const
{
  return _cuboid.getDeltaR();
}

template<typename T>
int BlockGeometry2D<T>::getNx() const
{
  return _cuboid.getNx();
}

template<typename T>
int BlockGeometry2D<T>::getNy() const
{
  return _cuboid.getNy();
}


template<typename T>
int& BlockGeometry2D<T>::get(int iX, int iY)
{
  resetStatistics();
  return BlockData2D<T,int>::get(iX,iY);
}

template<typename T>
int const& BlockGeometry2D<T>::get(int iX, int iY) const
{
  return BlockData2D<T,int>::get(iX,iY);
}

template<typename T>
int BlockGeometry2D<T>::getMaterial(int iX, int iY) const
{
  int material;
  if (iX < 0 || iX + 1 > getNx() || iY < 0 || iY + 1 > getNy() ) {
    material = 0;
  } else {
    material = BlockData2D<T,int>::get(iX,iY);
  }
  return material;
}

template<typename T>
void BlockGeometry2D<T>::getPhysR(T physR[2], const int& iX, const int& iY) const
{
  _cuboid.getPhysR(physR, iX, iY);
  return;
}

///TODO up, DONE down
/*
template<typename T>
olb::ScalarField3D<int>* BlockGeometry3D<T>::getRawData() {
  return &_geometryData;
}*/
/*
template<typename T>
void BlockGeometry3D<T>::resize(int X, int Y, int Z, int nX, int nY, int nZ) {
  olb::ScalarField3D<int> geometryDataTemp(nX, nY, nZ);
  geometryDataTemp.construct();
  geometryDataTemp.reset();
  for(int iX=0; iX<nX; iX++) {
    for(int iY=0; iY<nY; iY++) {
      for(int iZ=0; iZ<nZ; iZ++) {
        geometryDataTemp.get(iX, iY, iZ)
          = _geometryData.get(X + iX, Y + iY, Z + iZ);
      }
    }
  }
  _geometryData.swap(geometryDataTemp);
  geometryDataTemp.deConstruct();
}
*/
/*
template<typename T>
void BlockGeometry3D<T>::refineMesh(int level) {

  // new number of Voxels in X-,Y-,and Z-direction and new spacing

  int _nXnew = level * getNx();
  int _nYnew = level * getNy();
  int _nZnew = level * getNz();
  T _hnew = getDeltaR() / ((T) level);

  olb::ScalarField3D<int> _refinedGeometryData(_nXnew, _nYnew,
      _nZnew);
  _refinedGeometryData.construct();

  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        for (int li = 0; li < level; li++) {
          for (int lj = 0; lj < level; lj++) {
            for (int lk = 0; lk < level; lk++) {

              _refinedGeometryData.get(level * iX + li, level
                                       * iY + lj, level * iZ + lk) = getMaterial(
                                             iX, iY, iZ);

            }
          }
        }
      }
    }
  }

  for (int iX = 0; iX < _nXnew; iX++) {
    for (int iY = 0; iY < _nYnew; iY++) {
      for (int iZ = 0; iZ < _nZnew; iZ++) {

        if (iX == 0 || iY == 0 || iZ == 0 || iX == _nXnew - 1 || iY
            == _nYnew - 1 || iZ == _nZnew - 1)
          _refinedGeometryData.get(iX, iY, iZ) = 0;
      }
    }
  }
  reInit(getOrigin()[0], getOrigin()[1], getOrigin()[2], _hnew, _nXnew - 2, _nYnew - 2, _nZnew - 2, this->_iCglob,
         &_refinedGeometryData);
}*/


template<typename T>
void BlockGeometry2D<T>::addToStatisticsList(bool* statisticStatus)
{
  _statisticsUpdateNeeded.push_back(statisticStatus);
}

template<typename T>
void BlockGeometry2D<T>::removeFromStatisticsList(bool* statisticStatus)
{
  _statisticsUpdateNeeded.remove(statisticStatus);
}


template<typename T>
void BlockGeometry2D<T>::printLayer(int x0, int x1, int y0, int y1, bool linenumber)
{
  for (int x = x0; x <= x1; x++) {
    if (linenumber) {
      this->clout << x << ": ";
    }
    for (int y = y0; y <= y1; y++) {
      this->clout << getMaterial(x, y) << " ";
    }
    if (x1 - x0 != 0) {
      this->clout << std::endl;
    }
  }
  this->clout << std::endl;
}

template<typename T>
void BlockGeometry2D<T>::printLayer(int direction, int layer, bool linenumber)
{
  assert(direction >= 0 && direction <= 2);
  switch (direction) {
  case 0:
    printLayer(layer, layer, 0, getNy() - 1, linenumber);
    break;
  case 1:
    printLayer(0, getNx() - 1, layer, layer, linenumber);
    break;
  }
}

template<typename T>
void BlockGeometry2D<T>::printNode(int x0, int y0)
{
  for (int x = x0 - 1; x <= x0 + 1; x++) {
    this->clout << "x=" << x << std::endl;
    for (int y = y0 - 1; y <= y0 + 1; y++) {
      this->clout << getMaterial(x, y) << " ";
    }
    this->clout << std::endl;
  }
  this->clout << std::endl;
}

template<typename T>
void BlockGeometry2D<T>::resetStatistics()
{
  for (std::list<bool*>::iterator it = _statisticsUpdateNeeded.begin(); it != _statisticsUpdateNeeded.end(); it++) {
    **it = true;
  }
}

} // namespace olb

#endif
