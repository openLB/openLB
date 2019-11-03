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
 * Representation of a statistic for a 3D geometry -- generic implementation.
 */

#ifndef BLOCK_GEOMETRY_STATISTICS_3D_HH
#define BLOCK_GEOMETRY_STATISTICS_3D_HH

#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <cmath>

#include "geometry/blockGeometry3D.h"
#include "geometry/blockGeometryStatistics3D.h"

namespace olb {

template<typename T>
BlockGeometryStatistics3D<T>::BlockGeometryStatistics3D( BlockGeometryStructure3D<T>* blockGeometry)
  : _blockGeometry(blockGeometry), clout(std::cout,"BlockGeometryStatistics3D")
{
  _statisticsUpdateNeeded=true;
}


template<typename T>
bool& BlockGeometryStatistics3D<T>::getStatisticsStatus()
{
  return _statisticsUpdateNeeded;
}

template<typename T>
bool const & BlockGeometryStatistics3D<T>::getStatisticsStatus() const
{
  return _statisticsUpdateNeeded;
}


template<typename T>
void BlockGeometryStatistics3D<T>::update(bool verbose)
{

  if (getStatisticsStatus() ) {
    _material2n.clear();

    _nX = _blockGeometry->getNx();
    _nY = _blockGeometry->getNy();
    _nZ = _blockGeometry->getNz();
    _h = _blockGeometry->getDeltaR();

    for (int iX = 0; iX < _nX; ++iX) {
      for (int iY = 0; iY < _nY; ++iY) {
        for (int iZ = 0; iZ < _nZ; ++iZ) {
          takeStatistics(iX, iY, iZ);
        }
      }
    }

    _nMaterials=int();
    std::map<int, int>::iterator iter;
    for (iter = _material2n.begin(); iter != _material2n.end(); iter++) {
      _nMaterials++;
    }

    if (verbose) {
      clout << "updated" << std::endl;
    }
    getStatisticsStatus()=false;
  }
}


template<typename T>
int BlockGeometryStatistics3D<T>::getNmaterials()
{
  update();
  return _nMaterials;
}

template<typename T>
int BlockGeometryStatistics3D<T>::getNvoxel(int material)
{
  update();
  return _material2n[material];
}

template<typename T>
std::map<int, int> BlockGeometryStatistics3D<T>::getMaterial2n()
{
  update();
  return _material2n;
}

template<typename T>
int BlockGeometryStatistics3D<T>::getNvoxel()
{
  update();
  int total = 0;
  std::map<int, int>::iterator iter;
  for (iter = _material2n.begin(); iter != _material2n.end(); iter++) {
    if (iter->first!=0) {
      total+=iter->second;
    }
  }
  return total;
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::getMinLatticeR(int material)
{
  update();
  return _material2min[material];
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::getMaxLatticeR(int material)
{
  update();
  return _material2max[material];
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getMinPhysR(int material)
{
  std::vector<T> tmp(3,T());
  _blockGeometry->getPhysR(&(tmp[0]), &(getMinLatticeR(material)[0]));
  return tmp;
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getMaxPhysR(int material)
{
  std::vector<T> tmp(3,T());
  _blockGeometry->getPhysR(&(tmp[0]), &(getMaxLatticeR(material)[0]));
  return tmp;
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getLatticeExtend(int material)
{
  update();
  std::vector<T> extend;
  for (int iDim = 0; iDim < 3; iDim++) {
    extend.push_back(_material2max[material][iDim] - _material2min[material][iDim]);
  }
  return extend;
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getPhysExtend(int material)
{
  update();
  std::vector<T> extend;
  for (int iDim = 0; iDim < 3; iDim++) {
    extend.push_back(getMaxPhysR(material)[iDim] - getMinPhysR(material)[iDim]);
  }
  return extend;
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getPhysRadius(int material)
{
  update();
  std::vector<T> radius;
  for (int iDim=0; iDim<3; iDim++) {
    radius.push_back((getMaxPhysR(material)[iDim] - getMinPhysR(material)[iDim])/2.);
  }
  return radius;
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getCenterPhysR(int material)
{
  update();
  std::vector<T> center;
  for (int iDim=0; iDim<3; iDim++) {
    center.push_back(getMinPhysR(material)[iDim] + getPhysRadius(material)[iDim]);
  }
  return center;
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::getType(int iX, int iY, int iZ, bool anyNormal)
{

  update();
  std::vector<int> discreteNormal(4, 0);
  std::vector<int> discreteNormal2(4, 0);
  std::vector<int> nullVector(4, 0);

  if (_blockGeometry->getMaterial(iX, iY, iZ) != 1
      && _blockGeometry->getMaterial(iX, iY, iZ) != 0) {

    //boundary0N and boundary 0P
    if (_blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0) {

      if (_blockGeometry->getMaterial(iX + 1, iY, iZ) == 1) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 0;
          discreteNormal[1] = -1;
          discreteNormal[2] = 0;
          discreteNormal[3] = 0;
        } else {
          discreteNormal2[0] = 0;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = 0;
        }
      }

      if (_blockGeometry->getMaterial(iX - 1, iY, iZ) == 1) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 0;
          discreteNormal[1] = 1;
          discreteNormal[2] = 0;
          discreteNormal[3] = 0;
        } else {
          discreteNormal2[0] = 0;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = 0;
        }
      }
    }

    // boundary1N and boundary1P
    if (_blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0) {

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ) == 1) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 0;
          discreteNormal[1] = 0;
          discreteNormal[2] = -1;
          discreteNormal[3] = 0;
        } else {
          discreteNormal2[0] = 0;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 0;
        }
      }

      if (_blockGeometry->getMaterial(iX, iY - 1, iZ) == 1) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 0;
          discreteNormal[1] = 0;
          discreteNormal[2] = 1;
          discreteNormal[3] = 0;
        } else {
          discreteNormal2[0] = 0;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 0;
        }
      }
    }

    // boundary2N and boundary2P
    if (_blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0) {

      if (_blockGeometry->getMaterial(iX, iY, iZ + 1) == 1) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 0;
          discreteNormal[1] = 0;
          discreteNormal[2] = 0;
          discreteNormal[3] = -1;
        } else {
          discreteNormal2[0] = 0;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = -1;
        }
      }

      if (_blockGeometry->getMaterial(iX, iY, iZ - 1) == 1) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 0;
          discreteNormal[1] = 0;
          discreteNormal[2] = 0;
          discreteNormal[3] = 1;
        } else {
          discreteNormal2[0] = 0;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = 1;
        }
      }
    }

    // externalCornerNNN and externalCornerNPN
    if (_blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
        && _blockGeometry->getMaterial(iX + 1, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ + 1) != 0) {

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0
          && _blockGeometry->getMaterial(iX + 1, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX + 1, iY + 1, iZ) != 0
          && _blockGeometry->getMaterial(iX, iY + 1, iZ + 1) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ + 1) != 0
          && _blockGeometry->getMaterial(iX + 1, iY + 1, iZ + 1) == 1) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = -1;
          discreteNormal[2] = -1;
          discreteNormal[3] = -1;
        } else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = -1;
        }
      }

      if (_blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0
          && _blockGeometry->getMaterial(iX + 1, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX + 1, iY - 1, iZ) != 0
          && _blockGeometry->getMaterial(iX, iY - 1, iZ + 1) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ + 1) != 0
          && _blockGeometry->getMaterial(iX + 1, iY - 1, iZ + 1) == 1) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = -1;
          discreteNormal[2] = 1;
          discreteNormal[3] = -1;
        } else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = -1;
        }
      }
    }

    // externalCornerNPP and externalCornerNNP
    if (_blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
        && _blockGeometry->getMaterial(iX + 1, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ - 1) != 0) {

      if (_blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0
          && _blockGeometry->getMaterial(iX + 1, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX + 1, iY - 1, iZ) != 0
          && _blockGeometry->getMaterial(iX, iY - 1, iZ - 1) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ - 1) != 0
          && _blockGeometry->getMaterial(iX + 1, iY - 1, iZ - 1) == 1) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = -1;
          discreteNormal[2] = 1;
          discreteNormal[3] = 1;
        } else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 1;
        }
      }

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0
          && _blockGeometry->getMaterial(iX + 1, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX + 1, iY + 1, iZ) != 0
          && _blockGeometry->getMaterial(iX, iY + 1, iZ - 1) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ - 1) != 0
          && _blockGeometry->getMaterial(iX + 1, iY + 1, iZ - 1) == 1) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = -1;
          discreteNormal[2] = -1;
          discreteNormal[3] = 1;
        }

        else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 1;
        }
      }
    }

    // externalCornerPPP and externalCornerPNP
    if (_blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
        && _blockGeometry->getMaterial(iX - 1, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ - 1) != 0) {

      if (_blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0
          && _blockGeometry->getMaterial(iX, iY - 1, iZ - 1) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ - 1) != 0
          && _blockGeometry->getMaterial(iX - 1, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX - 1, iY - 1, iZ) != 0
          && _blockGeometry->getMaterial(iX - 1, iY - 1, iZ - 1) == 1) {

        if (discreteNormal == nullVector) {

          discreteNormal[0] = 1;
          discreteNormal[1] = 1;
          discreteNormal[2] = 1;
          discreteNormal[3] = 1;
        }

        else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 1;
        }
      }

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0
          && _blockGeometry->getMaterial(iX, iY + 1, iZ - 1) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ - 1) != 0
          && _blockGeometry->getMaterial(iX - 1, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX - 1, iY + 1, iZ) != 0
          && _blockGeometry->getMaterial(iX - 1, iY + 1, iZ - 1) == 1) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = 1;
          discreteNormal[2] = -1;
          discreteNormal[3] = 1;
        } else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 1;
        }
      }
    }

    // externalCornerPNN and externalCornerPPN
    if (_blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
        && _blockGeometry->getMaterial(iX - 1, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ + 1) != 0) {

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0
          && _blockGeometry->getMaterial(iX, iY + 1, iZ + 1) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ + 1) != 0
          && _blockGeometry->getMaterial(iX - 1, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX - 1, iY + 1, iZ) != 0
          && _blockGeometry->getMaterial(iX - 1, iY + 1, iZ + 1) == 1) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = 1;
          discreteNormal[2] = -1;
          discreteNormal[3] = -1;
        } else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = -1;
        }
      }

      if (_blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0
          && _blockGeometry->getMaterial(iX, iY - 1, iZ + 1) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ + 1) != 0
          && _blockGeometry->getMaterial(iX - 1, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX - 1, iY - 1, iZ) != 0
          && _blockGeometry->getMaterial(iX - 1, iY - 1, iZ + 1) == 1) {

        if (discreteNormal == nullVector) {

          discreteNormal[0] = 1;
          discreteNormal[1] = 1;
          discreteNormal[2] = 1;
          discreteNormal[3] = -1;
        }

        else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = -1;
        }
      }
    }

    // internalCornerPPP and internalCornerPNP
    if (_blockGeometry->getMaterial(iX - 1, iY, iZ) == 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) == 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0) {

      if (_blockGeometry->getMaterial(iX, iY - 1, iZ) == 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = 1;
          discreteNormal[2] = 1;
          discreteNormal[3] = 1;
        } else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 1;
        }
      }

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ) == 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = 1;
          discreteNormal[2] = -1;
          discreteNormal[3] = 1;
        }

        else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 1;
        }
      }
    }

    // internalCornerPNN and InternalCornerPPN
    if (_blockGeometry->getMaterial(iX - 1, iY, iZ) == 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) == 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0) {

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ) == 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = 1;
          discreteNormal[2] = -1;
          discreteNormal[3] = -1;
        } else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = -1;
        }
      }

      if (_blockGeometry->getMaterial(iX, iY - 1, iZ) == 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = 1;
          discreteNormal[2] = 1;
          discreteNormal[3] = -1;
        } else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = -1;
        }
      }
    }

    // internalCornerNPP and internalCornerNNP
    if (_blockGeometry->getMaterial(iX + 1, iY, iZ) == 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) == 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0) {

      if (_blockGeometry->getMaterial(iX, iY - 1, iZ) == 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = -1;
          discreteNormal[2] = 1;
          discreteNormal[3] = 1;
        } else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 1;
        }
      }

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ) == 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = -1;
          discreteNormal[2] = -1;
          discreteNormal[3] = 1;
        }

        else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 1;
        }
      }
    }

    // internalCornerNPN and internalCornerNNN
    if (_blockGeometry->getMaterial(iX + 1, iY, iZ) == 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) == 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0) {

      if (_blockGeometry->getMaterial(iX, iY - 1, iZ) == 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0) {

        if (discreteNormal == nullVector) {

          discreteNormal[0] = 2;
          discreteNormal[1] = -1;
          discreteNormal[2] = 1;
          discreteNormal[3] = -1;
        }

        else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = -1;
        }
      }

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ) == 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = -1;
          discreteNormal[2] = -1;
          discreteNormal[3] = -1;
        } else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = -1;
        }
      }
    }

    // externalEdge0PN and externalEdge0NN
    if (_blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
        && _blockGeometry->getMaterial(iX + 1, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ + 1) != 1) {

      if (_blockGeometry->getMaterial(iX, iY - 1, iZ + 1) == 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 0;
          discreteNormal[2] = 1;
          discreteNormal[3] = -1;
        } else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = -1;
        }
      }

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ + 1) == 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 0;
          discreteNormal[2] = -1;
          discreteNormal[3] = -1;
        } else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = -1;
        }
      }
    }

    // externalEdge0NP and externalEdge0PP
    if (_blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
        && _blockGeometry->getMaterial(iX + 1, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ - 1) != 1) {

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ - 1) == 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 0;
          discreteNormal[2] = -1;
          discreteNormal[3] = 1;
        } else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 1;
        }
      }

      if (_blockGeometry->getMaterial(iX, iY - 1, iZ - 1) == 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 0;
          discreteNormal[2] = 1;
          discreteNormal[3] = 1;
        } else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 1;
        }
      }
    }

    // externalEdge1NN and externalEdge1NP
    if (_blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0) {

      if (_blockGeometry->getMaterial(iX + 1, iY, iZ + 1) == 1
          && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
          && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0
          && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = -1;
          discreteNormal[2] = 0;
          discreteNormal[3] = -1;
        } else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = -1;
        }
      }

      if (_blockGeometry->getMaterial(iX - 1, iY, iZ + 1) == 1
          && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
          && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0
          && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 1;
          discreteNormal[2] = 0;
          discreteNormal[3] = -1;
        } else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = -1;
        }
      }
    }

    // externalEdge1PN and externalEdge1PP
    if (_blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0) {

      if (_blockGeometry->getMaterial(iX + 1, iY, iZ - 1) == 1
          && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
          && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0
          && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = -1;
          discreteNormal[2] = 0;
          discreteNormal[3] = 1;
        } else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = 1;
        }
      }

      if (_blockGeometry->getMaterial(iX - 1, iY, iZ - 1) == 1
          && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
          && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0
          && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 1;
          discreteNormal[2] = 0;
          discreteNormal[3] = 1;
        } else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = 1;
        }
      }
    }

    // externalEdge2NN and externalEdge2PN
    if (_blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0) {

      if (_blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
          && _blockGeometry->getMaterial(iX + 1, iY + 1, iZ) == 1) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = -1;
          discreteNormal[2] = -1;
          discreteNormal[3] = 0;
        } else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 0;
        }
      }

      if (_blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
          && _blockGeometry->getMaterial(iX - 1, iY + 1, iZ) == 1) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 1;
          discreteNormal[2] = -1;
          discreteNormal[3] = 0;
        } else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 0;
        }
      }
    }

    // externalEdge2PP and externalEdge2NP
    if (_blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0) {

      if (_blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
          && _blockGeometry->getMaterial(iX - 1, iY - 1, iZ) == 1) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 1;
          discreteNormal[2] = 1;
          discreteNormal[3] = 0;
        } else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 0;
        }
      }

      if (_blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
          && _blockGeometry->getMaterial(iX + 1, iY - 1, iZ) == 1) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = -1;
          discreteNormal[2] = 1;
          discreteNormal[3] = 0;
        } else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 0;
        }
      }
    }

    // internalEdge0NN and internalEdge0PN
    if (_blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) == 1) {

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ) == 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0
          && _blockGeometry->getMaterial(iX - 1, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX - 1, iY - 1, iZ) != 0
          && _blockGeometry->getMaterial(iX + 1, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX + 1, iY - 1, iZ) != 0) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 0;
        discreteNormal[2] = -1;
        discreteNormal[3] = -1;
      }
      if (_blockGeometry->getMaterial(iX, iY - 1, iZ) == 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0
          && _blockGeometry->getMaterial(iX - 1, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX - 1, iY + 1, iZ) != 0
          && _blockGeometry->getMaterial(iX + 1, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX + 1, iY + 1, iZ) != 0) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 0;
        discreteNormal[2] = 1;
        discreteNormal[3] = -1;
      }
    }

    // internalEdge0NP and internalEdge0PP
    if (_blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) == 1) {

      if (_blockGeometry->getMaterial(iX, iY + 1, iZ) == 1
          && _blockGeometry->getMaterial(iX - 1, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX - 1, iY - 1, iZ) != 0
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0
          && _blockGeometry->getMaterial(iX + 1, iY - 1, iZ) != 1
          && _blockGeometry->getMaterial(iX + 1, iY - 1, iZ) != 0) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 0;
        discreteNormal[2] = -1;
        discreteNormal[3] = 1;
      }

      if (_blockGeometry->getMaterial(iX, iY - 1, iZ) == 1
          && _blockGeometry->getMaterial(iX - 1, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX - 1, iY + 1, iZ) != 0
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0
          && _blockGeometry->getMaterial(iX + 1, iY + 1, iZ) != 1
          && _blockGeometry->getMaterial(iX + 1, iY + 1, iZ) != 0) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 0;
        discreteNormal[2] = 1;
        discreteNormal[3] = 1;
      }
    }

    // internalEdge1PP and internalEdge 1NP
    if (_blockGeometry->getMaterial(iX - 1, iY, iZ) == 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX + 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0) {

      if (_blockGeometry->getMaterial(iX, iY, iZ - 1) == 1
          && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
          && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
          && _blockGeometry->getMaterial(iX, iY + 1, iZ + 1) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ + 1) != 0
          && _blockGeometry->getMaterial(iX, iY - 1, iZ + 1) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ + 1) != 0) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 1;
        discreteNormal[2] = 0;
        discreteNormal[3] = 1;
      }

      if (_blockGeometry->getMaterial(iX, iY, iZ + 1) == 1
          && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
          && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
          && _blockGeometry->getMaterial(iX, iY + 1, iZ - 1) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ - 1) != 0
          && _blockGeometry->getMaterial(iX, iY - 1, iZ - 1) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ - 1) != 0) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 1;
        discreteNormal[2] = 0;
        discreteNormal[3] = -1;
      }
    }

    // internalEdge1PN and internalEdge1NN
    if (_blockGeometry->getMaterial(iX + 1, iY, iZ) == 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 1
        && _blockGeometry->getMaterial(iX - 1, iY, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0) {

      if (_blockGeometry->getMaterial(iX, iY, iZ - 1) == 1
          && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
          && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
          && _blockGeometry->getMaterial(iX, iY + 1, iZ + 1) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ + 1) != 0
          && _blockGeometry->getMaterial(iX, iY - 1, iZ + 1) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ + 1) != 0) {

        discreteNormal[0] = 4;
        discreteNormal[1] = -1;
        discreteNormal[2] = 0;
        discreteNormal[3] = 1;
      }

      if (_blockGeometry->getMaterial(iX, iY, iZ + 1) == 1
          && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
          && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
          && _blockGeometry->getMaterial(iX, iY + 1, iZ - 1) != 1
          && _blockGeometry->getMaterial(iX, iY + 1, iZ - 1) != 0
          && _blockGeometry->getMaterial(iX, iY - 1, iZ - 1) != 1
          && _blockGeometry->getMaterial(iX, iY - 1, iZ - 1) != 0) {

        discreteNormal[0] = 4;
        discreteNormal[1] = -1;
        discreteNormal[2] = 0;
        discreteNormal[3] = -1;
      }
    }

    // internalEdge2PP and internalEdge2NP
    if (_blockGeometry->getMaterial(iX, iY - 1, iZ) == 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY + 1, iZ) != 0) {

      if (_blockGeometry->getMaterial(iX - 1, iY, iZ) == 1
          && _blockGeometry->getMaterial(iX - 1, iY - 1, iZ) == 1) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 1;
        discreteNormal[2] = 1;
        discreteNormal[3] = 0;
      }

      if (_blockGeometry->getMaterial(iX + 1, iY, iZ) == 1
          && _blockGeometry->getMaterial(iX + 1, iY - 1, iZ) == 1) {

        discreteNormal[0] = 4;
        discreteNormal[1] = -1;
        discreteNormal[2] = 1;
        discreteNormal[3] = 0;
      }
    }

    // internalEdge2PN and internalEdge2NN
    if (_blockGeometry->getMaterial(iX, iY + 1, iZ) == 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ - 1) != 0
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 1
        && _blockGeometry->getMaterial(iX, iY, iZ + 1) != 0
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 1
        && _blockGeometry->getMaterial(iX, iY - 1, iZ) != 0) {

      if (_blockGeometry->getMaterial(iX - 1, iY, iZ) == 1
          && _blockGeometry->getMaterial(iX - 1, iY + 1, iZ) == 1) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 1;
        discreteNormal[2] = -1;
        discreteNormal[3] = 0;
      }

      if (_blockGeometry->getMaterial(iX + 1, iY, iZ) == 1
          && _blockGeometry->getMaterial(iX + 1, iY + 1, iZ) == 1) {

        discreteNormal[0] = 4;
        discreteNormal[1] = -1;
        discreteNormal[2] = -1;
        discreteNormal[3] = 0;
      }
    }

    // special boundary conditions
    if (discreteNormal2 != nullVector && anyNormal == false) {
      discreteNormal = checkExtraBoundary(discreteNormal, discreteNormal2);
    }
  }

  if (discreteNormal[1] == 0 && discreteNormal[2] == 0 && discreteNormal[3] == 0) {
    clout<<"WARNING: no discreteNormal is found"<<std::endl;
  } else if (_blockGeometry->getMaterial(iX-discreteNormal[1], iY-discreteNormal[2], iZ-discreteNormal[3]) != 1) {
#ifdef OLB_DEBUG
    clout<<"WARNING: discreteNormal is not pointing outside the fluid. Use option: anyNormal"<<std::endl;
#endif
  }

  return discreteNormal;
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::computeNormal(int iX, int iY, int iZ)
{
  update();
  std::vector<int> normal (3,int(0));

  if (iX != 0) {
    if (_blockGeometry->getMaterial(iX - 1, iY, iZ) == 1) {
      normal[0] = -1;
    }
  }
  if (iX != _nX - 1) {
    if (_blockGeometry->getMaterial(iX + 1, iY, iZ) == 1) {
      normal[0] = 1;
    }
  }
  if (iY != 0) {
    if (_blockGeometry->getMaterial(iX, iY - 1, iZ) == 1) {
      normal[1] = -1;
    }
  }
  if (iY != _nY - 1) {
    if (_blockGeometry->getMaterial(iX, iY + 1, iZ) == 1) {
      normal[1] = 1;
    }
  }
  if (iZ != 0) {
    if (_blockGeometry->getMaterial(iX, iY, iZ - 1) == 1) {
      normal[2] = -1;
    }
  }
  if (iZ != _nZ - 1) {
    if (_blockGeometry->getMaterial(iX, iY, iZ + 1) == 1) {
      normal[2] = 1;
    }
  }
  return normal;
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::computeNormal(int material)
{

  update();
  std::vector<T> normal (3,int(0));
  std::vector<int> minC = getMinLatticeR(material);
  std::vector<int> maxC = getMaxLatticeR(material);
  for (int iX = minC[0]; iX<=maxC[0]; iX++) {
    for (int iY = minC[1]; iY<=maxC[1]; iY++) {
      for (int iZ = minC[2]; iZ<=maxC[2]; iZ++) {
        if (_blockGeometry->getMaterial(iX,iY,iZ) == material) {
          normal[0]+=computeNormal(iX,iY,iZ)[0];
          normal[1]+=computeNormal(iX,iY,iZ)[1];
          normal[2]+=computeNormal(iX,iY,iZ)[2];
        }
      }
    }
  }
  T norm = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  if (norm>0.) {
    normal[0]/=norm;
    normal[1]/=norm;
    normal[2]/=norm;
  }
  return normal;
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::computeDiscreteNormal(int material, T maxNorm)
{

  update();
  std::vector<T> normal = computeNormal(material);
  std::vector<int> discreteNormal(3,int(0));

  T smallestAngle = T(0);
  for (int iX = -1; iX<=1; iX++) {
    for (int iY = -1; iY<=1; iY++) {
      for (int iZ = -1; iZ<=1; iZ++) {
        T norm = sqrt(iX*iX+iY*iY+iZ*iZ);
        if (norm>0.&& norm<maxNorm) {
          T angle = (iX*normal[0] + iY*normal[1] + iZ*normal[2])/norm;
          if (angle>=smallestAngle) {
            smallestAngle=angle;
            discreteNormal[0] = iX;
            discreteNormal[1] = iY;
            discreteNormal[2] = iZ;
          }
        }
      }
    }
  }
  return discreteNormal;
}


template<typename T>
bool BlockGeometryStatistics3D<T>::check(int material, int iX, int iY,
    int iZ, unsigned offsetX, unsigned offsetY, unsigned offsetZ)
{
  update();
  bool found = true;
  for (int iOffsetX = -offsetX; iOffsetX <= (int) offsetX; ++iOffsetX) {
    for (int iOffsetY = -offsetY; iOffsetY <= (int) offsetY; ++iOffsetY) {
      for (int iOffsetZ = -offsetZ; iOffsetZ <= (int) offsetZ; ++iOffsetZ) {
        if (_blockGeometry->getMaterial(iX + iOffsetX, iY + iOffsetY,
                                        iZ + iOffsetZ) != material) {
          found = false;
        }
      }
    }
  }
  return found;
}

template<typename T>
bool BlockGeometryStatistics3D<T>::find(int material, unsigned offsetX,
                                        unsigned offsetY, unsigned offsetZ, int& foundX, int& foundY,
                                        int& foundZ)
{
  update();
  bool found = false;
  for (foundX = 0; foundX < _nX; foundX++) {
    for (foundY = 0; foundY < _nY; foundY++) {
      for (foundZ = 0; foundZ < _nZ; foundZ++) {
        found = check(material, foundX, foundY, foundZ, offsetX,
                      offsetY, offsetZ);
        if (found) {
          return found;
        }
      }
    }
  }
  return found;
}


template<typename T>
void BlockGeometryStatistics3D<T>::print()
{

  update();
  std::map<int, int>::iterator iter;
  for (iter = _material2n.begin(); iter != _material2n.end(); iter++) {
    clout << "materialNumber=" << iter->first
          << "; count=" << iter->second
          << "; minLatticeR=(" << _material2min[iter->first][0] <<","<< _material2min[iter->first][1] <<","<< _material2min[iter->first][2] <<")"
          << "; maxLatticeR=(" << _material2max[iter->first][0] <<","<< _material2max[iter->first][1] <<","<< _material2max[iter->first][2] <<")"
          << std::endl;
  }
}


template<typename T>
void BlockGeometryStatistics3D<T>::takeStatistics(int iX, int iY, int iZ)
{

  int type = _blockGeometry->getMaterial(iX, iY, iZ);
  if (_material2n.count(type) == 0) {
    _material2n[type] = 1;
    std::vector<int> minCo;
    std::vector<int> maxCo;
    minCo.push_back(iX);
    minCo.push_back(iY);
    minCo.push_back(iZ);
    _material2min[type] = minCo;
    maxCo.push_back(iX);
    maxCo.push_back(iY);
    maxCo.push_back(iZ);
    _material2max[type] = maxCo;

  } else {
    _material2n[type]++;
    if (iX < _material2min[type][0]) {
      _material2min[type][0] = iX;
    }
    if (iY < _material2min[type][1]) {
      _material2min[type][1] = iY;
    }
    if (iZ < _material2min[type][2]) {
      _material2min[type][2] = iZ;
    }
    if (iX > _material2max[type][0]) {
      _material2max[type][0] = iX;
    }
    if (iY > _material2max[type][1]) {
      _material2max[type][1] = iY;
    }
    if (iZ > _material2max[type][2]) {
      _material2max[type][2] = iZ;
    }
  }
}

// This function compares two discrete normals (discreteNormal, discreteNormal2) in case of a duplicate assignment of boundary types.
// The goal of this function is to combine these special boundaryVoxels to an existing one (in this case boundary or externalEdge) according to
// the x-, y- and z-values of their discrete normals.
// In the following the algorithm is declared only for the x value, but it is also valid for the y and z values.
//
// for x1 = x2, the new value of x is x1 (1)
// for x1*x2 = -1, the new value of x is 0 (2)
// for x1*x2 = 0, the new value is 0   (3)
//
// It may be possible that all three values equal 0. To avoid that the values are tested again (x²+y²+z²==0) and the loosest assumption (3) is
// redefined to.
//
// If x,y and z == 0 --> find those x,y or z which are 0 because of (3) and redefine them to the value !=0
//
// Additionally the calculated entries are multiplied with (-1) to get the right existing boundary.

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::checkExtraBoundary(
  std::vector<int> discreteNormal, std::vector<int> discreteNormal2)
{

  update();
  std::vector<int> Data(6, 0);

  for (int i = 1; i < 4; i++) {
    if (discreteNormal[i] == discreteNormal2[i]) {
      Data[i - 1] = discreteNormal[i];
      Data[i + 2] = 1;
    } else if (discreteNormal[i] * discreteNormal2[i] == -1) {
      Data[i - 1] = 0;
      Data[i + 2] = 2;
    } else if (discreteNormal[i] * discreteNormal2[i] == 0) {
      Data[i - 1] = 0;
      Data[i + 2] = 3;
    }
  }

  std::vector<int> newDiscreteNormal(4, 0);

  for (int i = 1; i < 4; i++) {
    newDiscreteNormal[i] = Data[i - 1];
  }

  if (Data[0] * Data[0] + Data[1] * Data[1] + Data[2] * Data[2] == 0) {
    for (int i = 1; i < 4; i++) {
      if (Data[i + 2] == 3) {
        if (discreteNormal[i] == 0) {
          newDiscreteNormal[i] = (-1) * discreteNormal2[i];
        } else {
          newDiscreteNormal[i] = (-1) * discreteNormal[i];
        }
      }
    }
  }

  if (newDiscreteNormal[1] * newDiscreteNormal[1] + newDiscreteNormal[2]
      * newDiscreteNormal[2] + newDiscreteNormal[3]
      * newDiscreteNormal[3] == 1) {
    newDiscreteNormal[0] = 0;
  } else {
    newDiscreteNormal[0] = 3;
  }
  return newDiscreteNormal;
}


} // namespace olb

#endif
