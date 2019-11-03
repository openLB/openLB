/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2014 Mathias J. Krause
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
 * Representation of a statistic for a parallel 2D geometry -- generic implementation.
 */

#include <iostream>
#include <cmath>
#include <math.h>
#include <sstream>

#include "geometry/superGeometry2D.h"
#include "geometry/superGeometryStatistics2D.h"
#include "core/olbDebug.h"

namespace olb {

template<typename T>
SuperGeometryStatistics2D<T>::SuperGeometryStatistics2D(SuperGeometry2D<T>* superGeometry)
  : _superGeometry(superGeometry), _statisticsUpdateNeeded(true), _overlap(superGeometry->getOverlap()), clout(std::cout,"SuperGeometryStatistics2D")
{
}

template<typename T>
SuperGeometryStatistics2D<T>::SuperGeometryStatistics2D(SuperGeometryStatistics2D const& rhs)
  : _superGeometry(rhs._superGeometry), _statisticsUpdateNeeded(true), _overlap(rhs._superGeometry->getOverlap() ), clout(std::cout,"SuperGeometryStatistics2D")
{
}

template<typename T>
SuperGeometryStatistics2D<T>& SuperGeometryStatistics2D<T>::operator=(SuperGeometryStatistics2D const& rhs)
{
  _superGeometry = rhs._superGeometry;
  _statisticsUpdateNeeded = true;
  _overlap = rhs._overlap;
  return *this;
}


template<typename T>
bool& SuperGeometryStatistics2D<T>::getStatisticsStatus()
{
  return _statisticsUpdateNeeded;
}

template<typename T>
bool const & SuperGeometryStatistics2D<T>::getStatisticsStatus() const
{
  return _statisticsUpdateNeeded;
}


template<typename T>
void SuperGeometryStatistics2D<T>::update(bool verbose)
{

#ifdef PARALLEL_MODE_MPI
  int updateReallyNeededGlobal = 0;
  if (_statisticsUpdateNeeded ) {
    updateReallyNeededGlobal = 1;
  }
  singleton::mpi().reduceAndBcast(updateReallyNeededGlobal, MPI_SUM);
  if (updateReallyNeededGlobal>0) {
    _statisticsUpdateNeeded = true;
  }
  //singleton::mpi().reduceAndBcast(_statisticsUpdateNeeded, MPI_LOR);
#endif

  // check if update is really needed
  if (_statisticsUpdateNeeded ) {
    int updateReallyNeeded = 0;
    for (int iCloc=0; iCloc<_superGeometry->getLoadBalancer().size(); iCloc++) {
      if (_superGeometry->getBlockGeometry(iCloc).getStatistics().getStatisticsStatus() ) {
        _superGeometry->getBlockGeometry(iCloc).getStatistics().update(false);
        updateReallyNeeded++;
      }
    }


#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(updateReallyNeeded, MPI_SUM);
#endif

    if (updateReallyNeeded==0) {
      _statisticsUpdateNeeded = false;
      //      clout << "almost updated" << std::endl;
      return;
    }


    // get total number of different materials right
    _material2n = std::map<int, int>();
    _nMaterials = int();
    for (int iCloc=0; iCloc<_superGeometry->getLoadBalancer().size(); iCloc++) {
      if (_superGeometry->getBlockGeometry(iCloc).getStatistics(false).getNmaterials() > 0) {
        _nMaterials = _superGeometry->getBlockGeometry(iCloc).getStatistics(false).getNmaterials();
      }
    }


#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(_nMaterials, MPI_SUM);
#endif

    // store the number and min., max. possition for each rank
    for (int iCloc=0; iCloc<_superGeometry->getLoadBalancer().size(); iCloc++) {
      std::map<int, int> material2n = _superGeometry->getBlockGeometry(iCloc).getStatistics(false).getMaterial2n();
      std::map<int, int>::iterator iter;

      for (iter = material2n.begin(); iter != material2n.end(); iter++) {
        if (iter->second!=0) {
          std::vector<T> minPhysR = _superGeometry->getBlockGeometry(iCloc).getStatistics(false).getMinPhysR(iter->first);
          std::vector<T> maxPhysR = _superGeometry->getBlockGeometry(iCloc).getStatistics(false).getMaxPhysR(iter->first);
          if (_material2n.count(iter->first) == 0) {
            _material2n[iter->first] = iter->second;
            _material2min[iter->first] = minPhysR;
            _material2max[iter->first] = maxPhysR;
            //std::cout << iter->first<<":"<<_material2n[iter->first]<<std::endl;
          } else {
            _material2n[iter->first] += iter->second;
            for (int iDim=0; iDim<2; iDim++) {
              if (_material2min[iter->first][iDim] > minPhysR[iDim]) {
                _material2min[iter->first][iDim] = minPhysR[iDim];
              }
              if (_material2max[iter->first][iDim] < maxPhysR[iDim]) {
                _material2max[iter->first][iDim] = maxPhysR[iDim];
              }
            }
            //std::cout << iter->first<<":"<<_material2n[iter->first]<<std::endl;
          }
        }
      }
    }

    // store the number and min., max. possition for all ranks
#ifdef PARALLEL_MODE_MPI
    int materials[_nMaterials];
    int materialsInBuf[_nMaterials];
    int materialCount[_nMaterials];
    int materialCountInBuf[_nMaterials];
    T materialMinR[2*_nMaterials];
    T materialMaxR[2*_nMaterials];
    T materialMinRinBuf[2*_nMaterials];
    T materialMaxRinBuf[2*_nMaterials];

    for (int iM=0; iM<_nMaterials; iM++) {
      materials[iM]=-1;
      materialCount[iM]=0;
      for (int iDim=0; iDim<2; iDim++) {
        materialMinRinBuf[2*iM + iDim] = T();
        materialMaxRinBuf[2*iM + iDim] = T();
      }
    }
    int counter = 0;
    std::map<int, int>::iterator iMaterial;
    for (iMaterial = _material2n.begin(); iMaterial != _material2n.end(); iMaterial++) {
      materials[counter] = iMaterial->first;
      materialCount[counter] = iMaterial->second;
      for (int iDim=0; iDim<2; iDim++) {
        materialMinR[2*counter + iDim] = _material2min[iMaterial->first][iDim];
        materialMaxR[2*counter + iDim] = _material2max[iMaterial->first][iDim];
      }
      counter++;
    }

    for (int iRank=1; iRank<singleton::mpi().getSize(); iRank++) {
      int myRank = singleton::mpi().getRank();
      singleton::mpi().sendRecv(materials, materialsInBuf, _nMaterials, (myRank+iRank)%singleton::mpi().getSize(), (myRank-iRank+singleton::mpi().getSize())%singleton::mpi().getSize(), 0);
      singleton::mpi().sendRecv(materialCount, materialCountInBuf, _nMaterials, (myRank+iRank)%singleton::mpi().getSize(), (myRank-iRank+singleton::mpi().getSize())%singleton::mpi().getSize(), 1);
      singleton::mpi().sendRecv(materialMinR, materialMinRinBuf, 2*_nMaterials, (myRank+iRank)%singleton::mpi().getSize(), (myRank-iRank+singleton::mpi().getSize())%singleton::mpi().getSize(), 2);
      singleton::mpi().sendRecv(materialMaxR, materialMaxRinBuf, 2*_nMaterials, (myRank+iRank)%singleton::mpi().getSize(), (myRank-iRank+singleton::mpi().getSize())%singleton::mpi().getSize(), 2);
      for (int iM=0; iM<_nMaterials; iM++) {
        if (materialsInBuf[iM]!=-1) {
          std::vector<T> minPhysR(2,T());
          std::vector<T> maxPhysR(2,T());
          for (int iDim=0; iDim<2; iDim++) {
            minPhysR[iDim] = materialMinRinBuf[2*iM + iDim];
            maxPhysR[iDim] = materialMaxRinBuf[2*iM + iDim];
          }
          if (_material2n.count(materialsInBuf[iM]) == 0) {
            _material2n[materialsInBuf[iM]] = materialCountInBuf[iM];
            _material2min[materialsInBuf[iM]] = minPhysR;
            _material2max[materialsInBuf[iM]] = maxPhysR;
          } else {
            _material2n[materialsInBuf[iM]] += materialCountInBuf[iM];
            for (int iDim=0; iDim<2; iDim++) {
              if (_material2min[materialsInBuf[iM]][iDim] > minPhysR[iDim]) {
                _material2min[materialsInBuf[iM]][iDim] = minPhysR[iDim];
              }
              if (_material2max[materialsInBuf[iM]][iDim] < maxPhysR[iDim]) {
                _material2max[materialsInBuf[iM]][iDim] = maxPhysR[iDim];
              }
            }
          }
        }
      }
    }
#endif

    //clout.setMultiOutput(true);
    //print();
    //clout.setMultiOutput(false);

    if (verbose) {
      clout << "updated" << std::endl;
    }
    _statisticsUpdateNeeded = false;
  }
}

template<typename T>
int SuperGeometryStatistics2D<T>::getNmaterials()
{
  update();
  return _nMaterials;
}

template<typename T>
int SuperGeometryStatistics2D<T>::getNvoxel(int material)
{
  update(true);
  return _material2n[material];
}

template<typename T>
int SuperGeometryStatistics2D<T>::getNvoxel()
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
std::vector<T> SuperGeometryStatistics2D<T>::getMinPhysR(int material)
{
  update();
  return _material2min[material];
}

template<typename T>
std::vector<T> SuperGeometryStatistics2D<T>::getMaxPhysR(int material)
{
  update();
  return _material2max[material];
}

template<typename T>
std::vector<T> SuperGeometryStatistics2D<T>::getPhysExtend(int material)
{
  update();
  std::vector<T> extend;
  for (int iDim = 0; iDim < 2; iDim++) {
    extend.push_back(_material2max[material][iDim] - _material2min[material][iDim]);
  }
  return extend;
}

template<typename T>
std::vector<T> SuperGeometryStatistics2D<T>::getPhysRadius(int material)
{
  update();
  std::vector<T> radius;
  for (int iDim=0; iDim<2; iDim++) {
    radius.push_back((getMaxPhysR(material)[iDim] - getMinPhysR(material)[iDim])/2.);
  }
  return radius;
}

template<typename T>
std::vector<T> SuperGeometryStatistics2D<T>::getCenterPhysR(int material)
{
  update();
  std::vector<T> center;
  for (int iDim=0; iDim<2; iDim++) {
    center.push_back(getMinPhysR(material)[iDim] + getPhysRadius(material)[iDim]);
  }
  return center;
}

template<typename T>
std::vector<int> SuperGeometryStatistics2D<T>::getType(int iC, int iX, int iY)
{
  update();
  int iCloc=_superGeometry->getLoadBalancer().loc(iC);
  std::vector<int> discreteNormal = _superGeometry->getExtendedBlockGeometry(iCloc).getStatistics(false).getType(iX+_overlap, iY+_overlap);
  return discreteNormal;
}

template<typename T>
std::vector<T> SuperGeometryStatistics2D<T>::computeNormal(int material)
{
  update();

  std::vector<T> normal (2,int());
  for (int iCloc=0; iCloc<_superGeometry->getLoadBalancer().size(); iCloc++) {
    for (int iDim=0; iDim<2; iDim++) {
      if (_superGeometry->getBlockGeometry(iCloc).getStatistics(false).getNvoxel(material)!=0) {
        normal[iDim] += _superGeometry->getBlockGeometry(iCloc).getStatistics(false).computeNormal(material)[iDim]*_superGeometry->getBlockGeometry(iCloc).getStatistics(false).getNvoxel(material);
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (int iDim=0; iDim<2; iDim++) {
    singleton::mpi().reduceAndBcast((normal[iDim]), MPI_SUM);
  }
#endif

  int nVoxel = getNvoxel(material);
  if (nVoxel != 0) {
    for (int iDim=0; iDim<2; iDim++) {
      normal[iDim] /= nVoxel;
    }
  }
  OLB_ASSERT(nVoxel || (normal[0] == 0 && normal[1] == 0), "if no voxels found we expect the normal to be zero");

  T norm = sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
  if (norm>0.) {
    normal[0]/=norm;
    normal[1]/=norm;
  }
  return normal;
}

template<typename T>
std::vector<int> SuperGeometryStatistics2D<T>::computeDiscreteNormal(int material, T maxNorm)
{
  update();
  std::vector<T> normal = computeNormal(material);
  std::vector<int> discreteNormal(2,int(0));

  T smallestAngle = T(0);
  for (int iX = -1; iX<=1; iX++) {
    for (int iY = -1; iY<=1; iY++) {
      T norm = sqrt(iX*iX+iY*iY);
      if (norm>0.&& norm<maxNorm) {
        T angle = (iX*normal[0] + iY*normal[1])/norm;
        if (angle>=smallestAngle) {
          smallestAngle=angle;
          discreteNormal[0] = iX;
          discreteNormal[1] = iY;
        }
      }
    }
  }
  return discreteNormal;
}


template<typename T>
void SuperGeometryStatistics2D<T>::print()
{
  update();
  std::map<int, int>::iterator iter;
  for (iter = _material2n.begin(); iter != _material2n.end(); iter++) {
    clout << "materialNumber=" << iter->first
          << "; count=" << iter->second
          << "; minPhysR=(" << _material2min[iter->first][0] <<","<< _material2min[iter->first][1] <<")"
          << "; maxPhysR=(" << _material2max[iter->first][0] <<","<< _material2max[iter->first][1] <<")"
          << std::endl;
  }
}


} // namespace olb
