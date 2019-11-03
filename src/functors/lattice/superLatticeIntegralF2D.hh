/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause, Adrian Kummerlaender
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

#ifndef SUPER_LATTICE_INTEGRAL_F_2D_HH
#define SUPER_LATTICE_INTEGRAL_F_2D_HH

#include<vector>
#include<cmath>

#include "superLatticeIntegralF2D.h"
#include "blockLatticeIntegralF2D.h"
#include "indicator/superIndicatorF2D.h"
#include "utilities/vectorHelpers.h"
#include "core/vector.h"

using namespace olb::util;

namespace olb {


template <typename T>
SuperMax2D<T>::SuperMax2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Max("+_f.getName()+")";
}

template <typename T>
bool SuperMax2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i]=T();
    for (int iC = 0; iC < load.size(); ++iC) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          if (_superGeometry.get(load.glob(iC), iX, iY) == _material) {
            T outputTmp[_f.getTargetDim()];
            _f(outputTmp,load.glob(iC),iX,iY);
            if (fabs(outputTmp[i]) > output[i]) {
              output[i] = fabs(outputTmp[i]);
            }
          }
        }
      }
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(output[i], MPI_MAX);
#endif
  }
  return true;
}

template <typename T>
SuperMin2D<T>::SuperMin2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Min("+_f.getName()+")";
}

template <typename T>
bool SuperMin2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i]=T();
    for (int iC = 0; iC < load.size(); ++iC) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          if (_superGeometry.get(load.glob(iC), iX, iY) == _material) {
            T outputTmp[_f.getTargetDim()];
            _f(outputTmp,load.glob(iC),iX,iY);
            if (fabs(outputTmp[i]) < output[i]) {
              output[i] = fabs(outputTmp[i]);
            }
          }
        }
      }
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(output[i], MPI_MAX);
#endif
  }
  return true;
}


template <typename T>
SuperSum2D<T>::SuperSum2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()+1),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Sum("+_f.getName()+")";
}

template <typename T>
bool SuperSum2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  int numVoxels(0);
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i]=T();
  }
  for (int iC=0; iC<load.size(); ++iC) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    for (int iX = 0; iX < nX; iX++) {
      for (int iY = 0; iY < nY; iY++) {
        if (this->_superGeometry.get(load.glob(iC), iX, iY) == _material) {
          T outputTmp[_f.getTargetDim()];
          _f(outputTmp,load.glob(iC),iX,iY);
          for (int i = 0; i < this->getTargetDim()-1 /*f.getTargetDim()*/; ++i) {
            output[i] += outputTmp[i];
          }
          numVoxels++;
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim()-1; ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
  output[this->getTargetDim()-1] = numVoxels;
  return true;
}


template <typename T>
SuperIntegral2D<T>::SuperIntegral2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()), _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Integral("+_f.getName()+")";
}

template <typename T>
bool SuperIntegral2D<T>::operator() (T output[], const int input[])
{
  //  f.getSuperStructure().communicate();
  //  CuboidGeometry2D<T>& cGeometry = f.getSuperStructure().getCuboidGeometry();
  //  LoadBalancer<T>& load = f.getSuperStructure().getLoadBalancer();

  //  std::vector<T> tmp(this->_n, T() );
  //  for (int i=0; i<this->_n; ++i) {
  //    for (int iC=0; iC<load.size(); ++iC) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  ////      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      T weight = pow(this->superGeometry.getDeltaR(),3);
  //      for (int iX=0; iX<nX; iX++) {
  //        for (int iY=0; iY<nY; iY++) {
  ////          for (int iZ=0; iZ<nZ; iZ++) {
  //            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  ////            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  ////            if (this->superGeometry.getMaterial(globX, globY) == material) {
  ////              tmp[i]+=f(load.glob(iC),iX,iY)[i]*weight;
  ////            }
  ////            if (this->superGeometry.getMaterial(globX, globY, globZ) == material) {
  ////              tmp[i]+=f(load.glob(iC),iX,iY,iZ)[i]*weight;
  ////            }

  ////          }
  //        }
  //      }
  //    }
  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
  //#endif
  //  }
  //  return tmp;
  return false;
}


template <typename T>
template<typename DESCRIPTOR>
SuperGeometryFaces2D<T>::SuperGeometryFaces2D(SuperGeometry2D<T>& superGeometry,
    const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : GenericF<T,int>(7,3), _superGeometry(superGeometry), _material(material), _latticeL(converter.getConversionFactorLength())
{
  this->getName() = "superGeometryFaces";
}

template <typename T>
SuperGeometryFaces2D<T>::SuperGeometryFaces2D(SuperGeometry2D<T>& superGeometry,
    const int material, T latticeL)
  : GenericF<T,int>(7,3), _superGeometry(superGeometry), _material(material), _latticeL(latticeL)
{
  this->getName() = "superGeometryFaces";
}



template <typename T>
bool SuperGeometryFaces2D<T>::operator() (T output[], const int input[])
{
  for (int iDim = 0; iDim < 7; ++iDim) {
    output[iDim]=T();
  }
  _superGeometry.communicate();
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    BlockGeometryFaces2D<T> f(_superGeometry.getBlockGeometry(iC), _material, _latticeL);
    T outputTmp[f.getTargetDim()];
    f(outputTmp,input);
    for (int iDim = 0; iDim < 7; ++iDim) {
      output[iDim] += outputTmp[iDim];
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int iDim = 0; iDim < 7; ++iDim) {
    singleton::mpi().reduceAndBcast( output[iDim], MPI_SUM);
  }
#endif
  return true;
}

template <typename T, bool HLBM>
SuperGeometryFacesIndicator2D<T,HLBM>::SuperGeometryFacesIndicator2D(SuperGeometry2D<T>& superGeometry,
    SmoothIndicatorF2D<T,T,HLBM>& indicator, const int material, T deltaX)
  : GenericF<T,int>(7,0), _superGeometry(superGeometry), _indicator(indicator), _material(material), _latticeL(deltaX)
{
  this->getName() = "superGeometryFacesInd";
}

template <typename T, bool HLBM>
bool SuperGeometryFacesIndicator2D<T,HLBM>::operator() (T output[], const int input[])
{
  _superGeometry.communicate();
  for (int iDim = 0; iDim < 7; ++iDim) {
    output[iDim]=T();
  }
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    BlockGeometryFacesIndicator2D<T,HLBM> f(_superGeometry.getBlockGeometry(iC), _indicator, _material, _latticeL);
    T outputTmp[f.getTargetDim()];
    f(outputTmp,input);
    for (int iDim = 0; iDim < 7; ++iDim) {
      output[iDim] += outputTmp[iDim];
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int iDim = 0; iDim < 7; ++iDim) {
    singleton::mpi().reduceAndBcast( output[iDim], MPI_SUM);
  }
#endif
  return true;
}


template <typename T, typename DESCRIPTOR>
SuperLatticePhysDrag2D<T,DESCRIPTOR>::SuperLatticePhysDrag2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
 const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physDrag";
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePhysDrag2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  SuperGeometryFaces2D<T> faces(_superGeometry, _material, this->_converter.getConversionFactorLength());
  SuperLatticePhysBoundaryForce2D<T,DESCRIPTOR> f(this->_sLattice, _superGeometry, _material, this->_converter);
  SuperSum2D<T> sumF(f, _superGeometry, _material);

  T factor = 2. / (this->_converter.getPhysDensity() * this->_converter.getCharPhysVelocity() * this->_converter.getCharPhysVelocity());

  T facesTmp[faces.getTargetDim()], sumFTmp[sumF.getTargetDim()];
  sumF(sumFTmp, input);
  faces(facesTmp, input);
  output[0] = factor * sumFTmp[0] / facesTmp[0];
  output[1] = factor * sumFTmp[1] / facesTmp[1];

  return true;
}

template <typename T, typename DESCRIPTOR>
SuperLatticePhysCorrDrag2D<T,DESCRIPTOR>::SuperLatticePhysCorrDrag2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
 const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physCorrDrag";
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePhysCorrDrag2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  SuperGeometryFaces2D<T> faces(superGeometry, material, this->converter);

  //  SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR> f(this->sLattice, superGeometry, material, this->converter);
  //  SuperSum2D<T> sumF(f, superGeometry, material);

  //  T factor = 2. / (this->converter.getCharRho() * this->converter.getCharU() * this->converter.getCharU());

  //  std::vector<T> drag(2,T());
  //  drag[0] = factor * sumF(input)[0] / faces(input)[0];
  //  drag[1] = factor * sumF(input)[1] / faces(input)[1];
  ////  drag[2] = factor * sumF(input)[2] / faces(input)[2];
  //  return drag;
  return false;
}


} //end namespace olb

#endif
