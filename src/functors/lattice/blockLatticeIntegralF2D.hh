/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_LATTICE_INTEGRAL_F_2D_HH
#define BLOCK_LATTICE_INTEGRAL_F_2D_HH

#include<cmath>
#include<vector>

#include "blockLatticeIntegralF2D.h"
#include "blockLatticeLocalF2D.h"
#include "blockCalcF2D.h" // for IdentityF
#include "utilities/vectorHelpers.h"
#include "indicator/blockIndicatorBaseF2D.h"


namespace olb {


template <typename T, typename DESCRIPTOR>
BlockMax2D<T,DESCRIPTOR>::BlockMax2D(BlockLatticeF2D<T,DESCRIPTOR>& f,
                                     BlockGeometry2D<T>& blockGeometry, int material)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlockLattice(), f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "Max("+_f.getName()+")";
}

template <typename T, typename DESCRIPTOR>
bool BlockMax2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{

  //  f.getBlockLattice2D().communicate();
  //  CuboidGeometry2D<T><T>& cGeometry = f.getBlockLattice2D().get_cGeometry();
  //  loadBalancer& load = f.getBlockLattice2D().get_load();

  for (int i=0; i<this->getTargetDim(); ++i) {
    output[i]=T();
    //    for (int iC=0; iC<load.size(); iC++) {
    //      int nX = cGeometry.get(load.glob(iC)).getNx();
    //      int nY = cGeometry.get(load.glob(iC)).getNy();
    //      int nZ = cGeometry.get(load.glob(iC)).getNz();
    //      for (int iX=0; iX<nX; iX++) {
    //        for (int iY=0; iY<nY; iY++) {
    //          for (int iZ=0; iZ<nZ; iZ++) {
    //            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
    //            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
    //            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
    //            if (BlockGeometry.getMaterial(globX, globY, globZ) == material) {
    //              if (fabs(f(load.glob(iC),iX,iY,iZ)[i]) > tmp[i]) {
    //                tmp[i]=fabs(f(load.glob(iC),iX,iY,iZ)[i]);
    //              }
    //            }
    //          }
    //        }
    //      }
    //    }
    //#ifdef PARALLEL_MODE_MPI
    //    singleton::mpi().reduceAndBcast(tmp[i], MPI_MAX);
    //#endif
  }
  return true;


}

template <typename T, typename DESCRIPTOR>
BlockSum2D<T,DESCRIPTOR>::BlockSum2D(BlockLatticeF2D<T,DESCRIPTOR>& f,
                                     BlockGeometry2D<T>& blockGeometry, int material)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlockLattice(),f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "Sum("+_f.getName()+")";
}

template <typename T, typename DESCRIPTOR>
bool BlockSum2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  BlockIdentity2D<T> ff(_f); // exists only to prevent f from being deleted

  for (int i=0; i<this->getTargetDim(); ++i) {
    output[i]=T();
    //    int nX = f.getBlockLattice2D().getNx();
    //    int nY = f.getBlockLattice2D().getNy();
    //    int nZ = f.getBlockLattice2D().getNz();
    //    for (int iX=0; iX<nX; iX++) {
    //      for (int iY=0; iY<nY; iY++) {
    //        for (int iZ=0; iZ<nZ; iZ++) {
    //          if (this->BlockGeometry.getMaterial(iX, iY, iZ) == material) {
    //            tmp[i]+=f(iX,iY,iZ)[i];
    //          }
    //        }
    //      }
    //    }
  }
  return true;
}


template <typename T, typename DESCRIPTOR>
BlockIntegral2D<T,DESCRIPTOR>::BlockIntegral2D(BlockLatticeF2D<T,DESCRIPTOR>& f,
    BlockGeometry2D<T>& blockGeometry, int material)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlockLattice(),f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "Integral("+_f.getName()+")";
}

template <typename T, typename DESCRIPTOR>
bool BlockIntegral2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{

  //  f.getBlockLattice2D().communicate();
  //  CuboidGeometry2D<T><T>& cGeometry = f.getBlockLattice2D().get_cGeometry();
  //  loadBalancer& load = f.getBlockLattice2D().get_load();

  output[0]=T();
  //  for (int i=0; i<this->n; i++) {
  //    for (int iC=0; iC<load.size(); iC++) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  //      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      T weight = pow(this->BlockGeometry.getDeltaR(),3);
  //      for (int iX=0; iX<nX; iX++) {
  //        for (int iY=0; iY<nY; iY++) {
  //          for (int iZ=0; iZ<nZ; iZ++) {
  //            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  //            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  //            if (this->BlockGeometry.getMaterial(globX, globY, globZ) == material) {
  //              tmp[i]+=f(load.glob(iC),iX,iY,iZ)[i]*weight;
  //            }
  //          }
  //        }
  //      }
  //    }
  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
  //#endif
  //  }
  return true;
}


template <typename T, typename DESCRIPTOR>
BlockL1Norm2D<T,DESCRIPTOR>::BlockL1Norm2D(BlockLatticeF2D<T,DESCRIPTOR>& f,
    BlockGeometry2D<T>& blockGeometry, int material)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlockLattice(),f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "L1("+_f.getName()+")";
}

template <typename T, typename DESCRIPTOR>
bool BlockL1Norm2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{

  //  f.getBlockLattice2D().communicate();
  //  CuboidGeometry2D<T><T>& cGeometry = f.getBlockLattice2D().get_cGeometry();
  //  loadBalancer& load = f.getBlockLattice2D().get_load();

  //  int numVoxels(0);
  output[0]=T();
  //  for (int i=0; i<this->n; i++) {

  //    for (int iC=0; iC<load.size(); iC++) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  //      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      for (int iX=0; iX<nX; iX++) {
  //        for (int iY=0; iY<nY; iY++) {
  //          for (int iZ=0; iZ<nZ; iZ++) {
  //            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  //            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  //            if (this->BlockGeometry.getMaterial(globX, globY, globZ) == material) {
  //              tmp[i]+=fabs(f(load.glob(iC),iX,iY,iZ)[i]);
  //              numVoxels++;
  //            }
  //          }
  //        }
  //      }
  //    }
  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
  //#endif
  //  }
  //#ifdef PARALLEL_MODE_MPI
  //  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
  //#endif
  return true;
}


template <typename T, typename DESCRIPTOR>
BlockL222D<T,DESCRIPTOR>::BlockL222D(BlockLatticeF2D<T,DESCRIPTOR>& f,
                                     BlockGeometry2D<T>& blockGeometry, int material)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlockLattice(),f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "L22("+_f.getName()+")";
}

template <typename T, typename DESCRIPTOR>
bool BlockL222D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{

  //  f.getBlockLattice2D().communicate();
  //  CuboidGeometry2D<T><T>& cGeometry = f.getBlockLattice2D().get_cGeometry();
  //  loadBalancer& load = f.getBlockLattice2D().get_load();

  output[0]=T();
  //  for (int i=0; i<this->n; i++) {

  //    for (int iC=0; iC<load.size(); iC++) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  //      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      T weight = pow(this->BlockGeometry.getDeltaR(),3);
  //      for (int iX=0; iX<nX; iX++) {
  //        for (int iY=0; iY<nY; iY++) {
  //          for (int iZ=0; iZ<nZ; iZ++) {
  //            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  //            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  //            if (this->BlockGeometry.getMaterial(globX, globY, globZ) == material) {
  //              tmp[i]+=f(load.glob(iC),iX,iY,iZ)[i]*f(load.glob(iC),iX,iY,iZ)[i]*weight;
  //            }
  //          }
  //        }
  //      }
  //    }
  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
  //#endif
  //  }
  return true;
}


template <typename T>
template<typename DESCRIPTOR>
BlockGeometryFaces2D<T>::BlockGeometryFaces2D(BlockGeometryStructure2D<T>& blockGeometry, int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : GenericF<T,int>(7,4), _blockGeometry(blockGeometry), _material(material), _latticeL(converter.getConversionFactorLength())
{
  this->getName() = "blockGeometryFaces";
}

template <typename T>
BlockGeometryFaces2D<T>::BlockGeometryFaces2D(BlockGeometryStructure2D<T>& blockGeometry, int material, T latticeL)
  : GenericF<T,int>(7,4), _blockGeometry(blockGeometry), _material(material), _latticeL(latticeL)
{
  this->getName() = "blockGeometryFaces";
}

template <typename T>
bool BlockGeometryFaces2D<T>::operator() (T output[], const int input[])
{
  int counter[7] = {0,0,0,0,0,0,0};
  if (_blockGeometry.getStatistics().getNvoxel(_material)!=0) {
    const int x0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[0];
    const int y0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[1];
    const int x1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[0];
    const int y1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[1];

    // Iterate over all cells and count the cells of the face
    for (int iX = x0; iX <= x1; ++iX) {
      for (int iY = y0; iY <= y1; ++iY) {
        // Lock at solid nodes only
        if (_blockGeometry.getMaterial(iX, iY) == _material) {
          if (_blockGeometry.getMaterial(iX-1, iY) == 1) {
            counter[0]++;
          }
          if (_blockGeometry.getMaterial(iX, iY-1) == 1) {
            counter[1]++;
          }
          if (_blockGeometry.getMaterial(iX+1, iY) == 1) {
            counter[3]++;
          }
          if (_blockGeometry.getMaterial(iX, iY+1) == 1) {
            counter[4]++;
          }
        }
      }
    }

    T dx2 = _latticeL*_latticeL;
    T total = T();
    for (int i=0; i<6; ++i) {
      output[i]=(T) counter[i] * dx2;
      total+=(T) counter[i] * dx2;
    }
    output[6]=total;
    return true;
  } else {
    for (int i=0; i<7; ++i) {
      output[i]=T();
    }
    return true;
  }
  return false;
}


template <typename T, bool HLBM>
BlockGeometryFacesIndicator2D<T,HLBM>::BlockGeometryFacesIndicator2D(BlockGeometryStructure2D<T>& blockGeometry, SmoothIndicatorF2D<T,T,HLBM>& indicator, int material, T latticeL)
  : GenericF<T,int>(7,0), _blockGeometry(blockGeometry), _indicator(indicator), _material(material), _latticeL(latticeL)
{
  this->getName() = "facesInd";
}

template <typename T, bool HLBM>
bool BlockGeometryFacesIndicator2D<T,HLBM>::operator() (T output[], const int input[])
{
  int counter[7] = {0,0,0,0,0,0,0};
  T inside[1];
  T physR[2];
  if (_blockGeometry.getStatistics().getNvoxel(_material)!=0) {
    const int x0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[0];
    const int y0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[1];
    const int x1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[0];
    const int y1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[1];

    // Iterate over all cells and count the cells of the face
    for (int iX = x0; iX <= x1; ++iX) {
      for (int iY = y0; iY <= y1; ++iY) {
        // Look at solid nodes only
        _blockGeometry.getPhysR(physR, iX, iY);
        _indicator(inside, physR);
        if ( !util::nearZero(inside[0]) ) {
          _blockGeometry.getPhysR(physR, iX-1, iY);
          _indicator(inside, physR);
          if ( !util::nearZero(inside[0]) ) {
            counter[0]++;
          }
          _blockGeometry.getPhysR(physR, iX, iY-1);
          _indicator(inside, physR);
          if ( !util::nearZero(inside[0]) ) {
            counter[1]++;
          }
          _blockGeometry.getPhysR(physR, iX+1, iY);
          _indicator(inside, physR);
          if ( !util::nearZero(inside[0]) ) {
            counter[3]++;
          }
          _blockGeometry.getPhysR(physR, iX, iY+1);
          _indicator(inside, physR);
          if ( !util::nearZero(inside[0]) ) {
            counter[4]++;
          }
        }
      }
    }

    T total = T();
    for (int i=0; i<6; ++i) {
      output[i]=(T) counter[i] * _latticeL;
      total+=(T) counter[i] * _latticeL;
    }
    output[6]=total;
    return true;
  } else {
    for (int i=0; i<7; ++i) {
      output[i]=T();
    }
    return true;
  }
  return false;
}


template <typename T, typename DESCRIPTOR>
BlockLatticePhysDrag2D<T,DESCRIPTOR>::BlockLatticePhysDrag2D
(BlockLattice2D<T,DESCRIPTOR>& blockLattice, BlockGeometry2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "physDrag";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysDrag2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  BlockGeometryFaces2D<T> faces(_blockGeometry, _material, this->_converter);
  BlockIndicatorMaterial2D<T> indicator(_blockGeometry, _material);
  BlockLatticePhysBoundaryForce2D<T,DESCRIPTOR> fTemp(this->_blockLattice, indicator, this->_converter);
  BlockSum2D<T,DESCRIPTOR> sumF(fTemp, _blockGeometry, _material);

  T factor = 2. / (this->_converter.getPhysDensity() * this->_converter.getCharPhysVelocity() * this->_converter.getCharPhysVelocity());

  T outputSumF[2] = { T() };
  sumF(outputSumF,input);
  T outputFaces[2] = { T() };
  faces(outputFaces,input);

  output[0] = factor * outputSumF[0] / outputFaces[0];
  output[1] = factor * outputSumF[1] / outputFaces[1];

  return true;
}


template <typename T, typename DESCRIPTOR>
BlockLatticePhysCorrDrag2D<T,DESCRIPTOR>::BlockLatticePhysCorrDrag2D
(BlockLattice2D<T,DESCRIPTOR>& blockLattice, BlockGeometry2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "physCorrDrag";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysCorrDrag2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  BlockGeometryFaces2D<T> faces(_blockGeometry, _material, this->_converter);
  BlockLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR> tTemp(this->_blockLattice, _blockGeometry,
      _material, this->_converter);
  BlockSum2D<T,DESCRIPTOR> sumF(tTemp, _blockGeometry, _material);

  T factor = 2. / (this->_converter.getPhysDensity() * this->_converter.getCharPhysVelocity() * this->_converter.getCharPhysVelocity());

  T outputSumF[2] = { T() };
  sumF(outputSumF,input);
  T outputFaces[2] = { T() };
  faces(outputFaces,input);

  output[0] = factor * outputSumF[0] / outputFaces[0];
  output[1] = factor * outputSumF[1] / outputFaces[1];

  return true;
}


} // end namespace olb

#endif
