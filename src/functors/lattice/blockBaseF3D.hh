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

#ifndef BLOCK_BASE_F_3D_HH
#define BLOCK_BASE_F_3D_HH

#include "blockBaseF3D.h"

namespace olb {


template <typename T>
BlockF3D<T>::BlockF3D(BlockStructure3D& blockStructure, int targetDim)
  : GenericF<T,int>(targetDim,3), _blockStructure(blockStructure) { }

template <typename T>
BlockStructure3D& BlockF3D<T>::getBlockStructure() const
{
  return _blockStructure;
}

template <typename T>
std::vector<T> BlockF3D<T>::getMinValue()
{
  T min[this->getTargetDim()];
  T minTmp[this->getTargetDim()];
  this->operator()(min,0,0);
  for (int iX = 1; iX < _blockStructure.getNx(); ++iX) {
    for (int iY = 1; iY < _blockStructure.getNy(); ++iY) {
      for (int iZ = 1; iZ < _blockStructure.getNz(); ++iZ) {
        this->operator()(minTmp,iX,iY,iZ);
        for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
          if (min[iDim] > minTmp[iDim] ) {
            min[iDim] = minTmp[iDim];
          }
        }
      }
    }
  }
  std::vector<T> minV(min,min+this->getTargetDim());
  return minV;
}


template <typename T>
std::vector<T> BlockF3D<T>::getMaxValue()
{
  T max[this->getTargetDim()];
  T maxTmp[this->getTargetDim()];
  this->operator()(max,0,0);
  for (int iX = 1; iX < _blockStructure.getNx(); ++iX) {
    for (int iY = 1; iY < _blockStructure.getNy(); ++iY) {
      for (int iZ = 1; iZ < _blockStructure.getNz(); ++iZ) {
        this->operator()(maxTmp,iX,iY,iZ);
        for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
          if (max[iDim] > maxTmp[iDim] ) {
            max[iDim] = maxTmp[iDim];
          }
        }
      }
    }
  }
  std::vector<T> maxV(max,max+this->getTargetDim());
  return maxV;
}


template <typename T,typename BaseType>
BlockDataF3D<T,BaseType>::BlockDataF3D(BlockData3D<T,BaseType>& blockData)
  : BlockF3D<T>(blockData, blockData.getSize()),
    _blockData(blockData)
{ }

template <typename T,typename BaseType>
BlockDataF3D<T,BaseType>::BlockDataF3D(BlockF3D<BaseType>& f)
  : BlockF3D<T>(f.getBlockStructure(), f.getTargetDim()),
    _blockDataStorage(new BlockData3D<T,BaseType>(f)),
    _blockData(*_blockDataStorage)
{ }

template <typename T,typename BaseType>
BlockDataF3D<T,BaseType>::BlockDataF3D(int nx, int ny, int nz, int size)
// hacky solution to both managing BlockData3D using std::unique_ptr and
// passing it down the line to the base class
  : BlockF3D<T>(*(new BlockData3D<T,BaseType>(nx, ny, nz, size)), size),
    _blockDataStorage(static_cast<BlockData3D<T,BaseType>*>(&(this->getBlockStructure()))),
    _blockData(*_blockDataStorage)
{ }

template <typename T,typename BaseType>
BlockData3D<T,BaseType>& BlockDataF3D<T,BaseType>::getBlockData()
{
  return _blockData;
}

template <typename T, typename BaseType>
bool BlockDataF3D<T,BaseType>::operator() (BaseType output[], const int input[])
{
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = _blockData.get(input[0], input[1], input[2], iDim);
  }
  return true;
}


template <typename T,typename BaseType>
BlockDataViewF3D<T,BaseType>::BlockDataViewF3D(BlockData3D<T,BaseType>& blockData, int overlap)
  : BlockDataF3D<T,BaseType>(blockData),
    _overlap(overlap)
{ }

template <typename T, typename BaseType>
bool BlockDataViewF3D<T,BaseType>::operator() (BaseType output[], const int input[])
{
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = this->_blockData.get(input[0] + _overlap,
                                        input[1] + _overlap,
                                        input[2] + _overlap,
                                        iDim);
  }
  return true;
}


template <typename T>
BlockIdentity3D<T>::BlockIdentity3D(BlockF3D<T>& f)
  : BlockF3D<T>(f.getBlockStructure(),f.getTargetDim() ), _f(f)
{
  this->getName() = _f.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T>
bool BlockIdentity3D<T>::operator()(T output[], const int input[])
{
  return _f(output,input);
}


template <typename T>
BlockExtractComponentF3D<T>::BlockExtractComponentF3D(BlockF3D<T>& f, int extractDim)
  : BlockF3D<T>(f.getBlockStructure(),1 ), _f(f), _extractDim(extractDim)
{
  this->getName() = _f.getName();
}

template <typename T>
int BlockExtractComponentF3D<T>::getExtractDim()
{
  return _extractDim;
}

template <typename T>
bool BlockExtractComponentF3D<T>::operator()(T output[], const int input[])
{
  std::vector<T> outTmp(_f.getTargetDim(), T{});
  _f(outTmp.data(), input);
  output[0] = outTmp[_extractDim];
  return true;
}


template <typename T>
BlockExtractComponentIndicatorF3D<T>::BlockExtractComponentIndicatorF3D(
  BlockF3D<T>& f, int extractDim, BlockIndicatorF3D<T>& indicatorF)
  : BlockExtractComponentF3D<T>(f, extractDim),
    _indicatorF(indicatorF)
{
  this->getName() = f.getName();
}

template <typename T>
bool BlockExtractComponentIndicatorF3D<T>::operator()(T output[], const int input[])
{
  output[0] = T{};
  if (_indicatorF(input)) {
    return BlockExtractComponentF3D<T>::operator()(output, input);
  }
  return true;
}


template <typename T>
BlockExtractIndicatorF3D<T>::BlockExtractIndicatorF3D(
  BlockF3D<T>& f, BlockIndicatorF3D<T>& indicatorF)
  : BlockF3D<T>(f.getBlockStructure(), f.getTargetDim()),
    _f(f),
    _indicatorF(indicatorF)
{
  this->getName() = f.getName();
}

template <typename T>
bool BlockExtractIndicatorF3D<T>::operator()(T output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T{};
  }
  if (_indicatorF(input)) {
    _f(output, input);
  }
  return true;
}


template <typename T>
BlockDotProductF3D<T>::BlockDotProductF3D(BlockF3D<T>& f, T vector[])
  : BlockF3D<T>(f.getBlockStructure(),1 ), _f(f), _vector(vector)
{
  this->getName() = _f.getName();
  /*if ( (sizeof(_vector)/sizeof(T)) != _f.getTargetDim() ) {
    std::cout << "WARNING: dimension of vectors do not match!" << std::endl;
    exit(-1);
  }*/
}

template <typename T>
bool BlockDotProductF3D<T>::operator()(T output[], const int input[])
{
  T outTmp[3];
  output[0] = T();
  _f(outTmp, input);
  for (int iDim=0; iDim<_f.getTargetDim(); iDim++) {
    output[0] += outTmp[iDim]*_vector[iDim];
  }
  return true;
}


template <typename T, typename DESCRIPTOR>
BlockLatticeF3D<T,DESCRIPTOR>::BlockLatticeF3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockStructure, int targetDim)
  : BlockF3D<T>(blockStructure, targetDim), _blockLattice(blockStructure)
{ }
/*
template <typename T, typename DESCRIPTOR>
BlockLatticeF3D<T,DESCRIPTOR>::BlockLatticeF3D(BlockLatticeF3D<T,DESCRIPTOR> const& rhs)
  : BlockF3D<T>(rhs.getBlockStructure(), rhs.getTargetDim() ), _blockLattice(rhs.getBlockLattice())
{ }

template <typename T, typename DESCRIPTOR>
BlockLatticeF3D<T,DESCRIPTOR>& BlockLatticeF3D<T,DESCRIPTOR>::operator=(BlockLatticeF3D<T,DESCRIPTOR> const& rhs)
{
  BlockLatticeF3D<T,DESCRIPTOR> tmp(rhs);
  return tmp;
}
*/
template <typename T, typename DESCRIPTOR>
BlockLatticeStructure3D<T,DESCRIPTOR>& BlockLatticeF3D<T, DESCRIPTOR>::getBlockLattice()
{
  return _blockLattice;
}


template <typename T, typename DESCRIPTOR>
BlockLatticeIdentity3D<T,DESCRIPTOR>::BlockLatticeIdentity3D(
  BlockLatticeF3D<T,DESCRIPTOR>& f)
  : BlockLatticeF3D<T,DESCRIPTOR>(f.getBlockLattice(),f.getTargetDim()),
    _f(f)
{
  this->getName() = _f.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeIdentity3D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  return _f(output,input);
}


template <typename T, typename DESCRIPTOR>
BlockLatticePhysF3D<T,DESCRIPTOR>::BlockLatticePhysF3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, targetDim), _converter(converter)
{ }

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticeThermalPhysF3D
(BlockLatticeStructure3D<T,TDESCRIPTOR>& blockLattice, const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter, int targetDim)
  : BlockLatticeF3D<T,TDESCRIPTOR>(blockLattice, targetDim), _converter(converter)
{ }




} // end namespace olb

#endif
