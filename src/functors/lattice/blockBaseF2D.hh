/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Albert Mink, Lukas Baron, Mathias J. Krause
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

#ifndef BLOCK_BASE_F_2D_HH
#define BLOCK_BASE_F_2D_HH

#include "blockBaseF2D.h"

namespace olb {


// BlockF2D
template <typename T>
BlockF2D<T>::BlockF2D(BlockStructure2D& blockStructure, int targetDim)
  : GenericF<T,int>(targetDim,2), _blockStructure(&blockStructure) { }

template <typename T>
BlockF2D<T>::BlockF2D(int targetDim)
  : GenericF<T,int>(targetDim,2), _blockStructure(nullptr) { }

template <typename T>
BlockStructure2D& BlockF2D<T>::getBlockStructure()
{
  return *_blockStructure;
}

template <typename T>
void BlockF2D<T>::setBlockStructure(BlockStructure2D* blockStructure)
{
  _blockStructure = blockStructure;
}

//template <typename T>
//std::vector<T> BlockF2D<T>::getMinValue()
//{
//  T min[this->getTargetDim()];
//  T minTmp[this->getTargetDim()];
//  this->operator()(min,0,0);
//  for (int iX = 1; iX < _blockStructure.getNx(); ++iX) {
//    for (int iY = 1; iY < _blockStructure.getNy(); ++iY) {
//      this->operator()(minTmp,iX,iY);
//      for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
//        if (min[iDim] > minTmp[iDim] ) {
//          min[iDim] = minTmp[iDim];
//        }
//      }
//    }
//  }
//  std::vector<T> minV(min,min+this->getTargetDim());
//  return minV;
//}


//template <typename T>
//std::vector<T> BlockF2D<T>::getMaxValue()
//{
//  T max[this->getTargetDim()];
//  T maxTmp[this->getTargetDim()];
//  this->operator()(max,0,0);
//  for (int iX = 1; iX < _blockStructure.getNx(); ++iX) {
//    for (int iY = 1; iY < _blockStructure.getNy(); ++iY) {
//      this->operator()(maxTmp,iX,iY);
//      for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
//        if (max[iDim] > maxTmp[iDim] ) {
//          max[iDim] = maxTmp[iDim];
//        }
//      }
//    }
//  }
//  std::vector<T> maxV(max,max+this->getTargetDim());
//  return maxV;
//}


template <typename T,typename BaseType>
BlockDataF2D<T,BaseType>::BlockDataF2D(BlockData2D<T,BaseType>& blockData)
  : BlockF2D<T>(blockData, blockData.getSize()),
    _blockData(blockData)
{ }

template <typename T,typename BaseType>
BlockDataF2D<T,BaseType>::BlockDataF2D(BlockF2D<BaseType>& f)
  : BlockF2D<T>(f.getBlockStructure(), f.getTargetDim()),
    _blockDataStorage(new BlockData2D<T,BaseType>(f)),
    _blockData(*_blockDataStorage)
{ }

template <typename T,typename BaseType>
BlockDataF2D<T,BaseType>::BlockDataF2D(int nx, int ny, int size)
  // hacky solution to both managing BlockData2D using std::unique_ptr and
  // passing it down the line to the base class
  : BlockF2D<T>(*(new BlockData2D<T,BaseType>(nx, ny, size)), size),
    _blockDataStorage(static_cast<BlockData2D<T,BaseType>*>(&(this->getBlockStructure()))),
    _blockData(*_blockDataStorage)
{ }

template <typename T,typename BaseType>
BlockData2D<T,BaseType>& BlockDataF2D<T,BaseType>::getBlockData()
{
  return _blockData;
}

template <typename T, typename BaseType>
bool BlockDataF2D<T,BaseType>::operator()(T output[], const int input[])
{
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = static_cast<T>(_blockData.get(input[0], input[1], iDim));
  }
  return true;
}


template <typename T,typename BaseType>
BlockDataViewF2D<T,BaseType>::BlockDataViewF2D(BlockData2D<T,BaseType>& blockData, int overlap)
  : BlockDataF2D<T,BaseType>(blockData),
    _overlap(overlap)
{ }

template <typename T, typename BaseType>
bool BlockDataViewF2D<T,BaseType>::operator() (T output[], const int input[])
{
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = this->_blockData.get(input[0] + _overlap,
                                        input[1] + _overlap,
                                        iDim);
  }
  return true;
}


// BlockIdendity2D
template <typename T>
BlockIdentity2D<T>::BlockIdentity2D(BlockF2D<T>& f)
  : BlockF2D<T>(f.getBlockStructure() ,f.getTargetDim() ), _f(f)
{
  this->getName() = _f.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T>
bool BlockIdentity2D<T>::operator()(T output[], const int input[])
{
  _f(output,input);
  return true;
}


// BlockLatticeF2D
template <typename T, typename DESCRIPTOR>
BlockLatticeF2D<T,DESCRIPTOR>::BlockLatticeF2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockStructure, int targetDim)
  : BlockF2D<T>(blockStructure, targetDim), _blockLattice(blockStructure) { }

template <typename T, typename DESCRIPTOR>
BlockLatticeStructure2D<T,DESCRIPTOR>& BlockLatticeF2D<T,DESCRIPTOR>::getBlockLattice()
{
  return _blockLattice;
}


template <typename T, typename DESCRIPTOR>
BlockLatticePhysF2D<T,DESCRIPTOR>::BlockLatticePhysF2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice, targetDim), _converter(converter)
{ }

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticeThermalPhysF2D
(BlockLatticeStructure2D<T,TDESCRIPTOR>& blockLattice, const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter, int targetDim)
  : BlockLatticeF2D<T,TDESCRIPTOR>(blockLattice, targetDim), _converter(converter)
{ }

} // end namespace olb

#endif
