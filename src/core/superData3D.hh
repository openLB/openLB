/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Benjamin Förster
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
 * Dynamics for a generic 3D super data -- header file.
 */
#ifndef SUPER_DATA_3D_HH
#define SUPER_DATA_3D_HH

#include <numeric>
#include "superData3D.h"
#include "functors/lattice/superBaseF3D.h"
#include "geometry/cuboid3D.h"

namespace olb {

template<typename T, typename BaseType>
SuperData3D<T,BaseType>::~SuperData3D()
{
  deConstruct();
}

template<typename T, typename BaseType>
SuperData3D<T,BaseType>::SuperData3D(int size)
  : SuperStructure3D<T> (),
    _size(size)
{
}

template<typename T, typename BaseType>
SuperData3D<T,BaseType>::SuperData3D(CuboidGeometry3D<T>& cuboidGeometry,
                                     LoadBalancer<T>& loadBalancer, int overlap, int size)
  : SuperStructure3D<T>(cuboidGeometry, loadBalancer, overlap), _size(size)
{
  allocateMemory();
}

template<typename T, typename BaseType>
SuperData3D<T,BaseType>::SuperData3D(SuperData3D<T,BaseType>& rhs)
  : SuperStructure3D<T>(rhs.getCuboidGeometry(), rhs.getLoadBalancer(), rhs.getOverlap())
{
  // copy whole BlockData3D vector and size
  _extendedBlockData = rhs._extendedBlockData;
  _size = rhs._size;
}

template<typename T, typename BaseType>
SuperData3D<T,BaseType>::SuperData3D(SuperF3D<T,BaseType>& rhs)
  : SuperStructure3D<T>(rhs.getSuperStructure().getCuboidGeometry(),
                        rhs.getSuperStructure().getLoadBalancer(),
                        rhs.getSuperStructure().getOverlap()),
    _size(rhs.getTargetDim())
{
  allocateMemory();
  int i[4];
  int iX, iY, iZ;
  for (int iCloc=0; iCloc < this->getLoadBalancer().size(); ++iCloc) {
    i[0] = this->getLoadBalancer().glob(iCloc);
    for (iX=0; iX < get(iCloc).getNx(); iX++) {
      for (iY=0; iY < get(iCloc).getNy(); iY++) {
        for (iZ=0; iZ < get(iCloc).getNz(); iZ++) {
          i[1] = iX - this->_overlap;
          i[2] = iY - this->_overlap;
          i[3] = iZ - this->_overlap;
          BaseType tmp[rhs.getTargetDim()];
          rhs(tmp, i);
          for (int iDim=0; iDim<rhs.getTargetDim(); iDim++ ) {
            get(iCloc).get(iX, iY, iZ, iDim) = tmp[iDim];
          }
        }
      }
    }
  }
}

template<typename T, typename BaseType>
void SuperData3D<T,BaseType>::allocateMemory()
{
  for (int iCloc=0; iCloc < this->getLoadBalancer().size(); ++iCloc) {
    int iCglob = this->getLoadBalancer().glob(iCloc);
    Cuboid3D<T> extendedCuboid(this->getCuboidGeometry().get(iCglob), this->getOverlap());
    _extendedBlockData.emplace_back(extendedCuboid.getNx(),
                                    extendedCuboid.getNy(),
                                    extendedCuboid.getNz(),
                                    _size);
  }
};

template<typename T, typename BaseType>
SuperData3D<T,BaseType>& SuperData3D<T,BaseType>::operator=(SuperData3D<T,BaseType>& rhs)
{
  // Swap content of rhs with local content
  swap(rhs);
  return *this;
}

template<typename T, typename BaseType>
void SuperData3D<T,BaseType>::swap(std::vector< BlockData3D<T,BaseType> >& blockDatas)
{
  // Needed for vtiReader: After building geometry, vtiReader builds Block Data vector
  if ( blockDatas.size() == 0 ) {
    _size = 1;
  }
  else {
    _size = blockDatas[0].getSize();
  }
  std::swap(_extendedBlockData, blockDatas);
}

template<typename T, typename BaseType>
void SuperData3D<T,BaseType>::swap(SuperData3D<T,BaseType>& rhs)
{
  std::swap(_size, rhs._size);
  _extendedBlockData.swap(rhs._extendedBlockData);
}

template<typename T, typename BaseType>
bool SuperData3D<T,BaseType>::isConstructed() const
{
  // iterate through all blocks and check if isConstructed
  bool superIsConstructed = true;
  for (auto & blockData : _extendedBlockData) {
    superIsConstructed &= blockData.isConstructed();
  }
  return superIsConstructed;
}

template<typename T, typename BaseType>
void SuperData3D<T,BaseType>::deConstruct()
{
  // for all blocks: deConstruct
  for (auto & blockData : _extendedBlockData) {
    blockData.deConstruct();
  }
}

template<typename T, typename BaseType>
void SuperData3D<T,BaseType>::reset()
{
  // for all blocks: reset
  for (auto & blockData : _extendedBlockData) {
    blockData.reset();
  }
}

template<typename T, typename BaseType>
BlockData3D<T,BaseType>& SuperData3D<T,BaseType>::get(int iC)
{
  return _extendedBlockData[iC];
}

template<typename T, typename BaseType>
BlockData3D<T,BaseType> const& SuperData3D<T,BaseType>::get(int iC) const
{
  return _extendedBlockData[iC];
}

template<typename T, typename BaseType>
BaseType& SuperData3D<T,BaseType>::get(int iC, int iX, int iY, int iZ, int iData)
{
  return _extendedBlockData[iC].get(iX + this->_overlap,
                                    iY + this->_overlap,
                                    iZ + this->_overlap, iData);
}

template<typename T, typename BaseType>
BaseType const& SuperData3D<T,BaseType>::get(int iC, int iX, int iY, int iZ, int iData) const
{
  return _extendedBlockData[iC].get(iX + this->_overlap,
                                    iY + this->_overlap,
                                    iZ + this->_overlap, iData);
}

template<typename T, typename BaseType>
bool* SuperData3D<T,BaseType>::operator() (int iCloc, int iX, int iY, int iZ, int iData)
{
  return (bool*)&get(iCloc,iX,iY,iZ,iData);
}

template<typename T, typename BaseType>
int SuperData3D<T,BaseType>::getDataSize() const
{
  return _size;
}

template<typename T, typename BaseType>
int SuperData3D<T,BaseType>::getDataTypeSize() const
{
  return sizeof(BaseType);
}


template<typename T, typename BaseType>
size_t SuperData3D<T,BaseType>::getNblock() const
{
  return 1 + std::accumulate(_extendedBlockData.begin(), _extendedBlockData.end(), size_t(0),
                             Serializable::sumNblock());
}

template<typename T, typename BaseType>
size_t SuperData3D<T,BaseType>::getSerializableSize() const
{
  return   sizeof(int) // _size
           + std::accumulate(_extendedBlockData.begin(), _extendedBlockData.end(), size_t(0),
                             Serializable::sumSerializableSize());
};

template<typename T, typename BaseType>
bool* SuperData3D<T,BaseType>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  size_t sizeBufferIndex = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _size);
  registerStdVectorOfSerializablesOfConstSize(iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr,
      _extendedBlockData, loadingMode);
  return dataPtr;
}


}  // namespace olb

#endif
