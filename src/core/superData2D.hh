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
 * Dynamics for a generic 2D super data -- header file.
 */
#ifndef SUPER_DATA_2D_HH
#define SUPER_DATA_2D_HH

#include "superData2D.h"
#include "geometry/cuboid2D.h"
#include "geometry/cuboidGeometry2D.h"
#include "functors/lattice/superBaseF2D.h"

namespace olb {


template<typename T, typename BaseType>
SuperData2D<T,BaseType>::~SuperData2D()
{
  deConstruct();
}

template<typename T, typename BaseType>
SuperData2D<T,BaseType>::SuperData2D(int size) : SuperStructure2D<T> (), _size(size)
{
//  std::cout << "/// SuperData2D ctor" << std::endl;
}

template<typename T, typename BaseType>
SuperData2D<T,BaseType>::SuperData2D(CuboidGeometry2D<T>& cuboidGeometry,
                                     LoadBalancer<T>& loadBalancer, int overlap, int size)
  : SuperStructure2D<T>(cuboidGeometry, loadBalancer, overlap), _size(size)
{
//  std::cout << "/// SuperData2D ctor" << std::endl;
  allocateMemory();
}

template<typename T, typename BaseType>
SuperData2D<T,BaseType>::SuperData2D(SuperData2D<T,BaseType>& rhs)
  : SuperStructure2D<T>(rhs.getCuboidGeometry(), rhs.getLoadBalancer(), rhs.getOverlap())
{
//  std::cout << "/// SuperData2D copy ctor" << std::endl;
  // copy whole BlockData2D vector and size
  _extendedBlockData = rhs._extendedBlockData;
  _size = rhs._size;
}


template<typename T, typename BaseType>
SuperData2D<T,BaseType>::SuperData2D(SuperF2D<T,BaseType>& rhs)
  : SuperStructure2D<T>(rhs.getSuperStructure().getCuboidGeometry(),
                        rhs.getSuperStructure().getLoadBalancer(),
                        rhs.getSuperStructure().getOverlap()),
    _size(rhs.getTargetDim())
{
//  std::cout << "/// SuperData2D ctor" << std::endl;
  allocateMemory();

  int i[3];
  int iX, iY;
  for (int iCloc=0; iCloc < this->getLoadBalancer().size(); ++iCloc) {
    i[0] = this->getLoadBalancer().glob(iCloc);
    for (iX=0; iX < get(iCloc).getNx(); iX++) {
      for (iY=0; iY < get(iCloc).getNy(); iY++) {
        i[1] = iX - this->_overlap;
        i[2] = iY - this->_overlap;
        BaseType tmp[rhs.getTargetDim()];
        rhs(tmp, i);
        for (int iDim=0; iDim<rhs.getTargetDim(); iDim++) {
          get(iCloc).get(iX, iY, iDim) = (BaseType)(tmp[iDim]);
        }
      }
    }
  }
}

template<typename T, typename BaseType>
SuperData2D<T,BaseType>& SuperData2D<T,BaseType>::operator=(SuperData2D<T,BaseType>& rhs)
{
  // Swap content of rhs with local content
  swap(rhs);
  return *this;
}

template<typename T, typename BaseType>
void SuperData2D<T,BaseType>::allocateMemory()
{
  if (_extendedBlockData.size() == 0) {
    _extendedBlockData.reserve( this->getLoadBalancer().size() ); // speed up
    for (int iCloc = 0; iCloc < this->getLoadBalancer().size(); ++iCloc) {
      int iCglob = this->getLoadBalancer().glob(iCloc);
      Cuboid2D<T> extendedCuboid(this->getCuboidGeometry().get(iCglob), this->getOverlap());
      _extendedBlockData.emplace_back(extendedCuboid.getNx(),
                                      extendedCuboid.getNy(),
                                      _size);
    }
  }
}

template<typename T, typename BaseType>
void SuperData2D<T,BaseType>::swap(std::vector< BlockData2D<T,BaseType> > & blockDatas)
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
void SuperData2D<T,BaseType>::swap(SuperData2D<T,BaseType>& rhs)
{
  std::swap(_size, rhs._size);
  _extendedBlockData.swap(rhs._extendedBlockData);
}

template<typename T, typename BaseType>
bool SuperData2D<T,BaseType>::isConstructed() const
{
  // iterate through all blocks and check if isConstructed
  bool superIsConstructed = true;
  for (auto & blockData : _extendedBlockData) {
    superIsConstructed &= blockData.isConstructed();
  }
  return superIsConstructed;
}

template<typename T, typename BaseType>
void SuperData2D<T,BaseType>::deConstruct()
{
  // for all blocks: deConstruct
  for (auto & blockData : _extendedBlockData) {
    blockData.deConstruct();
  }
}

template<typename T, typename BaseType>
void SuperData2D<T,BaseType>::reset()
{
  // for all blocks: reset
  for (auto & blockData : _extendedBlockData) {
    blockData.reset();
  }
}

template<typename T, typename BaseType>
BlockData2D<T,BaseType>& SuperData2D<T,BaseType>::get(int iC)
{
  return _extendedBlockData[iC];
}

template<typename T, typename BaseType>
BlockData2D<T,BaseType> const& SuperData2D<T,BaseType>::get(int iC) const
{
  return _extendedBlockData[iC];
}

template<typename T, typename BaseType>
BaseType& SuperData2D<T,BaseType>::get(int iC, int iX, int iY, int iData)
{
  return _extendedBlockData[iC].get(iX + this->_overlap, iY + this->_overlap, iData);
}

template<typename T, typename BaseType>
BaseType const& SuperData2D<T,BaseType>::get(int iC, int iX, int iY, int iData) const
{
  return _extendedBlockData[iC].get(iX + this->_overlap, iY + this->_overlap, iData);
}

template<typename T, typename BaseType>
bool* SuperData2D<T,BaseType>::operator() (int iCloc, int iX, int iY, int iData)
{
  return (bool*)&get(iCloc,iX,iY,iData);
}

template<typename T, typename BaseType>
int SuperData2D<T,BaseType>::getDataSize() const
{
  return _size;
}

template<typename T, typename BaseType>
int SuperData2D<T,BaseType>::getDataTypeSize() const
{
  return sizeof(BaseType);
}


}  // namespace olb

#endif
