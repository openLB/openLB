/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause
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
 * Dynamics for a generic 2D block data -- header file.
 */
#ifndef BLOCK_DATA_2D_HH
#define BLOCK_DATA_2D_HH

#include <algorithm>
#include "olbDebug.h"
#include "blockData2D.h"
#include "geometry/cuboid2D.h"
#include "functors/lattice/blockBaseF2D.h"


namespace olb {


template<typename T, typename BaseType>
BlockData2D<T,BaseType>::BlockData2D() : BlockStructure2D(0,0), _size(0), _rawData(nullptr), _field(nullptr)
{
  construct();
}

template<typename T, typename BaseType>
BlockData2D<T,BaseType>::BlockData2D(Cuboid2D<T>& cuboid, int size)
  : BlockStructure2D(cuboid.getNx(), cuboid.getNy()), _size(size), _rawData(nullptr), _field(nullptr)
{
  construct();
}

template<typename T, typename BaseType>
BlockData2D<T,BaseType>::BlockData2D(int nx, int ny, int size)
  : BlockStructure2D(nx, ny), _size(size), _rawData(nullptr), _field(nullptr)
{
  construct();
}

template<typename T, typename BaseType>
BlockData2D<T,BaseType>::~BlockData2D()
{
  deConstruct();
}

template<typename T, typename BaseType>
BlockData2D<T,BaseType>::BlockData2D(BlockF2D<BaseType>& rhs)
  : BlockStructure2D(rhs.getBlockStructure().getNx(), rhs.getBlockStructure().getNy()),
    _size(rhs.getTargetDim())
{
  construct();
  int i[2];
  for (i[0] = 0; i[0] < this->_nx; ++i[0]) {
    for (i[1] = 0; i[1] < this->_ny; ++i[1]) {
      rhs(_field[i[0]][i[1]], i);
    }
  }
}

template<typename T, typename BaseType>
BlockData2D<T,BaseType>::BlockData2D(BlockData2D<T,BaseType> const& rhs)
  : BlockStructure2D(rhs._nx, rhs._ny), _size(rhs._size), _rawData(nullptr), _field(nullptr)
{
  if (rhs.isConstructed()) {
    construct();
    std::copy( rhs._rawData, rhs._rawData + getDataSize(), _rawData );
  }
}

template<typename T, typename BaseType>
BlockData2D<T,BaseType>& BlockData2D<T,BaseType>::operator=(BlockData2D<T,BaseType> const& rhs)
{
  BlockData2D<T,BaseType> tmp(rhs);
  swap(tmp);
  return *this;
}

// benefits of move operator: does not allocate memory or copys objects
template<typename T, typename BaseType>
BlockData2D<T,BaseType>& BlockData2D<T,BaseType>::operator=(BlockData2D<T,BaseType>&& rhs)
{
  if (this == &rhs) {
      return *this;
  }
//  this->releaseMemory();  // free data of object this

  _size = rhs._size;      // swap object data
  _rawData = rhs._rawData;
  _field = rhs._field;
  this->_nx = rhs._nx;
  this->_ny = rhs._ny;

  rhs._rawData = nullptr; // free data of object rhs
  rhs._field = nullptr;
  rhs._nx = 0;
  rhs._ny = 0;

  return *this;
}

// benefits of move operator: does not allocate memory
template<typename T, typename BaseType>
BlockData2D<T,BaseType>::BlockData2D(BlockData2D<T,BaseType>&& rhs)
  : BlockStructure2D(rhs._nx, rhs._ny), _size(0), _rawData(nullptr), _field(nullptr)
{
  *this = std::move(rhs); // https://msdn.microsoft.com/de-de/library/dd293665.aspx
}


template<typename T, typename BaseType>
bool BlockData2D<T,BaseType>::isConstructed() const
{
  return _rawData;
}

template<typename T, typename BaseType>
void BlockData2D<T,BaseType>::construct()
{
  if (!isConstructed()) {
    allocateMemory();
  }
}

template<typename T, typename BaseType>
void BlockData2D<T,BaseType>::deConstruct()
{
  if (isConstructed()) {
    releaseMemory();
  }
}

template<typename T, typename BaseType>
void BlockData2D<T,BaseType>::reset()
{
  OLB_PRECONDITION(isConstructed());
  for (size_t index = 0; index < getDataSize(); ++index) {
    (*this)[index] = BaseType();
  }
}

template<typename T, typename BaseType>
void BlockData2D<T,BaseType>::swap(BlockData2D<T,BaseType>& rhs)
{
  // Block2D
  std::swap(this->_nx, rhs._nx);
  std::swap(this->_ny, rhs._ny);
  // BlockData2D
  std::swap(_size, rhs._size);
  std::swap(_rawData, rhs._rawData);
  std::swap(_field, rhs._field);
}

template<typename T, typename BaseType>
void BlockData2D<T,BaseType>::allocateMemory()
{
  // The conversions to size_t ensure 64-bit compatibility. Note that
  //   nx and ny are of type int, which might by 32-bit types, even on
  //   64-bit platforms. Therefore, nx*ny may lead to a type overflow.
  _rawData = new BaseType[ getDataSize() ];
  _field   = new BaseType** [(size_t)(this->_nx)];
  for (int iX = 0; iX < this->_nx; ++iX) {
    _field[iX] = new BaseType* [(size_t)this->_ny];
    for (int iY = 0; iY < this->_ny; ++iY) {
      // connect matrix element to the corresponding array entry of _rawData
      _field[iX][iY] = _rawData + _size*( (size_t)iY + (size_t)(this->_ny)*(size_t)iX );
      for (int iDim = 0; iDim < _size; ++iDim) {
        // initialize data with zero
        _field[iX][iY][iDim] = BaseType();
      }
    }
  }
}

template<typename T, typename BaseType>
void BlockData2D<T,BaseType>::releaseMemory()
{
  delete [] _rawData;
  _rawData = nullptr;
  for (unsigned int iX = 0; iX < (size_t)(this->_nx); ++iX) {
    delete [] _field[iX];
  }
  delete [] _field;
}

template<typename T, typename BaseType>
BaseType& BlockData2D<T,BaseType>::get(int iX, int iY, int iSize)
{
  OLB_PRECONDITION(iX >= 0 && iX < this->_nx);
  OLB_PRECONDITION(iY >= 0 && iY < this->_ny);
  OLB_PRECONDITION(iSize >= 0 && iSize < _size);
  OLB_PRECONDITION(isConstructed());
  return _field[iX][iY][iSize];
}

template<typename T, typename BaseType>
BaseType const& BlockData2D<T,BaseType>::get(int iX, int iY, int iSize) const
{
  OLB_PRECONDITION(iX >= 0 && iX < this->_nx);
  OLB_PRECONDITION(iY >= 0 && iY < this->_ny);
  OLB_PRECONDITION(iSize >= 0 && iSize < _size);
  OLB_PRECONDITION(isConstructed());
  return _field[iX][iY][iSize];
}

template<typename T, typename BaseType>
BaseType& BlockData2D<T,BaseType>::operator[] (int ind)
{
  OLB_PRECONDITION(ind >= 0 && ind < this->_nx * this->_ny * this->_size);
  OLB_PRECONDITION(isConstructed());
  return _rawData[ind];
}

template<typename T, typename BaseType>
BaseType const& BlockData2D<T,BaseType>::operator[] (int ind) const
{
  OLB_PRECONDITION(ind >= 0 && ind < this->_nx * this->_ny * this->_size);
  OLB_PRECONDITION(isConstructed());
  return _rawData[ind];
}

template<typename T, typename BaseType>
bool* BlockData2D<T,BaseType>::operator() (int iX, int iY, int iData)
{
  return (bool*)&_field[iX][iY][iData];
}

template<typename T, typename BaseType>
BaseType BlockData2D<T,BaseType>::getMax()
{
  return ***std::max_element( _field, _field + getDataSize() ) ;
}

template<typename T, typename BaseType>
BaseType BlockData2D<T,BaseType>::getMin()
{
  return ***std::min_element( _field, _field + getDataSize() ) ;
}

template<typename T, typename BaseType>
BaseType* BlockData2D<T,BaseType>::getRawData() const
{
  return _rawData;
}

template<typename T, typename BaseType>
size_t BlockData2D<T,BaseType>::getDataSize() const
{
  return (size_t)(this->_nx) * (size_t)(this->_ny) * (size_t)(_size);
}

template<typename T, typename BaseType>
int BlockData2D<T,BaseType>::getSize() const
{
  return _size;
}

template<typename T, typename BaseType>
size_t BlockData2D<T,BaseType>::getSerializableSize() const
{
  return   3 * sizeof(int) // _size, _nX/Y
           + getDataSize() * sizeof(BaseType); // _rawData
};

template<typename T, typename BaseType>
bool* BlockData2D<T,BaseType>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _size);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_nx);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_ny);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *_rawData, getDataSize());

  return dataPtr;
}

}  // namespace olb

#endif
