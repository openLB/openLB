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
 * Dynamics for a generic 3D block data -- header file.
 */
#ifndef BLOCK_DATA_3D_HH
#define BLOCK_DATA_3D_HH

#include <algorithm>
#include "olbDebug.h"
#include "blockData3D.h"
#include "geometry/cuboid3D.h"
#include "functors/lattice/blockBaseF3D.h"

namespace olb {


template<typename T, typename BaseType>
BlockData3D<T,BaseType>::BlockData3D()
  : BlockStructure3D(0,0,0), _size(0), _rawData(nullptr), _field(nullptr)
{
  construct();
}

template<typename T, typename BaseType>
BlockData3D<T,BaseType>::BlockData3D(Cuboid3D<T>& cuboid, int size)
  : BlockStructure3D(cuboid.getNx(), cuboid.getNy(), cuboid.getNz()),
    _size(size), _rawData(nullptr), _field(nullptr)
{
  construct();
}

template<typename T, typename BaseType>
BlockData3D<T,BaseType>::BlockData3D(int nx, int ny, int nz, int size)
  : BlockStructure3D(nx, ny, nz), _size(size), _rawData(nullptr), _field(nullptr)
{
  construct();
}

template<typename T, typename BaseType>
BlockData3D<T,BaseType>::~BlockData3D()
{
  deConstruct();
}

template<typename T, typename BaseType>
BlockData3D<T,BaseType>::BlockData3D(BlockF3D<BaseType>& rhs)
  : BlockStructure3D(rhs.getBlockStructure().getNx(),
                     rhs.getBlockStructure().getNy(),
                     rhs.getBlockStructure().getNz()),
    _size(rhs.getTargetDim()), _rawData(nullptr), _field(nullptr)
{
  construct();
  int i[3];
  for (i[0] = 0; i[0] < this->_nx; ++i[0]) {
    for (i[1] = 0; i[1] < this->_ny; ++i[1]) {
      for (i[2] = 0; i[2] < this->_nz; ++i[2]) {
        rhs(_field[i[0]][i[1]][i[2]], i);
      }
    }
  }
}

template<typename T, typename BaseType>
BlockData3D<T,BaseType>::BlockData3D(BlockData3D<T,BaseType> const& rhs)
  : BlockStructure3D(rhs._nx, rhs._ny, rhs._nz), _size(rhs._size), _rawData(nullptr), _field(nullptr)
{
  if (rhs.isConstructed()) {
    construct();
//    for (size_t iData = 0; iData < getDataSize(); ++iData) {
//      (*this)[iData] = rhs[iData];
//    }
    std::copy( rhs._rawData, rhs._rawData + getDataSize(), _rawData );
  }
}

template<typename T, typename BaseType>
BlockData3D<T,BaseType>& BlockData3D<T,BaseType>::operator=(BlockData3D<T,BaseType> const& rhs)
{
  BlockData3D<T,BaseType> tmp(rhs);
  swap(tmp);
  return *this;
}

template<typename T, typename BaseType>
BlockData3D<T,BaseType>& BlockData3D<T,BaseType>::operator=(BlockData3D<T,BaseType>&& rhs)
{
  std::cout << "/// Move Operator BlockData3D" << std::endl;
  if (this != &rhs) {
  }
//  this->releaseMemory();  // free data of object this

  _size = rhs._size;      // swap object data
  _rawData = rhs._rawData;
  _field = rhs._field;
  this->_nx = rhs._nx;
  this->_ny = rhs._ny;
  this->_nz = rhs._nz;

  rhs._rawData = nullptr; // free data of object rhs
  rhs._field = nullptr;
  rhs._nx = 0;
  rhs._ny = 0;
  rhs._nz = 0;
  return *this;
}

template<typename T, typename BaseType>
BlockData3D<T,BaseType>::BlockData3D(BlockData3D<T,BaseType>&& rhs)
  : BlockStructure3D(rhs._nx, rhs._ny, rhs._nz), _size(rhs._size), _rawData(nullptr), _field(nullptr)
{
  std::cout << "/// Move Ctor BlockData3D" << std::endl;
  *this = std::move(rhs); // https://msdn.microsoft.com/de-de/library/dd293665.aspx
}


template<typename T, typename BaseType>
bool BlockData3D<T,BaseType>::isConstructed() const
{
  return _rawData;
}

template<typename T, typename BaseType>
void BlockData3D<T,BaseType>::construct()
{
  if (!isConstructed()) {
    allocateMemory();
  }
}

template<typename T, typename BaseType>
void BlockData3D<T,BaseType>::deConstruct()
{
  if (isConstructed()) {
    releaseMemory();
  }
}

template<typename T, typename BaseType>
void BlockData3D<T,BaseType>::reset()
{
  OLB_PRECONDITION(isConstructed());
  for (size_t index = 0; index < getDataSize(); ++index) {
    (*this)[index] = BaseType();
  }
}

template<typename T, typename BaseType>
void BlockData3D<T,BaseType>::swap(BlockData3D<T,BaseType>& rhs)
{
  // Block3D
  std::swap(this->_nx, rhs._nx);
  std::swap(this->_ny, rhs._ny);
  std::swap(this->_nz, rhs._nz);
  // BlockData3D
  std::swap(_size, rhs._size);
  std::swap(_rawData, rhs._rawData);
  std::swap(_field, rhs._field);
}

template<typename T, typename BaseType>
void BlockData3D<T,BaseType>::allocateMemory()
{
  // The conversions to size_t ensure 64-bit compatibility. Note that
  //   nx and ny are of type int, which might be 32-bit types, even on
  //   64-bit platforms. Therefore, nx*ny may lead to a type overflow.
  _rawData = new BaseType[ getDataSize() ];
  _field   = new BaseType*** [(size_t)(this->_nx)];
  for (int iX = 0; iX < this->_nx; ++iX) {
    _field[iX] = new BaseType** [(size_t)this->_ny];
    for (int iY = 0; iY < this->_ny; ++iY) {
      _field[iX][iY] = new BaseType* [(size_t)this->_nz];
      for (int iZ = 0; iZ < this->_nz; ++iZ) {
        _field[iX][iY][iZ] = _rawData + _size*( (size_t)iZ + (size_t)(this->_nz)*((size_t)iY + (size_t)(this->_ny)*(size_t)iX) );
        for (int iDim = 0; iDim < _size; ++iDim) {
          _field[iX][iY][iZ][iDim] = BaseType();
        }
      }
    }
  }
}

template<typename T, typename BaseType>
void BlockData3D<T,BaseType>::releaseMemory()
{
  delete [] _rawData;
  _rawData = nullptr;
  for (int iX = 0; iX < this->_nx; ++iX) {
    for (int iY = 0; iY < this->_ny; ++iY) {
      delete [] _field[iX][iY];
    }
    delete [] _field[iX];
  }
  delete [] _field;
}

template<typename T, typename BaseType>
BaseType& BlockData3D<T,BaseType>::operator[] (int ind)
{
  OLB_PRECONDITION(ind >= 0 && ind < this->_nx * this->_ny * this->_nz * this->_size);
  OLB_PRECONDITION(isConstructed());
  return _rawData[ind];
}

template<typename T, typename BaseType>
BaseType const& BlockData3D<T,BaseType>::operator[] (int ind) const
{
  OLB_PRECONDITION(ind >= 0 && ind < this->_nx * this->_ny * this->_nz * this->_size);
  OLB_PRECONDITION(isConstructed());
  return _rawData[ind];
}

template<typename T, typename BaseType>
bool* BlockData3D<T,BaseType>::operator() (int iX, int iY, int iZ, int iData)
{
  return (bool*)&_field[iX][iY][iZ][iData];
}

template<typename T, typename BaseType>
bool BlockData3D<T,BaseType>::operator() (T output[], const int input[])
{
  if ( input[0] >= 0 && input[1] >= 0 && input[2] >= 0 && input[0] < this->_nx && input[1] < this->_ny && input[2] < this->_nz ) {
    for (int i=0; i < _size; i++)
      output[i] = _field[input[0]][input[1]][input[2]][i];
    return true;
  } else {
    return false;
  }
}

template<typename T, typename BaseType>
BaseType& BlockData3D<T,BaseType>::get(int iX, int iY, int iZ, int iSize)
{
  OLB_PRECONDITION(iX >= 0 && iX < this->_nx);
  OLB_PRECONDITION(iY >= 0 && iY < this->_ny);
  OLB_PRECONDITION(iZ >= 0 && iZ < this->_nz);
  OLB_PRECONDITION(iSize >= 0 && iSize < _size);
  OLB_PRECONDITION(isConstructed());
  return _field[iX][iY][iZ][iSize];
}

template<typename T, typename BaseType>
BaseType const& BlockData3D<T,BaseType>::get(int iX, int iY, int iZ, int iSize) const
{
  OLB_PRECONDITION(iX >= 0 && iX < this->_nx);
  OLB_PRECONDITION(iY >= 0 && iY < this->_ny);
  OLB_PRECONDITION(iZ >= 0 && iZ < this->_nz);
  OLB_PRECONDITION(iSize >= 0 && iSize < _size);
  OLB_PRECONDITION(isConstructed());
  return _field[iX][iY][iZ][iSize];
}

template<typename T, typename BaseType>
BaseType BlockData3D<T,BaseType>::getMax()
{
  return ****std::max_element( _field, _field + getDataSize() ) ;
}

template<typename T, typename BaseType>
BaseType BlockData3D<T,BaseType>::getMin()
{
  return ****std::min_element( _field, _field + getDataSize() ) ;
}

template<typename T, typename BaseType>
BaseType* BlockData3D<T,BaseType>::getRawData() const
{
  return _rawData;
}

template<typename T, typename BaseType>
size_t BlockData3D<T,BaseType>::getDataSize() const
{
  return (size_t)(this->_nx) * (size_t)(this->_ny) * (size_t)(this->_nz) * (size_t)(_size);
}

template<typename T, typename BaseType>
int BlockData3D<T,BaseType>::getSize() const
{
  return _size;
}


template<typename T, typename BaseType>
size_t BlockData3D<T,BaseType>::getSerializableSize() const
{
  return 4 * sizeof(int) // _size, _nX/Y/Z
         + getDataSize() * sizeof(BaseType); // _rawData
};

template<typename T, typename BaseType>
bool* BlockData3D<T,BaseType>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _size);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_nx);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_ny);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_nz);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *_rawData, getDataSize());

  return dataPtr;
}

}  // namespace olb

#endif
