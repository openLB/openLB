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
#ifndef BLOCK_DATA_3D_H
#define BLOCK_DATA_3D_H


#include <stdio.h>
#include "blockStructure3D.h"
#include "serializer.h"


namespace olb {

// forward declaration is sufficient at this state, include is moved to .hh
template<typename BaseType> class BlockF3D;
template<typename T> class Cuboid3D;


/** A highly generic data structure class.
 *
 * Stored data is of type BaseType
 *
 * \param _size     data element size
 * \param _rawData  classic array representation of all elements
 * \param _filed    a matrix like representation of _rawData
 */
template<typename T,typename BaseType>
class BlockData3D : public BlockStructure3D, public Serializable {
protected:
  /// dimension of data element, vector, scalar, ...
  int _size;
  /// holds data as a 1D vector
  BaseType *_rawData;
  /// Pointer structure to [3D data]x[1D element]
  BaseType ****_field;
public:
  virtual ~BlockData3D();
  /// Construct empty cuboid
  BlockData3D();
  /// Construct from cuboid
  BlockData3D(Cuboid3D<T>& cuboid, int size=1);
  /// Construct from X-Y-Z node count
  BlockData3D(int nx, int ny, int nz, int size=1);
  /// Construct from Block Functor
  BlockData3D(BlockF3D<BaseType>& rhs);
  /// Copy Constructor
  BlockData3D(BlockData3D<T,BaseType> const& rhs);
  /// Assignment Operator
  BlockData3D<T,BaseType>& operator=(BlockData3D<T,BaseType> const& rhs);
  /// Move Operator
  BlockData3D<T,BaseType>& operator=(BlockData3D<T,BaseType>&& rhs);
  /// Move Constructor
  BlockData3D<T,BaseType>(BlockData3D<T,BaseType>&& rhs);
  /// Swap rhs Data into local fields
  void swap(BlockData3D<T,BaseType>& rhs);
  /// Memory Management
  bool isConstructed() const;
  void construct();
  void deConstruct();
  void reset();
  /// \return dataElement _rawData[ind], read and write
  BaseType& operator[] (int ind);
  /// \return dataElement _rawData[ind], read only
  BaseType const& operator[] (int ind) const;
  /// Write access to the memory of the data of the block data where (iX, iY, iZ) is the point providing the data iData
  bool* operator() (int iX, int iY, int iZ, int iData);
  bool operator() (T output[], const int input[]);
  /// read and write access to data element [iX][iY][iZ][iSize]
  BaseType& get(int iX, int iY, int iZ, int iSize=0);
  /// read only access to data element [iX][iY][iZ][iSize]
  BaseType const& get(int iX, int iY, int iZ, int iSize=0) const;
  /// \return max of data, for vector valued data it determines the max component
  BaseType getMax();
  /// \return min of data, for vector valued data it determines the max component
  BaseType getMin();
  /// \return _rawData array
  BaseType* getRawData() const;
  /// Number of all variables in the data field
  virtual size_t getDataSize() const;
  /// \return _size, the dimension of an data element
  int getSize() const;
  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override
  {
    return 5;
  };
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Returns a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;
private:
  /// Memory Management
  void allocateMemory();
  void releaseMemory();
};

}  // namespace olb

#endif
