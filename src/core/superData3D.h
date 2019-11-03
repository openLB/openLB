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

#ifndef SUPER_DATA_3D_H
#define SUPER_DATA_3D_H

#include <vector>
#include "blockData3D.h"
#include "serializer.h"
#include "communication/superStructure3D.h"


namespace olb {


template< typename T> class LoadBalancer;
template< typename T> class CuboidGeometry3D;
template< typename T> class SuperStructure3D;
template<typename T, typename W> class SuperF3D;
template<typename T,typename BaseType> class BlockData3D;


/** SuperData3D orchestrates BlockData3D.
 * 1. Specific BlockData3D can be accessed through get() methods.
 * 2. Specific BlockData3D elements can be accessed through get() methods.
 *
 * \param _size               corresponds to the size of an BlockData3D element, e.g. 3 for velocity
 * \param _extendedBlockData  vector of related BlockData3D
 */
template<typename T, typename BaseType>
class SuperData3D : public SuperStructure3D<T>, public BufferSerializable {
protected:
  /// Dimension of the data field
  int _size;
  /// Vector of BlockData
  std::vector< BlockData3D<T,BaseType> > _extendedBlockData;
public:
  /// dtor
  ~SuperData3D() override;
  /// Simple constructor with size
  SuperData3D(int size=1);
  /// This constructor should be commonly used
  SuperData3D(CuboidGeometry3D<T>& cuboidGeometry, LoadBalancer<T>& loadBalancer,
              int overlap=2, int size=1);
  /// Copy Constructor
  SuperData3D(SuperData3D<T,BaseType>& rhs);
  /// Construct `SuperData3D` out of `SuperF3D`
  SuperData3D(SuperF3D<T,BaseType>& rhs);
  /// Initialize `_extendedBlockData` from `CuboidGeometry3D` and `LoadBalancer` with `overlap`
  void allocateMemory();
  /// Assignment Operator (Swaps content of rhs with local content)
  SuperData3D<T,BaseType>& operator=(SuperData3D<T,BaseType>& rhs);

  /// Swap blockDatas into member variable and read size from first block in vector
  /// (Needed for vtiReader)
  // TODO: Move to vtiReader?
  void swap(std::vector< BlockData3D<T,BaseType> >& blockDatas);
  /// Swap SuperData object including size and blockData
  /// Careful: This method does not swap the complete SuperStructure,
  ///          only the SuperData-specific data
  void swap(SuperData3D<T,BaseType>& rhs);

  bool isConstructed() const;
  void deConstruct();
  void reset();

  /// \return BlockData object at position iC
  BlockData3D<T,BaseType>& get(int iC);
  /// \return BlockData object at position iC
  BlockData3D<T,BaseType> const& get(int iC) const;

  virtual BaseType& get(int iC, int iX, int iY, int nz, int iData=0);
  virtual BaseType const& get(int iC, int iX, int iY, int iZ, int iData=0) const;

  /// Write access to the memory of the data of the super structure where (iX, iY, iZ) is the point providing the data iData in the block iCloc
  bool* operator() (int iCloc, int iX, int iY, int iZ, int iData) override;
  /// Read only access to the dim of the data of the super structure
  int getDataSize() const override;
  /// Read only access to the data type dim of the data of the super structure
  int getDataTypeSize() const override;

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Returns a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;
};


}  // namespace olb

#endif
