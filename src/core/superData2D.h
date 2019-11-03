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

#ifndef SUPER_DATA_2D_H
#define SUPER_DATA_2D_H

#include <vector>
#include "blockData2D.h"
#include "communication/superStructure2D.h"


namespace olb {

template< typename T> class LoadBalancer;
template<typename T> class CuboidGeometry;
template<typename T, typename W> class SuperF2D;
template<typename T, typename BaseType> class BlockData2D;

/** SuperData2D orchestrates BlockData2D.
 * 1. Specific BlockData2D can be accessed through get() methods.
 * 2. Specific BlockData2D elements can be accessed through get() methods.
 *
 * \param _size               corresponds to the size of an BlockData2D element, e.g. 3 for velocity
 * \param _extendedBlockData  vector of related BlockData2D
 */
template<typename T, typename BaseType>
class SuperData2D : public SuperStructure2D<T> {
protected:
  /// Dimension of the data field
  int _size;
  /// Vector of BlockData
  std::vector< BlockData2D<T,BaseType> > _extendedBlockData;
public:
  virtual ~SuperData2D();
  /// Simple constructor with size
  SuperData2D(int size=1);
  /// This constructor should be commonly used
  SuperData2D(CuboidGeometry2D<T>& cuboidGeometry, LoadBalancer<T>& loadBalancer,
              int overlap=2, int size=1);
  /// Copy Constructor
  SuperData2D(SuperData2D<T,BaseType>& rhs);
  /// Construct `SuperData2D` out of `SuperF2D`
  SuperData2D(SuperF2D<T,BaseType>& rhs);
  /// Assignment Operator (Swaps content of rhs with local content)
  SuperData2D<T,BaseType>& operator=(SuperData2D<T,BaseType>& rhs);

  /// Initialize `_extendedBlockData` from `CuboidGeometry2D` and `LoadBalancer` with `overlap`
  void allocateMemory();

  /// Swap blockDatas into member variable and read size from first block in vector
  /// (Needed for vtiReader)
  // TODO: Move to vtiReader?
  void swap(std::vector< BlockData2D<T,BaseType> >& blockDatas);
  /// Swap SuperData object including size and blockData
  /// Careful: This method does not swap the complete SuperStructure,
  ///          only the SuperData-specific data
  void swap(SuperData2D<T,BaseType>& rhs);

  bool isConstructed() const;
  void deConstruct();
  void reset();

  /// Returns BlockData object iC
  BlockData2D<T,BaseType>& get(int iC);
  /// Returns BlockData object iC
  BlockData2D<T,BaseType> const& get(int iC) const;

  virtual BaseType& get(int iC, int iX, int iY, int iData=0);
  virtual BaseType const& get(int iC, int iX, int iY, int iData=0) const;

  /// Write access to the memory of the data of the super structure where (iX, iY, iZ) is the point providing the data iData in the block iCloc
  bool* operator() (int iCloc, int iX, int iY, int iData) override;
  /// Read only access to the dim of the data of the super structure
  int getDataSize() const override;
  /// Read only access to the data type dim of the data of the super structure
  int getDataTypeSize() const override;
};


}  // namespace olb

#endif
