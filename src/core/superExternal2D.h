/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Robin Trunk
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
 * The description of a 2D super external field -- header file.
 */


#ifndef SUPER_EXTERNAL_2D_H
#define SUPER_EXTERNAL_2D_H

#include <vector>
#include "superLattice2D.h"
#include "communication/superStructure2D.h"

/// All OpenLB code is contained in this namespace.
namespace olb {


template<typename T> class Communicator2D;
template<typename T> class SuperGeometry2D;
template<typename T, typename DESCRIPTOR> class SuperLattice2D;
template<typename T> class SuperStructure2D;

/// A super external field is needed to communicate values of the external field
template<typename T, typename DESCRIPTOR, typename FIELD>
class SuperExternal2D : public SuperStructure2D<T> {
private:
  int _overlap;
  SuperLattice2D<T, DESCRIPTOR>& _sLattice;
public:
  /// Construction of a super external field
  SuperExternal2D(SuperGeometry2D<T>& superGeometry,
                  SuperLattice2D<T,DESCRIPTOR>& sLattice,
                  int overlap);
  void communicate(bool verbose=true) override;
  /// Write access to the memory of the data of the super structure
  bool* operator() (int iCloc, int iX, int iY, int iData) override;
  /// Read only access to the dim of the data of the super structure
  int getDataSize() const override;
  /// Read only access to the data type dim of the data of the super structure
  int getDataTypeSize() const override;
};

} // namespace olb

#endif
