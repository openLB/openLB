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
 * The description of a 2D super external field -- generic implementation.
 */


#ifndef SUPER_EXTERNAL_2D_HH
#define SUPER_EXTERNAL_2D_HH

#include<limits>
#include "geometry/superGeometry2D.h"
#include "core/superExternal2D.h"


namespace olb {


////////////////////// Class SuperExternal2D /////////////////////////


template<typename T, typename DESCRIPTOR, typename FIELD>
SuperExternal2D<T,DESCRIPTOR,FIELD>::SuperExternal2D(SuperGeometry2D<T>& superGeometry,
    SuperLattice2D<T, DESCRIPTOR>& sLattice, int overlap)
  : SuperStructure2D<T>(superGeometry.getCuboidGeometry(), superGeometry.getLoadBalancer()),
    _overlap(overlap), _sLattice(sLattice)
{
  this->_communicator.init_nh();
  this->_communicator.add_cells(this->_overlap);
  this->_communicator.init();
}

template<typename T, typename DESCRIPTOR, typename FIELD>
void SuperExternal2D<T,DESCRIPTOR,FIELD>::communicate(bool verbose)
{
  this->_communicator.send();
  this->_communicator.receive();
  this->_communicator.write();
}

template<typename T, typename DESCRIPTOR, typename FIELD>
bool* SuperExternal2D<T,DESCRIPTOR,FIELD>::operator() (int iCloc, int iX, int iY, int iData)
{
  return (bool*) _sLattice.getExtendedBlockLattice(iCloc).get(iX+_overlap, iY+_overlap).template getFieldPointer<FIELD>();
}

template<typename T, typename DESCRIPTOR, typename FIELD>
int SuperExternal2D<T,DESCRIPTOR,FIELD>::getDataSize() const
{
  return DESCRIPTOR::template size<FIELD>();
}

template<typename T, typename DESCRIPTOR, typename FIELD>
int SuperExternal2D<T,DESCRIPTOR,FIELD>::getDataTypeSize() const
{
  return sizeof(T);
}

} // namespace olb

#endif
