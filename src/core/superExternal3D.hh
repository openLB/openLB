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
 * The description of a 3D super external field -- generic implementation.
 */


#ifndef SUPER_EXTERNAL_3D_HH
#define SUPER_EXTERNAL_3D_HH


#include "superExternal3D.h"
#include "geometry/superGeometry3D.h"



namespace olb {


template<typename T, typename DESCRIPTOR, typename FIELD>
SuperExternal3D<T,DESCRIPTOR,FIELD>::SuperExternal3D(SuperGeometry3D<T>& superGeometry,
    SuperLattice3D<T,DESCRIPTOR>& sLattice, int overlap)
  : SuperStructure3D<T>(superGeometry.getCuboidGeometry(), superGeometry.getLoadBalancer()),
    _overlap(overlap), _sLattice(sLattice)
{
  this->_communicator.init_nh();
  this->_communicator.add_cells(this->_overlap);
  this->_communicator.init();
}

template<typename T, typename DESCRIPTOR, typename FIELD>
void SuperExternal3D<T,DESCRIPTOR,FIELD>::communicate(bool verbose)
{
  this->_communicator.send();
  this->_communicator.receive();
  this->_communicator.write();
}

} // namespace olb

#endif
