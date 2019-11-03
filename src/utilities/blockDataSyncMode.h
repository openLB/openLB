/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef BLOCK_DATA_SYNC_MODE_H
#define BLOCK_DATA_SYNC_MODE_H

namespace olb {


/// Mode of synchronizing functor block data between processes
/**
* Required for optimizing functor operations to various usage patterns.
*
* i.e. the convention is for the full domain to be available on every rank
* but this is not ideal in most usage scenarios of reduction functors.
*
* e.g. the primary user of BlockReduction3D2D, BlockGifWriter, only requires
* full data to be available on the rank where its io is performed (rank 0).
*
* SuperLatticeFlux3D only requires rank-local data to be available. Any further
* synchronization would potentially impact performance in this critical area.
**/
enum class BlockDataSyncMode {
  /// default behavior, full block data available on all ranks after update
  ReduceAndBcast,
  /// optimize for usage in e.g. BlockGifWriter, full data only available on main rank
  ReduceOnly,
  /// optimize for usage in e.g. SuperLatticeFlux3D, only rank-local data available
  None
};


}

#endif
