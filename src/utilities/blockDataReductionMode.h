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

#ifndef BLOCK_DATA_REDUCTION_MODE_H
#define BLOCK_DATA_REDUCTION_MODE_H

namespace olb {


/// Mode of reducing block data from given, possibly higher dimensional data
/**
* Required for optimizing block reduction functors such as BlockReduction3D2D
* if hyperplane is axis-aligned i.e. trivially discretizable.
**/
enum class BlockDataReductionMode {
  /// Interpolate block data at exact physical locations
  Analytical,
  /// Read block data from discrete lattice locations
  Discrete
};


}

#endif
