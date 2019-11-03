/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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
 * Representation of a 2d block geometry structure -- template instantiation.
 */

#include "dynamics/latticeDescriptors.h"
 
#include "blockGeometryStructure2D.h"
#include "blockGeometryStructure2D.hh"

namespace olb {

template class BlockGeometryStructure2D<double>;

template bool BlockGeometryStructure2D<double>::findStreamDirections<double,descriptors::D2Q9<>>(int iX, int iY, BlockIndicatorF2D<double>& boundaryIndicator, BlockIndicatorF2D<double>& bulkIndicator, bool streamDirections[]);
template bool BlockGeometryStructure2D<double>::findStreamDirections<double,descriptors::D2Q9<>>(int iX, int iY, int material, std::list<int> bulkMaterials, bool streamDirections[]);

}
