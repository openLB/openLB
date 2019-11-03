/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Patrick Nathen
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

#ifndef WALE_LATTICE_DESCRIPTOR_H
#define WALE_LATTICE_DESCRIPTOR_H

#include "dynamics/latticeDescriptors.h"


namespace olb {

namespace descriptors {



/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Wall Adaptive Local Eddy Viscosity (WALE)

using WALED3Q19Descriptor = D3Q19<EFFECTIVE_OMEGA,VELO_GRAD>;

using WALED3Q27Descriptor = D3Q27<EFFECTIVE_OMEGA,VELO_GRAD>;

/// WALE 3D Forced

using WALEForcedD3Q19Descriptor = D3Q19<EFFECTIVE_OMEGA,VELO_GRAD,FORCE>;

using WALEForcedD3Q27Descriptor = D3Q27<EFFECTIVE_OMEGA,VELO_GRAD,FORCE>;


} // namespace descriptors

} // namespace olb

#endif
