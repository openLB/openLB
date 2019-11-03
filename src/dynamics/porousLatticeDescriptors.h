/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Mathias J. Krause, Jonas Latt
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

#ifndef POROUS_LATTICE_DESCRIPTOR_H
#define POROUS_LATTICE_DESCRIPTOR_H

#include "latticeDescriptors.h"


namespace olb {

namespace descriptors {

// 2D Descriptors for flow through porous media

using PorousD2Q9Descriptor = D2Q9<POROSITY>;

////////////////////////////////////////////////////////////////////////////////
// extended descriptor for drag computation - 2D

using ExtendedPorousD2Q9Descriptor = D2Q9<POROSITY,LOCAL_DRAG>;

////////////////////////////////////////////////////////////////////////////////
// extended descriptor for porous particles - 2D

using PorousParticleD2Q9Descriptor = D2Q9<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR>;

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow through porous media

using PorousD3Q19Descriptor = D3Q19<POROSITY>;

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow through porous media

using PorousForcedD3Q19Descriptor = D3Q19<POROSITY,FORCE>;

////////////////////////////////////////////////////////////////////////////////
// extended descriptor for drag computation - 3D

using ExtendedPorousD3Q19Descriptor = D3Q19<POROSITY,LOCAL_DRAG>;

////////////////////////////////////////////////////////////////////////////////
// extended descriptor for porous particles - 3D

using PorousParticleD3Q19Descriptor = D3Q19<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR>;

////////////////////////////////////////////////////////////////////////////////
// descriptor for PSM - 2D

using PSMD2Q9Descriptor = D2Q9<POROSITY,VELOCITY_SOLID>;

////////////////////////////////////////////////////////////////////////////////
// descriptor for PSM - 3D

using PSMD3Q19Descriptor = D3Q19<POROSITY,VELOCITY_SOLID>;

////////////////////////////////////////////////////////////////////////////////
// descriptor for ForcedPSM - 3D

using ForcedPSMD3Q19Descriptor = D3Q19<POROSITY,VELOCITY_SOLID,FORCE>;

} // namespace descriptors

} // namespace olb

#endif
