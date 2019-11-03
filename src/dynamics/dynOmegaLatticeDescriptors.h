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
 * Descriptor for 2D and 3D lattices with dynamic omega. In principle,
 * thanks to the fact that the OpenLB code is generic, it is sufficient
 * to write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */

#ifndef DYN_OMEGA_DESCRIPTOR_H
#define DYN_OMEGA_DESCRIPTOR_H

#include "latticeDescriptors.h"

namespace olb {

/// Descriptors for 2D and 3D lattices with variable omega.
/** The implementation is to be extended by combination with other
 * lattice descriptors.
 */

namespace descriptors {


/// 2D Descriptors for modells with variable omega

/// extended descriptor for porous particles and power-law rheology - 2D

using DynOmegaD2Q9Descriptor = D2Q9<OMEGA>;
using ForcedDynOmegaD2Q9Descriptor =  D2Q9<FORCE,OMEGA>;

using DynOmegaPorousParticleD2Q9Descriptor = D2Q9<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR,OMEGA>;

/// 3D Descriptors for modells with variable omega

/// extended descriptor for porous particles and power-law rheology - 3D

using DynOmegaD3Q19Descriptor = D3Q19<OMEGA>;
using DynOmegaD3Q27Descriptor = D3Q27<OMEGA>;

using ForcedDynOmegaD3Q19Descriptor = D3Q19<OMEGA,FORCE>;
using ForcedDynOmegaD3Q27Descriptor = D3Q27<OMEGA,FORCE>;

using DynOmegaPorousParticleD3Q19Descriptor = D3Q19<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR,OMEGA>;
using DynOmegaPorousParticleD3Q27Descriptor = D3Q27<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR,OMEGA>;


} // namespace descriptors

} // namespace olb

#endif
