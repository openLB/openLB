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

#ifndef DYN_SMAGORINSKY_LATTICE_DESCRIPTOR_H
#define DYN_SMAGORINSKY_LATTICE_DESCRIPTOR_H

#include "latticeDescriptors.h"

namespace olb {

namespace descriptors {

// 2D Descriptors for flow through DynSmagorinsky media

using DynSmagorinskyD2Q9Descriptor = D2Q9<SMAGO_CONST>;

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow through DynSmagorinsky media

using DynSmagorinskyD3Q19Descriptor = D3Q19<SMAGO_CONST>;


} // namespace descriptors

} // namespace olb

#endif
