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

#ifndef POROUS_ADVECTION_DIFFUSION_DESCRIPTORS_H
#define POROUS_ADVECTION_DIFFUSION_DESCRIPTORS_H

#include "latticeDescriptors.h"
#include "porousLatticeDescriptors.h"

namespace olb {

namespace descriptors {

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow through porous media with advection diffusion

using PorousAdvectionDiffusionD3Q7Descriptor = D3Q7<POROSITY,VELOCITY>;

} // namespace descriptors

} // namespace olb

#endif
