/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Marc Haussmann
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

#ifndef WALLFUNCTION_LATTICE_DESCRIPTOR_H
#define WALLFUNCTION_LATTICE_DESCRIPTOR_H

#include "dynamics/latticeDescriptors.h"


namespace olb {

namespace descriptors {

using WallFunctionForcedD3Q19Descriptor = D3Q19<AV_SHEAR,TAU_W,TAU_EFF,FORCE,V12>;

using WallFunctionForcedD3Q27Descriptor = D3Q27<AV_SHEAR,TAU_W,TAU_EFF,FORCE,V12>;

} // namespace descriptors

} // namespace olb

#endif
