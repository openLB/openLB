/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Andrea Parmigiani, Orestis Malaspinas,
 *  Jonas Latt
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
 * DESCRIPTORBASE for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the OpenLB code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */
#ifndef SHAN_CHEN_FORCED_LATTICE_DESCRIPTORS_H
#define SHAN_CHEN_FORCED_LATTICE_DESCRIPTORS_H

#include <vector>

namespace olb {

namespace descriptors {

/// 2D Descriptors

using ShanChenForcedD2Q9Descriptor = D2Q9<VELOCITY,FORCE>;

/// 2D Descriptors for modells with variable omega and Forced Shan Chen

using ShanChenDynOmegaD2Q9Descriptor = D2Q9<VELOCITY,FORCE,OMEGA>;

using ShanChenDynOmegaForcedD2Q9Descriptor = D2Q9<VELOCITY,FORCE,EXTERNAL_FORCE,OMEGA>;

/// 2D Descriptors for modells with variable G and Forced Shan Chen

using ShanChenDynGD2Q9Descriptor = D2Q9<VELOCITY,FORCE,G>;

using ShanChenDynGForcedD2Q9Descriptor = D2Q9<VELOCITY,FORCE,EXTERNAL_FORCE,G>;

/// 3D Descriptors

using ShanChenForcedD3Q19Descriptor = D3Q19<VELOCITY,FORCE>;

/// 3D Descriptors for modells with variable omega and Forced Shan Chen

using ShanChenDynOmegaD3Q19Descriptor = D3Q19<VELOCITY,FORCE,OMEGA>;

using ShanChenDynOmegaForcedD3Q19Descriptor = D3Q19<VELOCITY,FORCE,EXTERNAL_FORCE,OMEGA>;

}  // namespace descriptors

}  // namespace olb

#endif
