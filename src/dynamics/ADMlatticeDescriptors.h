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

#ifndef ADM_BGK_DYNAMICS_DESCRIPTOR_H
#define ADM_BGK_DYNAMICS_DESCRIPTOR_H

#include "dynamics/latticeDescriptors.h"
//#include <cmath>


namespace olb {

namespace descriptors {

// 2D Descriptors for ADM

using ADMD2Q9Descriptor = D2Q9<FIL_RHO,LOCAL_FIL_VEL_X,LOCAL_FIL_VEL_Y>;

////////////////////////////////////////////////////////////////////////////////
// extended descriptor for ADM

using ADMD3Q19Descriptor = D3Q19<FIL_RHO,LOCAL_FIL_VEL_X,LOCAL_FIL_VEL_Y,LOCAL_FIL_VEL_Z>;

////////////////////////////////////////
/// ADM Descriptors for forced fields

//// Forced 2D ADM scheme
using ForcedADMD2Q9Descriptor = D2Q9<FORCE,FIL_RHO,LOCAL_FIL_VEL_X,LOCAL_FIL_VEL_Y>;

//// Forced 3D ADM scheme

using ForcedADMD3Q19Descriptor = D3Q19<FORCE,FIL_RHO,LOCAL_FIL_VEL_X,LOCAL_FIL_VEL_Y,LOCAL_FIL_VEL_Z>;

//// Forced adapted 3D ADM scheme

using ForcedAdaptiveADMD3Q19Descriptor = D3Q19<FORCE,FIL_RHO,LOCAL_FIL_VEL_X,LOCAL_FIL_VEL_Y,LOCAL_FIL_VEL_Z,LOCAL_AV_DISS,LOCAL_AV_TKE,LOCAL_SIGMA_ADM,LOCAL_NU_EDDY,TAU_W>;


} // namespace descriptors

} // namespace olb

#endif
