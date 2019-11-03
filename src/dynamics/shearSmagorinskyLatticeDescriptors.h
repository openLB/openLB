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

#ifndef SHEAR_SMAGORINSKY_LATTICE_DESCRIPTOR_H
#define SHEAR_SMAGORINSKY_LATTICE_DESCRIPTOR_H

#include "dynamics/latticeDescriptors.h"


namespace olb {

namespace descriptors {

// 2D Descriptors for flow with Shear-Improved Smagorinsky

using ShearSmagorinskyD2Q9Descriptor = D2Q9<AV_SHEAR>;

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Shear-Improved Smagorinsky

using ShearSmagorinskyD3Q19Descriptor = D3Q19<AV_SHEAR>;

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Forced Shear-Improved Smagorinsky

using ShearSmagorinskyForcedD3Q19Descriptor = D3Q19<AV_SHEAR,FORCE>;

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Forced Shear-Improved Smagorinsky

using ForcedShearWallSmagorinskyD3Q19Descriptor = D3Q19<AV_SHEAR,FORCE,TAU_W>;

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Shear-Improved Kalman Finitie Difference Smagorinsky

using FDKalmanShearSmagorinskyD3Q19Descriptor = D3Q19<ERROR_COVARIANCE,VARIANCE,VELOCITY,FILTERED_VEL_GRAD,VELO_GRAD>;

/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Shear-Improved Kalman Finitie Difference Smagorinsky

using FDKalmanShearSmagorinskyForcedD3Q19Descriptor = D3Q19<ERROR_COVARIANCE,VARIANCE,VELOCITY,FILTERED_VEL_GRAD,VELO_GRAD,FORCE>;

////////////////////////////////////////////////////
// Kalman filter : Adaptive exponential smoothing //
////////////////////////////////////////////////////
// 3D Descriptors for flow with Shear-Improved Smagorinsky - Kalman Filter
// Boudet et al. (2016) A Kalman filter adapted of the estimation of mean gradients
//   in the a large-eddy simulation of unsteady turbulent flows.

using KalmanShearSmagorinskyD3Q19Descriptor = D3Q19<ERROR_COVARIANCE,VARIANCE,TAU_SGS,FILTERED_POPULATION>;


} // namespace descriptors

} // namespace olb

#endif
