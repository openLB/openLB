/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#include "blockLatticeLocalF3D.h"
#include "blockLatticeLocalF3D.hh"
#include "dynamics/latticeDescriptors.h"

namespace olb {

template class BlockLatticeFpop3D<double,descriptors::D3Q19<>>;
template class BlockLatticeDissipation3D<double,descriptors::D3Q19<>>;
template class BlockLatticePhysDissipation3D<double,descriptors::D3Q19<>>;
template class BlockLatticeEffevtiveDissipation3D<double,descriptors::D3Q19<>>;
template class BlockLatticePhysEffevtiveDissipation3D<double,descriptors::D3Q19<>>;
template class BlockLatticeDensity3D<double,descriptors::D3Q19<>>;
template class BlockLatticeVelocity3D<double,descriptors::D3Q19<>>;
template class BlockLatticeFlux3D<double,descriptors::D3Q19<>>;
template class BlockLatticeGeometry3D<double,descriptors::D3Q19<>>;
template class BlockLatticeRank3D<double,descriptors::D3Q19<>>;
template class BlockLatticeCuboid3D<double,descriptors::D3Q19<>>;
template class BlockLatticePhysPressure3D<double,descriptors::D3Q19<>>;
template class BlockLatticePhysVelocity3D<double,descriptors::D3Q19<>>;
template class BlockLatticeStrainRate3D<double,descriptors::D3Q19<>>;
template class BlockLatticePhysBoundaryForce3D<double,descriptors::D3Q19<>>;
template class BlockLatticePhysCorrBoundaryForce3D<double,descriptors::D3Q19<>>;
template class BlockEuklidNorm3D<double,descriptors::D3Q19<>>;
template class BlockLatticeInterpPhysVelocity3D<double,descriptors::D3Q19<>>;
template class BlockLatticeIndicatorSmoothIndicatorIntersection3D<double,descriptors::D3Q19<>,false>;

} // end namespace olb
