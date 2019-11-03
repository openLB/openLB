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

#include "blockLatticeLocalF2D.h"
#include "blockLatticeLocalF2D.hh"
#include "dynamics/latticeDescriptors.h"

namespace olb {

template class BlockLatticeDissipation2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysDissipation2D<double,descriptors::D2Q9<>>;
template class BlockLatticeDensity2D<double,descriptors::D2Q9<>>;
template class BlockLatticeVelocity2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysStrainRate2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysWallShearStress2D<double,descriptors::D2Q9<>>;
template class BlockLatticeGeometry2D<double,descriptors::D2Q9<>>;
template class BlockLatticeRank2D<double,descriptors::D2Q9<>>;
template class BlockLatticeCuboid2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysPressure2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysVelocity2D<double,descriptors::D2Q9<>>;
//template class BlockLatticePhysExternalPorosity2D<double,descriptors::D2Q9<>>;
//template class BlockLatticePhysExternalVelocity2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysBoundaryForce2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysCorrBoundaryForce2D<double,descriptors::D2Q9<>>;
template class BlockLatticePorosity2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysPermeability2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysDarcyForce2D<double,descriptors::D2Q9<>>;
template class BlockLatticeAverage2D<double,descriptors::D2Q9<>>;
template class BlockEuklidNorm2D<double,descriptors::D2Q9<>>;
template class BlockLatticeIndicatorSmoothIndicatorIntersection2D<double,descriptors::D2Q9<>,false>;
} // end namespace olb
