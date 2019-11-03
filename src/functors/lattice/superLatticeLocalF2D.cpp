/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause
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

#include "superLatticeLocalF2D.h"
#include "superLatticeLocalF2D.hh"
#include "dynamics/latticeDescriptors.h"

namespace olb {

template class SuperLatticeDissipation2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysDissipation2D<double,descriptors::D2Q9<>>;
template class SuperLatticeDensity2D<double,descriptors::D2Q9<>>;
template class SuperLatticeVelocity2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysStrainRate2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysWallShearStress2D<double,descriptors::D2Q9<>>;
template class SuperLatticeGeometry2D<double,descriptors::D2Q9<>>;
template class SuperLatticeRank2D<double,descriptors::D2Q9<>>;
template class SuperLatticeCuboid2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysPressure2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysVelocity2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysBoundaryForce2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysCorrBoundaryForce2D<double,descriptors::D2Q9<>>;
template class SuperLatticePorosity2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysPermeability2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysDarcyForce2D<double,descriptors::D2Q9<>>;
template class SuperEuklidNorm2D<double,descriptors::D2Q9<>>;
template class SuperLatticeIndicatorSmoothIndicatorIntersection2D<double,descriptors::D2Q9<>,false>;

} // end namespace olb
