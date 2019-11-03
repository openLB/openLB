/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2016 Jonas Kratzke, Mathias J. Krause
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

#include "offBoundaryCondition2D.h"
#include "offBoundaryCondition2D.hh"
#include "dynamics/latticeDescriptors.h"
 

namespace olb {

template class OffLatticeBoundaryCondition2D<double, descriptors::D2Q9<>>;

template class OffBoundaryConditionInstantiator2D
<double, descriptors::D2Q9<>, BouzidiBoundaryManager2D < double, descriptors::D2Q9<>, BGKdynamics<double,descriptors::D2Q9<>> > >;

template class OffBoundaryConditionInstantiator2D
<double, descriptors::D2Q9<>, BounceBackBoundaryManager2D < double, descriptors::D2Q9<>> >;

template OffLatticeBoundaryCondition2D<double,descriptors::D2Q9<>>*
createBouzidiBoundaryCondition2D < double,descriptors::D2Q9<>,
                                 BGKdynamics<double,descriptors::D2Q9<>> >
                                 (
                                   BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block
                                 );

template OffLatticeBoundaryCondition2D<double,descriptors::D2Q9<>>*
createBounceBackBoundaryCondition2D < double,descriptors::D2Q9<>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block);
}

