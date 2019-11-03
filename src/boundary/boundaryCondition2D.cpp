/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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

#include "boundaryCondition2D.h"
#include "boundaryCondition2D.hh"
#include "momentaOnBoundaries.h"
#include "momentaOnBoundaries.hh"
#include "momentaOnBoundaries2D.h"
#include "momentaOnBoundaries2D.hh"
#include "dynamics/latticeDescriptors.h"
 

namespace olb {

template class OnLatticeBoundaryCondition2D<double, descriptors::D2Q9<>>;

template class BoundaryConditionInstantiator2D<double, descriptors::D2Q9<>,
    RegularizedBoundaryManager2D < double, descriptors::D2Q9<>,
                                   RLBdynamics<double,descriptors::D2Q9<>> > >;

template OnLatticeBoundaryCondition2D<double,descriptors::D2Q9<>>*
createLocalBoundaryCondition2D < double,descriptors::D2Q9<>,
                               RLBdynamics<double,descriptors::D2Q9<>> >
                               (
                                 BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block
                               );

template OnLatticeBoundaryCondition2D<double,descriptors::D2Q9<>>*
createLocalBoundaryCondition2D < double,descriptors::D2Q9<>,
                               BGKdynamics<double,descriptors::D2Q9<>> >
                               (
                                 BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block
                               );

template class BoundaryConditionInstantiator2D<double, descriptors::D2Q9<>,
    InterpolationBoundaryManager2D < double, descriptors::D2Q9<>,
                                     BGKdynamics<double,descriptors::D2Q9<>> > >;

template OnLatticeBoundaryCondition2D<double,descriptors::D2Q9<>>*
createInterpBoundaryCondition2D < double,descriptors::D2Q9<>,
                                BGKdynamics<double,descriptors::D2Q9<>> >
                                (
                                  BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block
                                );

template OnLatticeBoundaryCondition2D<double,descriptors::D2Q9<>>*
createInterpBoundaryCondition2D < double,descriptors::D2Q9<>,
                                ConstRhoBGKdynamics<double,descriptors::D2Q9<>> >
                                (
                                  BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block
                                );
}
