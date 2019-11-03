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

/** \file
 * A helper for initialising 2D boundaries -- template instantiation.
 */


#include "communication/mpiManager.h"
#include "superOffBoundaryCondition2D.h"
#include "superOffBoundaryCondition2D.hh"
#include "dynamics/latticeDescriptors.h"
 
#include "dynamics/dynamics.h"
#include "dynamics/dynamics.hh"

namespace olb {

template class sOffLatticeBoundaryCondition2D
<double, descriptors::D2Q9<>>;

template void createBouzidiBoundaryCondition2D
<double,descriptors::D2Q9<>,
 BGKdynamics<double,descriptors::D2Q9<>> >
(sOffLatticeBoundaryCondition2D<double,descriptors::D2Q9<>>& sBC);


template void createBounceBackBoundaryCondition2D
<double,descriptors::D2Q9<>>
(sOffLatticeBoundaryCondition2D<double,descriptors::D2Q9<>>& sBC);
}
