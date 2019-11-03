/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2016 Robin Trunk
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

#include "advectionDiffusionBoundaryPostProcessor3D.h"
#include "advectionDiffusionBoundaryPostProcessor3D.hh"
#include "core/postProcessing.h"
#include "core/postProcessing.hh"
#include "latticeDescriptors.h"
 
//#include "advectionDiffusionLatticeDescriptors.h"
//#include "advectionDiffusionLatticeDescriptors.hh"


namespace olb {

template class ConvectionBoundaryProcessor3D<double, descriptors::D3Q7<descriptors::VELOCITY>>;
template class ConvectionBoundaryProcessorGenerator3D<double, descriptors::D3Q7<descriptors::VELOCITY>>;

}  // namespace olb
