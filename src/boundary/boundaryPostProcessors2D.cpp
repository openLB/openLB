/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

#include "momentaOnBoundaries.h"
#include "momentaOnBoundaries.hh"
#include "boundaryPostProcessors2D.h"
#include "boundaryPostProcessors2D.hh"
#include "core/postProcessing.h"
#include "core/postProcessing.hh"
#include "dynamics/latticeDescriptors.h"
 


namespace olb {

template class StraightFdBoundaryProcessor2D<double, descriptors::D2Q9<>,          0,1>;
template class StraightFdBoundaryProcessorGenerator2D<double, descriptors::D2Q9<>, 0,1>;
template class StraightFdBoundaryProcessor2D<double, descriptors::D2Q9<>,          0,-1>;
template class StraightFdBoundaryProcessorGenerator2D<double, descriptors::D2Q9<>, 0,-1>;
template class StraightFdBoundaryProcessor2D<double, descriptors::D2Q9<>,          1,1>;
template class StraightFdBoundaryProcessorGenerator2D<double, descriptors::D2Q9<>, 1,1>;
template class StraightFdBoundaryProcessor2D<double, descriptors::D2Q9<>,          1,-1>;
template class StraightFdBoundaryProcessorGenerator2D<double, descriptors::D2Q9<>, 1,-1>;

template class StraightConvectionBoundaryProcessor2D<double, descriptors::D2Q9<>,          0,1>;
template class StraightConvectionBoundaryProcessorGenerator2D<double, descriptors::D2Q9<>, 0,1>;
template class StraightConvectionBoundaryProcessor2D<double, descriptors::D2Q9<>,          0,-1>;
template class StraightConvectionBoundaryProcessorGenerator2D<double, descriptors::D2Q9<>, 0,-1>;
template class StraightConvectionBoundaryProcessor2D<double, descriptors::D2Q9<>,          1,1>;
template class StraightConvectionBoundaryProcessorGenerator2D<double, descriptors::D2Q9<>, 1,1>;
template class StraightConvectionBoundaryProcessor2D<double, descriptors::D2Q9<>,          1,-1>;
template class StraightConvectionBoundaryProcessorGenerator2D<double, descriptors::D2Q9<>, 1,-1>;

template class SlipBoundaryProcessor2D<double, descriptors::D2Q9<>>;
template class SlipBoundaryProcessorGenerator2D<double, descriptors::D2Q9<>>;

template class PartialSlipBoundaryProcessor2D<double, descriptors::D2Q9<>>;
template class PartialSlipBoundaryProcessorGenerator2D<double, descriptors::D2Q9<>>;

template class OuterVelocityCornerProcessor2D          <double, descriptors::D2Q9<>, 1,1>;
template class OuterVelocityCornerProcessorGenerator2D <double, descriptors::D2Q9<>, 1,1>;
template class OuterVelocityCornerProcessor2D          <double, descriptors::D2Q9<>, 1,-1>;
template class OuterVelocityCornerProcessorGenerator2D <double, descriptors::D2Q9<>, 1,-1>;
template class OuterVelocityCornerProcessor2D          <double, descriptors::D2Q9<>,-1,1>;
template class OuterVelocityCornerProcessorGenerator2D <double, descriptors::D2Q9<>,-1,1>;
template class OuterVelocityCornerProcessor2D          <double, descriptors::D2Q9<>,-1,-1>;
template class OuterVelocityCornerProcessorGenerator2D <double, descriptors::D2Q9<>,-1,-1>;

template class FreeEnergyWallProcessor2D<double, descriptors::D2Q9<>>;
template class FreeEnergyWallProcessorGenerator2D<double, descriptors::D2Q9<>>;

template class FreeEnergyChemPotBoundaryProcessor2D<double, descriptors::D2Q9<>>;
template class FreeEnergyChemPotBoundaryProcessorGenerator2D<double, descriptors::D2Q9<>>;

template class FreeEnergyConvectiveProcessor2D<double, descriptors::D2Q9<>>;
template class FreeEnergyConvectiveProcessorGenerator2D<double, descriptors::D2Q9<>>;

}  // namespace olb
