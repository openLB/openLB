/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Orestis Malaspinas, Jonas Latt
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

#include "extendedFiniteDifferenceBoundary3D.h"
#include "extendedFiniteDifferenceBoundary3D.hh"
#include "boundaryPostProcessors3D.h"
#include "boundaryPostProcessors3D.hh"
#include "core/postProcessing.h"
 


namespace olb {

template class ExtendedFdPlaneBoundaryPostProcessor3D<double, descriptors::D3Q19<>,          0,1>;
template class ExtendedFdPlaneBoundaryProcessorGenerator3D<double, descriptors::D3Q19<>, 0,1>;
template class ExtendedFdPlaneBoundaryPostProcessor3D<double, descriptors::D3Q19<>,          0,-1>;
template class ExtendedFdPlaneBoundaryProcessorGenerator3D<double, descriptors::D3Q19<>, 0,-1>;
template class ExtendedFdPlaneBoundaryPostProcessor3D<double, descriptors::D3Q19<>,          1,1>;
template class ExtendedFdPlaneBoundaryProcessorGenerator3D<double, descriptors::D3Q19<>, 1,1>;
template class ExtendedFdPlaneBoundaryPostProcessor3D<double, descriptors::D3Q19<>,          1,-1>;
template class ExtendedFdPlaneBoundaryProcessorGenerator3D<double, descriptors::D3Q19<>, 1,-1>;
template class ExtendedFdPlaneBoundaryPostProcessor3D<double, descriptors::D3Q19<>,          2,1>;
template class ExtendedFdPlaneBoundaryProcessorGenerator3D<double, descriptors::D3Q19<>, 2,1>;
template class ExtendedFdPlaneBoundaryPostProcessor3D<double, descriptors::D3Q19<>,          2,-1>;
template class ExtendedFdPlaneBoundaryProcessorGenerator3D<double, descriptors::D3Q19<>, 2,-1>;


template class BoundaryConditionInstantiator3D
<
  double, descriptors::D3Q19<>,
  ExtendedFdBoundaryManager3D < double, descriptors::D3Q19<>,
                                BGKdynamics<double,descriptors::D3Q19<>> >
  >;

template OnLatticeBoundaryCondition3D<double,descriptors::D3Q19<>>*
createExtendedFdBoundaryCondition3D < double,descriptors::D3Q19<>,
                                    BGKdynamics<double,descriptors::D3Q19<>> >
                                    (
                                      BlockLatticeStructure3D<double,descriptors::D3Q19<>>& block
                                    );

}  // namespace olb
