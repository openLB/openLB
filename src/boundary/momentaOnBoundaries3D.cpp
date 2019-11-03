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

/** \file
 * Implementation of boundary cell dynamics -- template instantiation.
 */
#include "momentaOnBoundaries.h"
#include "momentaOnBoundaries.hh"
#include "momentaOnBoundaries3D.h"
#include "momentaOnBoundaries3D.hh"
#include "dynamics/latticeDescriptors.h"
 

namespace olb {
template class EquilibriumBM<double, descriptors::D3Q19<>>;

template class FixedVelocityBM<double, descriptors::D3Q19<>>;

template class FreeStressBM<double, descriptors::D3Q19<>>;

template class VelocityBM<double, descriptors::D3Q19<>, 0, 1>;
template class VelocityBM<double, descriptors::D3Q19<>, 0, -1>;
template class VelocityBM<double, descriptors::D3Q19<>, 1, 1>;
template class VelocityBM<double, descriptors::D3Q19<>, 1, -1>;
template class VelocityBM<double, descriptors::D3Q19<>, 2, 1>;
template class VelocityBM<double, descriptors::D3Q19<>, 2, -1>;

template class RegularizedVelocityBM<double, descriptors::D3Q19<>, 0, 1>;
template class RegularizedVelocityBM<double, descriptors::D3Q19<>, 0, -1>;
template class RegularizedVelocityBM<double, descriptors::D3Q19<>, 1, 1>;
template class RegularizedVelocityBM<double, descriptors::D3Q19<>, 1, -1>;
template class RegularizedVelocityBM<double, descriptors::D3Q19<>, 2, 1>;
template class RegularizedVelocityBM<double, descriptors::D3Q19<>, 2, -1>;

template class PressureBM<double, descriptors::D3Q19<>, 0, 1>;
template class PressureBM<double, descriptors::D3Q19<>, 0, -1>;
template class PressureBM<double, descriptors::D3Q19<>, 1, 1>;
template class PressureBM<double, descriptors::D3Q19<>, 1, -1>;
template class PressureBM<double, descriptors::D3Q19<>, 2, 1>;
template class PressureBM<double, descriptors::D3Q19<>, 2, -1>;

template class RegularizedPressureBM<double, descriptors::D3Q19<>, 0, 1>;
template class RegularizedPressureBM<double, descriptors::D3Q19<>, 0, -1>;
template class RegularizedPressureBM<double, descriptors::D3Q19<>, 1, 1>;
template class RegularizedPressureBM<double, descriptors::D3Q19<>, 1, -1>;
template class RegularizedPressureBM<double, descriptors::D3Q19<>, 2, 1>;
template class RegularizedPressureBM<double, descriptors::D3Q19<>, 2, -1>;

template class BasicDirichletBM<double, descriptors::D3Q19<>, VelocityBM, 0, 1>;
template class BasicDirichletBM<double, descriptors::D3Q19<>, VelocityBM, 0, -1>;
template class BasicDirichletBM<double, descriptors::D3Q19<>, VelocityBM, 1, 1>;
template class BasicDirichletBM<double, descriptors::D3Q19<>, VelocityBM, 1, -1>;
template class BasicDirichletBM<double, descriptors::D3Q19<>, VelocityBM, 2, 1>;
template class BasicDirichletBM<double, descriptors::D3Q19<>, VelocityBM, 2, -1>;

template class BasicDirichletBM<double, descriptors::D3Q19<>, PressureBM, 0, 1>;
template class BasicDirichletBM<double, descriptors::D3Q19<>, PressureBM, 0, -1>;
template class BasicDirichletBM<double, descriptors::D3Q19<>, PressureBM, 1, 1>;
template class BasicDirichletBM<double, descriptors::D3Q19<>, PressureBM, 1, -1>;
template class BasicDirichletBM<double, descriptors::D3Q19<>, PressureBM, 2, 1>;
template class BasicDirichletBM<double, descriptors::D3Q19<>, PressureBM, 2, -1>;

template class InnerCornerVelBM3D<double, descriptors::D3Q19<>, 1, 1, 1>;
template class InnerCornerVelBM3D<double, descriptors::D3Q19<>, 1, 1,-1>;
template class InnerCornerVelBM3D<double, descriptors::D3Q19<>, 1,-1, 1>;
template class InnerCornerVelBM3D<double, descriptors::D3Q19<>, 1,-1,-1>;
template class InnerCornerVelBM3D<double, descriptors::D3Q19<>,-1, 1, 1>;
template class InnerCornerVelBM3D<double, descriptors::D3Q19<>,-1, 1,-1>;
template class InnerCornerVelBM3D<double, descriptors::D3Q19<>,-1,-1, 1>;
template class InnerCornerVelBM3D<double, descriptors::D3Q19<>,-1,-1,-1>;

template class InnerEdgeVelBM3D<double, descriptors::D3Q19<>, 0, 1, 1>;
template class InnerEdgeVelBM3D<double, descriptors::D3Q19<>, 0, 1,-1>;
template class InnerEdgeVelBM3D<double, descriptors::D3Q19<>, 0,-1, 1>;
template class InnerEdgeVelBM3D<double, descriptors::D3Q19<>, 0,-1,-1>;
template class InnerEdgeVelBM3D<double, descriptors::D3Q19<>, 1, 1, 1>;
template class InnerEdgeVelBM3D<double, descriptors::D3Q19<>, 1, 1,-1>;
template class InnerEdgeVelBM3D<double, descriptors::D3Q19<>, 1,-1, 1>;
template class InnerEdgeVelBM3D<double, descriptors::D3Q19<>, 1,-1,-1>;
template class InnerEdgeVelBM3D<double, descriptors::D3Q19<>, 2, 1, 1>;
template class InnerEdgeVelBM3D<double, descriptors::D3Q19<>, 2, 1,-1>;
template class InnerEdgeVelBM3D<double, descriptors::D3Q19<>, 2,-1, 1>;
template class InnerEdgeVelBM3D<double, descriptors::D3Q19<>, 2,-1,-1>;
}
