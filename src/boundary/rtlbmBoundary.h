/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Albert Mink
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

#ifndef RTLBM_BOUNDARY_H
#define RTLBM_BOUNDARY_H


#include "advectionDiffusionBoundaryCondition3D.h"

namespace olb {

template<typename T, typename DESCRIPTOR> class sOnLatticeBoundaryCondition3D;
template<typename T, typename DESCRIPTOR> class BlockLatticeStructure3D;

// blockLattice creator
template<typename T, typename DESCRIPTOR>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR>* createRtlbmDiffuseBoundaryCondition3D( BlockLatticeStructure3D<T,DESCRIPTOR>& block );

// superLattice creator
template<typename T, typename DESCRIPTOR>
void createRtlbmDiffuseBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sBC);


// blockLattice creator
template<typename T, typename DESCRIPTOR>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR>* createRtlbmDiffuseConstBoundaryCondition3D( BlockLatticeStructure3D<T,DESCRIPTOR>& block );

// superLattice creator
template<typename T, typename DESCRIPTOR>
void createRtlbmDiffuseConstBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sBC);

// blockLattice creator
template<typename T, typename DESCRIPTOR>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR>* createRtlbmDirectedBoundaryCondition3D( BlockLatticeStructure3D<T,DESCRIPTOR>& block );

// superLattice creator
template<typename T, typename DESCRIPTOR>
void createRtlbmDirectedBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sBC);
}  // namespace olb

#endif
