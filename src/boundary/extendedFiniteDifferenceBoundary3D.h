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

#ifndef EXTENDED_FINITE_DIFFERENCE_BOUNDARY_3D_H
#define EXTENDED_FINITE_DIFFERENCE_BOUNDARY_3D_H

#include "core/postProcessing.h"
#include "momentaOnBoundaries.h"
#include "core/blockLattice3D.h"
#include "boundaryCondition3D.h"


namespace olb {

/**
* This class computes the finite difference approximation to LB boundary conditions
* on a plane wall in 3D with all the terms of the CE expansion.
*/

template<typename T, typename DESCRIPTOR, int direction, int orientation>
class ExtendedFdPlaneBoundaryPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  ExtendedFdPlaneBoundaryPostProcessor3D (int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
private:
  template<int deriveDirection>
  void interpolateGradients( BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
                             T velDeriv[DESCRIPTOR::d], int iX, int iY, int iZ ) const;

  template<int deriveDirection>
  void interpolateGradients ( BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
                              T& rhoDeriv, int iX, int iY, int iZ) const;
private:
  int x0, x1, y0, y1, z0, z1;
};

template<typename T, typename DESCRIPTOR, int direction, int orientation>
class ExtendedFdPlaneBoundaryProcessorGenerator3D
  : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  ExtendedFdPlaneBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
};


////////// Factory function for Extended Finite Difference BC ///////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR> >
OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* createExtendedFdBoundaryCondition3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block);

}

#endif
