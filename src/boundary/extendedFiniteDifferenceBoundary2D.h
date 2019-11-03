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


#ifndef EXTENDED_FINITE_DIFFERENCE_BOUNDARY_2D_H
#define EXTENDED_FINITE_DIFFERENCE_BOUNDARY_2D_H

#include "core/postProcessing.h"
#include "momentaOnBoundaries.h"
#include "core/blockLattice2D.h"
#include "boundaryCondition2D.h"


namespace olb {

/**
* This class computes the finite difference approximation to LB boundary conditions
* on a flat wall in 2D with all the terms of the CE expansion. The details
* on the computations can be found in the thesis of
* Jonas Latt, "Hydrodynamic limit of lattice Boltzmann equations",
* University of Geneva, (2007).
*/
template<typename T, typename DESCRIPTOR, int direction, int orientation>
class ExtendedStraightFdBoundaryPostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  ExtendedStraightFdBoundaryPostProcessor2D(int x0_, int x1_, int y0_, int y1_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_ ) override;
private:
  template<int deriveDirection>
  void interpolateGradients(BlockLattice2D<T,DESCRIPTOR> const& blockLattice,
                            T velDeriv[DESCRIPTOR::d], int iX, int iY) const;
  template<int deriveDirection>
  void interpolateGradients (
    BlockLattice2D<T,DESCRIPTOR> const& blockLattice,T& rhoDeriv,
    int iX, int iY ) const;
private:
  int x0, x1, y0, y1;
};

template<typename T, typename DESCRIPTOR, int direction, int orientation>
class ExtendedStraightFdBoundaryProcessorGenerator2D : public PostProcessorGenerator2D<T,DESCRIPTOR> {
public:
  ExtendedStraightFdBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_);
  PostProcessor2D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator2D<T,DESCRIPTOR>*  clone() const override;
};


////////// Factory function for Extended Finite Difference BC ///////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR> >
OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
createExtendedFdBoundaryCondition2D(BlockLatticeStructure2D<T,DESCRIPTOR>& block);

}  // namespace olb

#endif
