/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
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

#ifndef SHAN_CHEN_FORCED_POST_PROCESSOR_3D_H
#define SHAN_CHEN_FORCED_POST_PROCESSOR_3D_H

#include "core/spatiallyExtendedObject3D.h"
#include "core/postProcessing.h"
#include "core/blockLattice3D.h"


namespace olb {

/**
* Multiphysics class for coupling between different lattices.
*/

// =========================================================================//
// ===========Shan Chen coupling without wall interaction===================//
// =========================================================================//

template<typename T, typename DESCRIPTOR>
class ShanChenForcedPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  ShanChenForcedPostProcessor3D (
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T G_, std::vector<T> rho0_,
    AnalyticalF1D<T,T>& iP_, std::vector<SpatiallyExtendedObject3D*> partners_);
  ShanChenForcedPostProcessor3D (
    T G_, std::vector<T> rho0_,
    AnalyticalF1D<T,T>& iP_, std::vector<SpatiallyExtendedObject3D*> partners_);
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
  int x0, x1, y0, y1, z0, z1;
  T G;
  std::vector<T> rho0;
  AnalyticalF1D<T,T>& interactionPotential;
  std::vector<SpatiallyExtendedObject3D*> partners;
};

template<typename T, typename DESCRIPTOR>
class ShanChenForcedGenerator3D : public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  ShanChenForcedGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_);
  ShanChenForcedGenerator3D(T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_);
  PostProcessor3D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject3D*> partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  T G;
  std::vector<T> rho0;
  AnalyticalF1D<T,T>& interactionPotential;
};

}

#endif
