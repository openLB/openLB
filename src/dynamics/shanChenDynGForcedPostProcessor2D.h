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

#ifndef SHAN_CHEN_DYN_G_FORCED_POST_PROCESSOR_2D_H
#define SHAN_CHEN_DYN_G_FORCED_POST_PROCESSOR_2D_H

#include "core/spatiallyExtendedObject2D.h"
#include "core/postProcessing.h"
#include "core/blockLattice2D.h"


namespace olb {

/**
* Multiphysics class for coupling between different lattices.
*/

// =========================================================================//
// ===========Shan Chen coupling without wall interaction===================//
// =========================================================================//

/*
template<typename T, typename DESCRIPTOR>
class ShanChenDynGForcedPostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  ShanChenDynGForcedPostProcessor2D(int x0_, int x1_, int y0_, int y1_, T G_,
                                    std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_,
                                    std::vector<SpatiallyExtendedObject2D*> partners_);
  ShanChenDynGForcedPostProcessor2D(T G_,
                                    std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_,
                                    std::vector<SpatiallyExtendedObject2D*> partners_);
  virtual int extent() const
  {
    return 1;
  }
  virtual int extent(int whichDirection) const
  {
    return 1;
  }
  virtual void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice);
  virtual void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_);
private:
  int x0, x1, y0, y1;
  T G;
  std::vector<T> rho0;
  AnalyticalF1D<T,T>& interactionPotential;
  std::vector<SpatiallyExtendedObject2D*> partners;
};

template<typename T, typename DESCRIPTOR>
class ShanChenDynGForcedGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  ShanChenDynGForcedGenerator2D(int x0_, int x1_, int y0_, int y1_, T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_);
  ShanChenDynGForcedGenerator2D(T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_);
  virtual PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject2D*> partners) const;
  virtual LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const;
private:
  T G;
  std::vector<T> rho0;
  AnalyticalF1D<T,T>& interactionPotential;
};
*/

}

#endif
