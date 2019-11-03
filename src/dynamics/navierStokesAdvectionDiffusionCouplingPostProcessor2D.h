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

#ifndef NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_2D_H
#define NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_2D_H

#include "core/spatiallyExtendedObject2D.h"
#include "core/postProcessing.h"
#include "core/blockLattice2D.h"
#include <cmath>


namespace olb {

 

/**
* Class for the coupling between a Navier-Stokes (NS) lattice and an
* Advection-Diffusion (AD) lattice.
*/

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
// This coupling must be necessarily be put on the Navier-Stokes lattice!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

//======================================================================
// ========Regularized NSDiffusion Coupling 2D ====================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class NavierStokesAdvectionDiffusionCouplingPostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  NavierStokesAdvectionDiffusionCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
      std::vector<SpatiallyExtendedObject2D* > partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_) override;
private:
  typedef DESCRIPTOR L;
  int x0, x1, y0, y1;
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  BlockLattice2D<T,descriptors::D2Q5<descriptors::VELOCITY>> *tPartner;
  T forcePrefactor[L::d];

  std::vector<SpatiallyExtendedObject2D*> partners;
};

template<typename T, typename DESCRIPTOR>
class NavierStokesAdvectionDiffusionCouplingGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  NavierStokesAdvectionDiffusionCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_);
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject2D* > partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
};

template<typename T, typename DESCRIPTOR>
class SmagorinskyBoussinesqCouplingPostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  SmagorinskyBoussinesqCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_,
      std::vector<SpatiallyExtendedObject2D* > partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_) override;
private:
  typedef DESCRIPTOR L;
  int x0, x1, y0, y1;
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  T PrTurb;
  BlockLattice2D<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF>> *tPartner;
  T forcePrefactor[L::d];
  T tauTurbADPrefactor;

  std::vector<SpatiallyExtendedObject2D*> partners;
  enum {
    velOffset = descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF>::template index<descriptors::VELOCITY>(),
    forceOffset = descriptors::D2Q9<descriptors::FORCE,descriptors::TAU_EFF>::template index<descriptors::FORCE>(),
    tauADoffset = descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF>::template index<descriptors::TAU_EFF>(),
    tauNSoffset = descriptors::D2Q9<descriptors::FORCE,descriptors::TAU_EFF>::template index<descriptors::TAU_EFF>()
  };
};

template<typename T, typename DESCRIPTOR>
class SmagorinskyBoussinesqCouplingGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  SmagorinskyBoussinesqCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_);
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject2D* > partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  T PrTurb;
};

template<typename T, typename DESCRIPTOR>
class MixedScaleBoussinesqCouplingPostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  MixedScaleBoussinesqCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_,
      std::vector<SpatiallyExtendedObject2D* > partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_) override;
private:
  typedef DESCRIPTOR L;
  int x0, x1, y0, y1;
  T gravity, T0, deltaTemp, PrTurb;
  std::vector<T> dir;
  BlockLattice2D<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF,descriptors::CUTOFF_HEAT_FLUX>> *tPartner;
  T forcePrefactor[L::d];
  T tauTurbADPrefactor;

  std::vector<SpatiallyExtendedObject2D*> partners;
};

template<typename T, typename DESCRIPTOR>
class MixedScaleBoussinesqCouplingGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  MixedScaleBoussinesqCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_);
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject2D* > partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp, PrTurb;
  std::vector<T> dir;
};

}

#endif
