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

#ifndef NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_3D_H
#define NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_3D_H

#include <cmath>

#include "core/spatiallyExtendedObject3D.h"
#include "core/postProcessing.h"
#include "core/blockLattice3D.h"
#include "advectionDiffusionForces.hh"
#include "advectionDiffusionForces.h"

namespace olb {

 

/**
* Multiphysics class for coupling between different lattices.
*/
//======================================================================
// ========Regularized NSDiffusion Coupling 3D ====================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class NavierStokesAdvectionDiffusionCouplingPostProcessor3D :
  public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  NavierStokesAdvectionDiffusionCouplingPostProcessor3D(
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
    std::vector<SpatiallyExtendedObject3D* > partners_);
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
                                int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_) override;
private:
  typedef DESCRIPTOR L;
  int x0, x1, y0, y1, z0, z1;
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  T forcePrefactor[L::d];

  std::vector<SpatiallyExtendedObject3D*> partners;
  BlockLattice3D<T,descriptors::D3Q7<descriptors::VELOCITY>> *tPartner;
};

template<typename T, typename DESCRIPTOR>
class NavierStokesAdvectionDiffusionCouplingGenerator3D :
  public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  NavierStokesAdvectionDiffusionCouplingGenerator3D(
    int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_);
  PostProcessor3D<T,DESCRIPTOR>* generate(
    std::vector<SpatiallyExtendedObject3D* > partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
};


//======================================================================
// =============  SmagorinskyBoussinesqCoupling 3D ===================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class SmagorinskyBoussinesqCouplingPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  SmagorinskyBoussinesqCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_,
      std::vector<SpatiallyExtendedObject3D* > partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
private:
  typedef DESCRIPTOR L;
  int x0, x1, y0, y1, z0, z1;
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  T PrTurb;
  BlockLattice3D<T,descriptors::D3Q7<descriptors::VELOCITY,descriptors::TAU_EFF>> *tPartner;
  T forcePrefactor[L::d];
  T tauTurbADPrefactor;

  std::vector<SpatiallyExtendedObject3D*> partners;
};

template<typename T, typename DESCRIPTOR>
class SmagorinskyBoussinesqCouplingGenerator3D : public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  SmagorinskyBoussinesqCouplingGenerator3D(
      int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_);
  PostProcessor3D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject3D* > partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  T PrTurb;
};

//==================================================================================================
// ========Coupling 3D of Navier-Stokes on Advection-Diffusion====================//
//==================================================================================================
template<typename T, typename DESCRIPTOR,
typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvectionDiffusionParticleCouplingPostProcessor3D :
  public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  AdvectionDiffusionParticleCouplingPostProcessor3D(
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int iC_, int offset_,
    std::vector<SpatiallyExtendedObject3D* > partners_,
    std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T, DESCRIPTOR,ADLattice> > > forces_);
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
                                int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_) override;
private:
  int x0, x1, y0, y1, z0, z1, iC;
  T dragCoeff;
  int offset;
  T *vel, *vel_new, *velXp, *velXn, *velYp, *velYn, *velZp, *velZn;
  bool par = true;
  Cell<T,ADLattice> *adCell;
  Cell<T,DESCRIPTOR> *nsCell;

protected:
  std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T, DESCRIPTOR, ADLattice> > > forces;
};

template<typename T, typename DESCRIPTOR,
typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvectionDiffusionParticleCouplingGenerator3D :
  public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  AdvectionDiffusionParticleCouplingGenerator3D(int offset_);
  PostProcessor3D<T,DESCRIPTOR>* generate(
    std::vector<SpatiallyExtendedObject3D* > partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;
  void addForce(AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> &force);

private:
  int offset;

protected:
  std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T, DESCRIPTOR, ADLattice> > > ADforces;
};


//======================================================================
// ======== Porous Regularized NSDiffusion Coupling 3D =================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D :
  public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D(
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
    std::vector<SpatiallyExtendedObject3D* > partners_);
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
                                int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_) override;
private:
  int x0, x1, y0, y1, z0, z1;
  T gravity, T0, deltaTemp;
  std::vector<T> dir;

  std::vector<SpatiallyExtendedObject3D*> partners;
};

template<typename T, typename DESCRIPTOR>
class PorousNavierStokesAdvectionDiffusionCouplingGenerator3D :
  public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  PorousNavierStokesAdvectionDiffusionCouplingGenerator3D(
    int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_);
  PostProcessor3D<T,DESCRIPTOR>* generate(
    std::vector<SpatiallyExtendedObject3D* > partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
};


}

#endif
