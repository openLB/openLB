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

#ifndef NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_MRT_COUPLING_POST_PROCESSOR_3D_H
#define NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_MRT_COUPLING_POST_PROCESSOR_3D_H

#include <cmath>

#include "core/spatiallyExtendedObject3D.h"
#include "core/postProcessing.h"
#include "core/blockLattice3D.h"
#include "advectionDiffusionForces.hh"
#include "advectionDiffusionForces.h"

namespace olb {

/**
 * @author Aleksandra Pachalieva
 * @date   17.02.2017
 *
 * Class for the coupling between a Navier-Stokes (NS) lattice and an
 * Advection-Diffusion MRT (AD) lattice.
 *
 * Multiphysics class for coupling between different lattices.
 *
 */
//======================================================================
// ======== Regularized NSDiffusion Coupling 3D ====================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class NavierStokesAdvectionDiffusionMRTCouplingPostProcessor3D :
  public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  NavierStokesAdvectionDiffusionMRTCouplingPostProcessor3D(
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
class NavierStokesAdvectionDiffusionMRTCouplingGenerator3D :
  public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  NavierStokesAdvectionDiffusionMRTCouplingGenerator3D(
    int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_);
  PostProcessor3D<T,DESCRIPTOR>* generate(
    std::vector<SpatiallyExtendedObject3D* > partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
};

//==================================================================================================
// ========Coupling 3D of Navier-Stokes on Advection-Diffusion with Stokes drag====================//
//==================================================================================================
template<typename T, typename DESCRIPTOR,
typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvectionDiffusionParticleMRTCouplingPostProcessor3D :
  public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  AdvectionDiffusionParticleMRTCouplingPostProcessor3D(
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
class AdvectionDiffusionParticleMRTCouplingGenerator3D :
  public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  AdvectionDiffusionParticleMRTCouplingGenerator3D(int offset_);
  PostProcessor3D<T,DESCRIPTOR>* generate(
    std::vector<SpatiallyExtendedObject3D* > partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;
  void addForce(AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> &force);

private:
  int offset;

protected:
  std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T, DESCRIPTOR, ADLattice> > > ADforces;
};

}

#endif
