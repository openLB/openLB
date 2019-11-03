/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software you can redistribute it and/or
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

#ifndef NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_MRT_COUPLING_POST_PROCESSOR_3D_HH
#define NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_MRT_COUPLING_POST_PROCESSOR_3D_HH

#include "latticeDescriptors.h"
#include "navierStokesAdvectionDiffusionMRTCouplingPostProcessor3D.h"
#include "core/blockLattice3D.h"
#include "core/util.h"
#include "core/finiteDifference3D.h"
#include "advectionDiffusionForces.hh"
#include "advectionDiffusionForces.h"

namespace olb {

//=====================================================================================
//==============  NavierStokesAdvectionDiffusionMRTCouplingPostProcessor3D ===========
//=====================================================================================
template<typename T, typename DESCRIPTOR>
NavierStokesAdvectionDiffusionMRTCouplingPostProcessor3D<T,DESCRIPTOR>::
NavierStokesAdvectionDiffusionMRTCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
    std::vector<SpatiallyExtendedObject3D* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
     gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
     dir(dir_), partners(partners_)
{
  // we normalize the direction of force vector
  T normDir = T();
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    normDir += dir[iD]*dir[iD];
  }
  normDir = sqrt(normDir);
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    dir[iD] /= normDir;
  }
}

template<typename T, typename DESCRIPTOR>
void NavierStokesAdvectionDiffusionMRTCouplingPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  typedef DESCRIPTOR L;
  enum {x,y,z};
  enum {
    velOffset = descriptors::AdvectionDiffusionMRTD3Q7Descriptor::template index<descriptors::VELOCITY>(),
    forceOffset = DESCRIPTOR::template index<descriptors::FORCE>()
  };

  BlockLattice3D<T,descriptors::AdvectionDiffusionMRTD3Q7Descriptor> *tPartner =
    dynamic_cast<BlockLattice3D<T,descriptors::AdvectionDiffusionMRTD3Q7Descriptor> *>(partners[0]);

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          //                  Velocity coupling
          T *u = tPartner->get(iX,iY,iZ).template getFieldPointer<descriptors::VELOCITY>();
          blockLattice.get(iX,iY,iZ).computeU(u);

          //coupling between the temperature and navier stokes.

          T *force = blockLattice.get(iX,iY,iZ).template getFieldPointer<descriptors::FORCE>();
          T temperature = tPartner->get(iX,iY,iZ).computeRho();
          T rho = blockLattice.get(iX,iY,iZ).computeRho();
          for (unsigned iD = 0; iD < L::d; ++iD) {
            force[iD] = gravity * rho * (temperature - T0) / deltaTemp * dir[iD];
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void NavierStokesAdvectionDiffusionMRTCouplingPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

//=====================================================================================
//==============  StokesDragCouplingPostProcessor3D ===========
//=====================================================================================

template<typename T, typename DESCRIPTOR,
typename ADLattice>
AdvectionDiffusionParticleMRTCouplingPostProcessor3D<T,DESCRIPTOR,ADLattice>::
AdvectionDiffusionParticleMRTCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int iC_,
    int offset_, std::vector<SpatiallyExtendedObject3D* > partners_,
    std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T, DESCRIPTOR, ADLattice> > > forces_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), iC(iC_), offset(offset_), forces(forces_)
{
  BlockLattice3D<T,ADLattice> *partnerLattice =
    dynamic_cast<BlockLattice3D<T,ADLattice> *>(partners_[0]);
  adCell = &partnerLattice->get(x0,y0,z0);
  vel = partnerLattice->get(x0,y0,z0)[offset];
  vel_new = partnerLattice->get(x0,y0,z0)[offset];
  velXp = partnerLattice->get(x0+1,y0,z0)[offset];
  velXn = partnerLattice->get(x0-1,y0,z0)[offset];
  velYp = partnerLattice->get(x0,y0+1,z0)[offset];
  velYn = partnerLattice->get(x0,y0-1,z0)[offset];
  velZp = partnerLattice->get(x0,y0,z0+1)[offset];
  velZn = partnerLattice->get(x0,y0,z0-1)[offset];
}

template<typename T, typename DESCRIPTOR,
typename ADLattice>
void AdvectionDiffusionParticleMRTCouplingPostProcessor3D<T,DESCRIPTOR,ADLattice>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;

  int off = (par) ? 3 : 0;
  int off2 = (par) ? 0 : 3;

  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {

          int latticeR[4] = {iC, iX, iY, iZ};
          T velGrad[3] = {0.,0.,0.};
          T forceValue[3] = {0.,0.,0.};
          T velF[3] = {0.,0.,0.};

          nsCell = &(blockLattice.get(iX,iY,iZ));
          //.computeU(velF);

          if (forces.begin() != forces.end()) {
            // calculating upwind Gradient
            if (vel[0+off]<0.) {
              velGrad[0] = vel[0+off]*(velXp[0+off]-vel[0+off]);
              velGrad[1] = vel[0+off]*(velXp[1+off]-vel[1+off]);
              velGrad[2] = vel[0+off]*(velXp[2+off]-vel[2+off]);
            } else {
              velGrad[0] = vel[0+off]*(vel[0+off]-velXn[0+off]);
              velGrad[1] = vel[0+off]*(vel[1+off]-velXn[1+off]);
              velGrad[2] = vel[0+off]*(vel[2+off]-velXn[2+off]);
            }
            if (vel[1+off]<0.) {
              velGrad[0] += vel[1+off]*(velYp[0+off]-vel[0+off]);
              velGrad[1] += vel[1+off]*(velYp[1+off]-vel[1+off]);
              velGrad[2] += vel[1+off]*(velYp[2+off]-vel[2+off]);
            } else {
              velGrad[0] += vel[1+off]*(vel[0+off]-velYn[0+off]);
              velGrad[1] += vel[1+off]*(vel[1+off]-velYn[1+off]);
              velGrad[2] += vel[1+off]*(vel[2+off]-velYn[2+off]);
            }
            if (vel[2+off]<0.) {
              velGrad[0] += vel[2+off]*(velZp[0+off]-vel[0+off]);
              velGrad[1] += vel[2+off]*(velZp[1+off]-vel[1+off]);
              velGrad[2] += vel[2+off]*(velZp[2+off]-vel[2+off]);
            } else {
              velGrad[0] += vel[2+off]*(vel[0+off]-velZn[0+off]);
              velGrad[1] += vel[2+off]*(vel[1+off]-velZn[1+off]);
              velGrad[2] += vel[2+off]*(vel[2+off]-velZn[2+off]);
            }

            for (AdvectionDiffusionForce3D<T, DESCRIPTOR,ADLattice>& f : forces) {
              f.applyForce(forceValue, nsCell, adCell, &vel[off], latticeR);
            }

            // compute new particle velocity
            for (int i=0; i < DESCRIPTOR::d; i++) {
              vel_new[i+off2] = vel[i+off] + forceValue[i] - velGrad[i];
            }
          } else {
            nsCell->computeU(velF);
            T scaleFactor = T(1);
            for (int i = 0; i < DESCRIPTOR::d; i++) {
              vel_new[i + off2] = velF[i] * scaleFactor;
            }
          }
        }
      }
    }
  }
  par = !par;
}

template<typename T, typename DESCRIPTOR,
typename ADLattice>
void AdvectionDiffusionParticleMRTCouplingPostProcessor3D<T,DESCRIPTOR,ADLattice>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
NavierStokesAdvectionDiffusionMRTCouplingGenerator3D<T,DESCRIPTOR>::
NavierStokesAdvectionDiffusionMRTCouplingGenerator3D(int x0_, int x1_, int y0_, int y1_,int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* NavierStokesAdvectionDiffusionMRTCouplingGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject3D* > partners) const
{
  return new NavierStokesAdvectionDiffusionMRTCouplingPostProcessor3D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, gravity, T0, deltaTemp, dir,partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* NavierStokesAdvectionDiffusionMRTCouplingGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new NavierStokesAdvectionDiffusionMRTCouplingGenerator3D<T,DESCRIPTOR>(*this);
}

// LatticeCouplingGenerator for one-way advectionDiffusion coupling with Stokes drag

template<typename T, typename DESCRIPTOR,
typename ADLattice>
AdvectionDiffusionParticleMRTCouplingGenerator3D<T,DESCRIPTOR,ADLattice>::
AdvectionDiffusionParticleMRTCouplingGenerator3D(int offset_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(0, 0, 0, 0, 0, 0), offset(offset_)
{ }

template<typename T, typename DESCRIPTOR,
typename ADLattice>
PostProcessor3D<T,DESCRIPTOR>* AdvectionDiffusionParticleMRTCouplingGenerator3D<T,DESCRIPTOR,ADLattice>::generate (
  std::vector<SpatiallyExtendedObject3D* > partners) const
{
  return new AdvectionDiffusionParticleMRTCouplingPostProcessor3D<T,DESCRIPTOR,ADLattice>(this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, this->iC, offset, partners, ADforces);
}

template<typename T, typename DESCRIPTOR,
typename ADLattice>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* AdvectionDiffusionParticleMRTCouplingGenerator3D<T,DESCRIPTOR,ADLattice>::clone() const
{
  return new AdvectionDiffusionParticleMRTCouplingGenerator3D<T,DESCRIPTOR,ADLattice>(*this);
}

template<typename T, typename DESCRIPTOR,
typename ADLattice>
void AdvectionDiffusionParticleMRTCouplingGenerator3D<T,DESCRIPTOR,ADLattice>::addForce(
    AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> &force)
{
  ADforces.push_back(force);
}


}  // namespace olb

#endif
