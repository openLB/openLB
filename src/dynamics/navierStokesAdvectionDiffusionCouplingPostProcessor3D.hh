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

#ifndef NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_3D_HH
#define NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_3D_HH

#include "latticeDescriptors.h"
#include "navierStokesAdvectionDiffusionCouplingPostProcessor3D.h"
#include "core/blockLattice3D.h"
#include "core/util.h"
#include "core/finiteDifference3D.h"
#include "advectionDiffusionForces.hh"
#include "advectionDiffusionForces.h"

#include "porousAdvectionDiffusionDescriptors.h"

namespace olb {

//=====================================================================================
//==============  NavierStokesAdvectionDiffusionCouplingPostProcessor3D ===========
//=====================================================================================

template<typename T, typename DESCRIPTOR>
NavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>::
NavierStokesAdvectionDiffusionCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
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

  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    forcePrefactor[iD] = gravity * dir[iD];
  }

  tPartner = dynamic_cast<BlockLattice3D<T,descriptors::D3Q7<descriptors::VELOCITY>> *>(partners[0]);
}

template<typename T, typename DESCRIPTOR>
void NavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          // computation of the bousinessq force
          T *force = blockLattice.get(iX,iY,iZ).template getFieldPointer<descriptors::FORCE>();
          T temperatureDifference = tPartner->get(iX,iY,iZ).computeRho() - T0;
          for (unsigned iD = 0; iD < L::d; ++iD) {
            force[iD] = forcePrefactor[iD] * temperatureDifference;
          }
          // Velocity coupling
          T *u = tPartner->get(iX,iY,iZ).template getFieldPointer<descriptors::VELOCITY>();
          blockLattice.get(iX,iY,iZ).computeU(u);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void NavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}



//=====================================================================================
//==============  SmagorinskyBoussinesqCouplingPostProcessor3D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR>
SmagorinskyBoussinesqCouplingPostProcessor3D<T,DESCRIPTOR>::
SmagorinskyBoussinesqCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_,
    std::vector<SpatiallyExtendedObject3D* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
     gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
     dir(dir_), PrTurb(PrTurb_), partners(partners_)
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

  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    forcePrefactor[iD] = gravity * dir[iD];
  }

  tauTurbADPrefactor = descriptors::invCs2<T,descriptors::D3Q7<descriptors::VELOCITY,descriptors::TAU_EFF>>() / descriptors::invCs2<T,DESCRIPTOR>() / PrTurb;
  tPartner = dynamic_cast<BlockLattice3D<T,descriptors::D3Q7<descriptors::VELOCITY,descriptors::TAU_EFF>> *>(partners[0]);
}

template<typename T, typename DESCRIPTOR>
void SmagorinskyBoussinesqCouplingPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {

          // computation of the bousinessq force
          T *force = blockLattice.get(iX,iY,iZ).template getFieldPointer<descriptors::FORCE>();
          T temperatureDifference = tPartner->get(iX,iY,iZ).computeRho() - T0;
          for (unsigned iD = 0; iD < L::d; ++iD) {
            force[iD] = forcePrefactor[iD] * temperatureDifference;
          }

          // Velocity coupling
          T *u = tPartner->get(iX,iY,iZ).template getFieldPointer<descriptors::VELOCITY>();

          // tau coupling
          T *tauNS = blockLattice.get(iX,iY,iZ).template getFieldPointer<descriptors::TAU_EFF>();
          T *tauAD = tPartner->get(iX,iY,iZ).template getFieldPointer<descriptors::TAU_EFF>();

          T rho, pi[util::TensorVal<DESCRIPTOR >::n];
          blockLattice.get(iX,iY,iZ).computeAllMomenta(rho, u, pi);
          T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
          if (util::TensorVal<DESCRIPTOR >::n == 6) {
            PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
          }
          T PiNeqNorm    = sqrt(PiNeqNormSqr);
          /// Molecular realaxation time
          T tau_mol_NS = 1. / blockLattice.get(iX,iY,iZ).getDynamics()->getOmega();
          T tau_mol_AD = 1. / tPartner->get(iX,iY,iZ).getDynamics()->getOmega();
          /// Turbulent realaxation time
          T tau_turb_NS = 0.5*(sqrt(tau_mol_NS*tau_mol_NS + dynamic_cast<SmagorinskyDynamics<T,DESCRIPTOR>*>(blockLattice.get(iX,iY,iZ).getDynamics())->getPreFactor()/rho*PiNeqNorm) - tau_mol_NS);
          /// Effective realaxation time
          tauNS[0] = tau_mol_NS+tau_turb_NS;

          T tau_turb_AD = tau_turb_NS * tauTurbADPrefactor;
          tauAD[0] = tau_mol_AD+tau_turb_AD;
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SmagorinskyBoussinesqCouplingPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
SmagorinskyBoussinesqCouplingGenerator3D<T,DESCRIPTOR>::
SmagorinskyBoussinesqCouplingGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_), PrTurb(PrTurb_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* SmagorinskyBoussinesqCouplingGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject3D* > partners) const
{
  return new SmagorinskyBoussinesqCouplingPostProcessor3D<T,DESCRIPTOR>(
       this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, gravity, T0, deltaTemp, dir, PrTurb, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* SmagorinskyBoussinesqCouplingGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new SmagorinskyBoussinesqCouplingGenerator3D<T,DESCRIPTOR>(*this);
}



//=====================================================================================
//==============  AdvectionDiffusionParticleCouplingPostProcessor3D ===========
//=====================================================================================

template<typename T, typename DESCRIPTOR, typename ADLattice>
AdvectionDiffusionParticleCouplingPostProcessor3D<T,DESCRIPTOR,ADLattice>::
AdvectionDiffusionParticleCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int iC_,
    int offset_, std::vector<SpatiallyExtendedObject3D* > partners_,
    std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T, DESCRIPTOR, ADLattice> > > forces_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), iC(iC_), offset(offset_), forces(forces_)
{
  BlockLattice3D<T,ADLattice> *partnerLattice =
    dynamic_cast<BlockLattice3D<T,ADLattice> *>(partners_[0]);
  adCell = &partnerLattice->get(x0,y0,z0);
  vel = &partnerLattice->get(x0,y0,z0)[offset];
  vel_new = &partnerLattice->get(x0,y0,z0)[offset];
  velXp = &partnerLattice->get(x0+1,y0,z0)[offset];
  velXn = &partnerLattice->get(x0-1,y0,z0)[offset];
  velYp = &partnerLattice->get(x0,y0+1,z0)[offset];
  velYn = &partnerLattice->get(x0,y0-1,z0)[offset];
  velZp = &partnerLattice->get(x0,y0,z0+1)[offset];
  velZn = &partnerLattice->get(x0,y0,z0-1)[offset];
}

template<typename T, typename DESCRIPTOR, typename ADLattice>
void AdvectionDiffusionParticleCouplingPostProcessor3D<T,DESCRIPTOR,ADLattice>::
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

          if (forces.begin() != forces.end()) {
            // calculating upwind Gradient
            // vel contains velocity information on ADlattice
            // velGrad contains upwind vel on ADlattice
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

            for (AdvectionDiffusionForce3D<T, DESCRIPTOR, ADLattice>& f : forces) {
              // writes force in forceValues, vel refers to ADlattice              
              f.applyForce(forceValue, nsCell, adCell, &vel[off], latticeR);
            }

            // compute new particle velocity
            for (int i=0; i < DESCRIPTOR::d; i++) {
              vel_new[i+off2] = vel[i+off] + forceValue[i] - velGrad[i];
            }
          } else { // set particle velocity to fluid velocity
            nsCell->computeU(velF);
            for (int i = 0; i < DESCRIPTOR::d; i++) {
              vel_new[i + off2] = velF[i];
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
void AdvectionDiffusionParticleCouplingPostProcessor3D<T,DESCRIPTOR,ADLattice>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
NavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>::
NavierStokesAdvectionDiffusionCouplingGenerator3D(int x0_, int x1_, int y0_, int y1_,int z0_, int z1_, T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* NavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject3D* > partners) const
{
  return new NavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, gravity, T0, deltaTemp, dir,partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* NavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new NavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>(*this);
}

// LatticeCouplingGenerator for one-way advectionDiffusion coupling with Stokes drag

template<typename T, typename DESCRIPTOR,
typename ADLattice>
AdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,ADLattice>::
AdvectionDiffusionParticleCouplingGenerator3D(int offset_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(0, 0, 0, 0, 0, 0), offset(offset_)
{ }

template<typename T, typename DESCRIPTOR,
typename ADLattice>
PostProcessor3D<T,DESCRIPTOR>* AdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,ADLattice>::generate (
  std::vector<SpatiallyExtendedObject3D* > partners) const
{
  return new AdvectionDiffusionParticleCouplingPostProcessor3D<T,DESCRIPTOR,ADLattice>(this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, this->iC, offset, partners, ADforces);
}

template<typename T, typename DESCRIPTOR,
typename ADLattice>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* AdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,ADLattice>::clone() const
{
  return new AdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,ADLattice>(*this);
}

template<typename T, typename DESCRIPTOR,
typename ADLattice>
void AdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,ADLattice>::addForce(
    AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> &force)
{
  ADforces.push_back(force);
}


//=====================================================================================
//==============  PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D ===========
//=====================================================================================

template<typename T, typename DESCRIPTOR>
PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>::
PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
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
void PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  typedef DESCRIPTOR L;
  enum {x,y,z};

  BlockLattice3D<T,descriptors::PorousAdvectionDiffusionD3Q7Descriptor> *tPartner =
    dynamic_cast<BlockLattice3D<T,descriptors::PorousAdvectionDiffusionD3Q7Descriptor> *>(partners[0]);

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
          // this should return the interpolated solid-fluid temperature
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
void PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

template<typename T, typename DESCRIPTOR>
PorousNavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>::
PorousNavierStokesAdvectionDiffusionCouplingGenerator3D(int x0_, int x1_, int y0_, int y1_,int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* PorousNavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject3D* > partners) const
{
  return new PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, gravity, T0, deltaTemp, dir,partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* PorousNavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new PorousNavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>(*this);
}


}  // namespace olb

#endif
