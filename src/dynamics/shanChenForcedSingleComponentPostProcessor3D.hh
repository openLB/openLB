/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani,
 *  Jonas Latt, 2013 Mathias J. Krause
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

#ifndef SHAN_CHEN_FORCED_SINGLE_COMPONENT_POST_PROCESSOR_3D_HH
#define SHAN_CHEN_FORCED_SINGLE_COMPONENT_POST_PROCESSOR_3D_HH

#include "shanChenForcedSingleComponentPostProcessor3D.h"
#include "interactionPotential.h"
#include "core/blockLattice3D.h"
#include "core/util.h"
#include "core/finiteDifference3D.h"

namespace olb {

////////  ShanChenForcedSingleComponentPostProcessor3D ///////////////////////////////////


template<typename T, typename DESCRIPTOR>
ShanChenForcedSingleComponentPostProcessor3D <T,DESCRIPTOR>::
ShanChenForcedSingleComponentPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_,
    std::vector<SpatiallyExtendedObject3D*> partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), G(G_), rho0(rho0_), interactionPotential(iP_), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
ShanChenForcedSingleComponentPostProcessor3D <T,DESCRIPTOR>::
ShanChenForcedSingleComponentPostProcessor3D(T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_,
    std::vector<SpatiallyExtendedObject3D*> partners_)
  :  x0(0), x1(0), y0(0), y1(0), z0(0), z1(0), G(G_), rho0(rho0_), interactionPotential(iP_), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
void ShanChenForcedSingleComponentPostProcessor3D<T,DESCRIPTOR>::
processSubDomain( BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ )
{
  typedef DESCRIPTOR L;

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect ( x0, x1, y0, y1, z0, z1,
                         x0_, x1_, y0_, y1_, z0_, z1_,
                         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    int nx = newX1-newX0+3; // include a one-cell boundary
    int ny = newY1-newY0+3; // include a one-cell boundary
    int nz = newZ1-newZ0+3; // include a one-cell boundary
    int offsetX = newX0-1;
    int offsetY = newY0-1;
    int offsetZ = newZ0-1;

    BlockData3D<T,T> rhoField1(nx,ny,nz);

    // Compute density and velocity on every site of first lattice, and store result
    //   in external scalars; envelope cells are included, because they are needed
    //   to compute the interaction potential in what follows.
    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {
        for (int iZ=newZ0-1; iZ<=newZ1+1; ++iZ) {
          Cell<T,DESCRIPTOR>& cell = blockLattice.get(iX,iY,iZ);
          rhoField1.get(iX-offsetX, iY-offsetY, iZ-offsetZ) = cell.computeRho()*rho0[0];
        }
      }
    }

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          Cell<T,DESCRIPTOR>& blockCell   = blockLattice.get(iX,iY,iZ);

          T* j = blockCell.template getFieldPointer<descriptors::VELOCITY>();
          lbHelpers<T,DESCRIPTOR>::computeJ(blockCell,j);

          T blockOmega   = blockLattice.getDynamics(iX, iY, iZ)->getOmega();

          // Computation of the common velocity, shared among the two populations
          T rhoTot = rhoField1.get(iX-offsetX, iY-offsetY, iZ-offsetZ)*blockOmega;

          T uTot[DESCRIPTOR::d];
          T *blockU = blockCell.template getFieldPointer<descriptors::VELOCITY>(); // contains precomputed value rho*u
          for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
            uTot[iD] = (blockU[iD]*rho0[0]*blockOmega) / rhoTot;
          }

          // Computation of the interaction potential
          T rhoBlockContribution[L::d]   = {T(), T(), T()};
          T psi;
          interactionPotential(&psi, &rhoField1.get(iX-offsetX, iY-offsetY, iZ-offsetZ));
          for (int iPop = 0; iPop < L::q; ++iPop) {
            int nextX = iX + descriptors::c<L>(iPop,0);
            int nextY = iY + descriptors::c<L>(iPop,1);
            int nextZ = iZ + descriptors::c<L>(iPop,2);
            T blockRho;
            interactionPotential(&blockRho, &rhoField1.get(nextX-offsetX, nextY-offsetY, nextZ-offsetZ));
            for (int iD = 0; iD < L::d; ++iD) {
              rhoBlockContribution[iD]   += psi * blockRho * descriptors::c<L>(iPop,iD)* descriptors::t<T,L>(iPop);
            }
          }

          // Computation and storage of the final velocity, consisting
          //   of u and the momentum difference due to interaction
          //   potential plus external force
          T *blockForce = blockCell.template getFieldPointer<descriptors::FORCE>();
          T *externalBlockForce = blockCell.template getFieldPointer<descriptors::EXTERNAL_FORCE>();

          for (int iD = 0; iD < L::d; ++iD) {
            blockU[iD] = uTot[iD];
            blockForce[iD] = externalBlockForce[iD] - G*rhoBlockContribution[iD]/rhoField1.get(iX-offsetX, iY-offsetY, iZ-offsetZ);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ShanChenForcedSingleComponentPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}


/// LatticeCouplingGenerator for NS coupling

template<typename T, typename DESCRIPTOR>
ShanChenForcedSingleComponentGenerator3D<T,DESCRIPTOR>::ShanChenForcedSingleComponentGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_), G(G_), rho0(rho0_), interactionPotential(iP_)
{ }

template<typename T, typename DESCRIPTOR>
ShanChenForcedSingleComponentGenerator3D<T,DESCRIPTOR>::ShanChenForcedSingleComponentGenerator3D (
  T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(0, 0, 0, 0, 0, 0), G(G_), rho0(rho0_), interactionPotential(iP_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* ShanChenForcedSingleComponentGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject3D*> partners) const
{
  return new ShanChenForcedSingleComponentPostProcessor3D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, G, rho0, interactionPotential, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* ShanChenForcedSingleComponentGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ShanChenForcedSingleComponentGenerator3D<T,DESCRIPTOR>(*this);
}



}  // namespace olb

#endif
