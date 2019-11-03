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

#ifndef SHAN_CHEN_FORCED_POST_PROCESSOR_2D_HH
#define SHAN_CHEN_FORCED_POST_PROCESSOR_2D_HH

#include "shanChenForcedPostProcessor2D.h"
#include "interactionPotential.h"
#include "core/blockLattice2D.h"
#include "core/util.h"
#include "core/finiteDifference2D.h"

namespace olb {

////////  ShanChenForcedPostProcessor2D ///////////////////////////////////


template<typename T, typename DESCRIPTOR>
ShanChenForcedPostProcessor2D <T,DESCRIPTOR>::
ShanChenForcedPostProcessor2D(int x0_, int x1_, int y0_, int y1_, T G_,
                              std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_,
                              std::vector<SpatiallyExtendedObject2D*> partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), G(G_), rho0(rho0_), interactionPotential(iP_), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
ShanChenForcedPostProcessor2D <T,DESCRIPTOR>::
ShanChenForcedPostProcessor2D(T G_,
                              std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_,
                              std::vector<SpatiallyExtendedObject2D*> partners_)
  :  x0(0), x1(0), y0(0), y1(0), G(G_), rho0(rho0_), interactionPotential(iP_), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
void ShanChenForcedPostProcessor2D<T,DESCRIPTOR>::
processSubDomain( BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                  int x0_, int x1_, int y0_, int y1_ )
{
  typedef DESCRIPTOR L;

  BlockLattice2D<T,DESCRIPTOR> *partnerLattice = dynamic_cast<BlockLattice2D<T,DESCRIPTOR> *>(partners[0]);

  int newX0, newX1, newY0, newY1;
  if ( util::intersect ( x0, x1, y0, y1,
                         x0_, x1_, y0_, y1_,
                         newX0, newX1, newY0, newY1 ) ) {
    int nx = newX1-newX0+3; // include a one-cell boundary
    int ny = newY1-newY0+3; // include a one-cell boundary
    int offsetX = newX0-1;
    int offsetY = newY0-1;

    BlockData2D<T,T> rhoField1(nx, ny);
    BlockData2D<T,T> rhoField2(nx, ny);

    // Compute density and velocity on every site of first lattice, and store result
    //   in external scalars; envelope cells are included, because they are needed
    //   to compute the interaction potential in what follows.
    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {
        Cell<T,DESCRIPTOR>& cell = blockLattice.get(iX,iY);
        rhoField1.get(iX-offsetX, iY-offsetY) = cell.computeRho()*rho0[0];
      }
    }

    // Compute density and velocity on every site of second lattice, and store result
    //   in external scalars; envelope cells are included, because they are needed
    //   to compute the interaction potential in what follows.
    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {
        Cell<T,DESCRIPTOR>& cell = partnerLattice->get(iX,iY);
        rhoField2.get(iX-offsetX, iY-offsetY) = cell.computeRho()*rho0[1];
      }
    }

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        Cell<T,DESCRIPTOR>& blockCell   = blockLattice.get(iX,iY);
        Cell<T,DESCRIPTOR>& partnerCell = partnerLattice->get(iX,iY);

        T* j = blockCell.template getFieldPointer<descriptors::VELOCITY>();
        lbHelpers<T,DESCRIPTOR>::computeJ(blockCell,j);
        j = partnerCell.template getFieldPointer<descriptors::VELOCITY>();
        lbHelpers<T,DESCRIPTOR>::computeJ(partnerCell,j);

        T blockOmega   = blockLattice.getDynamics(iX, iY)->getOmega();
        T partnerOmega = partnerLattice->getDynamics(iX, iY)->getOmega();
        // Computation of the common velocity, shared among the two populations
        T rhoTot = rhoField1.get(iX-offsetX, iY-offsetY)*blockOmega +
                   rhoField2.get(iX-offsetX, iY-offsetY)*partnerOmega;

        T uTot[DESCRIPTOR::d];
        T* blockU = blockCell.template getFieldPointer<descriptors::VELOCITY>();      // contains precomputed value rho*u
        T* partnerU = partnerCell.template getFieldPointer<descriptors::VELOCITY>();  // contains precomputed value rho*u
        for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
          uTot[iD] = (blockU[iD]*rho0[0]*blockOmega + partnerU[iD]*rho0[1]*partnerOmega) / rhoTot;
        }

        // Computation of the interaction potential
        T rhoBlockContribution[L::d]   = {T(), T()};
        T rhoPartnerContribution[L::d] = {T(), T()};
        T psi2;
        T psi1;
        interactionPotential(&psi2, &rhoField2.get(iX-offsetX, iY-offsetY));
        interactionPotential(&psi1, &rhoField1.get(iX-offsetX, iY-offsetY));
        for (int iPop = 0; iPop < L::q; ++iPop) {
          int nextX = iX + descriptors::c<L>(iPop,0);
          int nextY = iY + descriptors::c<L>(iPop,1);
          T blockRho;
          T partnerRho;
          interactionPotential(&blockRho, &rhoField1.get(nextX-offsetX, nextY-offsetY));//rho0[0];
          interactionPotential(&partnerRho, &rhoField2.get(nextX-offsetX, nextY-offsetY));///rho0[1];
          for (int iD = 0; iD < L::d; ++iD) {
            rhoBlockContribution[iD]   += psi2 * blockRho * descriptors::c<L>(iPop,iD)* descriptors::t<T,L>(iPop);
            rhoPartnerContribution[iD] += psi1 * partnerRho * descriptors::c<L>(iPop,iD)* descriptors::t<T,L>(iPop);
          }
        }

        // Computation and storage of the final velocity, consisting
        //   of u and the momentum difference due to interaction
        //   potential plus external force
        T *blockForce   = blockCell.template getFieldPointer<descriptors::FORCE>();
        T *partnerForce = partnerCell.template getFieldPointer<descriptors::FORCE>();
        T *externalBlockForce   = blockCell.template getFieldPointer<descriptors::EXTERNAL_FORCE>();
        T *externalPartnerForce = partnerCell.template getFieldPointer<descriptors::EXTERNAL_FORCE>();

        for (int iD = 0; iD < L::d; ++iD) {
          blockU[iD] = uTot[iD];
          blockForce[iD] = externalBlockForce[iD] - G*rhoPartnerContribution[iD]/rhoField1.get(iX-offsetX, iY-offsetY);
          partnerU[iD] = uTot[iD];
          partnerForce[iD] = externalPartnerForce[iD] - G*rhoBlockContribution[iD]/rhoField2.get(iX-offsetX, iY-offsetY);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ShanChenForcedPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}


/// LatticeCouplingGenerator for NS coupling

template<typename T, typename DESCRIPTOR>
ShanChenForcedGenerator2D<T,DESCRIPTOR>::ShanChenForcedGenerator2D (
  int x0_, int x1_, int y0_, int y1_, T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_ )
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), G(G_), rho0(rho0_), interactionPotential(iP_)
{ }

template<typename T, typename DESCRIPTOR>
ShanChenForcedGenerator2D<T,DESCRIPTOR>::ShanChenForcedGenerator2D (
  T G_, std::vector<T> rho0_, AnalyticalF1D<T,T>& iP_ )
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(0, 0, 0, 0), G(G_), rho0(rho0_), interactionPotential(iP_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* ShanChenForcedGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject2D*> partners) const
{
  return new ShanChenForcedPostProcessor2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1,G, rho0, interactionPotential, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* ShanChenForcedGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new ShanChenForcedGenerator2D<T,DESCRIPTOR>(*this);
}



}  // namespace olb

#endif
