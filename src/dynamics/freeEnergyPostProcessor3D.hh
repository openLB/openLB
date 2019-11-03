/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Sam Avis, Robin Trunk
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

#ifndef FREE_ENERGY_POST_PROCESSOR_3D_HH
#define FREE_ENERGY_POST_PROCESSOR_3D_HH

#include "freeEnergyPostProcessor3D.h"
#include "core/blockLattice3D.h"

namespace olb {

////////  FreeEnergyChemicalPotentialCoupling3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyChemicalPotentialCoupling3D <T,DESCRIPTOR>::FreeEnergyChemicalPotentialCoupling3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T alpha_, T kappa1_, T kappa2_, T kappa3_,
  std::vector<SpatiallyExtendedObject3D*> partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), alpha(alpha_),
     kappa1(kappa1_), kappa2(kappa2_), kappa3(kappa3_), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyChemicalPotentialCoupling3D <T,DESCRIPTOR>::FreeEnergyChemicalPotentialCoupling3D (
  T alpha_, T kappa1_, T kappa2_, T kappa3_,
  std::vector<SpatiallyExtendedObject3D*> partners_)
  :  x0(0), x1(0), y0(0), y1(0), z0(0), z1(0), alpha(alpha_),
     kappa1(kappa1_), kappa2(kappa2_), kappa3(kappa3_), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
void FreeEnergyChemicalPotentialCoupling3D<T,DESCRIPTOR>::processSubDomain (
  BlockLattice3D<T,DESCRIPTOR>& blockLattice,
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ )
{
  // If partners.size() == 1: two fluid components
  // If partners.size() == 2: three fluid components
  BlockLattice3D<T,DESCRIPTOR> *partnerLattice1 = dynamic_cast<BlockLattice3D<T,DESCRIPTOR> *>(partners[0]);
  BlockLattice3D<T,DESCRIPTOR> *partnerLattice2 = 0;
  if (partners.size() > 1) {
    partnerLattice2 = dynamic_cast<BlockLattice3D<T,DESCRIPTOR> *>(partners[1]);
  }

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect ( x0, x1, y0, y1, z0, z1,
                         x0_, x1_, y0_, y1_, z0, z1,
                         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    int nx = newX1-newX0+3; // include a one-cell boundary
    int ny = newY1-newY0+3; // include a one-cell boundary
    int nz = newZ1-newZ0+3; // include a one-cell boundary
    int offsetX = newX0-1;
    int offsetY = newY0-1;
    int offsetZ = newZ0-1;

    // compute the density fields for each lattice
    BlockData3D<T,T> rhoField1(nx, ny, nz);
    BlockData3D<T,T> rhoField2(nx, ny, nz);
    BlockData3D<T,T> rhoField3(nx, ny, nz);
    for (int iX=newX0-1; iX<=newX1+1; ++iX)
      for (int iY=newY0-1; iY<=newY1+1; ++iY)
        for (int iZ=newZ0-1; iZ<=newZ1+1; ++iZ) {
          rhoField1.get(iX-offsetX, iY-offsetY, iZ-offsetZ) = blockLattice.get(iX,iY,iZ).computeRho();
        }
    for (int iX=newX0-1; iX<=newX1+1; ++iX)
      for (int iY=newY0-1; iY<=newY1+1; ++iY)
        for (int iZ=newZ0-1; iZ<=newZ1+1; ++iZ) {
          rhoField2.get(iX-offsetX, iY-offsetY, iZ-offsetZ) = partnerLattice1->get(iX,iY,iZ).computeRho();
        }
    if (partners.size() > 1) {
      for (int iX=newX0-1; iX<=newX1+1; ++iX)
        for (int iY=newY0-1; iY<=newY1+1; ++iY)
          for (int iZ=newZ0-1; iZ<=newZ1+1; ++iZ) {
            rhoField3.get(iX-offsetX, iY-offsetY, iZ-offsetZ) = partnerLattice2->get(iX,iY,iZ).computeRho();
          }
    }

    // calculate chemical potential
    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          T densitySum = rhoField1.get(iX-offsetX, iY-offsetY, iZ-offsetZ)
                         + rhoField2.get(iX-offsetX, iY-offsetY, iZ-offsetZ);
          T densityDifference = rhoField1.get(iX-offsetX, iY-offsetY, iZ-offsetZ)
                                - rhoField2.get(iX-offsetX, iY-offsetY, iZ-offsetZ);
          if (partners.size() > 1) {
            densitySum -= rhoField3.get(iX-offsetX, iY-offsetY, iZ-offsetZ);
            densityDifference -= rhoField3.get(iX-offsetX, iY-offsetY, iZ-offsetZ);
          }
          T term1 = 0.125 * kappa1 * (densitySum)
                    * (densitySum-1.) * (densitySum-2.);
          T term2 = 0.125 * kappa2 * (densityDifference)
                    * (densityDifference-1.) * (densityDifference-2.);
          T term3 = 0.;
          if (partners.size() > 1) {
            T rho3 = rhoField3.get(iX-offsetX, iY-offsetY, iZ-offsetZ);
            term3 = kappa3 * rho3 * (rho3 - 1.) * (2.*rho3 - 1.);
          }

          T laplaceRho1 = 1.0 / 6.0 * (
                            rhoField1.get(iX-offsetX,   iY-offsetY-1, iZ-offsetZ-1)
                            +      rhoField1.get(iX-offsetX-1, iY-offsetY,   iZ-offsetZ-1)
                            + 2. * rhoField1.get(iX-offsetX,   iY-offsetY,   iZ-offsetZ-1)
                            +      rhoField1.get(iX-offsetX+1, iY-offsetY,   iZ-offsetZ-1)
                            +      rhoField1.get(iX-offsetX,   iY-offsetY+1, iZ-offsetZ-1)
                            +      rhoField1.get(iX-offsetX-1, iY-offsetY-1, iZ-offsetZ)
                            + 2. * rhoField1.get(iX-offsetX,   iY-offsetY-1, iZ-offsetZ)
                            +      rhoField1.get(iX-offsetX+1, iY-offsetY-1, iZ-offsetZ)
                            + 2. * rhoField1.get(iX-offsetX-1, iY-offsetY,   iZ-offsetZ)
                            -24. * rhoField1.get(iX-offsetX,   iY-offsetY,   iZ-offsetZ)
                            + 2. * rhoField1.get(iX-offsetX+1, iY-offsetY,   iZ-offsetZ)
                            +      rhoField1.get(iX-offsetX-1, iY-offsetY+1, iZ-offsetZ)
                            + 2. * rhoField1.get(iX-offsetX,   iY-offsetY+1, iZ-offsetZ)
                            +      rhoField1.get(iX-offsetX+1, iY-offsetY+1, iZ-offsetZ)
                            +      rhoField1.get(iX-offsetX,   iY-offsetY-1, iZ-offsetZ+1)
                            +      rhoField1.get(iX-offsetX-1, iY-offsetY,   iZ-offsetZ+1)
                            + 2. * rhoField1.get(iX-offsetX,   iY-offsetY,   iZ-offsetZ+1)
                            +      rhoField1.get(iX-offsetX+1, iY-offsetY,   iZ-offsetZ+1)
                            +      rhoField1.get(iX-offsetX,   iY-offsetY+1, iZ-offsetZ+1)
                          );

          T laplaceRho2 = 1.0 / 6.0 * (
                            rhoField2.get(iX-offsetX,   iY-offsetY-1, iZ-offsetZ-1)
                            +      rhoField2.get(iX-offsetX-1, iY-offsetY,   iZ-offsetZ-1)
                            + 2. * rhoField2.get(iX-offsetX,   iY-offsetY,   iZ-offsetZ-1)
                            +      rhoField2.get(iX-offsetX+1, iY-offsetY,   iZ-offsetZ-1)
                            +      rhoField2.get(iX-offsetX,   iY-offsetY+1, iZ-offsetZ-1)
                            +      rhoField2.get(iX-offsetX-1, iY-offsetY-1, iZ-offsetZ)
                            + 2. * rhoField2.get(iX-offsetX,   iY-offsetY-1, iZ-offsetZ)
                            +      rhoField2.get(iX-offsetX+1, iY-offsetY-1, iZ-offsetZ)
                            + 2. * rhoField2.get(iX-offsetX-1, iY-offsetY,   iZ-offsetZ)
                            -24. * rhoField2.get(iX-offsetX,   iY-offsetY,   iZ-offsetZ)
                            + 2. * rhoField2.get(iX-offsetX+1, iY-offsetY,   iZ-offsetZ)
                            +      rhoField2.get(iX-offsetX-1, iY-offsetY+1, iZ-offsetZ)
                            + 2. * rhoField2.get(iX-offsetX,   iY-offsetY+1, iZ-offsetZ)
                            +      rhoField2.get(iX-offsetX+1, iY-offsetY+1, iZ-offsetZ)
                            +      rhoField2.get(iX-offsetX,   iY-offsetY-1, iZ-offsetZ+1)
                            +      rhoField2.get(iX-offsetX-1, iY-offsetY,   iZ-offsetZ+1)
                            + 2. * rhoField2.get(iX-offsetX,   iY-offsetY,   iZ-offsetZ+1)
                            +      rhoField2.get(iX-offsetX+1, iY-offsetY,   iZ-offsetZ+1)
                            +      rhoField2.get(iX-offsetX,   iY-offsetY+1, iZ-offsetZ+1)
                          );

          T laplaceRho3 = 0.;
          if (partners.size() > 1) {
            laplaceRho3 = 1.0 / 6.0 * (
                            rhoField3.get(iX-offsetX,   iY-offsetY-1, iZ-offsetZ-1)
                            +      rhoField3.get(iX-offsetX-1, iY-offsetY,   iZ-offsetZ-1)
                            + 2. * rhoField3.get(iX-offsetX,   iY-offsetY,   iZ-offsetZ-1)
                            +      rhoField3.get(iX-offsetX+1, iY-offsetY,   iZ-offsetZ-1)
                            +      rhoField3.get(iX-offsetX,   iY-offsetY+1, iZ-offsetZ-1)
                            +      rhoField3.get(iX-offsetX-1, iY-offsetY-1, iZ-offsetZ)
                            + 2. * rhoField3.get(iX-offsetX,   iY-offsetY-1, iZ-offsetZ)
                            +      rhoField3.get(iX-offsetX+1, iY-offsetY-1, iZ-offsetZ)
                            + 2. * rhoField3.get(iX-offsetX-1, iY-offsetY,   iZ-offsetZ)
                            -24. * rhoField3.get(iX-offsetX,   iY-offsetY,   iZ-offsetZ)
                            + 2. * rhoField3.get(iX-offsetX+1, iY-offsetY,   iZ-offsetZ)
                            +      rhoField3.get(iX-offsetX-1, iY-offsetY+1, iZ-offsetZ)
                            + 2. * rhoField3.get(iX-offsetX,   iY-offsetY+1, iZ-offsetZ)
                            +      rhoField3.get(iX-offsetX+1, iY-offsetY+1, iZ-offsetZ)
                            +      rhoField3.get(iX-offsetX,   iY-offsetY-1, iZ-offsetZ+1)
                            +      rhoField3.get(iX-offsetX-1, iY-offsetY,   iZ-offsetZ+1)
                            + 2. * rhoField3.get(iX-offsetX,   iY-offsetY,   iZ-offsetZ+1)
                            +      rhoField3.get(iX-offsetX+1, iY-offsetY,   iZ-offsetZ+1)
                            +      rhoField3.get(iX-offsetX,   iY-offsetY+1, iZ-offsetZ+1)
                          );
          }

          // setting chemical potential to the respective lattices
          blockLattice.get(iX,iY,iZ).template setField<descriptors::CHEM_POTENTIAL>(
            term1 + term2
            + 0.25*alpha*alpha*( (kappa2 - kappa1) * laplaceRho2
                                 +(kappa2 + kappa1) * (laplaceRho3 - laplaceRho1) )
          );
          partnerLattice1->get(iX, iY, iZ).template setField<descriptors::CHEM_POTENTIAL>(
            term1 - term2
            + 0.25*alpha*alpha*( (kappa2 - kappa1) * (laplaceRho1 - laplaceRho3)
                                 -(kappa2 + kappa1) * laplaceRho2 )
          );
          if (partners.size() > 1) {
            partnerLattice2->get(iX, iY, iZ).template setField<descriptors::CHEM_POTENTIAL>(
              - term1 - term2 + term3
              + 0.25*alpha*alpha*( (kappa2 + kappa1) * laplaceRho1
                                   -(kappa2 - kappa1) * laplaceRho2
                                   -(kappa2 + kappa1 + 4.*kappa3) * laplaceRho3 )
            );
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyChemicalPotentialCoupling3D<T,DESCRIPTOR>::process (
  BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  FreeEnergyForceCoupling3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyForceCoupling3D <T,DESCRIPTOR>::FreeEnergyForceCoupling3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  std::vector<SpatiallyExtendedObject3D*> partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyForceCoupling3D <T,DESCRIPTOR>::FreeEnergyForceCoupling3D (
  std::vector<SpatiallyExtendedObject3D*> partners_)
  :  x0(0), x1(0), y0(0), y1(0), z0(0), z1(0), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
void FreeEnergyForceCoupling3D<T,DESCRIPTOR>::processSubDomain (
  BlockLattice3D<T,DESCRIPTOR>& blockLattice,
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ )
{
  // If partners.size() == 1: two fluid components
  // If partners.size() == 2: three fluid components
  BlockLattice3D<T,DESCRIPTOR> *partnerLattice1 = dynamic_cast<BlockLattice3D<T,DESCRIPTOR> *>(partners[0]);
  BlockLattice3D<T,DESCRIPTOR> *partnerLattice2 = 0;
  if (partners.size() > 1) {
    partnerLattice2 = dynamic_cast<BlockLattice3D<T,DESCRIPTOR> *>(partners[1]);
  }

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect ( x0, x1, y0, y1, z0, z1,
                         x0_, x1_, y0_, y1_, z0_, z1_,
                         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          T phi = blockLattice.get(iX,iY,iZ).computeRho();
          T rho = partnerLattice1->get(iX,iY,iZ).computeRho();

          T gradMuPhiX = 1./12. * ( -blockLattice.get(iX-1,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  blockLattice.get(iX-1,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  blockLattice.get(iX-1,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    - 2.* blockLattice.get(iX-1,iY,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  blockLattice.get(iX-1,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  blockLattice.get(iX+1,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    + 2.* blockLattice.get(iX+1,iY,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  blockLattice.get(iX+1,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  blockLattice.get(iX+1,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  blockLattice.get(iX+1,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>() );
          T gradMuPhiY = 1./12. * ( -blockLattice.get(iX-1,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  blockLattice.get(iX+1,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  blockLattice.get(iX,iY-1,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    - 2.* blockLattice.get(iX,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  blockLattice.get(iX,iY-1,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  blockLattice.get(iX,iY+1,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    + 2.* blockLattice.get(iX,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  blockLattice.get(iX,iY+1,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  blockLattice.get(iX-1,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  blockLattice.get(iX+1,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>() );
          T gradMuPhiZ = 1./12. * ( -blockLattice.get(iX,iY-1,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  blockLattice.get(iX,iY+1,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  blockLattice.get(iX-1,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    - 2.* blockLattice.get(iX,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  blockLattice.get(iX+1,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  blockLattice.get(iX-1,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    + 2.* blockLattice.get(iX,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  blockLattice.get(iX+1,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  blockLattice.get(iX,iY-1,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  blockLattice.get(iX,iY+1,iZ+1).template getField<descriptors::CHEM_POTENTIAL>() );

          T gradMuRhoX = 1./12. * ( -partnerLattice1->get(iX-1,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice1->get(iX-1,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice1->get(iX-1,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    - 2.* partnerLattice1->get(iX-1,iY,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice1->get(iX-1,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice1->get(iX+1,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    + 2.* partnerLattice1->get(iX+1,iY,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice1->get(iX+1,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice1->get(iX+1,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice1->get(iX+1,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>() );
          T gradMuRhoY = 1./12. * ( -partnerLattice1->get(iX-1,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice1->get(iX+1,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice1->get(iX,iY-1,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    - 2.* partnerLattice1->get(iX,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice1->get(iX,iY-1,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice1->get(iX,iY+1,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    + 2.* partnerLattice1->get(iX,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice1->get(iX,iY+1,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice1->get(iX-1,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice1->get(iX+1,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>() );
          T gradMuRhoZ = 1./12. * ( -partnerLattice1->get(iX,iY-1,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice1->get(iX,iY+1,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice1->get(iX-1,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    - 2.* partnerLattice1->get(iX,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice1->get(iX+1,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice1->get(iX-1,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    + 2.* partnerLattice1->get(iX,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice1->get(iX+1,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice1->get(iX,iY-1,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice1->get(iX,iY+1,iZ+1).template getField<descriptors::CHEM_POTENTIAL>() );

          T psi = 0.;
          T gradMuPsiX = 0.;
          T gradMuPsiY = 0.;
          T gradMuPsiZ = 0.;
          if (partners.size() > 1) {
            psi = partnerLattice2->get(iX,iY,iZ).computeRho();
            gradMuPsiX = 1./12. * ( -partnerLattice2->get(iX-1,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice2->get(iX-1,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice2->get(iX-1,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    - 2.* partnerLattice2->get(iX-1,iY,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice2->get(iX-1,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice2->get(iX+1,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    + 2.* partnerLattice2->get(iX+1,iY,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice2->get(iX+1,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice2->get(iX+1,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice2->get(iX+1,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>() );
            gradMuPsiY = 1./12. * ( -partnerLattice2->get(iX-1,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice2->get(iX+1,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice2->get(iX,iY-1,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    - 2.* partnerLattice2->get(iX,iY-1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice2->get(iX,iY-1,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice2->get(iX,iY+1,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    + 2.* partnerLattice2->get(iX,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice2->get(iX,iY+1,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice2->get(iX-1,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice2->get(iX+1,iY+1,iZ  ).template getField<descriptors::CHEM_POTENTIAL>() );
            gradMuPsiZ = 1./12. * ( -partnerLattice2->get(iX,iY-1,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice2->get(iX,iY+1,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice2->get(iX-1,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    - 2.* partnerLattice2->get(iX,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    -  partnerLattice2->get(iX+1,iY,iZ-1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice2->get(iX-1,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    + 2.* partnerLattice2->get(iX,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice2->get(iX+1,iY,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice2->get(iX,iY-1,iZ+1).template getField<descriptors::CHEM_POTENTIAL>()
                                    +  partnerLattice2->get(iX,iY+1,iZ+1).template getField<descriptors::CHEM_POTENTIAL>() );
          }

          T forceX = -rho*gradMuRhoX - phi*gradMuPhiX - psi*gradMuPsiX;
          T forceY = -rho*gradMuRhoY - phi*gradMuPhiY - psi*gradMuPsiY;
          T forceZ = -rho*gradMuRhoZ - phi*gradMuPhiZ - psi*gradMuPsiZ;
          partnerLattice1->get(iX,iY,iZ).template setField<descriptors::FORCE>({forceX, forceY, forceZ});
          T u[3];
          partnerLattice1->get(iX,iY,iZ).computeU(u);
          blockLattice.get(iX,iY,iZ).template setField<descriptors::FORCE>(u);
          if (partners.size() > 1) {
            partnerLattice2->get(iX,iY,iZ).template setField<descriptors::FORCE>(u);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyForceCoupling3D<T,DESCRIPTOR>::process (
  BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}


////////  FreeEnergyInletOutletCoupling3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyInletOutletCoupling3D <T,DESCRIPTOR>::FreeEnergyInletOutletCoupling3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  std::vector<SpatiallyExtendedObject3D*> partners_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyInletOutletCoupling3D <T,DESCRIPTOR>::FreeEnergyInletOutletCoupling3D (
  std::vector<SpatiallyExtendedObject3D*> partners_)
  : x0(0), x1(0), y0(0), y1(0), z0(0), z1(0), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
void FreeEnergyInletOutletCoupling3D<T,DESCRIPTOR>::processSubDomain (
  BlockLattice3D<T,DESCRIPTOR>& blockLattice,
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ )
{
  // If partners.size() == 1: two fluid components
  // If partners.size() == 2: three fluid components
  BlockLattice3D<T,DESCRIPTOR> *partnerLattice1 = dynamic_cast<BlockLattice3D<T,DESCRIPTOR> *>(partners[0]);
  BlockLattice3D<T,DESCRIPTOR> *partnerLattice2 = 0;
  if (partners.size() > 1) {
    partnerLattice2 = dynamic_cast<BlockLattice3D<T,DESCRIPTOR> *>(partners[1]);
  }

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect ( x0, x1, y0, y1, z0, z1,
                         x0_, x1_, y0_, y1_, z0_, z1_,
                         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          T u[DESCRIPTOR::d];
          partnerLattice1->get(iX,iY,iZ).computeU(u);
          blockLattice.get(iX,iY,iZ).defineU(u);
          if (partners.size() > 1) {
            partnerLattice2->get(iX,iY,iZ).defineU(u);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyInletOutletCoupling3D<T,DESCRIPTOR>::process(
  BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}


////////  FreeEnergyDensityOutletCoupling3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyDensityOutletCoupling3D <T,DESCRIPTOR>::FreeEnergyDensityOutletCoupling3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T rho_,
  std::vector<SpatiallyExtendedObject3D*> partners_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
    rho(rho_), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyDensityOutletCoupling3D <T,DESCRIPTOR>::FreeEnergyDensityOutletCoupling3D (
  T rho_, std::vector<SpatiallyExtendedObject3D*> partners_)
  : x0(0), x1(0), y0(0), y1(0), z0(0), z1(0), rho(rho_), partners(partners_)
{ }

template<typename T, typename DESCRIPTOR>
void FreeEnergyDensityOutletCoupling3D<T,DESCRIPTOR>::processSubDomain (
  BlockLattice3D<T,DESCRIPTOR>& blockLattice,
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ )
{
  // If partners.size() == 1: two fluid components
  // If partners.size() == 2: three fluid components
  BlockLattice3D<T,DESCRIPTOR> *partnerLattice1 = dynamic_cast<BlockLattice3D<T,DESCRIPTOR> *>(partners[0]);
  BlockLattice3D<T,DESCRIPTOR> *partnerLattice2 = 0;
  if (partners.size() > 1) {
    partnerLattice2 = dynamic_cast<BlockLattice3D<T,DESCRIPTOR> *>(partners[1]);
  }

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect ( x0, x1, y0, y1, z0, z1,
                         x0_, x1_, y0_, y1_, z0_, z1_,
                         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {

          T rho0, phi, psi;
          rho0 = blockLattice.get(iX,iY,iZ).computeRho();
          phi = partnerLattice1->get(iX,iY,iZ).computeRho();
          blockLattice.get(iX,iY,iZ).defineRho(rho);
          partnerLattice1->get(iX,iY,iZ).defineRho(phi * rho / rho0);
          if (partners.size() > 1) {
            psi = partnerLattice2->get(iX,iY,iZ).computeRho();
            partnerLattice2->get(iX,iY,iZ).defineRho(psi * rho / rho0);
          }

        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyDensityOutletCoupling3D<T,DESCRIPTOR>::process(
  BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}


////////  FreeEnergyChemicalPotentialGenerator3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyChemicalPotentialGenerator3D<T,DESCRIPTOR>::FreeEnergyChemicalPotentialGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  T alpha_, T kappa1_, T kappa2_ )
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_), alpha(alpha_),
    kappa1(kappa1_), kappa2(kappa2_), kappa3(0)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyChemicalPotentialGenerator3D<T,DESCRIPTOR>::FreeEnergyChemicalPotentialGenerator3D (
  T alpha_, T kappa1_, T kappa2_ )
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(0, 0, 0, 0, 0, 0), alpha(alpha_),
    kappa1(kappa1_), kappa2(kappa2_), kappa3(0)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyChemicalPotentialGenerator3D<T,DESCRIPTOR>::FreeEnergyChemicalPotentialGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  T alpha_, T kappa1_, T kappa2_, T kappa3_ )
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_), alpha(alpha_),
    kappa1(kappa1_), kappa2(kappa2_), kappa3(kappa3_)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyChemicalPotentialGenerator3D<T,DESCRIPTOR>::FreeEnergyChemicalPotentialGenerator3D (
  T alpha_, T kappa1_, T kappa2_, T kappa3_ )
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(0, 0, 0, 0, 0, 0), alpha(alpha_),
    kappa1(kappa1_), kappa2(kappa2_), kappa3(kappa3_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* FreeEnergyChemicalPotentialGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject3D*> partners) const
{
  return new FreeEnergyChemicalPotentialCoupling3D<T,DESCRIPTOR>(
           this->x0, this->x1, this->y0, this->y1, this->z0, this->z1,
           alpha, kappa1, kappa2, kappa3, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>*
FreeEnergyChemicalPotentialGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyChemicalPotentialGenerator3D<T,DESCRIPTOR>(*this);
}

////////  FreeEnergyForceGenerator3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyForceGenerator3D<T,DESCRIPTOR>::FreeEnergyForceGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyForceGenerator3D<T,DESCRIPTOR>::FreeEnergyForceGenerator3D ( )
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(0, 0, 0, 0, 0, 0)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* FreeEnergyForceGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject3D*> partners) const
{
  return new FreeEnergyForceCoupling3D<T,DESCRIPTOR>(
           this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* FreeEnergyForceGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyForceGenerator3D<T,DESCRIPTOR>(*this);
}

////////  FreeEnergyInletOutletGenerator3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyInletOutletGenerator3D<T,DESCRIPTOR>::FreeEnergyInletOutletGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyInletOutletGenerator3D<T,DESCRIPTOR>::FreeEnergyInletOutletGenerator3D ( )
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(0, 0, 0, 0, 0, 0)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* FreeEnergyInletOutletGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject3D*> partners) const
{
  return new FreeEnergyInletOutletCoupling3D<T,DESCRIPTOR>(
           this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* FreeEnergyInletOutletGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyInletOutletGenerator3D<T,DESCRIPTOR>(*this);
}

////////  FreeEnergyDensityOutletGenerator3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyDensityOutletGenerator3D<T,DESCRIPTOR>::FreeEnergyDensityOutletGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T rho_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), rho(rho_)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyDensityOutletGenerator3D<T,DESCRIPTOR>::FreeEnergyDensityOutletGenerator3D (
  T rho_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(0, 0, 0, 0, 0, 0), rho(rho_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* FreeEnergyDensityOutletGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject3D*> partners) const
{
  return new FreeEnergyDensityOutletCoupling3D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, rho, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* FreeEnergyDensityOutletGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyDensityOutletGenerator3D<T,DESCRIPTOR>(*this);
}


}  // namespace olb

#endif
