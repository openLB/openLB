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

#ifndef NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_2D_HH
#define NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_2D_HH

#include "latticeDescriptors.h"
#include "navierStokesAdvectionDiffusionCouplingPostProcessor2D.h"
#include "core/blockLattice2D.h"
#include "core/util.h"
#include "core/finiteDifference2D.h"

using namespace std;

namespace olb {

//=====================================================================================
//==============  NavierStokesAdvectionDiffusionCouplingPostProcessor2D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR>
NavierStokesAdvectionDiffusionCouplingPostProcessor2D<T,DESCRIPTOR>::
NavierStokesAdvectionDiffusionCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
    std::vector<SpatiallyExtendedObject2D* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_),
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

  tPartner = dynamic_cast<BlockLattice2D<T,descriptors::D2Q5<descriptors::VELOCITY>> *>(partners[0]);
}

template<typename T, typename DESCRIPTOR>
void NavierStokesAdvectionDiffusionCouplingPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_)
{

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        // computation of the bousinessq force
        T *force = blockLattice.get(iX,iY).template getFieldPointer<descriptors::FORCE>();
        T temperatureDifference = tPartner->get(iX,iY).computeRho() - T0;
        for (unsigned iD = 0; iD < L::d; ++iD) {
          force[iD] = forcePrefactor[iD] * temperatureDifference;
        }
        // Velocity coupling
        T *u = tPartner->get(iX,iY).template getFieldPointer<descriptors::VELOCITY>();
        blockLattice.get(iX,iY).computeU(u);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void NavierStokesAdvectionDiffusionCouplingPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
NavierStokesAdvectionDiffusionCouplingGenerator2D<T,DESCRIPTOR>::
NavierStokesAdvectionDiffusionCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* NavierStokesAdvectionDiffusionCouplingGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject2D* > partners) const
{
  return new NavierStokesAdvectionDiffusionCouplingPostProcessor2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, gravity, T0, deltaTemp, dir,partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* NavierStokesAdvectionDiffusionCouplingGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new NavierStokesAdvectionDiffusionCouplingGenerator2D<T,DESCRIPTOR>(*this);
}


//=====================================================================================
//==============  SmagorinskyBoussinesqCouplingPostProcessor2D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR>
SmagorinskyBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>::
SmagorinskyBoussinesqCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_,
    std::vector<SpatiallyExtendedObject2D* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_),
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

  tauTurbADPrefactor = descriptors::invCs2<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF>>() / descriptors::invCs2<T,DESCRIPTOR>() / PrTurb;
  tPartner = dynamic_cast<BlockLattice2D<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF>> *>(partners[0]);
}

template<typename T, typename DESCRIPTOR>
void SmagorinskyBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_)
{

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {

         // computation of the bousinessq force
        T *force = blockLattice.get(iX,iY).template getFieldPointer<descriptors::FORCE>();
        T temperatureDifference = tPartner->get(iX,iY).computeRho() - T0;
        for (unsigned iD = 0; iD < L::d; ++iD) {
          force[iD] = forcePrefactor[iD] * temperatureDifference;
        }

        // Velocity coupling
        T *u = tPartner->get(iX,iY).template getFieldPointer<descriptors::VELOCITY>();

        // tau coupling
        T *tauNS = blockLattice.get(iX,iY).template getFieldPointer<descriptors::TAU_EFF>();
        T *tauAD = tPartner->get(iX,iY).template getFieldPointer<descriptors::TAU_EFF>();

        T rho, pi[util::TensorVal<DESCRIPTOR >::n];
        blockLattice.get(iX,iY).computeAllMomenta(rho, u, pi);
        T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
        if (util::TensorVal<DESCRIPTOR >::n == 6) {
          PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
        }
        T PiNeqNorm    = sqrt(PiNeqNormSqr);
        /// Molecular realaxation time
        T tau_mol_NS = 1. / blockLattice.get(iX,iY).getDynamics()->getOmega();
        T tau_mol_AD = 1. / tPartner->get(iX,iY).getDynamics()->getOmega();
        /// Turbulent realaxation time
        T tau_turb_NS = 0.5*(sqrt(tau_mol_NS*tau_mol_NS + dynamic_cast<SmagorinskyDynamics<T,DESCRIPTOR>*>(blockLattice.get(iX,iY).getDynamics())->getPreFactor()/rho*PiNeqNorm) - tau_mol_NS);
        /// Effective realaxation time
        tauNS[0] = tau_mol_NS+tau_turb_NS;

        T tau_turb_AD = tau_turb_NS * tauTurbADPrefactor;
        tauAD[0] = tau_mol_AD+tau_turb_AD;
      }
    }
  }

}

template<typename T, typename DESCRIPTOR>
void SmagorinskyBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
SmagorinskyBoussinesqCouplingGenerator2D<T,DESCRIPTOR>::
SmagorinskyBoussinesqCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_), PrTurb(PrTurb_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* SmagorinskyBoussinesqCouplingGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject2D* > partners) const
{
  return new SmagorinskyBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, gravity, T0, deltaTemp, dir, PrTurb, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* SmagorinskyBoussinesqCouplingGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new SmagorinskyBoussinesqCouplingGenerator2D<T,DESCRIPTOR>(*this);
}


//=====================================================================================
//==============  MixedScaleBoussinesqCouplingPostProcessor2D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR>
MixedScaleBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>::
MixedScaleBoussinesqCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_,
    std::vector<SpatiallyExtendedObject2D* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_),
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

  tauTurbADPrefactor = descriptors::invCs2<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF,descriptors::CUTOFF_HEAT_FLUX>>() / descriptors::invCs2<T,DESCRIPTOR>() / PrTurb;
  tPartner = dynamic_cast<BlockLattice2D<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF,descriptors::CUTOFF_HEAT_FLUX>> *>(partners[0]);
}

template<typename T, typename DESCRIPTOR>
void MixedScaleBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_)
{

  const T C_nu = 0.04;
  const T C_alpha = 0.5;
  const T deltaT = 1.0;

  const T invCs2_g = descriptors::invCs2<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF,descriptors::CUTOFF_HEAT_FLUX>>();

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        T *cutoffKinEnergy_14 = blockLattice.get(iX, iY).template getFieldPointer<descriptors::CUTOFF_KIN_ENERGY>();
        T *cutoffHeatFlux_14 = tPartner->get(iX, iY).template getFieldPointer<descriptors::CUTOFF_HEAT_FLUX>();

        // Velocity coupling
        T *u = tPartner->get(iX,iY).template getFieldPointer<descriptors::VELOCITY>();

        // tau coupling
        T *tauNS = blockLattice.get(iX,iY).template getFieldPointer<descriptors::TAU_EFF>();
        T *tauAD = tPartner->get(iX,iY).template getFieldPointer<descriptors::TAU_EFF>();

        /// Molecular realaxation time
        T tau_mol_NS = 1. / blockLattice.get(iX,iY).getDynamics()->getOmega();
        T tau_mol_AD = 1. / tPartner->get(iX,iY).getDynamics()->getOmega();

        const T temperature = tPartner->get(iX,iY).computeRho();

        // computation of the bousinessq force
        T *force = blockLattice.get(iX,iY).template getFieldPointer<descriptors::FORCE>();
        T temperatureDifference = temperature - T0;
        for (unsigned iD = 0; iD < L::d; ++iD) {
          force[iD] = forcePrefactor[iD] * temperatureDifference;
        }

        T rho, pi[util::TensorVal<DESCRIPTOR>::n], j[DESCRIPTOR::d];
        blockLattice.get(iX,iY).computeAllMomenta(rho, u, pi);

        int iPi = 0;
        for (int Alpha=0; Alpha<DESCRIPTOR::d; ++Alpha) {
          for (int Beta=Alpha; Beta<DESCRIPTOR::d; ++Beta) {
            pi[iPi] += rho/2.*(force[Alpha]*u[Beta] + u[Alpha]*force[Beta]);
            ++iPi;
          }
        }
        const T piSqr[3] = {pi[0]*pi[0], pi[1]*pi[1], pi[2]*pi[2]};
        const T PiNeqNormSqr = piSqr[0] + 2.0*piSqr[1] + piSqr[2];
        const T PiNeqNorm    = sqrt(PiNeqNormSqr);

        tPartner->get(iX,iY).computeJ(j);
        const T tmp_preFactor = invCs2_g / rho / tauAD[0];
        const T jNeq[2] = {(j[0] - temperature * u[0]), (j[1] - temperature * u[1])};
        const T jNeqSqr[2] = {jNeq[0]*jNeq[0], jNeq[1]*jNeq[1]};
        const T jNeqSqr_prefacor = 2. * 0.25 * (jNeq[0] + jNeq[1]) * (jNeq[0] + jNeq[1]);

        const T TnormSqr = jNeqSqr_prefacor*PiNeqNormSqr;
        const T Tnorm    = sqrt(TnormSqr);

        /// Turbulent realaxation time
        // T tau_turb_NS = 0.5*(sqrt(tau_mol_NS*tau_mol_NS + dynamic_cast<SmagorinskyDynamics<T,DESCRIPTOR>*>(blockLattice.get(iX,iY).getDynamics())->getPreFactor()/rho*PiNeqNorm) - tau_mol_NS);

        // const T tmp_A = C_nu * sqrt(sqrt(2.)/2.) * descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>() * sqrt(PiNeqNorm / rho) * cutoffKinEnergy_14[0];
        // const T tmp_A_2 = tmp_A * tmp_A;
        // const T tmp_A_4 = tmp_A_2 * tmp_A_2;

        // const T tau_mol_NS_2 = tau_mol_NS * tau_mol_NS;
        // const T tau_mol_NS_3 = tau_mol_NS_2 * tau_mol_NS;

        // const T tmp_1_3 = 1./3.;
        // const T tmp_2_13 = pow(2., tmp_1_3);
        // const T tmp_3_3_12 = 3. * sqrt(3.);

        // const T tmp_sqrtA = sqrt(27.*tmp_A_4-4.*tmp_A_2*tau_mol_NS_3);

        //   // T tau_turb_NS = 1/3 ((27 A^2 + 3 sqrt(3) sqrt(27 A^4 - 4 A^2 b^3) - 2 b^3)^(1/3)/2^(1/3) + (2^(1/3) b^2)/(27 A^2 + 3 sqrt(3) sqrt(27 A^4 - 4 A^2 b^3) - 2 b^3)^(1/3) - b)
        // T tau_turb_NS = ( pow(27.*tmp_A_2 + tmp_3_3_12*sqrt(27.*tmp_A_4-4.*tmp_A_2*tau_mol_NS_3)-2.*tau_mol_NS_3, tmp_1_3) / tmp_2_13
        //                    + (tmp_2_13*tau_mol_NS_2) / pow(27.*tmp_A_2+tmp_3_3_12*sqrt(27.*tmp_A_4-4.*tmp_A_2*tau_mol_NS_3) - 2.*tau_mol_NS_3, tmp_1_3)
        //                    - tau_mol_NS
        //                    ) * tmp_1_3;

        // if ( tau_turb_NS != tau_turb_NS )
        //   tau_turb_NS = 0.;

        //cout << tau_turb_NS << " " << 27. * tmp_A_2 << " " << 4. * tau_mol_NS_3 << " " << PiNeqNorm << " " << " " << rho << endl;

        const T tmp_A = C_nu * sqrt(sqrt(2.)/2.) * descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>() * sqrt(PiNeqNorm / rho / tauNS[0]) * cutoffKinEnergy_14[0];
        const T tau_turb_NS = tmp_A;

        // T tau_turb_AD = tau_turb_NS * tauTurbADPrefactor;
        const T tmp_B = C_alpha * descriptors::invCs2<T,DESCRIPTOR>() / rho * sqrt(2.0 * Tnorm * invCs2_g / tauNS[0] / tauAD[0]) * cutoffHeatFlux_14[0];
        const T tau_turb_AD = tmp_B;
        // cout << jNeq[0] << " " << jNeq[1] << " " << sqrt(Tnorm * invCs2_g / tauNS[0] / tauAD[0]) << " " << TnormSqr << endl;

        /// Effective realaxation time
        tauNS[0] = tau_mol_NS+tau_turb_NS;
        tauAD[0] = tau_mol_AD+tau_turb_AD;

      }
    }
  }

}

template<typename T, typename DESCRIPTOR>
void MixedScaleBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
MixedScaleBoussinesqCouplingGenerator2D<T,DESCRIPTOR>::
MixedScaleBoussinesqCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_), PrTurb(PrTurb_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* MixedScaleBoussinesqCouplingGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<SpatiallyExtendedObject2D* > partners) const
{
  return new MixedScaleBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, gravity, T0, deltaTemp, dir, PrTurb, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* MixedScaleBoussinesqCouplingGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new MixedScaleBoussinesqCouplingGenerator2D<T,DESCRIPTOR>(*this);
}

}  // namespace olb

#endif
