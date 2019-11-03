/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2015 Mathias J. Krause, Vojtech Cvrcekt, Davide Dapelo
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

/** \file
 * Porous-particle BGK Dynamics with adjusted omega
 * and Smagorinsky turbulence model -- generic implementation.
 * Strain rate similar to "J.Boyd, J. Buick and S.Green: A second-order accurate lattice Boltzmann non-Newtonian flow model"
 * Power Law similar to "Huidan Yu, Sharath S. Girimaji, Li-Shi Luo - DNS and LES of decaying isotropic turbulence with and without frame rotation using lattice Boltzmann method"
 */
#ifndef SMAGORINSKY_POWER_LAW_BGK_DYNAMICS_HH
#define SMAGORINSKY_POWER_LAW_BGK_DYNAMICS_HH

#include "../dynamics/powerLawBGKdynamics.h"
#include "SmagorinskyPowerLawBGKdynamics.h"
#include "math.h"

namespace olb {

////////////////////// Class SmagorinskyPowerLawBGKdynamics //////////////////////////

/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
SmagorinskyPowerLawBGKdynamics<T,DESCRIPTOR>::SmagorinskyPowerLawBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta, T m, T n , T nuMin, T nuMax, T smagoConst)
  : SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega, momenta, smagoConst),
    PowerLawDynamics<T,DESCRIPTOR>(m, n, nuMin, nuMax)
{ }

template<typename T, typename DESCRIPTOR>
void SmagorinskyPowerLawBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);

  // Computation of the power-law omega.
  // An external is used in place of BGKdynamics::_omega to keep generality and flexibility.
  T oldOmega = cell.template getFieldPointer<descriptors::OMEGA>()[0];
  T intOmega = this->computeOmegaPL(cell, oldOmega, rho, pi);
  T newOmega = computeEffectiveOmega(cell, intOmega); // turbulent omega

  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  cell.template getFieldPointer<descriptors::OMEGA>()[0] = intOmega; // updating omega
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyPowerLawBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell, T omega0)
{
  T rho = this->_momenta.computeRho(cell);
  T PiNeqNorm    = sqrt(PiNeqNormSqr(cell));
  /// Molecular realaxation time
  T tau_mol = 1. /omega0;
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + this->getPreFactor()/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  T tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyPowerLawBGKdynamics<T,DESCRIPTOR>::PiNeqNormSqr(Cell<T,DESCRIPTOR>& cell )
{
  return lbHelpers<T,DESCRIPTOR>::computePiNeqNormSqr(cell);
}


////////////////////// Class SmagorinskyPowerLawForcedBGKdynamics //////////////////////////

/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
SmagorinskyPowerLawForcedBGKdynamics<T,DESCRIPTOR>::SmagorinskyPowerLawForcedBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta, T m, T n , T nuMin, T nuMax, T smagoConst)
  : SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>(omega, momenta, smagoConst),
    PowerLawDynamics<T,DESCRIPTOR>(m, n, nuMin, nuMax)
{ }

template<typename T, typename DESCRIPTOR>
void SmagorinskyPowerLawForcedBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);

  // Computation of the power-law omega.
  // An external is used in place of BGKdynamics::_omega to keep generality and flexibility.
  T oldOmega = cell.template getFieldPointer<descriptors::OMEGA>()[0];
  T intOmega = this->computeOmegaPL(cell, oldOmega, rho, pi);
  T newOmega = computeEffectiveOmega(cell, intOmega); // turbulent omega

  T* force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }

  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  cell.template getFieldPointer<descriptors::OMEGA>()[0] = intOmega; // updating omega
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyPowerLawForcedBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell, T omega0)
{
  T rho = this->_momenta.computeRho(cell);
  T PiNeqNorm    = sqrt(PiNeqNormSqr(cell));
  /// Molecular realaxation time
  T tau_mol = 1. /omega0;
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + this->getPreFactor()/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  T tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyPowerLawForcedBGKdynamics<T,DESCRIPTOR>::PiNeqNormSqr(Cell<T,DESCRIPTOR>& cell )
{
  return lbHelpers<T,DESCRIPTOR>::computeForcedPiNeqNormSqr(cell);
}

}

#endif
