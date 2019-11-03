/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Davide Dapelo
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
 * Smagorinsky BGK Dynamics for porous -- generic implementation.
 */
#ifndef SMAGORINSKY_POROUS_PARTICLE_BGK_DYNAMICS_HH
#define SMAGORINSKY_POROUS_PARTICLE_BGK_DYNAMICS_HH

#include "dynamics/porousBGKdynamics.hh"

namespace olb {

////////////////////// Class PorousParticleBGKdynamics //////////////////////////

template<typename T, typename DESCRIPTOR>
SmagorinskyPorousParticleBGKdynamics<T,DESCRIPTOR>::SmagorinskyPorousParticleBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_, T dx_, T dt_ )
  : PorousParticleBGKdynamics<T,DESCRIPTOR>(omega_,momenta_), smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_) )
{ }

template<typename T, typename DESCRIPTOR>
void SmagorinskyPorousParticleBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
/*
  T* porosity = cell.template getFieldPointer<descriptors::POROSITY>();
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= porosity[0];
  }
*/
  T* velDenominator = cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();
  T* velNumerator = cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
  T* porosity = cell.template getFieldPointer<descriptors::POROSITY>();
  if (*velDenominator > std::numeric_limits<T>::epsilon()) {
    *porosity = 1.-*porosity; // 1-prod(1-smoothInd)
    for (int i=0; i < DESCRIPTOR::d; i++)  {
      u[i] += *porosity * (*(velNumerator+i) / *velDenominator - u[i]);
    }
  }
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);

  cell.template defineField<descriptors::POROSITY>(1.0);
  cell.template defineField<descriptors::VELOCITY_DENOMINATOR>(0.0);
  cell.template defineField<descriptors::VELOCITY_NUMERATOR>(0.0);
}

template<typename T, typename DESCRIPTOR>
void SmagorinskyPorousParticleBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  PorousParticleBGKdynamics<T,DESCRIPTOR>::setOmega(omega_);
  preFactor = computePreFactor(omega_, smagoConst);
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyPorousParticleBGKdynamics<T,DESCRIPTOR>::getSmagorinskyOmega(Cell<T,DESCRIPTOR>& cell )
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  return newOmega;
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyPorousParticleBGKdynamics<T,DESCRIPTOR>::computePreFactor(T omega_, T smagoConst_)
{
  return (T)smagoConst_*smagoConst_*descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*2*sqrt(2);
}



template<typename T, typename DESCRIPTOR>
T SmagorinskyPorousParticleBGKdynamics<T,DESCRIPTOR>::computeOmega(T omega0, T preFactor_, T rho,
    T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /omega0;
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + preFactor_/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}

} // olb

#endif
