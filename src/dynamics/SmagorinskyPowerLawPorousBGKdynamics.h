/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2015 Mathias J. Krause, Vojtech Cvrcek, Davide Dapelo
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
 * and Smagorinsky turbulence model -- header file.
 * Strain rate similar to "J.Boyd, J. Buick and S.Green: A second-order accurate lattice Boltzmann non-Newtonian flow model"
 * Power Law similar to "Huidan Yu, Sharath S. Girimaji, Li-Shi Luo - DNS and LES of decaying isotropic turbulence with and without frame rotation using lattice Boltzmann method"
 *
 */
#ifndef SMAGORINSKY_POWER_LAW_POROUS_BGK_DYNAMICS_H
#define SMAGORINSKY_POWER_LAW_POROUS_BGK_DYNAMICS_H

#include "SmagorinskyPorousParticleBGKdynamics.h"
#include "core/cell.h"

namespace olb {

/// Implementation of the BGK collision step
template<typename T, typename DESCRIPTOR>
class SmagorinskyPowerLawPorousParticleBGKdynamics : public SmagorinskyPorousParticleBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  /// m,n...parameter in the power law model, dt...the explizit typed time step from LBM model
  SmagorinskyPowerLawPorousParticleBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T m_=0.1, T n_=.5, T dtPL_=T(0.0016), T nuMin=T(2.9686e-3), T nuMax=T(3.1667), T smagoConst_=T(0.14));

  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);

private:

  /// Computes the local powerLaw relaxation parameter
  T computeOmegaPL(T omega0_, T rho_, T pi_[util::TensorVal<DESCRIPTOR >::n] );

private:
  T m;
  T n;
  T dtPL;
  T omegaMin;
  T omegaMax;
};

}

#endif
