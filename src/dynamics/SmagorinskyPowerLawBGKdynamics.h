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
#ifndef SMAGORINSKY_POWER_LAW_BGK_DYNAMICS_H
#define SMAGORINSKY_POWER_LAW_BGK_DYNAMICS_H

#include "core/cell.h"

namespace olb {

/// Implementation of the BGK collision step
template<typename T, typename DESCRIPTOR>
class SmagorinskyPowerLawBGKdynamics : public SmagorinskyBGKdynamics<T,DESCRIPTOR>, PowerLawDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  /// m,n...parameter in the power law model
  SmagorinskyPowerLawBGKdynamics(T omega, Momenta<T,DESCRIPTOR>& momenta, T m=0.1, T n=.5, T nuMin=T(2.9686e-3), T nuMax=T(3.1667), T smagoConst=T(0.14));
  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// Computes the local smagorinsky relaxation parameter with omega0 as an input
  virtual T computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell, T omega0);
protected:
  /// Computes squared norm of non-equilibrium part of 2nd momentum (for power-law)
  virtual T PiNeqNormSqr(Cell<T,DESCRIPTOR>& cell ) override;
};

/// Implementation of the ForcedBGK collision step
template<typename T, typename DESCRIPTOR>
class SmagorinskyPowerLawForcedBGKdynamics : public SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>, public PowerLawDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  /// m,n...parameter in the power law model
  SmagorinskyPowerLawForcedBGKdynamics(T omega, Momenta<T,DESCRIPTOR>& momenta, T m=0.1, T n=.5, T nuMin=T(2.9686e-3), T nuMax=T(3.1667), T smagoConst=T(0.14));
  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;

protected:
  /// Computes the local smagorinsky relaxation parameter with omega0 as an input
  virtual T computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell, T omega0);
  /// Computes squared norm of non-equilibrium part of 2nd momentum (for power-law)
  virtual T PiNeqNormSqr(Cell<T,DESCRIPTOR>& cell ) override;
};

}

#endif
