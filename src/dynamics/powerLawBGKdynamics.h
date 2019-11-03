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
 * BGK Dynamics with adjusted omega -- header file.
 * Strain rate similar to "J.Boyd, J. Buick and S.Green: A second-order accurate lattice Boltzmann non-Newtonian flow model"
 * Power Law similar to "Huidan Yu, Sharath S. Girimaji, Li-Shi Luo - DNS and LES of decaying isotropic turbulence with and without frame rotation using lattice Boltzmann method"
 *
 */
#ifndef POWER_LAW_BGK_DYNAMICS_H
#define POWER_LAW_BGK_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

namespace olb {

/// Implementation of Power-Law Dynamics
template<typename T, typename DESCRIPTOR>
class PowerLawDynamics {
public:
  /// Constructor
  PowerLawDynamics(T m=0.1, T n=.5, T nuMin=T(2.9686e-3), T nuMax=T(3.1667));

protected:
  /// Computes the local power-Law relaxation parameter
  T computeOmegaPL ( Cell<T,DESCRIPTOR>& cell, T omega0,
           T rho, T pi[util::TensorVal<DESCRIPTOR >::n] );
  /// Computes squared norm of non-equilibrium part of 2nd momentum
  virtual T PiNeqNormSqr(Cell<T,DESCRIPTOR>& cell ) =0;

protected:
  T _m;
  T _n;
  T _omegaMin;
  T _omegaMax;
};

/// Implementation of the BGK collision step
template<typename T, typename DESCRIPTOR>
class PowerLawBGKdynamics : public BGKdynamics<T,DESCRIPTOR>, public PowerLawDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  /// m,n...parameter in the power law model
  PowerLawBGKdynamics(T omega, Momenta<T,DESCRIPTOR>& momenta, T m=0.1, T n=.5, T nuMin=T(2.9686e-3), T nuMax=T(3.1667));

  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;

protected:
  /// Computes squared norm of non-equilibrium part of 2nd momentum
  virtual T PiNeqNormSqr(Cell<T,DESCRIPTOR>& cell ) override;
};

/// Implementation of the forced BGK collision step
template<typename T, typename DESCRIPTOR>
class PowerLawForcedBGKdynamics : public ForcedBGKdynamics<T,DESCRIPTOR>, public PowerLawDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  /// m,n...parameter in the power law model
  PowerLawForcedBGKdynamics(T omega, Momenta<T,DESCRIPTOR>& momenta, T m=0.1, T n=.5, T nuMin=T(2.9686e-3), T nuMax=T(3.1667));

  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;

protected:
  /// Computes squared norm of non-equilibrium part of 2nd momentum
  virtual T PiNeqNormSqr(Cell<T,DESCRIPTOR>& cell ) override;
};

}

#endif
