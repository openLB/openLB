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
 * Smagorinsky BGK Dynamics for porous media -- header file.
 */
#ifndef SMAGORINSKY_POROUS_PARTICLE_BGK_DYNAMICS_H
#define SMAGORINSKY_POROUS_PARTICLE_BGK_DYNAMICS_H

#include "dynamics/porousBGKdynamics.h"

namespace olb {

/// Implementation of the BGK collision step for a porosity model
template<typename T, typename DESCRIPTOR>
class SmagorinskyPorousParticleBGKdynamics : public PorousParticleBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  SmagorinskyPorousParticleBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_,
                         T dx_ = 1, T dt_ = 1);

  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);

  /// set relaxation parameter
  void setOmega(T omega_);
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,DESCRIPTOR>& cell_);


protected:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);
  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0, T preFactor_, T rho,
                 T pi[util::TensorVal<DESCRIPTOR >::n] );

  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;
  T dx;
  T dt;
};

} // olb

#endif
